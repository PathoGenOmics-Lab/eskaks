use super::*;
use std::io::Write;

fn write_temp_vcf(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn sample_names_are_the_header_columns_in_column_order() {
    // The one place eskaks reads sample NAMES rather than column positions, so a tree's
    // tips can be matched to the cohort. Order matters absolutely: position i in this
    // list is the sample index every carrier set was built with.
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tERR_9\tERR_1\tERR_5
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\tGT\t1\t0\t1\n";
    let f = write_temp_vcf(vcf);
    assert_eq!(sample_names(f.path()).unwrap(), vec!["ERR_9", "ERR_1", "ERR_5"]);
    assert_eq!(sample_count(f.path()).unwrap(), 3, "count and names must agree");
    // And the carrier set is indexed by that same order: ERR_9 and ERR_5 carry the ALT.
    let snps = parse_vcf(f.path()).unwrap();
    let carriers = &snps[0].carriers.as_ref().expect("genotyped")[0];
    assert_eq!(carriers.samples().collect::<Vec<_>>(), vec![0, 2]);

    // An AF-only VCF has no sample columns, so there are no names to match at all.
    let af_only = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50;AF=0.5\n";
    let g = write_temp_vcf(af_only);
    assert!(sample_names(g.path()).unwrap().is_empty());
    assert_eq!(sample_count(g.path()).unwrap(), 0);
}

#[test]
fn parse_simple_snps() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50;AF=0.5
chr1\t200\t.\tC\tT\t30\tPASS\tDP=60;AF=0.3
chr1\t300\t.\tAT\tA\t30\tPASS\tDP=40\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps.len(), 2); // indel skipped
    assert_eq!(snps[0].pos, 100);
    assert_eq!(snps[0].ref_allele, b'A');
    assert_eq!(snps[0].alt_alleles, vec![b'G']);
    assert!((snps[0].alt_freqs[0] - 0.5).abs() < 1e-6);
    assert_eq!(snps[0].depth, Some(50));
}

#[test]
fn gt_counts_capture_genotype_truth_over_disagreeing_info_af() {
    // A record whose INFO/AF (0.50) disagrees with its genotypes: only 1 of 4 haploid
    // samples carries the ALT (true derived count = 1). The frequency path keeps the
    // declared INFO/AF (used for pN/pS), but the diversity path must see the exact GT
    // count, not round(0.50 * 4) = 2.
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t100\t.\tA\tG\t30\tPASS\tAF=0.50\tGT\t1\t0\t0\t0\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps.len(), 1);
    // Frequency path: unchanged, still the declared INFO/AF.
    assert!((snps[0].alt_freqs[0] - 0.50).abs() < 1e-9, "alt_freqs should keep INFO/AF");
    // Diversity path: the exact GT-derived count, independent of the disagreeing AF.
    let gc = snps[0].gt_counts.as_ref().expect("gt_counts present when GT columns exist");
    assert_eq!(gc.alt, vec![1], "GT-derived count must be 1, not the AF-implied 2");
    assert_eq!(gc.called, 4, "all four haploid samples were called");
    // And WHICH sample: the same-codon check intersects these sets, so a carrier credited
    // to the wrong column would pair this allele with the wrong neighbours.
    let carriers = snps[0].carriers.as_ref().expect("carriers present when GT columns exist");
    assert_eq!(carriers.len(), 1, "one carrier set per valid ALT");
    assert_eq!(carriers[0].len(), 1);
    assert!(carriers[0].contains(0), "S1 is the only carrier");
    assert!((1..4).all(|i| !carriers[0].contains(i)), "S2..S4 carry the reference");
}

#[test]
fn carriers_map_each_alt_to_the_samples_that_actually_carry_it() {
    // Two things a same-codon intersection depends on, and both are index arithmetic that
    // a plausible off-by-one would silently get wrong: the ALT index (`carriers[i]` must
    // be the samples carrying `alt_alleles[i]`, with REF at index 0 skipped) and the
    // sample index (bit `j` must be column `j`). A multiallelic record with a no-call and
    // a declared-but-invalid ALT pins both, since the invalid ALT shifts the ALT indices
    // apart from the GT allele numbers.
    //
    //   ALT list:  1 = G (valid), 2 = GG (not a SNP, dropped), 3 = C (valid)
    //   S1 = 3 (C), S2 = 0 (REF), S3 = 1 (G), S4 = . (no call), S5 = 3 (C)
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5
chr1\t100\t.\tA\tG,GG,C\t30\tPASS\tAF=0.2,0.1,0.4\tGT\t3\t0\t1\t.\t3\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps[0].alt_alleles, vec![b'G', b'C'], "the multi-base ALT is not a SNP");

    let carriers = snps[0].carriers.as_ref().expect("carriers");
    assert_eq!(carriers.len(), 2, "one set per VALID ALT, aligned with alt_alleles");
    // G is GT allele 1, carried by S3 alone.
    assert_eq!(carriers[0].len(), 1);
    assert!(carriers[0].contains(2), "G's carrier is column 2 (S3)");
    // C is GT allele 3, carried by S1 and S5. Its set must NOT have slid to allele 2's.
    assert_eq!(carriers[1].len(), 2);
    assert!(carriers[1].contains(0) && carriers[1].contains(4), "C's carriers are S1 and S5");
    assert!(!carriers[1].contains(1), "S2 carries the reference");
    assert!(!carriers[1].contains(3), "a no-call carries nothing");
    // The two ALTs exclude one another at one position, as haploid genotypes must.
    assert!(!carriers[0].intersects(&carriers[1]));
    // The counts they sit beside are unchanged and stay aligned with them.
    let gc = snps[0].gt_counts.as_ref().expect("gt_counts");
    assert_eq!(gc.alt, vec![1, 2]);
    assert_eq!(gc.called, 4, "the no-call is not a called allele");
}

#[test]
fn parse_multi_allelic() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG,C\t30\tPASS\tAF=0.3,0.2\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps.len(), 1);
    assert_eq!(snps[0].alt_alleles, vec![b'G', b'C']);
    assert_eq!(snps[0].alt_freqs.len(), 2);
}

#[test]
fn af_missing_token_keeps_positional_alignment() {
    // Regression: a missing AF token (".") must not shift the remaining
    // frequencies onto the wrong ALT. G has no AF (0.0), C=0.2, T=0.3.
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG,C,T\t30\tPASS\tAF=.,0.2,0.3\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps[0].alt_alleles, vec![b'G', b'C', b'T']);
    assert!((snps[0].alt_freqs[0] - 0.0).abs() < 1e-9, "G should be 0.0, got {}", snps[0].alt_freqs[0]);
    assert!((snps[0].alt_freqs[1] - 0.2).abs() < 1e-9, "C should be 0.2, got {}", snps[0].alt_freqs[1]);
    assert!((snps[0].alt_freqs[2] - 0.3).abs() < 1e-9, "T should be 0.3, got {}", snps[0].alt_freqs[2]);
}

#[test]
fn filter_pass_only() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tAF=0.5
chr1\t200\t.\tC\tT\t30\tLowQual\tAF=0.3\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    let filtered = filter_snps(snps, true, None, None, None);
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].pos, 100);
}

#[test]
fn filter_max_af_excludes_fixed() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tAF=0.3
chr1\t200\t.\tC\tT\t30\tPASS\tAF=1.0
chr1\t300\t.\tG\tA\t30\tPASS\tAF=0.95\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    let filtered = filter_snps(snps, false, None, Some(0.99), None);
    assert_eq!(filtered.len(), 2, "AF=1.0 should be excluded");
    assert_eq!(filtered[0].pos, 100);
    assert_eq!(filtered[1].pos, 300);
}

// ---- GT-based allele-frequency computation (no INFO/AF) ---------------

#[test]
fn af_computed_from_gt_when_info_af_absent() {
    // 4 diploid samples, biallelic A->G, no INFO/AF => AF is derived from GT.
    // Row 100 (0/0,0/1,1/1,1/1): ALT alleles = 1+2+2 = 5 of 8 total => 0.625.
    // Row 200 (0|1, . , 1|1, .|.): phased split, "." skipped => 3 of 4 => 0.75.
    // Row 300 (all missing): total_alleles == 0 => AF 0.0 (no divide-by-zero).
    // Row 400 (FORMAT has no GT key): falls back to assumed-fixed 1.0.
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\tGT\t0/0\t0/1\t1/1\t1/1
chr1\t200\t.\tC\tT\t30\tPASS\tDP=40\tGT\t0|1\t.\t1|1\t.|.
chr1\t300\t.\tG\tA\t30\tPASS\tDP=10\tGT\t.\t.\t./.\t.
chr1\t400\t.\tT\tC\t30\tPASS\tDP=10\tDP\t5\t6\t7\t8\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps.len(), 4);
    assert!((snps[0].alt_freqs[0] - 0.625).abs() < 1e-9, "row100 GT AF, got {}", snps[0].alt_freqs[0]);
    assert!((snps[1].alt_freqs[0] - 0.75).abs() < 1e-9, "row200 phased GT AF, got {}", snps[1].alt_freqs[0]);
    assert!((snps[2].alt_freqs[0] - 0.0).abs() < 1e-9, "row300 all-missing GT => 0.0, got {}", snps[2].alt_freqs[0]);
    assert!((snps[3].alt_freqs[0] - 1.0).abs() < 1e-9, "row400 no GT key => assumed fixed 1.0, got {}", snps[3].alt_freqs[0]);
}

// ---- parse_vcf robustness against malformed / non-SNP lines ----------

#[test]
fn parse_skips_non_snp_and_invalid_alleles() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tN\tG\t30\tPASS\tAF=0.5
chr1\t150\t.\tA\t<DEL>\t30\tPASS\tAF=0.5
chr1\t200\t.\tA\tA\t30\tPASS\tAF=0.5
chr1\t250\t.\tA\tG,<INS>\t30\tPASS\tAF=0.4,0.1
chr1\t300\t.\tA\tG
chr1\t350\t.\tC\tT\t30\tPASS\tAF=0.6\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    // Kept: 250 (G only, <INS> dropped) and 350. Skipped: N-ref, symbolic-only
    // ALT, ALT==REF, and the <8-field line.
    assert_eq!(snps.len(), 2, "got {:?}", snps.iter().map(|s| s.pos).collect::<Vec<_>>());
    assert_eq!(snps[0].pos, 250);
    assert_eq!(snps[0].alt_alleles, vec![b'G'], "symbolic ALT must be dropped, keeping G");
    assert_eq!(snps[1].pos, 350);
}

#[test]
fn parse_normalizes_lowercase_bases() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\ta\tg\t30\tPASS\tAF=0.5\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps.len(), 1);
    assert_eq!(snps[0].ref_allele, b'A');
    assert_eq!(snps[0].alt_alleles, vec![b'G']);
}

#[test]
fn parse_missing_dp_yields_none_depth() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tAF=0.5\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps[0].depth, None);
}

#[test]
fn parse_empty_vcf_is_ok_not_error() {
    // An empty (header-only) VCF is not fatal on its own — merge decides.
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert!(snps.is_empty());
}

// ---- filter_snps: depth + min_af paths -------------------------------

#[test]
fn filter_min_depth_and_min_af() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50;AF=0.3
chr1\t200\t.\tC\tT\t30\tPASS\tDP=50;AF=0.02
chr1\t300\t.\tG\tA\t30\tPASS\tAF=0.5\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    let filtered = filter_snps(snps, false, Some(0.05), None, Some(10));
    // 100 kept; 200 dropped by min_af; 300 dropped by missing DP under min_depth.
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].pos, 100);
}

#[test]
fn filter_dot_filter_treated_as_pass() {
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\t.\tAF=0.5\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    let filtered = filter_snps(snps, true, None, None, None);
    assert_eq!(filtered.len(), 1, "FILTER '.' must pass under pass_only");
}

// ---- merge_vcfs ------------------------------------------------------

fn merge_paths(handles: &[tempfile::NamedTempFile]) -> Vec<String> {
    handles.iter().map(|h| h.path().to_str().unwrap().to_string()).collect()
}

#[test]
fn merge_two_samples_computes_af_and_depth() {
    let a = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n");
    let b = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=30
chr1\t200\t.\tC\tT\t30\tPASS\tDP=40\n");
    let files = [a, b]; // keep the temp files alive for the duration of the merge
    let paths = merge_paths(&files);
    let merged = merge_vcfs(&paths, false, None).unwrap();
    assert_eq!(merged.len(), 2);
    // pos 100 carried by both samples => AF = 2/2 = 1.0, depth = (50+30)/2 = 40.
    assert_eq!(merged[0].pos, 100);
    assert_eq!(merged[0].ref_allele, b'A');
    assert_eq!(merged[0].alt_alleles, vec![b'G']);
    assert!((merged[0].alt_freqs[0] - 1.0).abs() < 1e-9);
    assert_eq!(merged[0].depth, Some(40));
    // pos 200 carried by one of two samples => AF = 1/2 = 0.5.
    assert_eq!(merged[1].pos, 200);
    assert!((merged[1].alt_freqs[0] - 0.5).abs() < 1e-9);
    assert_eq!(merged[1].depth, Some(40));
}

#[test]
fn merge_is_deterministic_across_runs() {
    let a = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr2\t500\t.\tG\tA\t30\tPASS\tDP=20
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n");
    let b = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t30\tPASS\tDP=30\n");
    let files = [a, b];
    let paths = merge_paths(&files);
    let key = |m: &[VcfSnp]| {
        m.iter()
            .map(|s| (s.chrom.clone(), s.pos, s.alt_alleles.clone(), s.alt_freqs.clone()))
            .collect::<Vec<_>>()
    };
    let r1 = merge_vcfs(&paths, false, None).unwrap();
    let r2 = merge_vcfs(&paths, false, None).unwrap();
    assert_eq!(key(&r1), key(&r2), "merge output must be byte-stable across runs");
    // Sorted by (chrom, pos): chr1:100 before chr2:500.
    assert_eq!(r1[0].chrom, "chr1");
    assert_eq!(r1[0].pos, 100);
    assert_eq!(r1[1].chrom, "chr2");
    // The multi-allelic chr1:100 (G from sample A, T from sample B) must emit its
    // ALTs in a stable, base-sorted order — not HashMap-iteration order.
    assert_eq!(r1[0].alt_alleles, vec![b'G', b'T'], "ALTs must be sorted deterministically");
}

#[test]
fn merge_drops_alt_equal_to_merged_ref() {
    // Samples disagree on REF at the same locus. The first-seen REF (A) wins;
    // sample B's ALT=A equals that merged REF and must be dropped, leaving G.
    let a = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\n");
    let b = write_temp_vcf("\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tG\tA\t30\tPASS\tDP=30\n");
    let files = [a, b];
    let paths = merge_paths(&files);
    let merged = merge_vcfs(&paths, false, None).unwrap();
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].ref_allele, b'A');
    assert_eq!(merged[0].alt_alleles, vec![b'G'], "ALT equal to merged REF must be dropped");
    assert!((merged[0].alt_freqs[0] - 0.5).abs() < 1e-9);
}

#[test]
fn merge_all_empty_is_error() {
    let a = write_temp_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    let b = write_temp_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    let files = [a, b];
    let paths = merge_paths(&files);
    assert!(merge_vcfs(&paths, false, None).is_err(), "merging only empty VCFs must error");
}


// ---- regression: malformed-input hardening ----

#[test]
fn af_nonfinite_or_out_of_range_is_rejected() {
    // "nan"/"2.0" AF tokens are invalid; they must become missing (0.0), not enter
    // alt_freqs where a NaN evades --min-af/--max-af and poisons genome-wide πN/πS.
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG,C,T\t30\tPASS\tAF=nan,2.0,0.3\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert_eq!(snps[0].alt_freqs, vec![0.0, 0.0, 0.3]);
    assert!(snps[0].alt_freqs.iter().all(|a| a.is_finite()));
}

#[test]
fn gt_out_of_range_allele_index_does_not_deflate_af() {
    // S2's GT references allele 9 (undeclared): a malformed index must be ignored,
    // not counted into total_alleles. AF(alt) = 1/2 = 0.5, not 1/4.
    let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tG\t30\tPASS\tDP=50\tGT\t0/1\t9/9\n";
    let f = write_temp_vcf(vcf);
    let snps = parse_vcf(f.path()).unwrap();
    assert!((snps[0].alt_freqs[0] - 0.5).abs() < 1e-9, "got {}", snps[0].alt_freqs[0]);
}

#[test]
fn merge_large_depth_does_not_overflow() {
    // Summed depth (4e9 + 4e9 = 8e9) exceeds u32::MAX; the u64 accumulator must not
    // panic (debug) or wrap (release). Average 4e9 fits back in u32.
    let big = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=4000000000\n";
    let files = [write_temp_vcf(big), write_temp_vcf(big)];
    let paths = merge_paths(&files);
    let merged = merge_vcfs(&paths, false, None).unwrap();
    assert_eq!(merged[0].depth, Some(4_000_000_000));
}

#[test]
fn merge_ignores_non_carrier_records_when_computing_frequency() {
    // A per-sample VCF that RECORDS a non-carrier site (GT 0, as gVCF / all-sites callers
    // do) must not be counted as a carrier. Sample 1 carries G; sample 2 is homozygous
    // REF at the same site: merged AF(G) = 1/2, not the inflated 2/2 that would drop this
    // real polymorphism as "fixed" in the diversity / MK / AF-filtered paths.
    let carrier = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t1\n";
    let non_carrier = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0\n";
    let files = [write_temp_vcf(carrier), write_temp_vcf(non_carrier)];
    let paths = merge_paths(&files);
    let merged = merge_vcfs(&paths, false, None).unwrap();
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].alt_alleles, vec![b'G']);
    assert!(
        (merged[0].alt_freqs[0] - 0.5).abs() < 1e-9,
        "AF should be 0.5 (1 of 2 carriers), got {}",
        merged[0].alt_freqs[0]
    );
}

#[test]
fn parse_rejects_empty_or_headerless_file_but_accepts_variant_free() {
    // A truly empty / header-less file is not a usable VCF and must error, so a merge does
    // not count it as a phantom sample that deflates every allele frequency. A valid
    // variant-free VCF (has the #CHROM header, no data rows) stays accepted as 0 SNPs.
    assert!(parse_vcf(write_temp_vcf("").path()).is_err(), "0-byte file should error");
    assert!(
        parse_vcf(write_temp_vcf("##fileformat=VCFv4.2\n").path()).is_err(),
        "meta-only file should error"
    );
    let header_only = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let snps = parse_vcf(write_temp_vcf(header_only).path()).unwrap();
    assert!(snps.is_empty(), "a variant-free VCF is a valid 0-SNP sample");
}


#[test]
fn af_filter_is_per_allele_at_multiallelic_sites() {
    // One multi-allelic record G@0.05, C@0.5, T@0.995. --min-af 0.1 --max-af 0.99
    // must keep only C (0.5), pruning the sub-threshold G and the fixed T, rather than
    // dropping the whole record or leaking the out-of-range alleles.
    let snp = VcfSnp {
        chrom: "chr1".into(), pos: 100, ref_allele: b'A',
        alt_alleles: vec![b'G', b'C', b'T'],
        alt_freqs: vec![0.05, 0.5, 0.995],
        gt_counts: None,
        carriers: None,
        filter: "PASS".into(), depth: None,
    };
    let out = filter_snps(vec![snp], false, Some(0.1), Some(0.99), None);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].alt_alleles, vec![b'C']);
    assert_eq!(out[0].alt_freqs, vec![0.5]);
}


#[test]
fn af_filter_keeps_gt_counts_aligned_at_multiallelic_sites() {
    // gt_counts.alt and carriers are both parallel to alt_alleles. Pruning a NON-last ALT
    // by frequency must prune them in lockstep; otherwise the diversity path (which reads
    // gt_counts by ALT index) reads the wrong survivor's derived-allele count and
    // silently corrupts piN/piS/theta_W/Tajima's D at multi-allelic sites, and the
    // same-codon check intersects the wrong allele's samples.
    let carriers = |samples: &[usize]| {
        let mut s = CarrierSet::new(9);
        for &i in samples {
            s.insert(i);
        }
        s
    };
    let snp = VcfSnp {
        chrom: "chr1".into(), pos: 100, ref_allele: b'A',
        alt_alleles: vec![b'G', b'C', b'T'],
        alt_freqs: vec![0.0, 0.5, 0.995], // G phantom (0.0), C kept (0.5), T fixed (0.995)
        gt_counts: Some(GtCounts { alt: vec![0, 2, 7], called: 9 }),
        carriers: Some(vec![
            carriers(&[]),
            carriers(&[3, 4]),
            carriers(&[0, 1, 2, 5, 6, 7, 8]),
        ]),
        filter: "PASS".into(), depth: None,
    };
    let out = filter_snps(vec![snp], false, Some(0.1), Some(0.99), None);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].alt_alleles, vec![b'C']);
    // The survivor C must carry ITS genotype count (2), not G's (0) or T's (7).
    assert_eq!(out[0].gt_counts.as_ref().unwrap().alt, vec![2]);
    // ...and ITS carriers (samples 3 and 4), not T's seven.
    let kept = out[0].carriers.as_ref().expect("carriers survive filtering");
    assert_eq!(kept.len(), 1, "one surviving ALT, one carrier set");
    assert_eq!(kept[0], carriers(&[3, 4]));
}


#[test]
fn fasta_content_as_vcf_reports_wrong_format() {
    let f = write_temp_vcf(">strain_A\nATGGCTGCT\n>strain_B\nATGGCTGCT\n");
    let err = parse_vcf(f.path()).unwrap_err().to_string();
    assert!(err.contains("valid VCF records"), "err: {}", err);
    assert!(err.contains("FASTA"), "should hint FASTA, err: {}", err);
}
