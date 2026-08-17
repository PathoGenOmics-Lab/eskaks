//! Same-codon detection, pinned against the bundled 20-sample toy VCF.
//!
//! eskaks scores every ALT against the REFERENCE codon, so two SNPs in one codon are
//! never evaluated jointly and the amino-acid change reported for such a codon can be
//! one no genome carries. Detecting those codons used to rely on an allele-frequency
//! pigeonhole bound, which is a floor: on `examples/toy_genome/variants_multisample.vcf`
//! it reached 8 of the 13 codons where a sample really does carry two alleles, missing
//! 38%. With the per-sample carriers retained it is a set intersection instead, and
//! exact in both directions.
//!
//! This file computes the truth from the VCF text by a route that shares no code with
//! eskaks (columns and coordinates, no parser, no codon table), then asserts eskaks
//! agrees. It also strips the genotype columns from the same file and checks the
//! frequency bound still behaves exactly as it did, since that remains the right
//! answer for an AF-only input.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use eskaks::vcf_analysis::{compute_pn_ps, parse_reference_fasta, ComputeDiagnostics};

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

const REF: &str = "examples/toy_genome/reference.fasta";
const GFF: &str = "examples/toy_genome/genes.gff3";
const VCF: &str = "examples/toy_genome/variants_multisample.vcf";

/// One coding SNP position as the oracle sees it: which samples carry a non-reference
/// base there, and the highest ALT frequency INFO/AF declares for it.
struct OraclePos {
    carriers: BTreeSet<usize>,
    max_af: f64,
}

/// Group every VCF position that falls in a CDS by (gene start, codon index), reading
/// the files as text.
///
/// This is only as simple as it is because the toy annotation is: twelve single-exon,
/// plus-strand, phase-0 CDS features, so the codon of a position is `(pos - start) / 3`
/// with no exon stitching and no reverse complement. `eskaks_diagnostics` asserts the
/// REF-allele check passed for every SNP, which is what makes skipping the reference
/// comparison here legitimate.
fn oracle_codons() -> BTreeMap<(usize, usize), Vec<OraclePos>> {
    let gff = std::fs::read_to_string(manifest(GFF)).expect("read gff");
    let mut cds: Vec<(usize, usize)> = Vec::new(); // (start, end), 1-based inclusive
    for line in gff.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() > 7 && f[2] == "CDS" {
            assert_eq!(f[6], "+", "the oracle assumes plus-strand toy genes");
            assert_eq!(f[7], "0", "the oracle assumes phase 0");
            cds.push((f[3].parse().unwrap(), f[4].parse().unwrap()));
        }
    }
    assert_eq!(cds.len(), 12, "the toy annotation has 12 CDS features");

    let vcf = std::fs::read_to_string(manifest(VCF)).expect("read vcf");
    let mut by_codon: BTreeMap<(usize, usize), Vec<OraclePos>> = BTreeMap::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let pos: usize = f[1].parse().expect("POS");
        // INFO/AF, one value per ALT; the maximum is the summary the frequency bound uses.
        let max_af = f[7]
            .split(';')
            .find_map(|e| e.strip_prefix("AF="))
            .map(|v| v.split(',').filter_map(|x| x.parse::<f64>().ok()).fold(0.0, f64::max))
            .expect("every toy record declares AF");
        // Haploid GT columns start at field 9; anything but "0" or "." is a carrier.
        let carriers: BTreeSet<usize> = f[9..]
            .iter()
            .enumerate()
            .filter(|(_, gt)| **gt != "0" && **gt != ".")
            .map(|(i, _)| i)
            .collect();
        for &(start, end) in &cds {
            if pos >= start && pos <= end {
                let codon = (pos - start) / 3;
                by_codon
                    .entry((start, codon))
                    .or_default()
                    .push(OraclePos { carriers: carriers.clone(), max_af });
            }
        }
    }
    by_codon
}

/// Codons carrying more than one SNP position, split into (a sample carries two of
/// them, the allele-frequency bound alone would have caught it).
fn oracle_counts() -> (usize, usize, usize) {
    let mut multi = 0usize;
    let mut shared = 0usize;
    let mut by_af_bound = 0usize;
    for group in oracle_codons().values() {
        if group.len() < 2 {
            continue;
        }
        multi += 1;
        let any_pair = (0..group.len()).any(|i| {
            ((i + 1)..group.len())
                .any(|j| group[i].carriers.intersection(&group[j].carriers).next().is_some())
        });
        if any_pair {
            shared += 1;
        }
        // The old test: sum of frequencies above k - 1 (see `forced_cooccurrence`).
        let sum: f64 = group.iter().map(|p| p.max_af).sum();
        if sum > group.len() as f64 - 1.0 {
            by_af_bound += 1;
        }
    }
    (multi, shared, by_af_bound)
}

/// Run the real pN/pS path over the toy genome with the given VCF.
fn eskaks_diagnostics(vcf: &Path) -> ComputeDiagnostics {
    let reference = parse_reference_fasta(&manifest(REF)).expect("reference");
    let genes = eskaks::gff::parse_gff3(&manifest(GFF)).expect("gff");
    let snps = eskaks::vcf::parse_vcf(vcf).expect("vcf");
    let gc = eskaks::genetic_code::get_table(11).expect("bacterial code");
    let (_results, diag) = compute_pn_ps(&reference, &genes, &snps, gc, false, 1.0, 0.99);
    assert_eq!(
        diag.ref_mismatch, 0,
        "every toy SNP must agree with the reference, or the oracle's coordinate \
         arithmetic is not comparable with eskaks's"
    );
    diag
}

/// The same VCF with its FORMAT and sample columns removed: an AF-only record set, the
/// input for which the frequency bound is the only sound inference.
fn af_only_copy() -> tempfile::NamedTempFile {
    let vcf = std::fs::read_to_string(manifest(VCF)).expect("read vcf");
    let mut out = tempfile::Builder::new().suffix(".vcf").tempfile().expect("tempfile");
    for line in vcf.lines() {
        let stripped = if line.starts_with("##") {
            line.to_string()
        } else {
            line.split('\t').take(8).collect::<Vec<_>>().join("\t")
        };
        writeln!(out, "{stripped}").expect("write");
    }
    out.flush().expect("flush");
    out
}

#[test]
fn every_codon_a_sample_carries_twice_is_found() {
    let (multi, shared, by_af_bound) = oracle_counts();

    // The numbers this whole change exists for, read off the file itself.
    assert_eq!(multi, 13, "codons of the toy VCF carrying more than one SNP position");
    assert_eq!(shared, 13, "codons where a sample really carries two of those alleles");
    assert_eq!(
        by_af_bound, 8,
        "the allele-frequency bound reaches only 8 of them, which is why carriers are kept"
    );

    let diag = eskaks_diagnostics(&manifest(VCF));
    assert_eq!(diag.multi_snp_codons, multi);
    assert_eq!(
        diag.cooccurring_codons, shared,
        "with genotype columns present the detector must find all {shared}, not the {by_af_bound} \
         the frequency bound reaches"
    );
    assert_eq!(diag.cooccurring_exact, shared, "all of them established from carriers");
    assert_eq!(
        diag.codons_with_carriers, multi,
        "every multi-SNP codon had carriers, so its answer is exact in both directions"
    );
}

#[test]
fn stripping_the_genotypes_falls_back_to_the_unweakened_frequency_bound() {
    // Removing the sample columns removes the only thing that made the answer exact.
    // The bound must then behave exactly as it always did (8 codons on this file):
    // for an AF-only VCF it is the best available inference and must not be softened.
    let af_only = af_only_copy();
    let diag = eskaks_diagnostics(af_only.path());
    let (multi, _shared, by_af_bound) = oracle_counts();

    assert_eq!(diag.multi_snp_codons, multi, "the same codons still carry two SNPs");
    assert_eq!(diag.codons_with_carriers, 0, "no genotypes, so nothing is checked exactly");
    assert_eq!(diag.cooccurring_exact, 0);
    assert_eq!(
        diag.cooccurring_codons, by_af_bound,
        "the AF-only path must keep reporting exactly the {by_af_bound} the bound proves"
    );
}

/// Explode the bundled multi-sample VCF into one single-sample file per column, the
/// `--vcf-list` shape: same data, a completely different code path into the analysis.
fn per_sample_copies(dir: &Path) -> Vec<String> {
    let vcf = std::fs::read_to_string(manifest(VCF)).expect("read vcf");
    let lines: Vec<&str> = vcf.lines().collect();
    let meta: Vec<&str> = lines.iter().copied().filter(|l| l.starts_with("##")).collect();
    let header: Vec<&str> =
        lines.iter().find(|l| l.starts_with("#CHROM")).expect("header").split('\t').collect();
    let records: Vec<Vec<&str>> =
        lines.iter().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();

    (9..header.len())
        .map(|col| {
            let path = dir.join(format!("{}.vcf", header[col]));
            let mut out = String::new();
            for m in &meta {
                out.push_str(m);
                out.push('\n');
            }
            out.push_str(&header[..9].join("\t"));
            out.push('\t');
            out.push_str(header[col]);
            out.push('\n');
            for r in &records {
                if r[col] != "0" && r[col] != "." {
                    out.push_str(&r[..9].join("\t"));
                    out.push_str("\t1\n");
                }
            }
            std::fs::write(&path, out).expect("write per-sample vcf");
            path.to_str().expect("utf-8 path").to_string()
        })
        .collect()
}

#[test]
fn one_vcf_per_sample_sees_the_same_co_occurrences() {
    // `--vcf-list` reduces each file to a count in the merge, but the file index IS the
    // sample index, so co-occurrence is directly observable there too. Splitting the
    // bundled VCF one column per file must therefore reach the same 13 codons, or the
    // two supported inputs disagree about the same data.
    let dir = tempfile::tempdir().expect("tempdir");
    let paths = per_sample_copies(dir.path());
    assert_eq!(paths.len(), 20, "the toy VCF has 20 sample columns");

    let reference = parse_reference_fasta(&manifest(REF)).expect("reference");
    let genes = eskaks::gff::parse_gff3(&manifest(GFF)).expect("gff");
    let merged = eskaks::vcf::merge_vcfs(&paths, false, None).expect("merge");
    let gc = eskaks::genetic_code::get_table(11).expect("bacterial code");
    let (results, diag) = compute_pn_ps(&reference, &genes, &merged, gc, false, 1.0, 0.99);

    let (multi, shared, by_af_bound) = oracle_counts();
    assert_eq!(diag.multi_snp_codons, multi);
    assert_eq!(diag.codons_with_carriers, multi, "the merge keeps per-sample identity");
    assert_eq!(diag.cooccurring_exact, shared);
    assert_eq!(
        diag.cooccurring_codons, shared,
        "the merged path must find all {shared}, not the {by_af_bound} the bound reaches"
    );

    // And row by row: the marks must match the single multi-sample VCF exactly.
    let single = {
        let snps = eskaks::vcf::parse_vcf(&manifest(VCF)).expect("vcf");
        let (r, _) = compute_pn_ps(&reference, &genes, &snps, gc, false, 1.0, 0.99);
        r
    };
    let marks = |rs: &[eskaks::vcf_analysis::GenePnPs]| -> Vec<(usize, u8, Option<bool>)> {
        let mut v: Vec<_> = rs
            .iter()
            .flat_map(|g| g.variants.iter().map(|x| (x.pos, x.alt_allele, x.codon_shared)))
            .collect();
        v.sort();
        v
    };
    assert_eq!(marks(&results), marks(&single), "the two input paths must agree row for row");
}

/// A single-sample VCF carrying the given (POS, REF, ALT) records, all at AF 1.0.
fn one_sample_vcf(dir: &Path, name: &str, records: &[(usize, char, char)]) -> String {
    let path = dir.join(format!("{name}.vcf"));
    let mut out = String::from("##fileformat=VCFv4.2\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (haploid)\">\n");
    out.push_str(&format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n"));
    for (pos, r, a) in records {
        out.push_str(&format!("chr1\t{pos}\t.\t{r}\t{a}\t60\tPASS\tDP=30;AF=1.0\tGT\t1\n"));
    }
    std::fs::write(&path, out).expect("write vcf");
    path.to_str().expect("utf-8 path").to_string()
}

#[test]
fn the_merge_credits_each_allele_to_the_file_it_came_from() {
    // The merge reduces each file to a count, and a carrier set recorded against the
    // wrong index would still intersect and still look right on a cohort where every
    // sample carries something. So: two codons of gene01, both with two SNPs, both at
    // merged AF 0.5 + 0.5 = 1.0 so the frequency bound stays silent for both. The only
    // thing separating them is WHICH file each allele came from.
    //
    //   codon 2  (chr1:4 G>A in file A, chr1:6 C>T in file B): different genomes, so
    //            scoring each against the reference codon GCC is correct.
    //   codon 10 (chr1:28 C>T and chr1:30 T>A, BOTH in file A): one genome carries CTT
    //            to TTA, Leu to Leu, one synonymous change, reported as missense +
    //            synonymous.
    let dir = tempfile::tempdir().expect("tempdir");
    let paths = vec![
        one_sample_vcf(dir.path(), "sampleA", &[(4, 'G', 'A'), (28, 'C', 'T'), (30, 'T', 'A')]),
        one_sample_vcf(dir.path(), "sampleB", &[(6, 'C', 'T')]),
    ];

    let reference = parse_reference_fasta(&manifest(REF)).expect("reference");
    let genes = eskaks::gff::parse_gff3(&manifest(GFF)).expect("gff");
    let merged = eskaks::vcf::merge_vcfs(&paths, false, None).expect("merge");
    let gc = eskaks::genetic_code::get_table(11).expect("bacterial code");
    let (results, diag) = compute_pn_ps(&reference, &genes, &merged, gc, false, 1.0, 0.99);

    assert_eq!(diag.ref_mismatch, 0, "all four REF alleles must match the toy reference");
    assert_eq!(diag.multi_snp_codons, 2, "codon 2 and codon 10 each carry two SNPs");
    assert_eq!(diag.codons_with_carriers, 2, "the merge keeps per-sample identity for both");
    assert_eq!(
        diag.cooccurring_codons, 1,
        "only codon 10 has both alleles in one file; crediting an allele to the wrong \
         sample index would flag codon 2 as well"
    );
    assert_eq!(diag.cooccurring_exact, 1);

    // Row by row: the two alleles of codon 10 are marked, the two of codon 2 are not.
    let mut marks: Vec<(usize, Option<bool>)> = results
        .iter()
        .flat_map(|g| g.variants.iter().map(|v| (v.pos, v.codon_shared)))
        .collect();
    marks.sort();
    assert_eq!(
        marks,
        vec![
            (4, Some(false)),
            (6, Some(false)),
            (28, Some(true)),
            (30, Some(true)),
        ]
    );
}

#[test]
fn the_variants_table_marks_the_rows_that_share_a_codon() {
    // The mark is per ALLELE and lands on real rows: `--variants --shared-codons` must
    // flag every row of the 13 shared codons and nothing else.
    let reference = parse_reference_fasta(&manifest(REF)).expect("reference");
    let genes = eskaks::gff::parse_gff3(&manifest(GFF)).expect("gff");
    let snps = eskaks::vcf::parse_vcf(&manifest(VCF)).expect("vcf");
    let gc = eskaks::genetic_code::get_table(11).expect("bacterial code");
    let (results, diag) = compute_pn_ps(&reference, &genes, &snps, gc, false, 1.0, 0.99);

    let rows: Vec<_> = results.iter().flat_map(|g| g.variants.iter()).collect();
    assert_eq!(rows.len(), 264, "the toy variants table has 264 rows");
    assert!(
        rows.iter().all(|v| v.codon_shared.is_some()),
        "with genotypes present no row may be left unknowable"
    );

    // Every marked row must sit on a codon carrying more than one SNP position, and the
    // marked rows must cover exactly the codons the oracle found.
    let flagged: BTreeSet<(usize, usize)> = results
        .iter()
        .flat_map(|g| g.variants.iter().map(move |v| (g.genome_start, v.aa_pos)))
        .zip(rows.iter())
        .filter(|(_, v)| v.codon_shared == Some(true))
        .map(|(key, _)| key)
        .collect();
    assert_eq!(
        flagged.len(),
        diag.cooccurring_codons,
        "one distinct codon per co-occurring codon, no more and no fewer"
    );
    assert_eq!(flagged.len(), 13);
}
