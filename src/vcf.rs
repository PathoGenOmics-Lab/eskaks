//! VCF file parser for SNP extraction.
//!
//! Parses VCF files to extract single-nucleotide polymorphisms (SNPs),
//! skipping indels and structural variants. Supports allele frequency
//! extraction from INFO/AF or calculation from GT fields.

use anyhow::Context;
use log::warn;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single SNP record parsed from a VCF file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VcfSnp {
    /// Chromosome/contig name
    pub chrom: String,
    /// 1-based position on the chromosome
    pub pos: usize,
    /// Reference allele (single base for SNPs)
    pub ref_allele: u8,
    /// Alternative alleles (single base each for SNPs)
    pub alt_alleles: Vec<u8>,
    /// Allele frequencies for each ALT allele (if available)
    pub alt_freqs: Vec<f64>,
    /// FILTER field value
    pub filter: String,
    /// Read depth from INFO/DP (if available)
    pub depth: Option<u32>,
}

/// Parse a VCF file and return all SNP records.
///
/// Skips header lines, indels, and multi-base variants.
/// Extracts allele frequencies from INFO/AF or calculates from GT fields.
pub fn parse_vcf(path: &Path) -> anyhow::Result<Vec<VcfSnp>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open VCF file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut snps = Vec::new();
    let mut sample_count = 0usize;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {} of VCF", line_no + 1))?;
        let line = line.trim_end();

        // Skip meta-information lines
        if line.starts_with("##") {
            continue;
        }

        // Parse header line to determine sample columns
        if line.starts_with('#') {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 9 {
                sample_count = fields.len() - 9;
            }
            continue;
        }

        // Parse data line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            warn!("VCF line {} has fewer than 8 fields, skipping", line_no + 1);
            continue;
        }

        let chrom = fields[0].to_string();
        let pos: usize = fields[1]
            .parse()
            .with_context(|| format!("Invalid POS at VCF line {}", line_no + 1))?;
        let ref_allele = fields[3];
        let alt_field = fields[4];
        let filter = fields[6].to_string();
        let info = fields[7];

        // Skip non-SNPs: ref must be single base
        if ref_allele.len() != 1 {
            continue;
        }
        let ref_base = ref_allele.as_bytes()[0].to_ascii_uppercase();
        if !matches!(ref_base, b'A' | b'C' | b'G' | b'T') {
            continue;
        }

        // Parse ALT alleles, keeping only single-base SNPs
        let alt_alleles_raw: Vec<&str> = alt_field.split(',').collect();
        let mut alt_alleles = Vec::new();
        let mut valid_alt_indices = Vec::new();

        for (i, alt) in alt_alleles_raw.iter().enumerate() {
            if alt.len() == 1 {
                let alt_base = alt.as_bytes()[0].to_ascii_uppercase();
                if matches!(alt_base, b'A' | b'C' | b'G' | b'T') && alt_base != ref_base {
                    alt_alleles.push(alt_base);
                    valid_alt_indices.push(i);
                }
            }
        }

        if alt_alleles.is_empty() {
            continue;
        }

        // Parse depth from INFO/DP
        let depth = parse_info_field(info, "DP")
            .and_then(|v| v.parse::<u32>().ok());

        // Parse allele frequencies
        let alt_freqs = if let Some(af_str) = parse_info_field(info, "AF") {
            // AF from INFO field (Number=A: one per ALT). Keep positional alignment —
            // an unparseable token (e.g. ".") becomes None at ITS position so the
            // remaining ALTs keep the right frequency, instead of shifting left.
            let all_afs: Vec<Option<f64>> = af_str
                .split(',')
                .map(|v| v.parse::<f64>().ok())
                .collect();
            valid_alt_indices
                .iter()
                .map(|&i| all_afs.get(i).copied().flatten().unwrap_or(0.0))
                .collect()
        } else if sample_count > 0 && fields.len() > 9 {
            // The GT position is read from THIS record's FORMAT column — VCF allows the
            // FORMAT order to differ between records, so it must not be cached file-wide.
            let gt_idx = fields.get(8).and_then(|f| f.split(':').position(|k| k == "GT"));
            if let Some(gt_idx) = gt_idx {
                calculate_af_from_gt(&fields[9..], gt_idx, alt_alleles_raw.len(), &valid_alt_indices)
            } else {
                vec![1.0; alt_alleles.len()]
            }
        } else {
            // No frequency info available, assume fixed
            vec![1.0; alt_alleles.len()]
        };

        snps.push(VcfSnp {
            chrom,
            pos,
            ref_allele: ref_base,
            alt_alleles,
            alt_freqs,
            filter,
            depth,
        });
    }

    // An empty VCF is not fatal on its own: in a multi-sample merge one
    // over-filtered or variant-free isolate must not abort the whole run.
    // Callers decide whether zero total SNPs is an error.
    if snps.is_empty() {
        warn!("No SNPs found in VCF file: {}", path.display());
    }

    Ok(snps)
}

/// Extract a value from the INFO field by key.
fn parse_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    for entry in info.split(';') {
        if let Some(rest) = entry.strip_prefix(key) {
            if let Some(val) = rest.strip_prefix('=') {
                return Some(val);
            }
        }
    }
    None
}

/// Calculate allele frequencies from GT (genotype) fields.
fn calculate_af_from_gt(
    samples: &[&str],
    gt_idx: usize,
    total_alt_count: usize,
    valid_indices: &[usize],
) -> Vec<f64> {
    let mut allele_counts = vec![0u32; total_alt_count + 1]; // index 0 = REF
    let mut total_alleles = 0u32;

    for sample in samples {
        let fields: Vec<&str> = sample.split(':').collect();
        if let Some(gt) = fields.get(gt_idx) {
            // Handle both phased (|) and unphased (/) genotypes
            for allele_str in gt.split(['/', '|']) {
                if allele_str == "." {
                    continue;
                }
                if let Ok(allele_idx) = allele_str.parse::<usize>() {
                    if allele_idx < allele_counts.len() {
                        allele_counts[allele_idx] += 1;
                    }
                    total_alleles += 1;
                }
            }
        }
    }

    if total_alleles == 0 {
        return vec![0.0; valid_indices.len()];
    }

    valid_indices
        .iter()
        .map(|&i| allele_counts.get(i + 1).copied().unwrap_or(0) as f64 / total_alleles as f64)
        .collect()
}

/// Apply filters to SNPs based on user criteria.
pub fn filter_snps(
    snps: Vec<VcfSnp>,
    pass_only: bool,
    min_af: Option<f64>,
    max_af: Option<f64>,
    min_depth: Option<u32>,
) -> Vec<VcfSnp> {
    snps.into_iter()
        .filter(|snp| {
            if pass_only && snp.filter != "PASS" && snp.filter != "." {
                return false;
            }
            if let Some(min_dp) = min_depth {
                match snp.depth {
                    Some(dp) if dp < min_dp => return false,
                    None => return false, // No DP info → filter out when min_depth requested
                    _ => {}
                }
            }
            if let Some(min) = min_af {
                if snp.alt_freqs.iter().all(|&af| af < min) {
                    return false;
                }
            }
            if let Some(max) = max_af {
                if snp.alt_freqs.iter().any(|&af| af > max) {
                    return false;
                }
            }
            true
        })
        .collect()
}

/// Merge SNPs from multiple single-sample VCF files.
///
/// For each unique (CHROM, POS, ALT) combination, the allele frequency is
/// computed as the fraction of samples carrying that variant.
/// Pre-filters each VCF by PASS and min_depth before merging.
pub fn merge_vcfs(
    vcf_paths: &[String],
    pass_only: bool,
    min_depth: Option<u32>,
) -> anyhow::Result<Vec<VcfSnp>> {
    use std::collections::HashMap;

    let n_samples = vcf_paths.len() as f64;

    // Key: (chrom, pos, alt_base) → count of samples with this variant
    let mut variant_counts: HashMap<(String, usize, u8), (u8, f64)> = HashMap::new();
    // Also store ref_allele per (chrom, pos)
    let mut ref_alleles: HashMap<(String, usize), u8> = HashMap::new();
    // Track depth sums for averaging
    let mut depth_sums: HashMap<(String, usize), (u32, u32)> = HashMap::new(); // (sum, count)

    // Parse + filter every VCF in parallel (the expensive, I/O-bound step). The
    // map reduction below stays serial and is order-independent (AF = count /
    // n_samples), so the merged result is deterministic regardless of threads.
    let per_file: Vec<Vec<VcfSnp>> = vcf_paths
        .par_iter()
        .map(|path| {
            let snps = parse_vcf(std::path::Path::new(path))?;
            Ok(filter_snps(snps, pass_only, None, None, min_depth))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    for (i, snps) in per_file.into_iter().enumerate() {
        log::info!("  Sample {}/{}: {} — {} SNPs", i + 1, vcf_paths.len(), vcf_paths[i], snps.len());
        if snps.is_empty() {
            warn!("Sample {} contributed 0 SNPs (empty or fully filtered); skipping", vcf_paths[i]);
        }

        // AF = fraction of SAMPLES carrying the variant, so each sample contributes at
        // most 1 per (chrom, pos, alt) — dedupe within the file to guard against a
        // caller emitting the same allele on more than one record.
        let mut seen_this_sample: std::collections::HashSet<(String, usize, u8)> =
            std::collections::HashSet::new();
        for snp in snps {
            ref_alleles
                .entry((snp.chrom.clone(), snp.pos))
                .or_insert(snp.ref_allele);

            if let Some(dp) = snp.depth {
                let entry = depth_sums
                    .entry((snp.chrom.clone(), snp.pos))
                    .or_insert((0, 0));
                entry.0 += dp;
                entry.1 += 1;
            }

            for alt_base in &snp.alt_alleles {
                let key = (snp.chrom.clone(), snp.pos, *alt_base);
                if seen_this_sample.insert(key.clone()) {
                    variant_counts.entry(key).or_insert((snp.ref_allele, 0.0)).1 += 1.0;
                }
            }
        }
    }

    // Convert to VcfSnp records with AF = count / n_samples
    // Group by (chrom, pos) to handle multi-allelic
    let mut pos_map: HashMap<(String, usize), Vec<(u8, f64)>> = HashMap::new();
    for ((chrom, pos, alt), (_ref_base, count)) in &variant_counts {
        pos_map
            .entry((chrom.clone(), *pos))
            .or_default()
            .push((*alt, *count / n_samples));
    }

    let mut merged: Vec<VcfSnp> = pos_map
        .into_iter()
        .filter_map(|((chrom, pos), mut alts)| {
            let ref_base = ref_alleles
                .get(&(chrom.clone(), pos))
                .copied()
                .unwrap_or(b'N');
            let avg_depth = depth_sums.get(&(chrom.clone(), pos)).map(|(sum, cnt)| {
                if *cnt > 0 { sum / cnt } else { 0 }
            });
            // `alts` came out of a HashMap, whose iteration order is randomized per
            // run. Sort by ALT base so a multi-allelic position emits its ALTs (and
            // their aligned frequencies) in a stable order — otherwise the merged
            // output is not reproducible across runs.
            alts.sort_by_key(|(alt, _)| *alt);
            // Samples can disagree on REF; drop any ALT that equals the merged REF (it
            // is not a variant against it) and skip a position left with no real ALT.
            let (alt_alleles, alt_freqs): (Vec<u8>, Vec<f64>) =
                alts.into_iter().filter(|(alt, _)| *alt != ref_base).unzip();
            if alt_alleles.is_empty() {
                return None;
            }
            Some(VcfSnp {
                chrom,
                pos,
                ref_allele: ref_base,
                alt_alleles,
                alt_freqs,
                filter: "PASS".to_string(),
                depth: avg_depth,
            })
        })
        .collect();

    merged.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)));

    if merged.is_empty() {
        anyhow::bail!("No SNPs found after merging {} VCF files", vcf_paths.len());
    }

    Ok(merged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_vcf(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
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
}
