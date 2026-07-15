//! VCF file parser for SNP extraction.
//!
//! Parses VCF files to extract single-nucleotide polymorphisms (SNPs),
//! skipping indels and structural variants. Supports allele frequency
//! extraction from INFO/AF or calculation from GT fields.

use anyhow::Context;
use log::warn;
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
    let mut gt_field_idx: Option<usize> = None;

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
            // AF from INFO field
            let all_afs: Vec<f64> = af_str
                .split(',')
                .filter_map(|v| v.parse::<f64>().ok())
                .collect();
            valid_alt_indices
                .iter()
                .map(|&i| all_afs.get(i).copied().unwrap_or(0.0))
                .collect()
        } else if sample_count > 0 && fields.len() > 9 {
            // Calculate from GT fields
            if gt_field_idx.is_none() {
                if let Some(format_field) = fields.get(8) {
                    gt_field_idx = format_field
                        .split(':')
                        .position(|f| f == "GT");
                }
            }
            if let Some(gt_idx) = gt_field_idx {
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

    for (i, path) in vcf_paths.iter().enumerate() {
        let snps = parse_vcf(std::path::Path::new(path))?;
        let snps = filter_snps(snps, pass_only, None, None, min_depth);
        log::info!("  Sample {}/{}: {} — {} SNPs", i + 1, vcf_paths.len(), path, snps.len());
        if snps.is_empty() {
            warn!("Sample {} contributed 0 SNPs (empty or fully filtered); skipping", path);
        }

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
                let entry = variant_counts.entry(key).or_insert((snp.ref_allele, 0.0));
                entry.1 += 1.0;
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
        .map(|((chrom, pos), alts)| {
            let ref_base = ref_alleles
                .get(&(chrom.clone(), pos))
                .copied()
                .unwrap_or(b'N');
            let avg_depth = depth_sums.get(&(chrom.clone(), pos)).map(|(sum, cnt)| {
                if *cnt > 0 { sum / cnt } else { 0 }
            });
            let (alt_alleles, alt_freqs): (Vec<u8>, Vec<f64>) = alts.into_iter().unzip();
            VcfSnp {
                chrom,
                pos,
                ref_allele: ref_base,
                alt_alleles,
                alt_freqs,
                filter: "PASS".to_string(),
                depth: avg_depth,
            }
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
}
