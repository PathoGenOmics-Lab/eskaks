//! VCF parsing: SNP extraction and allele-frequency computation.

use super::*;

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
    // Format diagnostics: distinguish "valid VCF, no SNPs" from "this isn't a VCF".
    let mut data_lines = 0usize;
    let mut malformed_lines = 0usize;
    let mut first_bad: Option<String> = None;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {} of VCF", line_no + 1))?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }

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
        data_lines += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            malformed_lines += 1;
            if first_bad.is_none() {
                first_bad = Some(line.chars().take(60).collect());
            }
            // Cap the per-line spam so a wrong-format file doesn't scroll the terminal.
            if malformed_lines <= 3 {
                warn!("VCF line {} has fewer than 8 tab-separated fields, skipping", line_no + 1);
            } else if malformed_lines == 4 {
                warn!("… further malformed-line warnings suppressed (see the format check below).");
            }
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
            // Reject non-finite or out-of-[0,1] AF tokens (e.g. "nan", "inf", "2.0"):
            // one malformed token must not poison the genome-wide πN/πS nor slip past
            // the --min-af/--max-af filters (every NaN comparison is false).
            let all_afs: Vec<Option<f64>> = af_str
                .split(',')
                .map(|v| v.parse::<f64>().ok().filter(|x| x.is_finite() && (0.0..=1.0).contains(x)))
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

    // If EVERY data line was malformed, this is almost certainly not a VCF — fail
    // loudly with a format hint rather than returning an empty, "successful" result.
    if snps.is_empty() && data_lines > 0 && malformed_lines == data_lines {
        let hint = match first_bad.as_deref() {
            Some(l) if l.starts_with('>') => " (the first data line starts with '>': is this a FASTA?)",
            _ => "",
        };
        anyhow::bail!(
            "{}: {} data line(s) but 0 valid VCF records — none had 8+ tab-separated columns. \
             Is this really a tab-delimited VCF?{}",
            path.display(), data_lines, hint
        );
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
                    // An out-of-range allele index (a malformed GT referencing an
                    // undeclared ALT) is treated as missing: it must not inflate
                    // total_alleles and deflate every real frequency at the site.
                    if allele_idx < allele_counts.len() {
                        allele_counts[allele_idx] += 1;
                        total_alleles += 1;
                    }
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

