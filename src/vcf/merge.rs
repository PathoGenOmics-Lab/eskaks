//! Merging per-sample VCFs into one variant set.

use super::*;

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
    let mut depth_sums: HashMap<(String, usize), (u64, u32)> = HashMap::new(); // (sum, count)

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
                entry.0 += dp as u64; // u64 sum: many samples × high DP can exceed u32
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
                if *cnt > 0 { (sum / *cnt as u64) as u32 } else { 0 }
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

