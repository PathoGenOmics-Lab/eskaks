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

    /// One merged ALT at a position: the base, its frequency, and the samples carrying it.
    type MergedAlt = (u8, f64, CarrierSet);

    let n_samples = vcf_paths.len() as f64;

    // Key: (chrom, pos, alt_base) → (ref base, count of samples, which samples).
    // The merge model is one sample per file, so the file index IS the sample index:
    // co-occurrence between two alleles is directly observable here, with no phasing.
    let mut variant_counts: HashMap<(String, usize, u8), (u8, f64, CarrierSet)> = HashMap::new();
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
            let p = std::path::Path::new(path);
            // The merge model is one sample per file (AF = carriers / n_files). A file
            // with multiple sample columns would be miscounted as a single presence
            // unit, so warn rather than silently produce a wrong frequency.
            if sample_count(p).unwrap_or(0) > 1 {
                warn!(
                    "{}: multiple sample columns; the merge treats each --vcf file as ONE \
                     sample, so per-sample genotypes are not used. Split to one sample per file.",
                    path
                );
            }
            let snps = parse_vcf(p)?;
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

            // `alt_idx`, not `i`: `i` is the FILE index above, which is also this
            // record's sample index, and shadowing it here would record every carrier
            // as the ALT's position in the record instead of the sample it came from.
            for (alt_idx, alt_base) in snp.alt_alleles.iter().enumerate() {
                // Count a sample as a carrier only when it actually carries the ALT. A
                // per-sample VCF that records non-carrier sites (GT 0/0, ./., or INFO
                // AF=0.0, as gVCF / all-sites callers do) has alt_freqs[i] == 0.0 there;
                // counting those as carriers would inflate the merged allele frequency
                // and, downstream, drop a real polymorphism as "fixed".
                if snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0) <= 0.0 {
                    continue;
                }
                let key = (snp.chrom.clone(), snp.pos, *alt_base);
                if seen_this_sample.insert(key.clone()) {
                    let entry = variant_counts
                        .entry(key)
                        .or_insert_with(|| (snp.ref_allele, 0.0, CarrierSet::new(vcf_paths.len())));
                    entry.1 += 1.0;
                    entry.2.insert(i);
                }
            }
        }
    }

    // Convert to VcfSnp records with AF = count / n_samples
    // Group by (chrom, pos) to handle multi-allelic
    let mut pos_map: HashMap<(String, usize), Vec<MergedAlt>> = HashMap::new();
    for ((chrom, pos, alt), (_ref_base, count, carriers)) in &variant_counts {
        pos_map
            .entry((chrom.clone(), *pos))
            .or_default()
            .push((*alt, *count / n_samples, carriers.clone()));
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
            alts.sort_by_key(|(alt, _, _)| *alt);
            // Samples can disagree on REF; drop any ALT that equals the merged REF (it
            // is not a variant against it) and skip a position left with no real ALT.
            // The three per-ALT vectors are built in one pass so they cannot drift apart.
            let (mut alt_alleles, mut alt_freqs, mut alt_carriers) =
                (Vec::new(), Vec::new(), Vec::new());
            for (alt, af, carriers) in alts {
                if alt != ref_base {
                    alt_alleles.push(alt);
                    alt_freqs.push(af);
                    alt_carriers.push(carriers);
                }
            }
            if alt_alleles.is_empty() {
                return None;
            }
            Some(VcfSnp {
                chrom,
                pos,
                ref_allele: ref_base,
                alt_alleles,
                alt_freqs,
                // A merged record has no per-sample GENOTYPE COLUMNS; its af is already
                // the exact carriers / n_samples, so round(af * n) recovers the count
                // and the diversity path needs no GtCounts here.
                gt_counts: None,
                // Per-sample IDENTITY, however, survives the merge: each input file is
                // one haploid sample, so the file index is the sample index and two
                // alleles sharing one index are two alleles in one genome.
                carriers: Some(alt_carriers),
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

