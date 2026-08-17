//! SNP filtering by FILTER/PASS, allele frequency, and depth.

use super::*;

/// Apply filters to SNPs based on user criteria.
pub fn filter_snps(
    snps: Vec<VcfSnp>,
    pass_only: bool,
    min_af: Option<f64>,
    max_af: Option<f64>,
    min_depth: Option<u32>,
) -> Vec<VcfSnp> {
    snps.into_iter()
        .filter_map(|mut snp| {
            if pass_only && snp.filter != "PASS" && snp.filter != "." {
                return None;
            }
            if let Some(min_dp) = min_depth {
                match snp.depth {
                    Some(dp) if dp < min_dp => return None,
                    None => return None, // No DP info → filter out when min_depth requested
                    _ => {}
                }
            }
            // Allele-frequency filtering is PER-ALLELE: at a multi-allelic site, prune
            // only the ALTs whose frequency is outside [min_af, max_af] and keep the
            // rest. This stops a sub-threshold allele from leaking through and a
            // co-located in-range allele from being dropped with it. The record is
            // dropped only when no ALT survives.
            if min_af.is_some() || max_af.is_some() {
                let lo = min_af.unwrap_or(f64::NEG_INFINITY);
                let hi = max_af.unwrap_or(f64::INFINITY);
                // Which ALT positions survive the frequency window.
                let keep: Vec<bool> =
                    snp.alt_freqs.iter().map(|&af| af >= lo && af <= hi).collect();
                if !keep.iter().any(|&k| k) {
                    return None;
                }
                // Prune EVERY per-ALT vector in lockstep. gt_counts.alt is documented to
                // be parallel to alt_alleles; leaving it stale shifts the surviving ALTs
                // out of sync with their genotype counts and silently corrupts the
                // diversity statistics (piN/piS, theta_W, Tajima's D) at multi-allelic
                // sites, since the diversity path reads gt_counts by ALT index.
                let mut it = keep.iter();
                snp.alt_alleles.retain(|_| *it.next().unwrap());
                let mut it = keep.iter();
                snp.alt_freqs.retain(|_| *it.next().unwrap());
                if let Some(gc) = snp.gt_counts.as_mut() {
                    if gc.alt.len() == keep.len() {
                        let mut it = keep.iter();
                        gc.alt.retain(|_| *it.next().unwrap());
                    }
                }
                // Carriers are per-ALT too, and the same-codon check reads them by ALT
                // index: a stale vector here would intersect the wrong allele's samples
                // and report co-occurrence for a variant that was filtered out.
                if let Some(cs) = snp.carriers.as_mut() {
                    if cs.len() == keep.len() {
                        let mut it = keep.iter();
                        cs.retain(|_| *it.next().unwrap());
                    }
                }
            }
            Some(snp)
        })
        .collect()
}

