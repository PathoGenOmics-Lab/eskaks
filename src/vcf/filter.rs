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
                let (alts, freqs): (Vec<u8>, Vec<f64>) = snp
                    .alt_alleles
                    .iter()
                    .zip(&snp.alt_freqs)
                    .filter(|(_, &af)| af >= lo && af <= hi)
                    .map(|(&a, &f)| (a, f))
                    .unzip();
                if alts.is_empty() {
                    return None;
                }
                snp.alt_alleles = alts;
                snp.alt_freqs = freqs;
            }
            Some(snp)
        })
        .collect()
}

