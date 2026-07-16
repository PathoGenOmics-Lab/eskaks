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

