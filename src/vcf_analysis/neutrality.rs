//! Multiple-testing correction and genomic-control (lambda).

use super::*;

/// Fill in Benjamini-Hochberg FDR q-values and Bonferroni-corrected p-values
/// for the per-gene neutrality test, across every gene with a finite p-value.
/// Genes not tested (no SNPs / --af-weighted) keep NaN.
pub fn apply_multiple_testing(results: &mut [GenePnPs], exclude_repetitive: bool) {
    // Excluding repetitive genes removes them from the test family entirely, so
    // they neither shrink other genes' q-values nor appear as hits.
    let pvals: Vec<f64> = results
        .iter()
        .map(|r| if exclude_repetitive && r.repetitive { f64::NAN } else { r.p_value })
        .collect();
    let qvals = crate::stats::benjamini_hochberg(&pvals);
    let bonf = crate::stats::bonferroni(&pvals);
    for (r, (q, b)) in results.iter_mut().zip(qvals.into_iter().zip(bonf)) {
        r.q_value = q;
        r.p_bonferroni = b;
    }
}

/// Genomic-control inflation factor λ = median(χ²) / median(χ²₁) over the tested
/// genes (finite q-value = in the correction family). Per gene, χ² =
/// (Φ⁻¹(1 − p/2))². λ ≈ 1 means well-calibrated; λ ≫ 1 flags inflation, expected
/// in clonal organisms where genes are not independent. NaN with < 2 tested genes.
pub fn genomic_inflation_lambda(results: &[GenePnPs]) -> f64 {
    let mut chi2: Vec<f64> = results
        .iter()
        .filter(|r| r.q_value.is_finite())
        .map(gene_chi2)
        .filter(|c| c.is_finite())
        .collect();
    if chi2.len() < 2 {
        return f64::NAN;
    }
    chi2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    crate::stats::percentile_sorted(&chi2, 50.0) / 0.454_936_4
}

/// The χ²₁ statistic for a gene's neutrality test, computed from the log-space
/// −log10(p) so it stays finite even where the raw p underflowed (falls back to
/// the raw p-value when `neglog10p` is unavailable, e.g. in synthetic tests).
fn gene_chi2(r: &GenePnPs) -> f64 {
    let nlp = if r.neglog10p.is_finite() {
        r.neglog10p
    } else if r.p_value.is_finite() && r.p_value > 0.0 {
        -r.p_value.log10()
    } else {
        f64::NAN
    };
    crate::stats::chi2_from_two_sided_neglog10p(nlp)
}

/// Apply genomic control: divide every tested gene's χ² by λ (floored at 1, so
/// the correction only ever deflates), then recompute a corrected p-value and a
/// BH q-value into `p_gc` / `q_gc`. Untested genes stay NaN.
pub fn apply_genomic_control(results: &mut [GenePnPs], lambda: f64) {
    let lam = if lambda.is_finite() && lambda > 1.0 { lambda } else { 1.0 };
    let pgc: Vec<f64> = results
        .iter()
        .map(|r| {
            if r.q_value.is_finite() {
                let chi2 = gene_chi2(r);
                if chi2.is_finite() {
                    crate::stats::normal_two_sided_p((chi2 / lam).sqrt())
                } else {
                    f64::NAN
                }
            } else {
                f64::NAN
            }
        })
        .collect();
    let qgc = crate::stats::benjamini_hochberg(&pgc);
    for (r, (p, q)) in results.iter_mut().zip(pgc.into_iter().zip(qgc)) {
        r.p_gc = p;
        r.q_gc = q;
    }
}

