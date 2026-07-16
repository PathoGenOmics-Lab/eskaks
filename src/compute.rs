//! Unified compute engine for dN/dS calculation across models.
//!
//! Encapsulates the precomputed lookup tables and provides a single
//! interface for pairwise computation, eliminating duplicated closures.

use log::info;

use crate::genetic_code::GeneticCode;
use crate::input::SequenceData;
use crate::models::{DsDn, Model};

/// Unified computation engine that holds precomputed tables for the selected model.
pub enum ComputeEngine {
    Nei(Box<crate::models::nei::NeiTables>),
    Li(Box<crate::models::li::LiTables>),
}

impl ComputeEngine {
    /// Build the compute engine for the given model and genetic code.
    pub fn new(model: Model, gc: &GeneticCode) -> Self {
        match model {
            Model::Nei => {
                info!("Precomputing lookup tables for Nei-Gojobori (1986) model...");
                let tables = crate::models::nei::NeiTables::with_genetic_code(gc);
                info!("Nei precomputation finished.");
                ComputeEngine::Nei(tables)
            }
            Model::Li => {
                info!("Precomputing lookup tables for Li (1993) model...");
                let tables = crate::models::li::LiTables::with_genetic_code(&gc.aa_table);
                info!("Li precomputation finished.");
                ComputeEngine::Li(tables)
            }
        }
    }

    /// Compute dN/dS for a pair of unique sequence indices.
    #[inline]
    pub fn compute_pair(&self, data: &SequenceData, u_i: usize, u_j: usize) -> DsDn {
        if u_i == u_j {
            // A sequence versus itself has 0 divergence — but only if it has at least
            // one valid codon. An all-N / all-gap sequence has no comparable codons, so
            // its dN/dS is undefined (NaN), not a spurious 0.0 that reads as strong
            // purifying selection.
            return if data.unique_codon_indices[u_i]
                .iter()
                .all(|&c| c == crate::codon::INVALID_CODON)
            {
                DsDn { dn: f64::NAN, ds: f64::NAN }
            } else {
                DsDn { dn: 0.0, ds: 0.0 }
            };
        }
        let (dn, ds) = self.compute_slices(
            &data.unique_codon_indices[u_i],
            &data.unique_codon_indices[u_j],
        );
        DsDn { dn, ds }
    }

    /// Compute dN/dS for raw codon slices (used by sliding window mode).
    #[inline]
    pub fn compute_slices(&self, s1: &[u8], s2: &[u8]) -> (f64, f64) {
        match self {
            ComputeEngine::Nei(tables) => tables.compute_pair(s1, s2),
            ComputeEngine::Li(tables) => tables.compute_pair(s1, s2),
        }
    }

    /// Percentile bootstrap 95% CIs for a pair's dN, dS, and dN/dS by resampling
    /// codon columns with replacement `n_boot` times (seeded). Model-agnostic,
    /// so it works for both Nei and Li. Returns
    /// `(dn_lo, dn_hi, ds_lo, ds_hi, ratio_lo, ratio_hi)`.
    pub fn bootstrap_ci(
        &self,
        s1: &[u8],
        s2: &[u8],
        n_boot: usize,
        seed: u64,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let l = s1.len().min(s2.len());
        let nan6 = (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN);
        if l == 0 || n_boot == 0 {
            return nan6;
        }
        let mut rng = crate::stats::SplitMix64::new(seed);
        let (mut dns, mut dss, mut ratios) = (Vec::new(), Vec::new(), Vec::new());
        let mut b1 = vec![0u8; l];
        let mut b2 = vec![0u8; l];
        for _ in 0..n_boot {
            for k in 0..l {
                let idx = rng.below(l);
                b1[k] = s1[idx];
                b2[k] = s2[idx];
            }
            let (dn, ds) = self.compute_slices(&b1, &b2);
            if dn.is_finite() {
                dns.push(dn);
            }
            if ds.is_finite() {
                dss.push(ds);
            }
            if dn.is_finite() && ds.is_finite() && ds > 0.0 {
                ratios.push(dn / ds);
            }
        }
        let ci = |mut v: Vec<f64>| -> (f64, f64) {
            if v.is_empty() {
                return (f64::NAN, f64::NAN);
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            (
                crate::stats::percentile_sorted(&v, 2.5),
                crate::stats::percentile_sorted(&v, 97.5),
            )
        };
        let (dn_lo, dn_hi) = ci(dns);
        let (ds_lo, ds_hi) = ci(dss);
        let (r_lo, r_hi) = ci(ratios);
        (dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi)
    }

    /// Compute `(dN, dS, var_dN, var_dS)` for a pair of unique indices. The
    /// Nei-Gojobori analytic variances are only defined for the Nei model; the
    /// Li model returns NaN variances (use bootstrap for Li instead).
    #[inline]
    pub fn compute_pair_stats(&self, data: &SequenceData, u_i: usize, u_j: usize) -> (f64, f64, f64, f64) {
        if u_i == u_j {
            // See compute_pair: an all-invalid self-comparison is undefined, not 0.
            return if data.unique_codon_indices[u_i]
                .iter()
                .all(|&c| c == crate::codon::INVALID_CODON)
            {
                (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
        }
        let s1 = &data.unique_codon_indices[u_i];
        let s2 = &data.unique_codon_indices[u_j];
        match self {
            ComputeEngine::Nei(tables) => tables.compute_pair_stats(s1, s2),
            ComputeEngine::Li(tables) => {
                let (dn, ds) = tables.compute_pair(s1, s2);
                (dn, ds, f64::NAN, f64::NAN)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::fasta_to_codon_indices;

    fn make_data(seqs: &[&[u8]], model: Model) -> SequenceData {
        let indices: Vec<Vec<u8>> = seqs.iter().map(|s| fasta_to_codon_indices(s, model)).collect();
        let ids: Vec<String> = (0..seqs.len()).map(|i| format!("s{}", i)).collect();
        let uidx_by_id: Vec<usize> = (0..seqs.len()).collect();
        let n_unique = seqs.len();
        SequenceData {
            ids,
            uidx_by_id,
            unique_codon_indices: indices,
            n_unique,
        }
    }

    #[test]
    fn nei_engine_identical_gives_zero() {
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let data = make_data(&[b"ATGGCTGCT", b"ATGGCTGCT"], Model::Nei);
        let result = engine.compute_pair(&data, 0, 0);
        assert_eq!(result.dn, 0.0);
        assert_eq!(result.ds, 0.0);
    }

    #[test]
    fn li_engine_identical_gives_zero() {
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Li, gc);
        let data = make_data(&[b"ATGGCTGCTGCT", b"ATGGCTGCTGCT"], Model::Li);
        let result = engine.compute_pair(&data, 0, 0);
        assert_eq!(result.dn, 0.0);
        assert_eq!(result.ds, 0.0);
    }

    #[test]
    fn all_invalid_self_pair_is_nan_not_zero() {
        // A sequence with no comparable codons (all-N / all-gap) compared with
        // itself carries no information: dN/dS must be NaN, not a spurious 0.0
        // that would read as perfect purifying selection. Regression for all-N
        // pairs rendering as 0.0.
        let gc = crate::genetic_code::get_table(1).unwrap();
        for model in [Model::Nei, Model::Li] {
            let engine = ComputeEngine::new(model, gc);
            let data = make_data(&[b"NNNNNNNNN", b"---------"], model);
            for u in [0usize, 1usize] {
                let d = engine.compute_pair(&data, u, u);
                assert!(d.dn.is_nan(), "{:?} all-invalid self dN should be NaN", model);
                assert!(d.ds.is_nan(), "{:?} all-invalid self dS should be NaN", model);
                let (dn, ds, vn, vs) = engine.compute_pair_stats(&data, u, u);
                assert!(dn.is_nan() && ds.is_nan() && vn.is_nan() && vs.is_nan());
            }
        }
    }

    #[test]
    fn identical_with_one_valid_codon_is_zero_not_nan() {
        // A single comparable codon is enough information to report zero
        // divergence; only truly all-invalid sequences collapse to NaN.
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let data = make_data(&[b"ATGNNNNNN", b"ATGNNNNNN"], Model::Nei);
        let d = engine.compute_pair(&data, 0, 0);
        assert_eq!(d.dn, 0.0);
        assert_eq!(d.ds, 0.0);
        assert_eq!(engine.compute_pair_stats(&data, 0, 0), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn nei_engine_synonymous_change() {
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let data = make_data(&[b"ATGGCTGCT", b"ATGGCCGCT"], Model::Nei);
        let result = engine.compute_pair(&data, 0, 1);
        assert_eq!(result.dn, 0.0);
        assert!(result.ds > 0.0, "dS should be > 0 for synonymous change");
    }

    #[test]
    fn compute_slices_matches_compute_pair() {
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let s1 = fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        let s2 = fasta_to_codon_indices(b"ATGATTGCT", Model::Nei);
        let data = make_data(&[b"ATGGCTGCT", b"ATGATTGCT"], Model::Nei);
        let pair = engine.compute_pair(&data, 0, 1);
        let (dn, ds) = engine.compute_slices(&s1, &s2);
        assert!((pair.dn - dn).abs() < 1e-10);
        assert!((pair.ds - ds).abs() < 1e-10);
    }

    #[test]
    fn mito_genetic_code_changes_results() {
        let std_gc = crate::genetic_code::get_table(1).unwrap();
        let mito_gc = crate::genetic_code::get_table(2).unwrap();
        let engine_std = ComputeEngine::new(Model::Nei, std_gc);
        let engine_mito = ComputeEngine::new(Model::Nei, mito_gc);
        // AGA codes Arg in standard, stop in vert mito
        let data_std = make_data(&[b"ATGGCTAGA", b"ATGATTAGA"], Model::Nei);
        let data_mito = make_data(&[b"ATGGCTAGA", b"ATGATTAGA"], Model::Nei);
        let r_std = engine_std.compute_pair(&data_std, 0, 1);
        let r_mito = engine_mito.compute_pair(&data_mito, 0, 1);
        // Results should differ because AGA changes meaning
        assert!(
            (r_std.dn - r_mito.dn).abs() > 1e-6 || (r_std.ds - r_mito.ds).abs() > 1e-6,
            "Standard and mito should give different results for AGA-containing sequences"
        );
    }

    // cov_: NaN-tolerant float comparison (bootstrap CIs can be NaN by design).
    fn cov_eq_or_nan(a: f64, b: f64) -> bool {
        (a.is_nan() && b.is_nan()) || a == b
    }

    #[test]
    fn cov_bootstrap_ci_empty_or_zero_boot_is_all_nan() {
        // bootstrap_ci returns six NaNs when there is nothing to resample:
        // l == 0 (empty overlap) or n_boot == 0 (early `return nan6`).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let s = fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        // First slice empty => l = 0.min(3) = 0 => early return.
        let re = engine.bootstrap_ci(&s[..0], &s, 100, 1);
        assert!(re.0.is_nan() && re.1.is_nan() && re.2.is_nan()
            && re.3.is_nan() && re.4.is_nan() && re.5.is_nan());
        // n_boot == 0 => same early return.
        let rz = engine.bootstrap_ci(&s, &s, 0, 1);
        assert!(rz.0.is_nan() && rz.1.is_nan() && rz.2.is_nan()
            && rz.3.is_nan() && rz.4.is_nan() && rz.5.is_nan());
    }

    #[test]
    fn cov_bootstrap_ci_identical_zero_and_undefined_ratio() {
        // Identical sequences => every resample yields (dN, dS) = (0, 0):
        //  * dN and dS are finite each iteration => their 2.5/97.5 percentiles are 0.
        //  * dS is never > 0 => the ratio vector stays empty => its CI is (NaN, NaN),
        //    exercising the empty-vector branch of the inner `ci` closure.
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let s = fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        let (dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi) = engine.bootstrap_ci(&s, &s, 64, 12345);
        assert_eq!(dn_lo, 0.0);
        assert_eq!(dn_hi, 0.0);
        assert_eq!(ds_lo, 0.0);
        assert_eq!(ds_hi, 0.0);
        assert!(r_lo.is_nan());
        assert!(r_hi.is_nan());
    }

    #[test]
    fn cov_bootstrap_ci_deterministic_and_ordered() {
        // Same seed => identical output (SplitMix64 is deterministic); and each
        // finite interval is ordered lo <= hi (2.5th percentile <= 97.5th).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        // ATG (Met, unchanged) . GCT->GCC (Ala, synonymous => dS) . AAA->AAT (Lys->Asn => dN)
        let s1 = fasta_to_codon_indices(b"ATGGCTAAA", Model::Nei);
        let s2 = fasta_to_codon_indices(b"ATGGCCAAT", Model::Nei);
        let a = engine.bootstrap_ci(&s1, &s2, 200, 999);
        let b = engine.bootstrap_ci(&s1, &s2, 200, 999);
        assert!(cov_eq_or_nan(a.0, b.0) && cov_eq_or_nan(a.1, b.1) && cov_eq_or_nan(a.2, b.2)
            && cov_eq_or_nan(a.3, b.3) && cov_eq_or_nan(a.4, b.4) && cov_eq_or_nan(a.5, b.5));
        if a.0.is_finite() && a.1.is_finite() { assert!(a.0 <= a.1); }
        if a.2.is_finite() && a.3.is_finite() { assert!(a.2 <= a.3); }
        if a.4.is_finite() && a.5.is_finite() { assert!(a.4 <= a.5); }
    }

    #[test]
    fn cov_compute_pair_stats_identity_is_zero() {
        // u_i == u_j short-circuits to (dN, dS, var_dN, var_dS) = (0, 0, 0, 0).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let data = make_data(&[b"ATGGCTGCT"], Model::Nei);
        assert_eq!(engine.compute_pair_stats(&data, 0, 0), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn cov_compute_pair_stats_li_variances_are_nan() {
        // Li has no analytic variance => compute_pair_stats returns (dN, dS, NaN, NaN),
        // with dN/dS matching compute_pair (the Li match arm).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Li, gc);
        let data = make_data(&[b"ATGGCTGCTGCT", b"ATGGCCGCTGCT"], Model::Li);
        let (dn, ds, var_dn, var_ds) = engine.compute_pair_stats(&data, 0, 1);
        assert!(var_dn.is_nan(), "Li var_dN must be NaN");
        assert!(var_ds.is_nan(), "Li var_dS must be NaN");
        let pair = engine.compute_pair(&data, 0, 1);
        assert!(cov_eq_or_nan(dn, pair.dn));
        assert!(cov_eq_or_nan(ds, pair.ds));
    }

    #[test]
    fn cov_compute_pair_stats_nei_matches_pair() {
        // Nei arm: dN/dS equal compute_pair; analytic variances, when finite, are
        // non-negative (they are variances).
        let gc = crate::genetic_code::get_table(1).unwrap();
        let engine = ComputeEngine::new(Model::Nei, gc);
        let data = make_data(&[b"ATGGCTAAA", b"ATGGCCAAT"], Model::Nei);
        let (dn, ds, var_dn, var_ds) = engine.compute_pair_stats(&data, 0, 1);
        let pair = engine.compute_pair(&data, 0, 1);
        assert!(cov_eq_or_nan(dn, pair.dn));
        assert!(cov_eq_or_nan(ds, pair.ds));
        if var_dn.is_finite() { assert!(var_dn >= 0.0); }
        if var_ds.is_finite() { assert!(var_ds >= 0.0); }
    }
}
