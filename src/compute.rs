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
            return DsDn { dn: 0.0, ds: 0.0 };
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
            return (0.0, 0.0, 0.0, 0.0);
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
}
