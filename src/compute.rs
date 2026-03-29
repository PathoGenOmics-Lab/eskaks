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
}
