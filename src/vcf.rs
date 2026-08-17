//! VCF file parser for SNP extraction.
//!
//! Parses VCF files to extract single-nucleotide polymorphisms (SNPs),
//! skipping indels and structural variants. Supports allele frequency
//! extraction from INFO/AF or calculation from GT fields.

use anyhow::Context;
use log::warn;
use rayon::prelude::*;
use std::io::BufRead;
use std::path::Path;

mod carriers;
mod filter;
mod merge;
mod parse;
#[cfg(test)]
mod tests;

pub use carriers::CarrierSet;
pub use filter::filter_snps;
pub use merge::merge_vcfs;
pub use parse::{parse_vcf, sample_count, sample_names};

/// Exact derived-allele counts read from the per-sample GT columns, when present.
/// Used by the diversity statistics, which need the true per-site allele counts
/// rather than a count reconstructed from a (possibly imprecise) INFO/AF frequency.
#[derive(Debug, Clone)]
pub struct GtCounts {
    /// Derived-allele count per *valid* ALT (same order as `VcfSnp::alt_alleles`).
    pub alt: Vec<usize>,
    /// Total called alleles at the site (the haploid sample size actually observed);
    /// the diversity path uses it as the per-site denominator so no-calls do not
    /// dilute the frequency, matching the AF path.
    pub called: usize,
}

/// A single SNP record parsed from a VCF file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VcfSnp {
    /// Chromosome/contig name
    pub chrom: String,
    /// 1-based position on the chromosome
    pub pos: usize,
    /// Reference allele (single base for SNPs)
    pub ref_allele: u8,
    /// Alternative alleles (single base each for SNPs)
    pub alt_alleles: Vec<u8>,
    /// Allele frequencies for each ALT allele (if available)
    pub alt_freqs: Vec<f64>,
    /// Exact GT-derived allele counts, `Some` only when the record carries genotype
    /// columns. The diversity path prefers these over `round(AF * n)`.
    pub gt_counts: Option<GtCounts>,
    /// Which samples carry each *valid* ALT (same order as `alt_alleles`), `Some`
    /// whenever per-sample identity survived parsing: genotype columns in a single
    /// multi-sample VCF, or the file index in a `--vcf-list` merge (one haploid sample
    /// per file). Kept alongside `gt_counts`, never instead of it, so every existing
    /// count stays byte-identical.
    ///
    /// This is what turns the same-codon check from an allele-frequency bound into an
    /// observation: two alleles listed for the same sample are on the same molecule,
    /// because the genotypes are haploid.
    pub carriers: Option<Vec<CarrierSet>>,
    /// FILTER field value
    pub filter: String,
    /// Read depth from INFO/DP (if available)
    pub depth: Option<u32>,
}

