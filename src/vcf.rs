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

mod filter;
mod merge;
mod parse;
#[cfg(test)]
mod tests;

pub use filter::filter_snps;
pub use merge::merge_vcfs;
pub use parse::{parse_vcf, sample_count};

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
    /// FILTER field value
    pub filter: String,
    /// Read depth from INFO/DP (if available)
    pub depth: Option<u32>,
}

