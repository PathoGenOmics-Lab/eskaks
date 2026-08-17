//! Core pN/pS analysis from VCF data.
//!
//! For each gene, reconstructs reference and alternate codons from SNPs,
//! classifies mutations as synonymous or nonsynonymous using the genetic
//! code table, and computes pN/pS ratios.

use crate::gff::{Gene, Strand};
use crate::genetic_code::GeneticCode;
use crate::vcf::{CarrierSet, VcfSnp};
use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

mod codons;
mod mnv;
mod neutrality;
mod origins;
mod output;
mod plots;
mod pnps;
mod sites;
#[cfg(test)]
mod tests;

// Public API (re-exported so callers keep using `vcf_analysis::…`).
pub use codons::{compute_codon_scan, CodonRow, CodonScan};
pub use origins::{cohort_sample_names, AlleleOrigins};
pub use neutrality::{apply_genomic_control, apply_multiple_testing, genomic_inflation_lambda};
pub use output::{
    genome_wide_diversity, write_codon_scan, write_diversity, write_mk_results, write_results,
    write_variants, GenomeDiversity,
};
pub use plots::{write_pnps_plot, write_pvalue_manhattan};
pub use pnps::{
    bootstrap_genome_wide_ci, compute_pn_ps, genome_wide_core_repetitive, genome_wide_pn_ps,
    ComputeDiagnostics,
};
pub use sites::parse_reference_fasta;

// Internal helpers shared across submodules.
pub(crate) use sites::{
    codon_change_split, codon_space, codon_to_aa, complement, count_sites, count_sites_weighted,
    extract_cds_sequence, genomic_to_cds_offset,
};
// Helpers exercised only by the test suite.
#[cfg(test)]
pub(crate) use sites::{base_to_li, is_transition, reverse_complement};

/// Result of pN/pS analysis for a single gene.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GenePnPs {
    /// Gene name
    pub name: String,
    /// Total CDS length in bp
    pub length_bp: usize,
    /// Number of nonsynonymous sites (fractional)
    pub n_sites: f64,
    /// Number of synonymous sites (fractional)
    pub s_sites: f64,
    /// Proportion of nonsynonymous substitutions (nonsyn_snps / N_sites)
    pub pn: f64,
    /// Proportion of synonymous substitutions (syn_snps / S_sites)
    pub ps: f64,
    /// pN/pS ratio
    pub pn_ps: f64,
    /// Count of nonsynonymous SNPs (AF-weighted if --af-weighted)
    pub nonsyn_snps: f64,
    /// Count of synonymous SNPs (AF-weighted if --af-weighted)
    pub syn_snps: f64,
    /// Total SNPs in this gene (AF-weighted if --af-weighted)
    pub total_snps: f64,
    /// Raw (unweighted) count of classified SNP alleles in this gene. Unlike
    /// `total_snps`, this is never AF-weighted, so `--min-snps` filters on the real
    /// SNP count even under `--af-weighted`.
    pub n_snps: usize,
    /// Genomic start position (1-based, for plotting/output)
    pub genome_start: usize,
    /// Genomic end position (1-based, max exon end)
    pub genome_end: usize,
    /// Strand ('+' or '-')
    pub strand: char,
    /// Chromosome
    pub chrom: String,
    /// Two-sided mid-p binomial p-value for H0: pN/pS = 1 (NaN if untested:
    /// no SNPs, --af-weighted, or a degenerate expected fraction). Mid-p, not the
    /// exact test: see `stats::binomial_two_sided_p` for the calibration reason.
    pub p_value: f64,
    /// Benjamini-Hochberg FDR q-value across all tested genes (NaN if untested)
    pub q_value: f64,
    /// Bonferroni-corrected p-value across all tested genes (NaN if untested)
    pub p_bonferroni: f64,
    /// McDonald-Kreitman fixed nonsynonymous count (Dn): ALTs with AF >= threshold
    pub mk_dn: u32,
    /// McDonald-Kreitman fixed synonymous count (Ds)
    pub mk_ds: u32,
    /// McDonald-Kreitman polymorphic nonsynonymous count (Pn): ALTs with AF < threshold
    pub mk_pn: u32,
    /// McDonald-Kreitman polymorphic synonymous count (Ps)
    pub mk_ps: u32,
    /// Classified ALT alleles this gene kept OUT of the McDonald-Kreitman 2x2 and the
    /// SFS bins because they sit inside an observed multi-nucleotide codon change.
    /// Fisher's exact test and the SFS need whole alleles, and a multi-nucleotide change
    /// splits between the synonymous and nonsynonymous classes (see
    /// [`crate::vcf_analysis::mnv`]), so those alleles are excluded and counted here
    /// rather than rounded into a cell they do not belong in. Always 0 for an input
    /// without per-sample genotypes.
    pub mk_mnv_excluded: u32,
    /// Two-sided binomial −log10(p) in log space — finite even where `p_value`
    /// underflows to 0 (NaN when untested).
    pub neglog10p: f64,
    /// Lower/upper bound of a 95% Wilson confidence interval on pN/pS (NaN when
    /// untested or --af-weighted). `pn_ps_hi` is +∞ when the CI reaches pS = 0.
    pub pn_ps_lo: f64,
    pub pn_ps_hi: f64,
    /// Genomic-control-corrected p-value and BH q-value (NaN unless
    /// `--genomic-control`; the χ² statistic is divided by the inflation λ).
    pub p_gc: f64,
    pub q_gc: f64,
    /// True for repetitive / hard-to-map genes (PE/PPE/PGRS, IS elements, …).
    pub repetitive: bool,
    /// Site-frequency-spectrum counts of nonsynonymous / synonymous SNPs binned
    /// by allele frequency (see [`SFS_EDGES`]). Pooled genome-wide for the SFS
    /// panel; not written to the per-gene table. Like the McDonald-Kreitman 2x2 these
    /// need whole alleles, so they count only alleles that stand alone in their codon;
    /// `mk_mnv_excluded` is how many were left out.
    pub sfs_nonsyn: [u32; SFS_NBINS],
    pub sfs_syn: [u32; SFS_NBINS],
    /// Per-coding-SNP records behind this gene's counts, in genomic order. Written
    /// to `<output>_variants.tsv` under `--variants`. Includes nonsense/stop-loss
    /// changes (excluded from the pN/pS counts but important for resistance triage).
    pub variants: Vec<Variant>,
    /// Codons of this gene's reference CDS with at least one possible nonsynonymous
    /// single-nucleotide change (`A_c > 0`), and the sum of `A_c` over them. These are
    /// the gene's contribution to the per-codon recurrence scan's multiple-testing
    /// family and to the denominator of its plug-in rate; see `sites::codon_space`.
    /// Both are independent of `--kappa` and of how many SNPs landed in the gene.
    pub scan_codons: usize,
    pub scan_poss_nonsyn: u64,
}

/// The coding effect of a single ALT allele on its codon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnpEffect {
    /// Same amino acid (counts toward pS).
    Synonymous,
    /// Different amino acid, neither a stop (counts toward pN).
    Missense,
    /// A coding codon → stop (excluded from pN/pS site & SNP counts).
    Nonsense,
    /// A reference stop → a coding amino acid (excluded from pN/pS counts).
    StopLoss,
}

impl SnpEffect {
    pub fn label(self) -> &'static str {
        match self {
            SnpEffect::Synonymous => "synonymous",
            SnpEffect::Missense => "missense",
            SnpEffect::Nonsense => "nonsense",
            SnpEffect::StopLoss => "stop_loss",
        }
    }
}

/// One coding single-nucleotide variant located within a gene's CDS. Gene name,
/// chromosome and strand come from the parent [`GenePnPs`], so this stays lean.
#[derive(Debug, Clone)]
pub struct Variant {
    /// 1-based genomic position.
    pub pos: usize,
    /// Reference and alternate base on the genomic (VCF) strand.
    pub ref_allele: u8,
    pub alt_allele: u8,
    /// 1-based residue (codon) number within the protein.
    pub aa_pos: usize,
    /// The reference codon this ALT sits in, on the CODING strand (so it reads in the
    /// same orientation as `ref_aa`, not as the VCF's `ref_allele`). Carried here
    /// because the per-codon scan needs the codon's own mutational opportunity, which
    /// is a property of these three bases; recomputing it would mean reconstructing
    /// the CDS a second time.
    pub ref_codon: [u8; 3],
    /// The codon this ALT's carriers actually hold, on the coding strand.
    ///
    /// Equal to `ref_codon` with this ALT's base substituted, **except** where a sample
    /// carrying this ALT also carries another SNP in the same codon: then it is the
    /// joint (multi-nucleotide) codon that sample encodes, which is what `alt_aa`,
    /// `Change` and `effect` describe. Where the carriers disagree (some hold the ALT
    /// alone, some with a neighbour) it is the codon most of them hold. How many bases
    /// it changes is its distance from `ref_codon`, and `codon_shared` says the same
    /// thing as a flag.
    pub alt_codon: [u8; 3],
    /// Reference and alternate amino acid (one-letter; `*` = stop).
    pub ref_aa: u8,
    pub alt_aa: u8,
    /// Allele frequency of this ALT.
    pub af: f64,
    /// Exact derived-allele count for this ALT from the GT columns, when available.
    /// The diversity path uses this over `round(af * n)`, which loses precision when
    /// `af` came from an INFO/AF that disagrees with the genotypes.
    pub gt_derived: Option<usize>,
    /// Total called alleles at this site from the GT columns (the observed sample
    /// size), when available. The diversity path uses it as the per-site denominator
    /// so a no-call-heavy site is scored on the samples actually genotyped (matching
    /// the AF path) instead of being diluted by the full column count.
    pub gt_called: Option<usize>,
    pub effect: SnpEffect,
    /// Does a sample carrying this ALT also carry another SNP in the same codon?
    ///
    /// `Some(true)` means the amino-acid change on this row is the JOINT
    /// (multi-nucleotide) change that sample's genome encodes, not the effect of this base
    /// on its own; `alt_codon` spells it out. `Some(false)` means the row is an ordinary
    /// single substitution. `None` means unknowable from the input (an AF-only VCF whose
    /// frequencies also fall short of the forcing bound), in which case the row was scored
    /// against the reference codon and may be wrong.
    ///
    /// Written as the `Codon_Shared` column under `--variants --shared-codons`. Kept out
    /// of the default table so the existing one stays byte-identical for its parsers.
    pub codon_shared: Option<bool>,
    /// How many times this ALT arose independently on the phylogeny supplied with
    /// `--tree`, by Fitch parsimony over the samples carrying it, counting only origins
    /// whose clade carries it in at least `--min-origin-support` samples (see
    /// [`crate::tree`]).
    ///
    /// `None` without `--tree`, and also for an allele with no per-sample carriers,
    /// which is `NA` rather than zero. This is the one column in which *rpoB* S450L is
    /// visible for what it is: one allele that arose many times over.
    ///
    /// Written as the `Origins` column under `--variants --tree`.
    pub origins: Option<u32>,
}

/// Number of allele-frequency bins for the site-frequency-spectrum panel.
pub const SFS_NBINS: usize = 6;
/// Upper edges of the AF bins (a variant with AF ≤ edge falls in that bin).
pub const SFS_EDGES: [f64; SFS_NBINS] = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0];

/// Bin an allele frequency into `0..SFS_NBINS`.
fn sfs_bin(af: f64) -> usize {
    for (i, &e) in SFS_EDGES.iter().enumerate() {
        if af <= e {
            return i;
        }
    }
    SFS_NBINS - 1
}

/// Heuristic flag for repetitive / hard-to-map *M. tuberculosis* genes whose
/// SNP calls are frequently mapping artefacts (PE/PPE/PGRS, IS6110, maturases).
/// Mirrors the report's badge regex so the two never disagree.
pub fn is_repetitive(name: &str) -> bool {
    let n = name.to_ascii_uppercase();
    // PE/PPE family: the prefix must be followed by a digit, '_', or end of name
    // (PE13, PE_PGRS, PPE18) — so real genes like pepN/pepA are NOT flagged.
    let pe_family = |pre: &str| {
        n.strip_prefix(pre).is_some_and(|rest| {
            rest.is_empty() || rest.starts_with('_') || rest.starts_with(|c: char| c.is_ascii_digit())
        })
    };
    pe_family("PE") || pe_family("PPE") || ["PGRS", "MATURASE", "TRANSPOS", "IS6110"].iter().any(|m| n.contains(m))
}

/// Genome-wide (pooled) pN/pS aggregated across all analyzed genes.
///
/// Unlike averaging per-gene pN/pS ratios — which over-weights genes carrying
/// only a handful of sites — this pools SNP counts and site counts across the
/// whole coding genome *before* taking the ratio:
///
/// ```text
/// pN = Σ nonsyn_snps / Σ N_sites
/// pS = Σ syn_snps    / Σ S_sites
/// ```
///
/// This is the standard way to summarise the overall strength and direction of
/// selection over a set of genes. Counts are AF-weighted when the per-gene
/// results were computed with `--af-weighted` (making this πN/πS).
#[derive(Debug, Clone)]
pub struct GenomeWidePnPs {
    /// Total nonsynonymous sites summed over all genes
    pub n_sites: f64,
    /// Total synonymous sites summed over all genes
    pub s_sites: f64,
    /// Total nonsynonymous SNPs (AF-weighted if applicable)
    pub nonsyn_snps: f64,
    /// Total synonymous SNPs (AF-weighted if applicable)
    pub syn_snps: f64,
    /// Pooled pN (Σ nonsyn_snps / Σ N_sites)
    pub pn: f64,
    /// Pooled pS (Σ syn_snps / Σ S_sites)
    pub ps: f64,
    /// Pooled pN/pS ratio
    pub pn_ps: f64,
}

/// A short qualitative interpretation of a pN/pS ratio.
///
/// The 0.9–1.1 neutral band is a coarse convenience heuristic, not a
/// statistical test — genome-wide pN/pS is routinely below 1 for real
/// populations, and formal inference needs a null model.
pub fn selection_label(pn_ps: f64) -> &'static str {
    if pn_ps.is_nan() {
        "undetermined (no coding variation)"
    } else if pn_ps.is_infinite() {
        "no synonymous variation (ratio undefined)"
    } else if pn_ps < 0.9 {
        "purifying selection (pN/pS < 1)"
    } else if pn_ps > 1.1 {
        "positive/diversifying selection (pN/pS > 1)"
    } else {
        "near-neutral (pN/pS ~ 1)"
    }
}

/// Expected nonsynonymous fraction of mutations under neutrality: N/(N+S).
fn exp_n_frac(r: &GenePnPs) -> f64 {
    let sites = r.n_sites + r.s_sites;
    if sites > 0.0 {
        r.n_sites / sites
    } else {
        f64::NAN
    }
}

/// Format a p-value / probability for a text column: "NA" for NaN, scientific
/// for very small values (so significance is not rounded to 0), else 6 dp.
fn format_pval(v: f64) -> String {
    let v = if v == 0.0 { 0.0 } else { v }; // normalise -0.0 → 0.0
    if v.is_nan() {
        "NA".to_string()
    } else if v != 0.0 && v.abs() < 1e-3 {
        format!("{:.3e}", v)
    } else {
        format!("{:.6}", v)
    }
}

/// Format a number for JSON: `null` for non-finite, else a round-tripping
/// literal (preserves small p-values that {:.6} would flatten to 0).
fn format_json_num(v: f64) -> String {
    if v.is_finite() {
        // Normalise IEEE -0.0 → 0.0 (an empty diversity class sums to -0.0), so the
        // JSON matches the delimited writers and never emits a stray "-0".
        let v = if v == 0.0 { 0.0 } else { v };
        format!("{}", v)
    } else {
        "null".to_string()
    }
}

/// Format a pN/pS ratio for human-readable output, mapping NaN/Infinity to
/// stable textual tokens instead of Rust's default `NaN`/`inf` Display.
pub fn format_ratio(v: f64) -> String {
    let v = if v == 0.0 { 0.0 } else { v }; // normalise -0.0 → 0.0
    if v.is_nan() {
        "NaN".to_string()
    } else if v.is_infinite() {
        "inf".to_string()
    } else {
        format!("{:.6}", v)
    }
}

fn format_json_f64(v: f64) -> String {
    if v.is_nan() || v.is_infinite() {
        "null".to_string()
    } else if v == 0.0 {
        "0.000000".to_string()
    } else {
        format!("{:.6}", v)
    }
}

