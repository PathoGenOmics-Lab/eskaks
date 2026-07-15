//! Core pN/pS analysis from VCF data.
//!
//! For each gene, reconstructs reference and alternate codons from SNPs,
//! classifies mutations as synonymous or nonsynonymous using the genetic
//! code table, and computes pN/pS ratios.

use crate::gff::{Gene, Strand};
use crate::genetic_code::GeneticCode;
use crate::vcf::VcfSnp;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

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
    /// Genomic start position (1-based, for plotting/output)
    pub genome_start: usize,
    /// Genomic end position (1-based, max exon end)
    pub genome_end: usize,
    /// Strand ('+' or '-')
    pub strand: char,
    /// Chromosome
    pub chrom: String,
    /// Two-sided exact-binomial p-value for H0: pN/pS = 1 (NaN if untested:
    /// no SNPs, --af-weighted, or a degenerate expected fraction)
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

/// Compute pN/pS for all genes given a reference sequence, gene annotations, and SNPs.
///
/// If `af_weighted` is true, each SNP contributes its allele frequency to the
/// syn/nonsyn count instead of 1.0 (πN/πS instead of pN/pS).
///
/// `kappa` is the transition/transversion rate ratio used when counting N and S
/// SITES. At `kappa == 1.0` the classic equal-rates Nei-Gojobori counting is
/// used unchanged; any other value activates mutation-spectrum-weighted counting
/// (see [`count_sites_weighted`]). Only the site denominators are affected —
/// the observed-SNP syn/nonsyn classification (the numerators) is empirical and
/// never rate-weighted.
pub fn compute_pn_ps(
    reference: &HashMap<String, Vec<u8>>,
    genes: &[Gene],
    snps: &[VcfSnp],
    gc: &GeneticCode,
    af_weighted: bool,
    kappa: f64,
    mk_fixed_af: f64,
) -> Vec<GenePnPs> {
    // kappa == 1 reproduces the original counting byte-for-byte; only a
    // non-neutral kappa switches to the rate-weighted path. (count_sites_weighted
    // is provably identical at kappa == 1, so this gate is purely defensive.)
    let spectrum_weighted = kappa != 1.0;
    // Index SNPs by chromosome and position for fast lookup
    let mut snp_map: HashMap<&str, Vec<&VcfSnp>> = HashMap::new();
    for snp in snps {
        snp_map.entry(snp.chrom.as_str()).or_default().push(snp);
    }
    // Sort SNPs by position within each chromosome
    for snps_list in snp_map.values_mut() {
        snps_list.sort_by_key(|s| s.pos);
    }

    // Aggregate diagnostics so a wrong reference doesn't spam one line per SNP.
    // Atomics let the per-gene work run in parallel; par_iter().collect()
    // preserves gene order, so results are deterministic regardless of threads.
    let ref_checked = AtomicUsize::new(0);
    let ref_mismatch = AtomicUsize::new(0);

    let results: Vec<GenePnPs> = genes
        .par_iter()
        .filter_map(|gene| {
            let ref_seq = match reference.get(&gene.seqid) {
                Some(seq) => seq,
                None => {
                    warn!("Reference sequence not found for {}, skipping gene {}", gene.seqid, gene.name);
                    return None;
                }
            };

            let chrom_snps = snp_map.get(gene.seqid.as_str());

            // Extract the full CDS sequence from the reference
            let cds_seq = extract_cds_sequence(gene, ref_seq);
            if cds_seq.len() < 3 {
                warn!("Gene {} CDS too short ({} bp), skipping", gene.name, cds_seq.len());
                return None;
            }

            // Count S and N sites from reference codons
            let (n_sites, s_sites) = if spectrum_weighted {
                count_sites_weighted(&cds_seq, gc, kappa)
            } else {
                count_sites(&cds_seq, gc)
            };

            // Find SNPs that fall within this gene's CDS regions
            // Counts are f64 to support AF-weighted mode (πN/πS)
            let mut nonsyn_count = 0.0f64;
            let mut syn_count = 0.0f64;
            let mut local_checked = 0usize;
            let mut local_mismatch = 0usize;
            // McDonald-Kreitman: fixed (AF >= threshold) vs polymorphic counts.
            let (mut mk_dn, mut mk_ds, mut mk_pn, mut mk_ps) = (0u32, 0u32, 0u32, 0u32);

            if let Some(snps_list) = chrom_snps {
                // SNPs are position-sorted, so binary-search the gene's genomic
                // span instead of scanning the whole chromosome for every gene
                // (O(S + G·log S) rather than O(G·S)). genomic_to_cds_offset
                // still filters precisely per exon, so results are byte-identical.
                let lo = gene.exons.iter().map(|e| e.start).min().unwrap_or(gene.start);
                let hi = gene.exons.iter().map(|e| e.end).max().unwrap_or(gene.start);
                let from = snps_list.partition_point(|s| s.pos < lo);
                let to = snps_list.partition_point(|s| s.pos <= hi);
                for snp in &snps_list[from..to] {
                    // Check if this SNP falls within any exon of this gene
                    if let Some(cds_offset) = genomic_to_cds_offset(gene, snp.pos) {
                        let codon_idx = cds_offset / 3;
                        let pos_in_codon = cds_offset % 3;
                        let codon_start = codon_idx * 3;

                        if codon_start + 3 > cds_seq.len() {
                            continue;
                        }

                        // Get reference codon from the extracted CDS
                        let ref_codon = [
                            cds_seq[codon_start],
                            cds_seq[codon_start + 1],
                            cds_seq[codon_start + 2],
                        ];

                        // Verify VCF REF allele matches the reference sequence
                        let expected_ref = if gene.strand == Strand::Minus {
                            complement(cds_seq[codon_start + pos_in_codon])
                        } else {
                            cds_seq[codon_start + pos_in_codon]
                        };
                        local_checked += 1;
                        if snp.ref_allele != expected_ref {
                            local_mismatch += 1;
                            debug!(
                                "VCF REF mismatch at {}:{} — VCF says {}, reference has {}. Skipping.",
                                snp.chrom, snp.pos,
                                snp.ref_allele as char, expected_ref as char
                            );
                            continue;
                        }

                        // For each ALT allele at this position
                        for (alt_idx, alt_base) in snp.alt_alleles.iter().enumerate() {
                            // Weight: AF if weighted mode, 1.0 otherwise
                            let weight = if af_weighted {
                                snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0)
                            } else {
                                1.0
                            };

                            // Build alternate codon
                            let mut alt_codon = ref_codon;
                            let alt_in_cds = if gene.strand == Strand::Minus {
                                complement(*alt_base)
                            } else {
                                *alt_base
                            };
                            alt_codon[pos_in_codon] = alt_in_cds;

                            // Look up amino acids
                            let ref_aa = codon_to_aa(&ref_codon, gc);
                            let alt_aa = codon_to_aa(&alt_codon, gc);

                            match (ref_aa, alt_aa) {
                                (Some(r), Some(a)) if r != b'*' && a != b'*' => {
                                    // MK split uses this ALT's allele frequency,
                                    // as raw counts (not AF-weighted).
                                    let af = snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0);
                                    let fixed = af >= mk_fixed_af;
                                    if r == a {
                                        syn_count += weight;
                                        if fixed { mk_ds += 1 } else { mk_ps += 1 }
                                    } else {
                                        nonsyn_count += weight;
                                        if fixed { mk_dn += 1 } else { mk_pn += 1 }
                                    }
                                }
                                // Skip: ambiguous codons or mutations to/from stop
                                _ => continue,
                            }
                        }
                    }
                }
            }

            ref_checked.fetch_add(local_checked, Ordering::Relaxed);
            ref_mismatch.fetch_add(local_mismatch, Ordering::Relaxed);

            let total_snps = nonsyn_count + syn_count;
            let pn = if n_sites > 0.0 {
                nonsyn_count / n_sites
            } else {
                0.0
            };
            let ps = if s_sites > 0.0 {
                syn_count / s_sites
            } else {
                0.0
            };
            let pn_ps = if ps > 0.0 {
                pn / ps
            } else if pn > 0.0 {
                f64::INFINITY
            } else {
                f64::NAN
            };

            // Per-gene neutrality test: under H0 (pN/pS = 1) the nonsynonymous
            // fraction of SNPs equals the mutational opportunity N/(N+S). Only
            // valid with integer counts, so skip it under --af-weighted.
            let sites = n_sites + s_sites;
            let p_value = if !af_weighted && total_snps > 0.0 && sites > 0.0 {
                crate::stats::binomial_two_sided_p(
                    nonsyn_count.round() as u64,
                    total_snps.round() as u64,
                    n_sites / sites,
                )
            } else {
                f64::NAN
            };

            let genome_end = gene.exons.iter().map(|e| e.end).max().unwrap_or(gene.start);
            let strand = if gene.strand == Strand::Minus { '-' } else { '+' };

            Some(GenePnPs {
                name: gene.name.clone(),
                length_bp: gene.length_bp,
                n_sites,
                s_sites,
                pn,
                ps,
                pn_ps,
                nonsyn_snps: nonsyn_count,
                syn_snps: syn_count,
                total_snps,
                genome_start: gene.start,
                genome_end,
                strand,
                chrom: gene.seqid.clone(),
                p_value,
                q_value: f64::NAN,
                p_bonferroni: f64::NAN,
                mk_dn,
                mk_ds,
                mk_pn,
                mk_ps,
            })
        })
        .collect();

    let ref_checked = ref_checked.load(Ordering::Relaxed);
    let ref_mismatch = ref_mismatch.load(Ordering::Relaxed);

    if ref_mismatch > 0 {
        let frac = ref_mismatch as f64 / ref_checked.max(1) as f64;
        if frac > 0.10 {
            warn!(
                "{}/{} in-CDS SNPs ({:.1}%) had a REF allele disagreeing with the reference and were skipped — \
                 this usually means a wrong reference build, an off-by-one, or a strand mix-up. Run with -vv for positions.",
                ref_mismatch, ref_checked, 100.0 * frac
            );
        } else {
            info!(
                "{}/{} in-CDS SNPs had a REF-allele mismatch and were skipped (run with -vv for positions).",
                ref_mismatch, ref_checked
            );
        }
    }

    info!("Computed pN/pS for {} genes", results.len());
    results
}

/// Percentile bootstrap confidence interval for the genome-wide pooled pN/pS,
/// resampling genes with replacement `n_boot` times (seeded for reproducibility).
/// Returns `(lo, hi)` at the given two-sided confidence, or None if there is no
/// data / no finite replicate. Pools counts and sites within each replicate,
/// exactly like [`genome_wide_pn_ps`].
pub fn bootstrap_genome_wide_ci(
    results: &[GenePnPs],
    n_boot: usize,
    seed: u64,
    confidence: f64,
) -> Option<(f64, f64)> {
    if results.is_empty() || n_boot == 0 {
        return None;
    }
    let n = results.len();
    let mut rng = crate::stats::SplitMix64::new(seed);
    let mut ratios: Vec<f64> = Vec::with_capacity(n_boot);
    let mut undefined = 0usize;
    for _ in 0..n_boot {
        let (mut n_sites, mut s_sites, mut nonsyn, mut syn) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for _ in 0..n {
            let r = &results[rng.below(n)];
            n_sites += r.n_sites;
            s_sites += r.s_sites;
            nonsyn += r.nonsyn_snps;
            syn += r.syn_snps;
        }
        let pn = if n_sites > 0.0 { nonsyn / n_sites } else { 0.0 };
        let ps = if s_sites > 0.0 { syn / s_sites } else { 0.0 };
        let ratio = if ps > 0.0 { pn / ps } else { f64::NAN };
        if ratio.is_finite() {
            ratios.push(ratio);
        } else {
            undefined += 1;
        }
    }
    if ratios.is_empty() {
        return None;
    }
    // Replicates with no synonymous variation give an undefined ratio and are
    // excluded from the percentiles, which truncates the upper tail. Warn so a
    // biased-low CI on sparse-synonymous data isn't read as precise.
    if undefined > 0 {
        warn!(
            "{}/{} bootstrap replicates had no synonymous variation (undefined pN/pS) and were \
             excluded; the CI may be biased low when synonymous SNPs are sparse.",
            undefined, n_boot
        );
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tail = (1.0 - confidence) / 2.0 * 100.0;
    Some((
        crate::stats::percentile_sorted(&ratios, tail),
        crate::stats::percentile_sorted(&ratios, 100.0 - tail),
    ))
}

/// Fill in Benjamini-Hochberg FDR q-values and Bonferroni-corrected p-values
/// for the per-gene neutrality test, across every gene with a finite p-value.
/// Genes not tested (no SNPs / --af-weighted) keep NaN.
pub fn apply_multiple_testing(results: &mut [GenePnPs]) {
    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let qvals = crate::stats::benjamini_hochberg(&pvals);
    let bonf = crate::stats::bonferroni(&pvals);
    for (r, (q, b)) in results.iter_mut().zip(qvals.into_iter().zip(bonf)) {
        r.q_value = q;
        r.p_bonferroni = b;
    }
}

/// Aggregate per-gene results into a single genome-wide (pooled) pN/pS estimate.
///
/// Returns `None` when there are no genes to aggregate. The pooled ratio uses
/// the same NaN/Infinity conventions as the per-gene computation: `NaN` when
/// there is no variation at all, `+Infinity` when there are nonsynonymous but
/// no synonymous changes.
pub fn genome_wide_pn_ps(results: &[GenePnPs]) -> Option<GenomeWidePnPs> {
    if results.is_empty() {
        return None;
    }

    let n_sites: f64 = results.iter().map(|r| r.n_sites).sum();
    let s_sites: f64 = results.iter().map(|r| r.s_sites).sum();
    let nonsyn_snps: f64 = results.iter().map(|r| r.nonsyn_snps).sum();
    let syn_snps: f64 = results.iter().map(|r| r.syn_snps).sum();

    let pn = if n_sites > 0.0 { nonsyn_snps / n_sites } else { 0.0 };
    let ps = if s_sites > 0.0 { syn_snps / s_sites } else { 0.0 };
    let pn_ps = if ps > 0.0 {
        pn / ps
    } else if pn > 0.0 {
        f64::INFINITY
    } else {
        f64::NAN
    };

    Some(GenomeWidePnPs {
        n_sites,
        s_sites,
        nonsyn_snps,
        syn_snps,
        pn,
        ps,
        pn_ps,
    })
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

/// Extract the CDS sequence from the reference, handling strand and multi-exon genes.
fn extract_cds_sequence(gene: &Gene, ref_seq: &[u8]) -> Vec<u8> {
    let mut cds = Vec::with_capacity(gene.length_bp);

    // Exons are already sorted in coding order (ascending for +, descending for -)
    for exon in &gene.exons {
        let start_idx = exon.start.saturating_sub(1); // Convert 1-based to 0-based
        let end_idx = exon.end.min(ref_seq.len());
        if start_idx >= ref_seq.len() || start_idx >= end_idx {
            continue;
        }
        cds.extend_from_slice(&ref_seq[start_idx..end_idx]);
    }

    if gene.strand == Strand::Minus {
        // Reverse complement
        cds = reverse_complement(&cds);
    }

    // Apply phase: skip leading bases
    if !gene.exons.is_empty() && gene.exons[0].phase > 0 {
        let skip = gene.exons[0].phase as usize;
        if skip < cds.len() {
            cds = cds[skip..].to_vec();
        }
    }

    // Uppercase
    for b in cds.iter_mut() {
        *b = b.to_ascii_uppercase();
    }

    // A CDS whose length is not a multiple of 3 signals a frame/annotation
    // problem (bad phase, wrong exon bounds, or an origin-spanning gene on a
    // circular genome). The trailing 1-2 bases are dropped by the len/3 loops;
    // flag it rather than silently skewing the site counts.
    if !cds.is_empty() && cds.len() % 3 != 0 {
        warn!(
            "Gene {} CDS length {} bp is not a multiple of 3 (phase/annotation issue?); trailing {} base(s) ignored",
            gene.name, cds.len(), cds.len() % 3
        );
    }

    cds
}

/// Convert a genomic position (1-based) to a CDS offset (0-based).
/// Returns None if the position is not within any exon of the gene.
fn genomic_to_cds_offset(gene: &Gene, pos: usize) -> Option<usize> {
    let mut cds_offset = 0usize;

    // Handle phase offset from first exon
    let phase_offset = if !gene.exons.is_empty() {
        gene.exons[0].phase as usize
    } else {
        0
    };

    match gene.strand {
        Strand::Plus => {
            for exon in &gene.exons {
                if pos >= exon.start && pos <= exon.end {
                    let offset = cds_offset + (pos - exon.start);
                    if offset >= phase_offset {
                        return Some(offset - phase_offset);
                    } else {
                        return None; // Within phase region
                    }
                }
                cds_offset += exon.end - exon.start + 1;
            }
        }
        Strand::Minus => {
            for exon in &gene.exons {
                // For minus strand, exons are sorted descending by start
                if pos >= exon.start && pos <= exon.end {
                    let offset = cds_offset + (exon.end - pos);
                    if offset >= phase_offset {
                        return Some(offset - phase_offset);
                    } else {
                        return None;
                    }
                }
                cds_offset += exon.end - exon.start + 1;
            }
        }
    }

    None
}

/// Count nonsynonymous (N) and synonymous (S) sites for a CDS sequence.
///
/// For each codon, enumerates all 9 possible single-nucleotide changes
/// and classifies each as synonymous or nonsynonymous. Fractional sites
/// are assigned as: S_sites = syn_changes/3, N_sites = nonsyn_changes/3.
fn count_sites(cds: &[u8], gc: &GeneticCode) -> (f64, f64) {
    let mut n_sites = 0.0f64;
    let mut s_sites = 0.0f64;

    let codons = cds.len() / 3;
    for i in 0..codons {
        let codon = [cds[i * 3], cds[i * 3 + 1], cds[i * 3 + 2]];
        let ref_aa = codon_to_aa(&codon, gc);
        if ref_aa.is_none() {
            continue; // Skip ambiguous codons
        }
        let ref_aa = ref_aa.unwrap();

        if ref_aa == b'*' {
            continue; // Skip stop codons
        }

        let mut syn = 0;
        let mut nonsyn = 0;

        for pos in 0..3 {
            for &alt_base in b"ACGT" {
                if alt_base == codon[pos] {
                    continue;
                }
                let mut alt_codon = codon;
                alt_codon[pos] = alt_base;
                if let Some(alt_aa) = codon_to_aa(&alt_codon, gc) {
                    if alt_aa == b'*' {
                        continue; // Exclude changes to stop codons from site counts
                    } else if alt_aa == ref_aa {
                        syn += 1;
                    } else {
                        nonsyn += 1;
                    }
                }
                // Skip ambiguous codons entirely (don't count as either)
            }
        }

        // Each codon contributes exactly 3 sites. With stop-codon changes
        // excluded, redistribute proportionally: S = 3 × syn/(syn+nonsyn)
        let valid_changes = (syn + nonsyn) as f64;
        if valid_changes > 0.0 {
            s_sites += 3.0 * syn as f64 / valid_changes;
            n_sites += 3.0 * nonsyn as f64 / valid_changes;
        }
    }

    (n_sites, s_sites)
}

/// Is the change `from`→`to` a transition (A↔G or C↔T)? Purine↔purine or
/// pyrimidine↔pyrimidine. Strand-symmetric, so it holds on either strand.
#[inline]
fn is_transition(from: u8, to: u8) -> bool {
    matches!(
        (from, to),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// Count N and S sites with mutation-spectrum weighting by a transition/
/// transversion rate ratio `kappa`.
///
/// This is the mutation-spectrum-aware generalisation of [`count_sites`]: each
/// candidate single-nucleotide change is weighted by its relative mutation rate
/// (`kappa` for a transition, `1.0` for a transversion) instead of counting 1.
/// It uses the *same* codon-level normalisation as [`count_sites`] — each codon
/// contributes exactly 3 sites, split between S and N by the rate-weighted
/// synonymous fraction — so at `kappa == 1.0` it reduces to [`count_sites`]
/// bit-for-bit (the weights become all-1 and the sums equal the plain counts).
/// Keeping the normalisation identical means `kappa` is the *only* thing that
/// changes, so the correction is not confounded with a normalisation switch.
///
/// Under a transition-biased spectrum (`kappa > 1`), synonymous changes at
/// 2-fold degenerate sites — reached almost exclusively by transitions — get
/// up-weighted, so a 2-fold synonymous-transition site moves from `1/3` toward
/// `kappa/(kappa+2)`. Fully (4-fold) degenerate positions are synonymous for
/// every change, so they stay `kappa`-invariant. Across the coding genome this
/// generally raises total S and lowers total N (individual codons can move
/// either way depending on whether their transition-reachable changes are
/// mostly synonymous or nonsynonymous).
fn count_sites_weighted(cds: &[u8], gc: &GeneticCode, kappa: f64) -> (f64, f64) {
    let mut n_sites = 0.0f64;
    let mut s_sites = 0.0f64;

    let codons = cds.len() / 3;
    for i in 0..codons {
        let codon = [cds[i * 3], cds[i * 3 + 1], cds[i * 3 + 2]];
        let ref_aa = match codon_to_aa(&codon, gc) {
            Some(aa) if aa != b'*' => aa,
            _ => continue, // Skip ambiguous and stop codons
        };

        // Pool rate-weighted synonymous / total over all three positions of the
        // codon, exactly as count_sites pools raw counts.
        let mut syn_rate = 0.0f64;
        let mut tot_rate = 0.0f64;
        for pos in 0..3 {
            for &alt_base in b"ACGT" {
                if alt_base == codon[pos] {
                    continue;
                }
                let mut alt_codon = codon;
                alt_codon[pos] = alt_base;
                if let Some(alt_aa) = codon_to_aa(&alt_codon, gc) {
                    if alt_aa == b'*' {
                        continue; // Exclude changes to stop codons
                    }
                    let w = if is_transition(codon[pos], alt_base) { kappa } else { 1.0 };
                    tot_rate += w;
                    if alt_aa == ref_aa {
                        syn_rate += w;
                    }
                }
                // Ambiguous alternates are skipped (contribute to neither)
            }
        }

        // Each codon contributes exactly 3 sites, split by the rate-weighted
        // synonymous fraction (reduces to count_sites at kappa == 1.0).
        if tot_rate > 0.0 {
            s_sites += 3.0 * syn_rate / tot_rate;
            n_sites += 3.0 * (tot_rate - syn_rate) / tot_rate;
        }
    }

    (n_sites, s_sites)
}

/// Look up amino acid for a codon using the genetic code table (Li indexing).
/// Returns None for codons with ambiguous bases.
fn codon_to_aa(codon: &[u8; 3], gc: &GeneticCode) -> Option<u8> {
    let b1 = base_to_li(codon[0])?;
    let b2 = base_to_li(codon[1])?;
    let b3 = base_to_li(codon[2])?;
    let idx = 16 * b1 + 4 * b2 + b3;
    Some(gc.aa_table[idx])
}

/// Convert a base to Li index: A=0, C=1, G=2, T=3.
#[inline]
fn base_to_li(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Complement a DNA base.
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'U' | b'u' => b'A',
        _ => b'N',
    }
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Parse a reference FASTA file into a map of sequence_name -> sequence.
/// Uses simple manual parsing (no needletail dependency for this).
pub fn parse_reference_fasta(path: &std::path::Path) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path)
        .with_context(|| format!("Failed to open reference FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                seqs.insert(name, std::mem::take(&mut current_seq));
            }
            // Take the first word as the sequence name
            let name = line.strip_prefix('>').unwrap_or("").split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
            current_seq = Vec::new();
        } else {
            current_seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(name) = current_name {
        seqs.insert(name, current_seq);
    }

    if seqs.is_empty() {
        anyhow::bail!("No sequences found in reference FASTA: {}", path.display());
    }

    info!("Loaded {} reference sequences", seqs.len());
    Ok(seqs)
}

/// Write pN/pS results to a file.
pub fn write_results(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_pnps.{}", prefix, ext);

    match format {
        crate::models::OutputFormat::Json => {
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(file, "[")?;
            for (i, r) in results.iter().enumerate() {
                let comma = if i + 1 < results.len() { "," } else { "" };
                writeln!(
                    file,
                    "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"length_bp\":{},\"N_sites\":{:.4},\"S_sites\":{:.4},\"exp_N_frac\":{},\"pN\":{},\"pS\":{},\"pN_pS\":{},\"nonsyn_snps\":{:.4},\"syn_snps\":{:.4},\"total_snps\":{:.4},\"p_value\":{},\"q_value_bh\":{},\"p_bonferroni\":{}}}{}",
                    r.name, r.chrom, r.genome_start, r.genome_end, r.strand, r.length_bp,
                    r.n_sites, r.s_sites, format_json_num(exp_n_frac(r)),
                    format_json_f64(r.pn), format_json_f64(r.ps), format_json_f64(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    format_json_num(r.p_value), format_json_num(r.q_value), format_json_num(r.p_bonferroni),
                    comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(
                file,
                "Gene{s}Length_bp{s}N_sites{s}S_sites{s}pN{s}pS{s}pN/pS{s}Nonsyn_SNPs{s}Syn_SNPs{s}Total_SNPs{s}Chrom{s}Start{s}End{s}Strand{s}Exp_N_frac{s}P_value{s}Q_value_BH{s}P_Bonferroni",
                s = sep
            )?;
            for r in results {
                writeln!(
                    file,
                    "{}{s}{}{s}{:.4}{s}{:.4}{s}{:.6}{s}{:.6}{s}{}{s}{:.4}{s}{:.4}{s}{:.4}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    r.name, r.length_bp, r.n_sites, r.s_sites,
                    r.pn, r.ps, format_ratio(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    r.chrom, r.genome_start, r.genome_end, r.strand,
                    format_pval(exp_n_frac(r)), format_pval(r.p_value),
                    format_pval(r.q_value), format_pval(r.p_bonferroni),
                    s = sep
                )?;
            }
        }
    }

    Ok(output_path)
}

/// Write per-gene McDonald-Kreitman results: the 2×2 fixed/polymorphic table
/// (Dn, Ds, Pn, Ps), Neutrality Index, alpha (proportion of adaptive
/// substitutions), a two-sided Fisher exact p-value, and its BH-FDR q-value.
///
/// "Fixed" means AF >= the chosen threshold *within the sample* (reference-
/// polarized MK); this conflates high-frequency derived alleles with true
/// between-species divergence, so it is a screen, not a substitute for an
/// outgroup-based MK test.
pub fn write_mk_results(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Only genes with at least one classified difference are informative.
    let genes: Vec<&GenePnPs> = results
        .iter()
        .filter(|r| r.mk_dn + r.mk_ds + r.mk_pn + r.mk_ps > 0)
        .collect();

    // Two-sided Fisher p per gene, then Benjamini-Hochberg across tested genes.
    let pvals: Vec<f64> = genes
        .iter()
        .map(|r| {
            crate::stats::fisher_exact_two_sided(
                r.mk_dn as u64,
                r.mk_ds as u64,
                r.mk_pn as u64,
                r.mk_ps as u64,
            )
        })
        .collect();
    let qvals = crate::stats::benjamini_hochberg(&pvals);

    let ext = format.extension();
    let output_path = format!("{}_mk.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);

    let ni = |r: &GenePnPs| {
        let (dn, ds, pn, ps) = (r.mk_dn as f64, r.mk_ds as f64, r.mk_pn as f64, r.mk_ps as f64);
        if ps > 0.0 && dn > 0.0 {
            (pn * ds) / (ps * dn)
        } else {
            f64::NAN
        }
    };
    let alpha = |r: &GenePnPs| {
        let (dn, ds, pn, ps) = (r.mk_dn as f64, r.mk_ds as f64, r.mk_pn as f64, r.mk_ps as f64);
        if dn > 0.0 && ps > 0.0 {
            1.0 - (ds * pn) / (dn * ps)
        } else {
            f64::NAN
        }
    };

    if let crate::models::OutputFormat::Json = format {
        writeln!(file, "[")?;
        for (i, r) in genes.iter().enumerate() {
            let comma = if i + 1 < genes.len() { "," } else { "" };
            writeln!(
                file,
                "  {{\"gene\":\"{}\",\"chrom\":\"{}\",\"start\":{},\"end\":{},\"strand\":\"{}\",\"Dn\":{},\"Ds\":{},\"Pn\":{},\"Ps\":{},\"NI\":{},\"alpha\":{},\"fisher_p\":{},\"fisher_q_bh\":{}}}{}",
                r.name, r.chrom, r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_json_num(ni(r)), format_json_num(alpha(r)),
                format_json_num(pvals[i]), format_json_num(qvals[i]), comma
            )?;
        }
        writeln!(file, "]")?;
    } else {
        let sep = format.separator();
        writeln!(
            file,
            "Gene{s}Chrom{s}Start{s}End{s}Strand{s}Dn{s}Ds{s}Pn{s}Ps{s}NI{s}alpha{s}Fisher_p{s}Fisher_q_BH",
            s = sep
        )?;
        for (i, r) in genes.iter().enumerate() {
            writeln!(
                file,
                "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                r.name, r.chrom, r.genome_start, r.genome_end, r.strand,
                r.mk_dn, r.mk_ds, r.mk_pn, r.mk_ps,
                format_pval(ni(r)), format_pval(alpha(r)),
                format_pval(pvals[i]), format_pval(qvals[i]),
                s = sep
            )?;
        }
    }

    Ok(output_path)
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
        format!("{}", v)
    } else {
        "null".to_string()
    }
}

/// Format a pN/pS ratio for human-readable output, mapping NaN/Infinity to
/// stable textual tokens instead of Rust's default `NaN`/`inf` Display.
pub fn format_ratio(v: f64) -> String {
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

/// Generate an SVG Manhattan-style plot of pN/pS per gene along the genome.
pub fn write_pnps_plot(results: &[GenePnPs], prefix: &str, fdr: f64) -> anyhow::Result<String> {
    use std::fmt::Write as FmtWrite;
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Color constants — passed as format arguments so they never appear
    // literally inside r#"..."# delimiters (which `"#` would terminate).
    const C_GRID: &str = "#e0e0e0";
    const C_POSITIVE: &str = "#d94a4a";
    const C_PURIFYING: &str = "#4a90d9";
    const C_AXIS: &str = "#333333";

    let plot_path = format!("{}_pnps_manhattan.svg", prefix);

    // Filter out genes with no SNPs or NaN/infinite pN/pS
    let plot_data: Vec<&GenePnPs> = results
        .iter()
        .filter(|r| r.total_snps > 0.0 && r.pn_ps.is_finite())
        .collect();

    // Genes with variation but an undefined ratio (pS = 0, i.e. only
    // nonsynonymous SNPs — often the strongest positive-selection candidates)
    // can't be placed on the ratio axis. They stay in the TSV, but say so
    // rather than dropping them from the plot without notice.
    let undefined = results
        .iter()
        .filter(|r| r.total_snps > 0.0 && !r.pn_ps.is_finite())
        .count();
    if undefined > 0 {
        warn!(
            "{} gene(s) with only nonsynonymous variation (pN/pS = ∞) are not shown on the Manhattan plot; they are in the {}_pnps table.",
            undefined, prefix
        );
    }

    if plot_data.is_empty() {
        info!("No valid data points for pN/pS plot");
        return Ok(plot_path);
    }

    let width = 900.0f64;
    let height = 500.0f64;
    let margin_top = 50.0f64;
    let margin_right = 40.0f64;
    let margin_bottom = 80.0f64;
    let margin_left = 80.0f64;
    let plot_w = width - margin_left - margin_right;
    let plot_h = height - margin_top - margin_bottom;

    let max_pos = plot_data.iter().map(|r| r.genome_start).max().unwrap_or(1) as f64;
    let min_pos = plot_data.iter().map(|r| r.genome_start).min().unwrap_or(0) as f64;
    let pos_range = (max_pos - min_pos).max(1.0);

    let max_y = plot_data
        .iter()
        .map(|r| r.pn_ps)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.5)
        * 1.1;

    let to_x = |pos: usize| margin_left + ((pos as f64 - min_pos) / pos_range) * plot_w;
    let to_y = |v: f64| margin_top + plot_h * (1.0 - v / max_y);

    let mut svg = String::with_capacity(4096);

    // Header
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" width=\"{w}\" height=\"{h}\">",
        w = width, h = height
    );
    let _ = writeln!(svg, "<style>");
    let _ = writeln!(svg, "  text {{ font-family: sans-serif; fill: {c}; }}", c = C_AXIS);
    let _ = writeln!(svg, "  .title {{ font-size: 16px; font-weight: bold; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .axis-label {{ font-size: 12px; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .tick-label {{ font-size: 10px; }}");
    let _ = writeln!(svg, "</style>");
    let _ = writeln!(svg, "<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>", w = width, h = height);
    let _ = writeln!(
        svg,
        "<text x=\"{cx}\" y=\"30\" class=\"title\">pN/pS per Gene (Manhattan Plot)</text>",
        cx = width / 2.0
    );

    // Y-axis grid lines and labels
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = max_y * frac;
        let y = to_y(val);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"0.5\"/>",
            margin_left, y, margin_left + plot_w, y, C_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x}\" y=\"{y:.1}\" class=\"tick-label\" text-anchor=\"end\" dominant-baseline=\"middle\">{val:.2}</text>",
            x = margin_left - 8.0, y = y, val = val
        );
    }

    // Neutral line at pN/pS = 1.0
    if 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\" stroke-dasharray=\"6,3\"/>",
            margin_left, y1, margin_left + plot_w, y1, C_POSITIVE
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"tick-label\" fill=\"{c}\">pN/pS = 1</text>",
            x = margin_left + plot_w + 3.0, y = y1 + 3.0, c = C_POSITIVE
        );
    }

    // Data points
    for r in &plot_data {
        let x = to_x(r.genome_start);
        let y = to_y(r.pn_ps);
        let color = if r.pn_ps < 1.0 { C_PURIFYING } else { C_POSITIVE };
        let radius = r.total_snps.sqrt().clamp(2.0, 8.0);
        // Outline genes significant in the per-gene neutrality test (BH-FDR).
        let sig = r.q_value.is_finite() && r.q_value < fdr;
        let stroke = if sig {
            " stroke=\"#000000\" stroke-width=\"1.5\""
        } else {
            ""
        };
        let _ = writeln!(
            svg,
            "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{r:.1}\" fill=\"{c}\" opacity=\"0.7\"{stroke}>",
            x = x, y = y, r = radius, c = color, stroke = stroke
        );
        let _ = writeln!(
            svg,
            "  <title>{name}: pN/pS={ratio:.4} ({s}S/{n}N SNPs), q={q}</title>",
            name = r.name, ratio = r.pn_ps, s = r.syn_snps, n = r.nonsyn_snps,
            q = format_pval(r.q_value)
        );
        svg.push_str("</circle>\n");
    }

    // Axes
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top, margin_left, margin_top + plot_h, C_AXIS
    );
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top + plot_h, margin_left + plot_w, margin_top + plot_h, C_AXIS
    );

    // Axis labels
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"axis-label\">Genome Position</text>",
        x = margin_left + plot_w / 2.0, y = height - 10.0
    );
    let _ = writeln!(
        svg,
        "<text x=\"15\" y=\"{y:.1}\" class=\"axis-label\" transform=\"rotate(-90,15,{y:.1})\">pN/pS</text>",
        y = margin_top + plot_h / 2.0
    );

    // Legend
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"40\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_PURIFYING
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"50\" class=\"tick-label\">Purifying (pN/pS &lt; 1)</text>",
        x = width - 184.0
    );
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"56\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_POSITIVE
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"66\" class=\"tick-label\">Positive (pN/pS &ge; 1)</text>",
        x = width - 184.0
    );

    svg.push_str("</svg>\n");

    let mut file = BufWriter::new(File::create(&plot_path)?);
    file.write_all(svg.as_bytes())?;

    Ok(plot_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_gc() -> &'static GeneticCode {
        crate::genetic_code::get_table(1).unwrap()
    }

    #[test]
    fn test_codon_to_aa() {
        let gc = make_gc();
        // ATG = Met
        assert_eq!(codon_to_aa(&[b'A', b'T', b'G'], gc), Some(b'M'));
        // TAA = Stop
        assert_eq!(codon_to_aa(&[b'T', b'A', b'A'], gc), Some(b'*'));
        // GCT = Ala
        assert_eq!(codon_to_aa(&[b'G', b'C', b'T'], gc), Some(b'A'));
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"AATTCC"), b"GGAATT");
    }

    #[test]
    fn test_count_sites() {
        let gc = make_gc();
        // ATG GCT = Met Ala (2 codons)
        let cds = b"ATGGCT";
        let (n, s) = count_sites(cds, gc);
        // Each codon contributes exactly 3 sites total (S + N = 3 per codon, 6 total)
        let total = n + s;
        assert!((total - 6.0).abs() < 0.01, "Total sites should be 6.0, got {}", total);
        // ATG (Met): 0 syn, 8 nonsyn (1 change → stop TAA excluded) → N=3*8/8=3, S=0
        // GCT (Ala): 3 syn (GCA,GCC,GCG), 6 nonsyn, 0 stops → N=3*6/9=2, S=3*3/9=1
        assert!(n > 4.5 && n < 5.5, "N_sites: expected ~5.0, got {}", n);
        assert!(s > 0.5 && s < 1.5, "S_sites: expected ~1.0, got {}", s);
    }

    #[test]
    fn test_is_transition() {
        assert!(is_transition(b'A', b'G'));
        assert!(is_transition(b'G', b'A'));
        assert!(is_transition(b'C', b'T'));
        assert!(is_transition(b'T', b'C'));
        assert!(!is_transition(b'A', b'C')); // transversion
        assert!(!is_transition(b'A', b'T')); // transversion
        assert!(!is_transition(b'G', b'T')); // transversion
    }

    #[test]
    fn test_weighted_reduces_to_unweighted_at_kappa_1() {
        let gc = make_gc();
        // At kappa == 1 the weighted count must equal the unweighted count for
        // EVERY codon — including stop-adjacent ones (TAC, TCA), where changes
        // to stop codons are excluded — since both use codon-level pooling.
        for cds in [
            &b"ATGGCT"[..],       // Met Ala — no stop-adjacent changes
            &b"TACTGCTCACAA"[..], // Tyr Cys Ser Gln — TAC/TCA are stop-adjacent
            &b"TTTTATCATAAT"[..], // Phe Tyr His Asn — 2-fold codons
        ] {
            let (n0, s0) = count_sites(cds, gc);
            let (n1, s1) = count_sites_weighted(cds, gc, 1.0);
            assert!((n0 - n1).abs() < 1e-9, "N: unweighted {} vs weighted@1 {} for {:?}", n0, n1, cds);
            assert!((s0 - s1).abs() < 1e-9, "S: unweighted {} vs weighted@1 {} for {:?}", s0, s1, cds);
        }
    }

    #[test]
    fn test_weighted_two_fold_site_follows_kappa_formula() {
        let gc = make_gc();
        // TTT = Phe. Only the 3rd position is (partly) synonymous: TTC (Phe) is a
        // transition (T->C); TTA/TTG (Leu) are transversions. Positions 1 and 2
        // are fully nonsynonymous with no stop changes.
        // So S = kappa/(kappa+2), N = 3 - S, for any kappa.
        for &k in &[0.5f64, 1.0, 2.0, 5.0] {
            let (n, s) = count_sites_weighted(b"TTT", gc, k);
            let expected_s = k / (k + 2.0);
            assert!((s - expected_s).abs() < 1e-9, "kappa={}: S expected {}, got {}", k, expected_s, s);
            assert!((n - (3.0 - expected_s)).abs() < 1e-9, "kappa={}: N expected {}, got {}", k, 3.0 - expected_s, n);
        }
    }

    #[test]
    fn test_weighted_four_fold_site_is_kappa_invariant() {
        let gc = make_gc();
        // GCT = Ala, GCN all Ala: the 3rd position is 4-fold degenerate, so it is
        // fully synonymous regardless of ts/tv weighting. Its synonymous-site
        // contribution must be exactly 1.0 for every kappa.
        let s_at = |k: f64| count_sites_weighted(b"GCT", gc, k).1;
        let base = s_at(1.0);
        for &k in &[0.5f64, 2.0, 10.0] {
            assert!((s_at(k) - base).abs() < 1e-9, "4-fold S must be kappa-invariant (k={})", k);
        }
        // GCT: pos1 & pos2 fully nonsyn, pos3 fully syn → S == 1.0 exactly.
        assert!((base - 1.0).abs() < 1e-9, "GCT S should be 1.0, got {}", base);
    }

    #[test]
    fn test_weighted_transition_bias_raises_s_lowers_n() {
        let gc = make_gc();
        // Across a mixed CDS, kappa>1 must raise total S and lower total N
        // relative to the equal-rates count (the whole point of the correction).
        let cds = b"TTTGCTAAACGT"; // Phe Ala Lys Arg
        let (n1, s1) = count_sites_weighted(cds, gc, 1.0);
        let (n5, s5) = count_sites_weighted(cds, gc, 5.0);
        assert!(s5 > s1, "S should rise with kappa: {} -> {}", s1, s5);
        assert!(n5 < n1, "N should fall with kappa: {} -> {}", n1, n5);
        // Sites are conserved: each position is exactly one site.
        assert!((n1 + s1 - (n5 + s5)).abs() < 1e-9, "total sites must be kappa-invariant");
    }

    #[test]
    fn test_genomic_to_cds_offset_plus_strand() {
        let gene = Gene {
            name: "test".to_string(),
            seqid: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![crate::gff::CdsExon {
                seqid: "chr1".to_string(),
                start: 100,
                end: 109,
                strand: Strand::Plus,
                phase: 0,
            }],
            length_bp: 10,
            start: 100,
        };

        assert_eq!(genomic_to_cds_offset(&gene, 100), Some(0));
        assert_eq!(genomic_to_cds_offset(&gene, 105), Some(5));
        assert_eq!(genomic_to_cds_offset(&gene, 109), Some(9));
        assert_eq!(genomic_to_cds_offset(&gene, 110), None);
        assert_eq!(genomic_to_cds_offset(&gene, 99), None);
    }

    /// Build a minimal `GenePnPs` for aggregation tests. Only the fields read
    /// by `genome_wide_pn_ps` need to be meaningful.
    fn gene_result(n_sites: f64, s_sites: f64, nonsyn: f64, syn: f64) -> GenePnPs {
        let pn = if n_sites > 0.0 { nonsyn / n_sites } else { 0.0 };
        let ps = if s_sites > 0.0 { syn / s_sites } else { 0.0 };
        let pn_ps = if ps > 0.0 {
            pn / ps
        } else if pn > 0.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };
        GenePnPs {
            name: "g".to_string(),
            length_bp: 0,
            n_sites,
            s_sites,
            pn,
            ps,
            pn_ps,
            nonsyn_snps: nonsyn,
            syn_snps: syn,
            total_snps: nonsyn + syn,
            genome_start: 0,
            genome_end: 0,
            strand: '+',
            chrom: "chr1".to_string(),
            p_value: f64::NAN,
            q_value: f64::NAN,
            p_bonferroni: f64::NAN,
            mk_dn: 0,
            mk_ds: 0,
            mk_pn: 0,
            mk_ps: 0,
        }
    }

    #[test]
    fn test_genome_wide_pools_counts_not_ratios() {
        // Gene 1: many sites, ratio 0.5. Gene 2: few sites, ratio 4.0.
        // Averaging ratios would give 2.25; pooling counts must not.
        let results = vec![
            gene_result(100.0, 100.0, 10.0, 20.0), // pN=0.10, pS=0.20
            gene_result(2.0, 2.0, 4.0, 1.0),       // pN=2.00, pS=0.50
        ];
        let gw = genome_wide_pn_ps(&results).expect("should aggregate");
        assert!((gw.n_sites - 102.0).abs() < 1e-9);
        assert!((gw.s_sites - 102.0).abs() < 1e-9);
        assert!((gw.nonsyn_snps - 14.0).abs() < 1e-9);
        assert!((gw.syn_snps - 21.0).abs() < 1e-9);
        // pooled pN = 14/102, pS = 21/102, ratio = 14/21
        assert!((gw.pn - 14.0 / 102.0).abs() < 1e-9);
        assert!((gw.ps - 21.0 / 102.0).abs() < 1e-9);
        assert!((gw.pn_ps - 14.0 / 21.0).abs() < 1e-9, "got {}", gw.pn_ps);
    }

    #[test]
    fn test_genome_wide_empty_is_none() {
        assert!(genome_wide_pn_ps(&[]).is_none());
    }

    #[test]
    fn test_genome_wide_no_variation_is_nan() {
        // Genes exist but carry no SNPs → ratio is NaN, not 0 or Inf.
        let results = vec![gene_result(50.0, 25.0, 0.0, 0.0)];
        let gw = genome_wide_pn_ps(&results).unwrap();
        assert_eq!(gw.pn, 0.0);
        assert_eq!(gw.ps, 0.0);
        assert!(gw.pn_ps.is_nan(), "expected NaN, got {}", gw.pn_ps);
    }

    #[test]
    fn test_genome_wide_no_syn_is_infinite() {
        // Nonsynonymous variation but zero synonymous → +Infinity.
        let results = vec![gene_result(50.0, 25.0, 3.0, 0.0)];
        let gw = genome_wide_pn_ps(&results).unwrap();
        assert!(gw.pn_ps.is_infinite() && gw.pn_ps > 0.0);
    }

    #[test]
    fn test_selection_label_bands() {
        assert!(selection_label(0.3).contains("purifying"));
        assert!(selection_label(1.0).contains("neutral"));
        assert!(selection_label(2.0).contains("diversifying"));
        assert!(selection_label(f64::NAN).contains("undetermined"));
        assert!(selection_label(f64::INFINITY).contains("synonymous"));
    }

    #[test]
    fn test_format_ratio_special_values() {
        assert_eq!(format_ratio(f64::NAN), "NaN");
        assert_eq!(format_ratio(f64::INFINITY), "inf");
        assert_eq!(format_ratio(0.5), "0.500000");
    }

    #[test]
    fn test_genomic_to_cds_offset_minus_strand() {
        let gene = Gene {
            name: "test".to_string(),
            seqid: "chr1".to_string(),
            strand: Strand::Minus,
            exons: vec![crate::gff::CdsExon {
                seqid: "chr1".to_string(),
                start: 100,
                end: 109,
                strand: Strand::Minus,
                phase: 0,
            }],
            length_bp: 10,
            start: 100,
        };

        // Minus strand: position 109 maps to CDS offset 0
        assert_eq!(genomic_to_cds_offset(&gene, 109), Some(0));
        assert_eq!(genomic_to_cds_offset(&gene, 100), Some(9));
    }
}
