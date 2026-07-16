//! Per-gene pN/pS computation and genome-wide pooled estimates.

use super::*;

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

            // A CDS exon extending past the reference end would desync
            // genomic_to_cds_offset (which counts full exon spans) from the
            // reference-clamped CDS, silently mis-mapping or dropping SNPs. Skip such a
            // gene loudly rather than analyse a truncated / mismatched annotation.
            if gene.exons.iter().any(|e| e.end > ref_seq.len()) {
                warn!(
                    "Gene {} has a CDS exon extending past the end of {} ({} bp), skipping",
                    gene.name, gene.seqid, ref_seq.len()
                );
                return None;
            }

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
            // Raw (unweighted) count of classified SNP alleles, for --min-snps and the
            // "genes with SNPs" tally — independent of AF weighting.
            let mut n_snps_raw = 0usize;
            let mut local_checked = 0usize;
            let mut local_mismatch = 0usize;
            // McDonald-Kreitman: fixed (AF >= threshold) vs polymorphic counts.
            let (mut mk_dn, mut mk_ds, mut mk_pn, mut mk_ps) = (0u32, 0u32, 0u32, 0u32);
            // Site frequency spectrum: nonsyn/syn SNP counts per allele-frequency bin.
            let mut sfs_nonsyn = [0u32; SFS_NBINS];
            let mut sfs_syn = [0u32; SFS_NBINS];

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
                                    let bin = sfs_bin(af);
                                    n_snps_raw += 1;
                                    if r == a {
                                        syn_count += weight;
                                        sfs_syn[bin] += 1;
                                        if fixed { mk_ds += 1 } else { mk_ps += 1 }
                                    } else {
                                        nonsyn_count += weight;
                                        sfs_nonsyn[bin] += 1;
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
            let (p_value, neglog10p) = if !af_weighted && total_snps > 0.0 && sites > 0.0 {
                let k = nonsyn_count.round() as u64;
                let n = total_snps.round() as u64;
                let p0 = n_sites / sites;
                (
                    crate::stats::binomial_two_sided_p(k, n, p0),
                    crate::stats::binomial_two_sided_neglog10p(k, n, p0),
                )
            } else {
                (f64::NAN, f64::NAN)
            };

            // 95% Wilson CI on pN/pS: bound the nonsynonymous fraction of SNPs,
            // then map q → (q/(1−q))·(S_sites/N_sites). Undefined under AF weighting.
            let (pn_ps_lo, pn_ps_hi) = if !af_weighted && total_snps > 0.0 && n_sites > 0.0 && s_sites > 0.0 {
                let k = nonsyn_count.round() as u64;
                let n = total_snps.round() as u64;
                let (qlo, qhi) = crate::stats::wilson_interval(k, n, 0.95);
                let scale = s_sites / n_sites;
                let map = |q: f64| if q >= 1.0 { f64::INFINITY } else { (q / (1.0 - q)) * scale };
                (map(qlo), map(qhi))
            } else {
                (f64::NAN, f64::NAN)
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
                n_snps: n_snps_raw,
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
                neglog10p,
                pn_ps_lo,
                pn_ps_hi,
                p_gc: f64::NAN,
                q_gc: f64::NAN,
                repetitive: is_repetitive(&gene.name),
                sfs_nonsyn,
                sfs_syn,
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
        // ps == 0 with pn > 0 is a genuine upper-tail extreme (pN/pS = +∞), so keep it
        // toward the upper percentile instead of discarding it and biasing the CI low.
        // Only 0/0 (no variation at all) is genuinely undefined and excluded.
        let ratio = if ps > 0.0 {
            pn / ps
        } else if pn > 0.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };
        if ratio.is_nan() {
            undefined += 1;
        } else {
            ratios.push(ratio);
        }
    }
    // Warn BEFORE the empty-set early return so an all-undefined result is never silent.
    if undefined > 0 {
        warn!(
            "{}/{} bootstrap replicates had no variation at all (0/0) and were excluded \
             from the genome-wide pN/pS CI.",
            undefined, n_boot
        );
    }
    if ratios.is_empty() {
        return None;
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tail = (1.0 - confidence) / 2.0 * 100.0;
    Some((
        crate::stats::percentile_sorted(&ratios, tail),
        crate::stats::percentile_sorted(&ratios, 100.0 - tail),
    ))
}

/// Aggregate per-gene results into a single genome-wide (pooled) pN/pS estimate.
///
/// Returns `None` when there are no genes to aggregate. The pooled ratio uses
/// the same NaN/Infinity conventions as the per-gene computation: `NaN` when
/// there is no variation at all, `+Infinity` when there are nonsynonymous but
/// no synonymous changes.
pub fn genome_wide_pn_ps(results: &[GenePnPs]) -> Option<GenomeWidePnPs> {
    pool_pn_ps(results.iter())
}

/// Pool a (possibly filtered) set of genes into one pN/pS estimate. Shared by
/// the whole-genome and the core-vs-repetitive stratified estimates.
fn pool_pn_ps<'a>(it: impl Iterator<Item = &'a GenePnPs>) -> Option<GenomeWidePnPs> {
    let (mut n_sites, mut s_sites, mut nonsyn_snps, mut syn_snps, mut any) =
        (0.0, 0.0, 0.0, 0.0, false);
    for r in it {
        any = true;
        n_sites += r.n_sites;
        s_sites += r.s_sites;
        nonsyn_snps += r.nonsyn_snps;
        syn_snps += r.syn_snps;
    }
    if !any {
        return None;
    }
    let pn = if n_sites > 0.0 { nonsyn_snps / n_sites } else { 0.0 };
    let ps = if s_sites > 0.0 { syn_snps / s_sites } else { 0.0 };
    let pn_ps = if ps > 0.0 {
        pn / ps
    } else if pn > 0.0 {
        f64::INFINITY
    } else {
        f64::NAN
    };
    Some(GenomeWidePnPs { n_sites, s_sites, nonsyn_snps, syn_snps, pn, ps, pn_ps })
}

/// Pooled pN/pS split into (core, repetitive) — repetitive = PE/PPE/PGRS/IS/etc.
/// A big gap between the two flags mapping-artefact inflation in the repeats.
pub fn genome_wide_core_repetitive(
    results: &[GenePnPs],
) -> (Option<GenomeWidePnPs>, Option<GenomeWidePnPs>) {
    (
        pool_pn_ps(results.iter().filter(|r| !r.repetitive)),
        pool_pn_ps(results.iter().filter(|r| r.repetitive)),
    )
}

