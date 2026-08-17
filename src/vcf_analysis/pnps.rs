//! Per-gene pN/pS computation and genome-wide pooled estimates.

use super::*;
use indicatif::{ProgressBar, ProgressStyle};

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
) -> (Vec<GenePnPs>, ComputeDiagnostics) {
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
    // Genes whose reference CDS translates to an internal stop codon (a strong signal
    // of a wrong --genetic-code, wrong frame/phase, or a wrong reference build).
    let internal_stops = AtomicUsize::new(0);
    // Codons carrying two or more retained SNPs. Every SNP is scored against the
    // REFERENCE codon (see the alternate-codon construction below), so such codons are
    // classified one SNP at a time and their joint amino-acid change is never evaluated.
    let multi_snp_codons = AtomicUsize::new(0);
    // The subset of those codons whose allele frequencies force the SNPs onto the same
    // haplotype (see `forced_cooccurrence`), i.e. where the independent scoring is
    // provably answering the wrong question.
    let cooccurring_codons = AtomicUsize::new(0);
    let genes_with_cooccurring = AtomicUsize::new(0);

    // Per-gene progress bar for liveness on genome-scale runs (auto-hidden when stderr
    // is not a terminal, so tests and pipelines stay quiet).
    let pb = ProgressBar::new(genes.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Genes: {pos}/{len} ({eta})")
            .expect("progress template is a literal, so it always parses")
            .progress_chars("#>-"),
    );

    let results: Vec<GenePnPs> = genes
        .par_iter()
        .filter_map(|gene| {
            pb.inc(1);
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

            // Diagnostic: a stop codon anywhere but the last codon of the reference CDS
            // means the gene does not translate cleanly — usually the wrong genetic
            // code, reading frame/phase, or reference build.
            let ncodons = cds_seq.len() / 3;
            if (0..ncodons.saturating_sub(1)).any(|c| {
                codon_to_aa(&[cds_seq[c * 3], cds_seq[c * 3 + 1], cds_seq[c * 3 + 2]], gc) == Some(b'*')
            }) {
                internal_stops.fetch_add(1, Ordering::Relaxed);
            }

            // Count S and N sites from reference codons
            let (n_sites, s_sites) = if spectrum_weighted {
                count_sites_weighted(&cds_seq, gc, kappa)
            } else {
                count_sites(&cds_seq, gc)
            };

            // A CDS with neither synonymous nor nonsynonymous sites has no translatable
            // codons at all (every codon is ambiguous/all-N in the reference). pN and pS
            // would both be 0/0, so analysing it yields only spurious zeros that read as a
            // real gene under total constraint. Skip it loudly, like the other
            // untranslatable-CDS guards, rather than emit a misleading all-zero row.
            if n_sites == 0.0 && s_sites == 0.0 {
                warn!(
                    "Gene {} CDS has no usable (translatable) codons — an all-ambiguous / all-N \
                     reference region, skipping",
                    gene.name
                );
                return None;
            }

            // The gene's per-codon mutational opportunity, for the optional recurrence
            // scan's family size and plug-in rate. Counted here, where the CDS is
            // already reconstructed, and never rate-weighted: kappa moves the pN/pS
            // site denominators, not the number of changes a codon can undergo.
            let (scan_codons, scan_poss_nonsyn) = codon_space(&cds_seq, gc);

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
            // Per-coding-SNP records for the optional variants table (all effects,
            // including nonsense/stop-loss that the pN/pS counts below exclude).
            let mut variants: Vec<Variant> = Vec::new();
            // (codon index, genomic position, highest ALT frequency) for every SNP that
            // received a coding classification, so the same-codon (MNV) check below can
            // group them. Recorded per SNP, not per ALT: alternative alleles of one site
            // are mutually exclusive, never a haplotype.
            let mut codon_hits: Vec<(usize, usize, f64)> = Vec::new();

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

                        // Highest ALT frequency among this position's classified alleles,
                        // for the same-codon check. The maximum (not the sum) is the right
                        // summary: the alleles exclude one another, so co-occurrence with a
                        // neighbouring SNP is driven by whichever ALT is commonest here.
                        let mut coding_af: Option<f64> = None;

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

                            // Look up amino acids; an ambiguous codon (None) is skipped.
                            let (Some(r), Some(a)) = (codon_to_aa(&ref_codon, gc), codon_to_aa(&alt_codon, gc))
                            else {
                                continue;
                            };
                            let af = snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0);

                            // Record the variant for the optional table, classifying every
                            // coding effect — including nonsense/stop-loss.
                            let effect = if r == a {
                                SnpEffect::Synonymous
                            } else if a == b'*' {
                                SnpEffect::Nonsense
                            } else if r == b'*' {
                                SnpEffect::StopLoss
                            } else {
                                SnpEffect::Missense
                            };
                            variants.push(Variant {
                                pos: snp.pos,
                                ref_allele: snp.ref_allele,
                                alt_allele: *alt_base,
                                aa_pos: codon_idx + 1,
                                ref_codon,
                                ref_aa: r,
                                alt_aa: a,
                                af,
                                gt_derived: snp
                                    .gt_counts
                                    .as_ref()
                                    .and_then(|gc| gc.alt.get(alt_idx).copied()),
                                gt_called: snp.gt_counts.as_ref().map(|gc| gc.called),
                                effect,
                            });
                            coding_af = Some(coding_af.map_or(af, |best: f64| best.max(af)));

                            // pN/pS + MK + SFS counting: synonymous and missense only
                            // (changes to/from a stop are excluded, as before). Uses this
                            // ALT's allele frequency as raw counts (not AF-weighted).
                            if r != b'*' && a != b'*' {
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
                        }

                        if let Some(af) = coding_af {
                            codon_hits.push((codon_idx, snp.pos, af));
                        }
                    }
                }
            }

            // Same-codon (MNV) check. eskaks builds every alternate codon by copying the
            // REFERENCE codon and mutating one position, so two SNPs in one codon are each
            // scored against the reference and never against each other. Group the retained
            // SNPs by codon and record how many codons that affects, so the run can say so
            // instead of silently reporting a change no haplotype carries.
            let (mut gene_multi, mut gene_cooccurring) = (0usize, 0usize);
            if codon_hits.len() > 1 {
                // Sort by codon, then position, then descending AF, so collapsing repeated
                // positions (a VCF may split one site across records) keeps the highest AF.
                codon_hits.sort_by(|a, b| {
                    a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(b.2.total_cmp(&a.2))
                });
                codon_hits.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);
                for group in codon_hits.chunk_by(|a, b| a.0 == b.0) {
                    if group.len() < 2 {
                        continue;
                    }
                    gene_multi += 1;
                    let afs: Vec<f64> = group.iter().map(|h| h.2).collect();
                    let forced = forced_cooccurrence(&afs);
                    if forced {
                        gene_cooccurring += 1;
                    }
                    let detail: Vec<String> = group
                        .iter()
                        .map(|h| format!("{} (AF {:.3})", h.1, h.2))
                        .collect();
                    debug!(
                        "Gene {} residue {}: {} SNPs share one codon at {}:{}. {}",
                        gene.name,
                        group[0].0 + 1,
                        group.len(),
                        gene.seqid,
                        detail.join(", "),
                        if forced {
                            "Their allele frequencies put them on the same haplotype, so the joint \
                             codon change is the real one and the per-SNP calls reported for this \
                             residue are not it."
                        } else {
                            "Phase is unknown at these frequencies, so each SNP was scored against \
                             the reference codon on its own."
                        }
                    );
                }
            }
            multi_snp_codons.fetch_add(gene_multi, Ordering::Relaxed);
            cooccurring_codons.fetch_add(gene_cooccurring, Ordering::Relaxed);
            if gene_cooccurring > 0 {
                genes_with_cooccurring.fetch_add(1, Ordering::Relaxed);
            }

            ref_checked.fetch_add(local_checked, Ordering::Relaxed);
            ref_mismatch.fetch_add(local_mismatch, Ordering::Relaxed);

            let total_snps = nonsyn_count + syn_count;
            // A zero site denominator makes the density undefined (0/0), not 0.0: a gene of
            // only Met/Trp codons has s_sites == 0, so pS cannot be estimated. Report NaN so
            // it is not read as an observed synonymous rate of exactly zero. The pn_ps guard
            // below already treats a NaN/absent pS the same as the old 0.0 (→ +inf when
            // pn > 0, NaN when both are 0), so the headline ratio is unchanged.
            let pn = if n_sites > 0.0 {
                nonsyn_count / n_sites
            } else {
                f64::NAN
            };
            let ps = if s_sites > 0.0 {
                syn_count / s_sites
            } else {
                f64::NAN
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
                variants,
                scan_codons,
                scan_poss_nonsyn,
            })
        })
        .collect();
    pb.finish_and_clear();

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

    let multi_snp_codons = multi_snp_codons.load(Ordering::Relaxed);
    let cooccurring_codons = cooccurring_codons.load(Ordering::Relaxed);
    let genes_with_cooccurring = genes_with_cooccurring.load(Ordering::Relaxed);

    // Warn only where co-occurrence is a fact of the allele frequencies. In a
    // multi-sample population VCF two SNPs in one codon usually sit on different
    // genomic backgrounds, and scoring them separately is then the correct thing to
    // do; a warning on every such codon would cry wolf on every population run.
    if cooccurring_codons > 0 {
        warn!(
            "{} codon(s) in {} gene(s) carry two or more SNPs that the allele frequencies place on \
             the same haplotype in at least some samples. Each SNP is classified against the \
             REFERENCE codon on its own, so the joint (multi-nucleotide) codon change is never \
             evaluated and the synonymous / nonsynonymous calls for those codons, and the pN/pS \
             they feed, can be wrong. Haplotype-aware calling needs phase and is out of scope; run \
             with -vv to list the codons.",
            cooccurring_codons, genes_with_cooccurring
        );
    }
    // The rest are reported, not warned about: their frequencies leave phase open, so
    // the independent per-SNP classification may well be the right answer.
    let phase_unknown = multi_snp_codons - cooccurring_codons;
    if phase_unknown > 0 {
        info!(
            "{} further codon(s) carry more than one SNP, but their allele frequencies do not force \
             the SNPs onto one haplotype (they may sit on different genomic backgrounds), so each \
             was classified against the reference codon (run with -vv to list them).",
            phase_unknown
        );
    }

    info!("Computed pN/pS for {} genes", results.len());
    let diagnostics = ComputeDiagnostics {
        snps_in_cds: ref_checked,
        ref_mismatch,
        genes_with_internal_stops: internal_stops.load(Ordering::Relaxed),
        multi_snp_codons,
        cooccurring_codons,
        genes_with_cooccurring,
    };
    (results, diagnostics)
}

/// Do the allele frequencies of the SNPs sharing one codon force them onto the same
/// haplotype?
///
/// With `k` variant positions at frequencies `p_1..p_k`, the fraction of sampled
/// genomes carrying all of them at once is at least `Σ p_i − (k − 1)` (each variant is
/// absent from `1 − p_i` of the genomes; the absences can cover everything only while
/// their total is at least 1). When that bound is positive the co-occurrence is a
/// property of the data rather than a guess, so scoring each SNP against the reference
/// codon alone is provably answering the wrong question: at AF 1.0 and 1.0 every sample
/// carries both, and a single-sample VCF (no AF, so 1.0 by convention) is the same case.
///
/// Below the bound the variants may sit on different genomic backgrounds, where the
/// independent classification is correct, so this returns false and nothing is warned.
///
/// The per-codon recurrence scan reuses this over the same codon's distinct SNP
/// positions: when the frequencies force the changes onto one haplotype they are one
/// multi-nucleotide event, not several independent origins, so its independence
/// assumption is provably false and the codon is suppressed rather than tested.
pub(crate) fn forced_cooccurrence(afs: &[f64]) -> bool {
    if afs.len() < 2 {
        return false;
    }
    let sum: f64 = afs.iter().sum();
    sum > afs.len() as f64 - 1.0
}

/// Whole-run diagnostics from [`compute_pn_ps`], surfaced in the CLI summary so an
/// empty or garbage result is never mistaken for a clean run.
#[derive(Debug, Clone, Copy, Default)]
pub struct ComputeDiagnostics {
    /// SNPs whose position fell inside a CDS (the REF-allele check ran on these).
    pub snps_in_cds: usize,
    /// Of the in-CDS SNPs, how many had a REF allele disagreeing with the reference.
    pub ref_mismatch: usize,
    /// Genes whose reference CDS contains an internal stop codon.
    pub genes_with_internal_stops: usize,
    /// Codons carrying two or more retained SNPs at distinct positions. Each SNP is
    /// classified against the reference codon on its own, so the joint codon change is
    /// never evaluated for these.
    pub multi_snp_codons: usize,
    /// The subset of `multi_snp_codons` whose allele frequencies force the SNPs onto the
    /// same haplotype (see `forced_cooccurrence`), where that independent classification
    /// is provably reporting a change no genome carries.
    pub cooccurring_codons: usize,
    /// Genes contributing at least one `cooccurring_codons` codon.
    pub genes_with_cooccurring: usize,
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
        let pn = if n_sites > 0.0 { nonsyn / n_sites } else { f64::NAN };
        let ps = if s_sites > 0.0 { syn / s_sites } else { f64::NAN };
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
    // Pooled density is undefined (NaN), not 0.0, when its site denominator is 0 — the
    // whole retained set has no synonymous (or nonsynonymous) sites. The pn_ps guard below
    // keeps the same +inf / NaN outcome, so only the surfaced pN / pS components change.
    let pn = if n_sites > 0.0 { nonsyn_snps / n_sites } else { f64::NAN };
    let ps = if s_sites > 0.0 { syn_snps / s_sites } else { f64::NAN };
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


#[cfg(test)]
mod same_codon_tests {
    use super::*;
    use crate::gff::CdsExon;
    use std::sync::{Mutex, Once};

    /// Log sink for the warning assertions. The `log` facade allows exactly one logger
    /// per process, so it is installed once and every test reads the shared buffer,
    /// taking a snapshot index first so a concurrent test's lines are never read as ours.
    struct Capture;
    static CAPTURE: Capture = Capture;
    static LINES: Mutex<Vec<String>> = Mutex::new(Vec::new());
    /// Serialises the tests that read the buffer: the sink is process-wide, so two of
    /// them running at once would each see the other's lines. Poisoning is ignored on
    /// purpose, so one failing test does not cascade into the rest.
    static LOG_TESTS: Mutex<()> = Mutex::new(());

    impl log::Log for Capture {
        fn enabled(&self, _: &log::Metadata) -> bool {
            true
        }
        fn log(&self, record: &log::Record) {
            LINES
                .lock()
                .expect("log buffer mutex")
                .push(format!("{}|{}", record.level(), record.args()));
        }
        fn flush(&self) {}
    }

    /// Install the capture logger and return the current end of the buffer.
    fn start_capture() -> usize {
        static ONCE: Once = Once::new();
        ONCE.call_once(|| {
            log::set_logger(&CAPTURE).expect("no other logger in the test binary");
            log::set_max_level(log::LevelFilter::Debug);
        });
        LINES.lock().expect("log buffer mutex").len()
    }

    /// Lines logged at `level` since `from` that contain `needle`.
    fn captured(from: usize, level: &str, needle: &str) -> Vec<String> {
        LINES
            .lock()
            .expect("log buffer mutex")
            .iter()
            .skip(from)
            .filter(|l| l.starts_with(level) && l.contains(needle))
            .cloned()
            .collect()
    }

    /// The reference of the bug report: ATG CTT TAA, one plus-strand gene on chr1.
    fn leu_gene(name: &str) -> (HashMap<String, Vec<u8>>, Gene) {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ATGCTTTAA".to_vec());
        let gene = Gene {
            name: name.to_string(),
            seqid: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![CdsExon {
                seqid: "chr1".to_string(),
                start: 1,
                end: 9,
                strand: Strand::Plus,
                phase: 0,
            }],
            length_bp: 9,
            start: 1,
        };
        (reference, gene)
    }

    /// Single-ALT SNP builder.
    fn snp(pos: usize, ref_allele: u8, alt: u8, af: f64) -> VcfSnp {
        VcfSnp {
            chrom: "chr1".to_string(),
            pos,
            ref_allele,
            alt_alleles: vec![alt],
            alt_freqs: vec![af],
            gt_counts: None,
            filter: "PASS".to_string(),
            depth: None,
        }
    }

    #[test]
    fn fixed_snps_in_one_codon_warn_and_are_still_scored_apart() {
        // Reference CDS ATG CTT TAA with chr1:4 C>T and chr1:6 T>A, both at AF 1.0.
        // Scored one at a time against the reference codon CTT: C>T gives TTT (Phe,
        // missense L2F) and T>A gives CTA (Leu, synonymous). The haplotype every sample
        // actually carries is TTA, which is Leu: one synonymous multi-nucleotide change.
        // Calling it needs phase and is out of scope, so the run must SAY so.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_fixed");
        let snps = vec![snp(4, b'C', b'T', 1.0), snp(6, b'T', b'A', 1.0)];

        let _serialised = LOG_TESTS.lock().unwrap_or_else(|e| e.into_inner());
        let from = start_capture();
        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        // The classification itself is unchanged (this is a detection change, not a fix).
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].nonsyn_snps, 1.0, "the invented L2F still counts as one nonsyn");
        assert_eq!(res[0].syn_snps, 1.0);

        // Detection: one codon, forced onto one haplotype, in one gene.
        assert_eq!(diag.multi_snp_codons, 1, "codon 2 carries both SNPs");
        assert_eq!(diag.cooccurring_codons, 1, "AF 1.0 + 1.0 forces co-occurrence");
        assert_eq!(diag.genes_with_cooccurring, 1);

        // The warning fires, at warn level, naming the count.
        let warns = captured(from, "WARN", "same haplotype");
        assert_eq!(
            warns.len(),
            1,
            "expected exactly one same-codon warning, got {:?}",
            warns
        );
        assert!(warns[0].contains("1 codon(s) in 1 gene(s)"), "warning text: {}", warns[0]);

        // -vv detail names the residue and both positions, so the user can act on it.
        let details = captured(from, "DEBUG", "share one codon");
        assert_eq!(details.len(), 1, "expected one debug detail line, got {:?}", details);
        assert!(
            details[0].contains("residue 2")
                && details[0].contains("chr1:4 (AF 1.000)")
                && details[0].contains("6 (AF 1.000)"),
            "debug detail: {}",
            details[0]
        );
    }

    #[test]
    fn rare_snps_in_one_codon_are_counted_but_not_warned_about() {
        // Same codon, but at AF 0.3 and 0.2 the two SNPs can sit on entirely different
        // genomic backgrounds, where scoring them apart is the CORRECT thing to do.
        // Warning here would cry wolf on every multi-sample population VCF.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_rare");
        let snps = vec![snp(4, b'C', b'T', 0.3), snp(6, b'T', b'A', 0.2)];

        let _serialised = LOG_TESTS.lock().unwrap_or_else(|e| e.into_inner());
        let from = start_capture();
        let (_res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        assert_eq!(diag.multi_snp_codons, 1, "the codon is still counted");
        assert_eq!(diag.cooccurring_codons, 0, "0.3 + 0.2 leaves phase open");
        assert_eq!(diag.genes_with_cooccurring, 0);
        assert!(
            captured(from, "WARN", "same haplotype").is_empty(),
            "an unphased pair must not warn"
        );
    }

    #[test]
    fn snps_in_neighbouring_codons_are_not_flagged() {
        // chr1:3 (codon 1) and chr1:4 (codon 2) are adjacent bases in DIFFERENT codons,
        // so the reference-codon scoring is exact and nothing should be reported.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_split");
        let snps = vec![snp(3, b'G', b'A', 1.0), snp(4, b'C', b'T', 1.0)];
        let (_res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);
        assert_eq!(diag.multi_snp_codons, 0);
        assert_eq!(diag.cooccurring_codons, 0);
    }

    #[test]
    fn two_alts_at_one_position_are_not_a_haplotype() {
        // A multi-allelic site is ONE position: its ALTs exclude one another, so they can
        // never form a multi-nucleotide codon change however high their frequencies are.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_multiallelic");
        let multi = VcfSnp {
            chrom: "chr1".to_string(),
            pos: 4,
            ref_allele: b'C',
            alt_alleles: vec![b'T', b'A'],
            alt_freqs: vec![0.6, 0.4],
            gt_counts: None,
            filter: "PASS".to_string(),
            depth: None,
        };
        let (_res, diag) = compute_pn_ps(&reference, &[gene], &[multi], gc, false, 1.0, 0.95);
        assert_eq!(diag.multi_snp_codons, 0, "one position is not two positions");
        assert_eq!(diag.cooccurring_codons, 0);
    }

    #[test]
    fn forced_cooccurrence_follows_the_pigeonhole_bound() {
        // Sum of frequencies above k-1 guarantees an overlap; at or below it, none.
        assert!(forced_cooccurrence(&[1.0, 1.0]));
        assert!(forced_cooccurrence(&[0.6, 0.6]));
        assert!(!forced_cooccurrence(&[0.5, 0.5]), "sum exactly 1 guarantees nothing");
        assert!(!forced_cooccurrence(&[0.9, 0.05]));
        assert!(forced_cooccurrence(&[0.9, 0.9, 0.9]), "3 SNPs: 2.7 > 2");
        assert!(!forced_cooccurrence(&[0.9, 0.9, 0.1]), "3 SNPs: 1.9 <= 2");
        // Fewer than two positions is never a multi-nucleotide change.
        assert!(!forced_cooccurrence(&[1.0]));
        assert!(!forced_cooccurrence(&[]));
    }
}
