//! Per-gene pN/pS computation and genome-wide pooled estimates.

use super::mnv::{call_codon, AlleleCall, CodonAllele};
use super::*;
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashSet;

/// One coding ALT allele, collected on a first pass over the gene's SNPs and classified
/// on a second, once every SNP sharing its codon is known.
///
/// The two passes are the whole point: an ALT cannot be classified until the codon its
/// carriers actually hold is known, and that needs its neighbours. They are kept in the
/// order their records were read so the variants table's row order is exactly what a
/// single pass produced.
struct PendingAllele<'a> {
    /// 0-based codon (residue) index in the CDS, and which of its bases this ALT changes.
    codon_idx: usize,
    pos_in_codon: usize,
    /// 1-based genomic position, and the REF/ALT bases on the VCF (genomic) strand.
    pos: usize,
    ref_allele: u8,
    alt_allele: u8,
    /// The alternate base on the CODING strand (complemented for a minus-strand gene).
    alt_in_cds: u8,
    af: f64,
    /// This allele's contribution to the counts: its AF under `--af-weighted`, else 1.
    weight: f64,
    gt_derived: Option<usize>,
    gt_called: Option<usize>,
    carriers: Option<&'a CarrierSet>,
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
    // Codons carrying two or more retained SNPs at distinct positions.
    let multi_snp_codons = AtomicUsize::new(0);
    // The subset of those codons where the SNPs provably occur together in one genome.
    // With genotypes those codons are scored jointly; without them, they are the codons
    // whose joint change eskaks cannot evaluate.
    let cooccurring_codons = AtomicUsize::new(0);
    let genes_with_cooccurring = AtomicUsize::new(0);
    // Of `cooccurring_codons`, those established by intersecting the per-sample carrier
    // sets (a fact about a named sample) rather than by the allele-frequency bound.
    let cooccurring_exact = AtomicUsize::new(0);
    // Multi-SNP codons whose every position carried a usable carrier set, so BOTH the
    // positive and the negative answer are exact rather than "phase unknown".
    let codons_with_carriers = AtomicUsize::new(0);
    // ALT alleles scored as part of a multi-nucleotide change (their carriers hold a
    // codon differing from the reference at two or three bases), how many of those left
    // the McDonald-Kreitman 2x2 and the SFS, and the codons whose joint change is a
    // premature stop that no single SNP of theirs creates.
    let mnv_alleles = AtomicUsize::new(0);
    let mnv_excluded_from_mk = AtomicUsize::new(0);
    let mnv_new_stops = AtomicUsize::new(0);

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
            // Counts are f64 to support AF-weighted mode (πN/πS), and because
            // Nei-Gojobori pathway averaging splits a multi-nucleotide change between the
            // two classes: each classified allele still contributes exactly 1.0 in total.
            let mut nonsyn_count = 0.0f64;
            let mut syn_count = 0.0f64;
            // Raw (unweighted) count of classified SNP alleles, for --min-snps and the
            // "genes with SNPs" tally — independent of AF weighting.
            let mut n_snps_raw = 0usize;
            let mut local_checked = 0usize;
            let mut local_mismatch = 0usize;
            // McDonald-Kreitman: fixed (AF >= threshold) vs polymorphic counts, restricted
            // to alleles that stand alone in their codon (Fisher's exact needs whole
            // alleles, and a multi-nucleotide change splits between the two classes).
            let (mut mk_dn, mut mk_ds, mut mk_pn, mut mk_ps) = (0u32, 0u32, 0u32, 0u32);
            let mut mk_mnv_excluded = 0u32;
            // Alleles scored inside a multi-nucleotide change, and the codons whose joint
            // change is a premature stop that no single SNP of theirs creates.
            let mut gene_mnv_alleles = 0usize;
            let mut gene_new_stops: HashSet<usize> = HashSet::new();
            // Site frequency spectrum: nonsyn/syn SNP counts per allele-frequency bin.
            let mut sfs_nonsyn = [0u32; SFS_NBINS];
            let mut sfs_syn = [0u32; SFS_NBINS];
            // Per-coding-SNP records for the optional variants table (all effects,
            // including nonsense/stop-loss that the pN/pS counts below exclude).
            let mut variants: Vec<Variant> = Vec::new();
            // (codon index, genomic position, highest ALT frequency, carriers) for every
            // SNP that received a coding classification, so the same-codon (MNV) check
            // below can group them. Recorded per SNP, not per ALT: alternative alleles of
            // one site are mutually exclusive, never a haplotype. `carriers` is therefore
            // the UNION over this position's classified ALTs, i.e. the samples carrying
            // any non-reference base here, and `None` when the input has no genotypes.
            let mut codon_hits: Vec<(usize, usize, f64, Option<CarrierSet>)> = Vec::new();
            // Carriers of each pushed `variants` row, parallel to it, so the shared-codon
            // mark can be decided per ALLELE and not merely per codon. Borrowed from the
            // SNP records, never cloned, so this costs one pointer per row.
            let mut variant_carriers: Vec<Option<&CarrierSet>> = Vec::new();
            // Pass 1 output: every coding ALT of this gene, in record order, plus the
            // reference codon of each residue they landed on. An ALT cannot be classified
            // until its codon's other SNPs are known, so nothing is scored here.
            let mut pending: Vec<PendingAllele> = Vec::new();
            let mut ref_codons: HashMap<usize, [u8; 3]> = HashMap::new();

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
                        // Samples carrying ANY classified ALT at this position, and whether
                        // every classified ALT contributed one. A haploid sample carries at
                        // most one ALT here, so the union is exactly "has a variant base at
                        // this position", which is what a neighbouring SNP must intersect.
                        let mut coding_carriers: Option<CarrierSet> = None;
                        let mut carriers_complete = true;

                        // An untranslatable reference codon classifies nothing, exactly as
                        // when every ALT of it was skipped one at a time.
                        if codon_to_aa(&ref_codon, gc).is_none() {
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

                            let alt_in_cds = if gene.strand == Strand::Minus {
                                complement(*alt_base)
                            } else {
                                *alt_base
                            };
                            // An ALT that cannot be translated in its own codon is skipped,
                            // as it always was, and since every retained base is then a
                            // plain A/C/G/T, no codon built from them can be ambiguous.
                            let mut solo_codon = ref_codon;
                            solo_codon[pos_in_codon] = alt_in_cds;
                            if codon_to_aa(&solo_codon, gc).is_none() {
                                continue;
                            }
                            let af = snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0);

                            let alt_carriers = snp.carriers.as_ref().and_then(|c| c.get(alt_idx));
                            pending.push(PendingAllele {
                                codon_idx,
                                pos_in_codon,
                                pos: snp.pos,
                                ref_allele: snp.ref_allele,
                                alt_allele: *alt_base,
                                alt_in_cds,
                                af,
                                weight,
                                gt_derived: snp
                                    .gt_counts
                                    .as_ref()
                                    .and_then(|gc| gc.alt.get(alt_idx).copied()),
                                gt_called: snp.gt_counts.as_ref().map(|gc| gc.called),
                                carriers: alt_carriers,
                            });
                            ref_codons.insert(codon_idx, ref_codon);
                            match alt_carriers {
                                Some(cs) => match coding_carriers.as_mut() {
                                    Some(acc) => acc.union_with(cs),
                                    None => coding_carriers = Some(cs.clone()),
                                },
                                None => carriers_complete = false,
                            }
                            coding_af = Some(coding_af.map_or(af, |best: f64| best.max(af)));
                        }

                        if let Some(af) = coding_af {
                            let carriers = if carriers_complete { coding_carriers } else { None };
                            codon_hits.push((codon_idx, snp.pos, af, carriers));
                        }
                    }
                }
            }

            // Pass 2: classify codon by codon, scoring each ALT against the codon its
            // carriers actually hold rather than against the reference codon. Codons whose
            // SNPs sit at a single base (every codon of an ordinary run) and codons from
            // an input without genotypes take the historic per-SNP path unchanged.
            let mut by_codon: HashMap<usize, Vec<usize>> = HashMap::new();
            for (i, p) in pending.iter().enumerate() {
                by_codon.entry(p.codon_idx).or_default().push(i);
            }
            let mut calls: Vec<Option<AlleleCall>> = vec![None; pending.len()];
            for (codon_idx, slots) in &by_codon {
                let ref_codon = ref_codons[codon_idx];
                let ref_aa = codon_to_aa(&ref_codon, gc)
                    .expect("pass 1 only records alleles on translatable reference codons");
                let group: Vec<CodonAllele> = slots
                    .iter()
                    .map(|&i| CodonAllele {
                        pos_in_codon: pending[i].pos_in_codon,
                        alt_in_cds: pending[i].alt_in_cds,
                        carriers: pending[i].carriers,
                    })
                    .collect();
                for (&slot, call) in slots.iter().zip(call_codon(&ref_codon, ref_aa, &group, gc)) {
                    calls[slot] = call;
                }
            }

            // Pass 3: emit the rows and accumulate the counts, in record order.
            for (p, call) in pending.iter().zip(&calls) {
                let Some(call) = call else { continue };
                let ref_codon = ref_codons[&p.codon_idx];
                let ref_aa = codon_to_aa(&ref_codon, gc).expect("checked in pass 1");
                variants.push(Variant {
                    pos: p.pos,
                    ref_allele: p.ref_allele,
                    alt_allele: p.alt_allele,
                    aa_pos: p.codon_idx + 1,
                    ref_codon,
                    alt_codon: call.alt_codon,
                    ref_aa,
                    alt_aa: call.alt_aa,
                    af: p.af,
                    gt_derived: p.gt_derived,
                    gt_called: p.gt_called,
                    effect: call.effect,
                    // Filled in below, once every SNP of the gene is known and
                    // the codons can be grouped.
                    codon_shared: None,
                    // Filled in by `AlleleOrigins::attach` when a tree was supplied.
                    // Keyed back to the SNP records rather than plumbed through here, so
                    // an allele shared by two overlapping genes costs one carrier bitset
                    // and not one per row.
                    origins: None,
                });
                variant_carriers.push(p.carriers);

                if call.diffs > 1 {
                    gene_mnv_alleles += 1;
                    // A joint change to a stop that this ALT does not create on its own is
                    // the case the per-SNP scoring reported as an ordinary missense.
                    if call.alt_aa == b'*' {
                        let mut solo = ref_codon;
                        solo[p.pos_in_codon] = p.alt_in_cds;
                        if codon_to_aa(&solo, gc) != Some(b'*') {
                            gene_new_stops.insert(p.codon_idx);
                        }
                    }
                }
                // Alleles whose count is a blend of a solo and a joint background have
                // `diffs == 1` but still split between the two classes, so the tally of
                // what the joint scoring touched is the wider `joint`, not `diffs`.
                if call.joint && call.diffs == 1 {
                    gene_mnv_alleles += 1;
                }

                // pN/pS counting: synonymous and missense only (changes to/from a stop are
                // excluded, as before). `syn + nonsyn == 1` per classified allele, so a
                // multi-nucleotide change still contributes its nucleotide count.
                if !call.counted {
                    continue;
                }
                n_snps_raw += 1;
                syn_count += p.weight * call.syn;
                nonsyn_count += p.weight * call.nonsyn;

                // The McDonald-Kreitman 2x2 and the SFS bins need whole alleles in one
                // class, so they take only alleles that stand alone in their codon for
                // EVERY sample carrying them. An allele carried on both a solo and a joint
                // background splits between the two classes just as surely as one carried
                // only jointly, so it leaves too, and is counted out loud rather than
                // rounded into a cell it only mostly belongs in.
                if call.joint {
                    mk_mnv_excluded += 1;
                    continue;
                }
                let fixed = p.af >= mk_fixed_af;
                let bin = sfs_bin(p.af);
                if call.effect == SnpEffect::Synonymous {
                    sfs_syn[bin] += 1;
                    if fixed { mk_ds += 1 } else { mk_ps += 1 }
                } else {
                    sfs_nonsyn[bin] += 1;
                    if fixed { mk_dn += 1 } else { mk_pn += 1 }
                }
            }

            // Same-codon (MNV) bookkeeping. Pass 2 above already scored each ALT against
            // the codon its carriers hold, so this no longer decides what is counted; it
            // decides what is REPORTED: which residues carry a multi-nucleotide change,
            // and, for an input with no genotypes, which ones eskaks could not check.
            //
            // Where carriers survived parsing this is an OBSERVATION, not an inference:
            // the genotypes are haploid, so two alleles listed for one sample are two
            // alleles in one molecule. The allele-frequency bound is used only where no
            // carriers exist, which is the one input for which it is the best available
            // answer, and the one input whose joint changes stay unevaluated.
            let (mut gene_multi, mut gene_cooccurring) = (0usize, 0usize);
            let (mut gene_exact, mut gene_known) = (0usize, 0usize);
            // Codons carrying more than one SNP position, and what is known about each.
            let mut codon_state: HashMap<usize, CodonShare<'_>> = HashMap::new();
            if codon_hits.len() > 1 {
                // Sort by codon, then position, then descending AF, so collapsing repeated
                // positions (a VCF may split one site across records) keeps the highest AF.
                codon_hits.sort_by(|a, b| {
                    a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(b.2.total_cmp(&a.2))
                });
                // Duplicated records describe ONE position, so their carriers are unioned
                // into the survivor rather than discarded with the losing record: dropping
                // them would lose samples that only the second record listed.
                codon_hits.dedup_by(|dropped, kept| {
                    if dropped.0 != kept.0 || dropped.1 != kept.1 {
                        return false;
                    }
                    match (kept.3.as_mut(), dropped.3.as_ref()) {
                        (Some(k), Some(d)) => k.union_with(d),
                        // One record without carriers makes the merged position unknown.
                        _ => kept.3 = None,
                    }
                    true
                });
                for group in codon_hits.chunk_by(|a, b| a.0 == b.0) {
                    if group.len() < 2 {
                        continue;
                    }
                    gene_multi += 1;
                    let (state, cooccurring) = if group.iter().all(|h| h.3.is_some()) {
                        gene_known += 1;
                        let positions: Vec<(usize, &CarrierSet)> = group
                            .iter()
                            .map(|h| (h.1, h.3.as_ref().expect("checked by the all() above")))
                            .collect();
                        // Exact: is there a sample carrying variants at two of these
                        // positions? k is 2 or 3 (a codon has three bases), so the
                        // pairwise scan is at most three intersections.
                        let shared = (0..positions.len()).any(|i| {
                            ((i + 1)..positions.len())
                                .any(|j| positions[i].1.intersects(positions[j].1))
                        });
                        if shared {
                            gene_exact += 1;
                        }
                        (CodonShare::Observed { positions, shared }, shared)
                    } else {
                        let afs: Vec<f64> = group.iter().map(|h| h.2).collect();
                        if forced_cooccurrence(&afs) {
                            (CodonShare::ForcedByFrequency, true)
                        } else {
                            (CodonShare::Unknown, false)
                        }
                    };
                    if cooccurring {
                        gene_cooccurring += 1;
                    }
                    let detail: Vec<String> = group
                        .iter()
                        .map(|h| match &h.3 {
                            Some(cs) => format!("{} (AF {:.3}, {} carrier(s))", h.1, h.2, cs.len()),
                            None => format!("{} (AF {:.3})", h.1, h.2),
                        })
                        .collect();
                    debug!(
                        "Gene {} residue {}: {} SNPs share one codon at {}:{}. {}",
                        gene.name,
                        group[0].0 + 1,
                        group.len(),
                        gene.seqid,
                        detail.join(", "),
                        state.explanation()
                    );
                    codon_state.insert(group[0].0, state);
                }
            }
            // Mark the variant rows a reader must not trust: this ALT is carried by a
            // sample that also carries another SNP in the same codon, so the amino-acid
            // change on the row is not the one that genome encodes.
            for (v, alt_carriers) in variants.iter_mut().zip(&variant_carriers) {
                v.codon_shared = match codon_state.get(&(v.aa_pos - 1)) {
                    // The codon carries no other SNP position: scoring against the
                    // reference codon is exact, whatever the input offered.
                    None => Some(false),
                    Some(CodonShare::Observed { positions, .. }) => alt_carriers.map(|cs| {
                        positions.iter().any(|(p, other)| *p != v.pos && cs.intersects(other))
                    }),
                    // No carriers to attribute the overlap to a particular ALT, so every
                    // row of a forced codon is marked: the column exists to warn, and
                    // under-warning here would be the expensive mistake.
                    Some(CodonShare::ForcedByFrequency) => Some(true),
                    Some(CodonShare::Unknown) => None,
                };
            }
            multi_snp_codons.fetch_add(gene_multi, Ordering::Relaxed);
            cooccurring_codons.fetch_add(gene_cooccurring, Ordering::Relaxed);
            cooccurring_exact.fetch_add(gene_exact, Ordering::Relaxed);
            codons_with_carriers.fetch_add(gene_known, Ordering::Relaxed);
            mnv_alleles.fetch_add(gene_mnv_alleles, Ordering::Relaxed);
            mnv_excluded_from_mk.fetch_add(mk_mnv_excluded as usize, Ordering::Relaxed);
            mnv_new_stops.fetch_add(gene_new_stops.len(), Ordering::Relaxed);
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
            //
            // `total_snps` is the number of classified alleles exactly (each contributes
            // 1.0, however its codon was scored), so `n` is exact. `k` can carry a half
            // or a quarter when a multi-nucleotide change split between the two classes,
            // and is rounded: unlike the McDonald-Kreitman table, this test is a statement
            // about the same counts the reported ratio is built from, so restricting it to
            // distance-1 alleles would test a different gene from the one tabulated. The
            // rounding moves `k` by at most 0.5 of one SNP, and only in a gene that has a
            // multi-nucleotide change at all.
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
                mk_mnv_excluded,
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
    let cooccurring_exact = cooccurring_exact.load(Ordering::Relaxed);
    let codons_with_carriers = codons_with_carriers.load(Ordering::Relaxed);

    let diagnostics = ComputeDiagnostics {
        snps_in_cds: ref_checked,
        ref_mismatch,
        genes_with_internal_stops: internal_stops.load(Ordering::Relaxed),
        multi_snp_codons,
        cooccurring_codons,
        genes_with_cooccurring,
        cooccurring_exact,
        codons_with_carriers,
        mnv_alleles: mnv_alleles.load(Ordering::Relaxed),
        mnv_excluded_from_mk: mnv_excluded_from_mk.load(Ordering::Relaxed),
        mnv_new_stops: mnv_new_stops.load(Ordering::Relaxed),
    };
    for (level, message) in same_codon_messages(&diagnostics) {
        log::log!(level, "{}", message);
    }

    info!("Computed pN/pS for {} genes", results.len());
    (results, diagnostics)
}

/// The run-level same-codon report, as `(level, message)` pairs.
///
/// Built as data rather than written straight to the log so the wording can be asserted
/// on directly: the `log` facade allows one sink per process, so a test that greps the
/// captured lines is really asserting about every other test running beside it.
///
/// Four disjoint groups, and the split is the substance rather than presentation. Where
/// carriers survived parsing, both the positive and the negative are facts about named
/// samples, and the codon those samples carry is what gets scored; where they did not,
/// the positive is a frequency bound (a floor, so its count understates), the negative is
/// not a negative at all but an absence of proof, and nothing can be scored jointly.
///
/// That is also why only the second group warns: a warning should fire where the answer
/// may be wrong, and after this stage that is exactly the input eskaks cannot check.
pub(crate) fn same_codon_messages(d: &ComputeDiagnostics) -> Vec<(log::Level, String)> {
    let mut out = Vec::new();
    // Co-occurrences the carrier sets proved, and those only the frequency bound reached.
    let cooccurring_bound = d.cooccurring_codons.saturating_sub(d.cooccurring_exact);
    // Observed co-occurrence: a named sample carries two SNPs of one codon, so the codon
    // it carries is the one that was scored. Reported, not warned about: it is handled.
    if d.cooccurring_exact > 0 {
        out.push((
            log::Level::Info,
            format!(
                "{} codon(s) in {} gene(s) carry two or more SNPs that at least one sample carries \
                 TOGETHER, read from the per-sample genotypes (eskaks assumes haploid genotypes, so \
                 two alleles in one sample lie on one molecule). Each was scored as the codon those \
                 samples actually carry, by Nei-Gojobori pathway averaging over the orders the \
                 changes could have happened in, so the joint (multi-nucleotide) change is what the \
                 counts and the `--variants` rows report. {} allele(s) were classified that way, of \
                 which {} are excluded from the McDonald-Kreitman table and the site-frequency \
                 spectrum, which need whole alleles. Run with -vv to list the codons.",
                d.cooccurring_exact, d.genes_with_cooccurring, d.mnv_alleles, d.mnv_excluded_from_mk
            ),
        ));
    }
    // A joint stop that no single SNP creates: the case where the per-SNP classification
    // was not merely imprecise but reported a functioning protein where there is none.
    if d.mnv_new_stops > 0 {
        out.push((
            log::Level::Warn,
            format!(
                "{} codon(s) encode a premature STOP only when their SNPs are taken together: \
                 each SNP alone leaves a coding amino acid. They are now classified as nonsense \
                 and, like any change to a stop, excluded from the pN/pS counts; `--variants` \
                 reports them as `nonsense`. Scored one SNP at a time they would read as ordinary \
                 missense or synonymous changes, so check these residues before trusting any \
                 earlier annotation of them.",
                d.mnv_new_stops
            ),
        ));
    }
    // No genotypes to check: the allele-frequency pigeonhole bound is the only sound
    // inference left, and it is stated as the bound it is rather than as an observation.
    if cooccurring_bound > 0 {
        out.push((
            log::Level::Warn,
            format!(
                "{} codon(s) carry two or more SNPs with no per-sample genotypes to check, but \
                 their allele frequencies leave no room for the SNPs to avoid one another (sum of \
                 AF > k - 1), so some sample must carry them together. With no genotypes there is \
                 no codon to score jointly, so each was classified against the REFERENCE codon on \
                 its own and the joint change was never evaluated: those calls, and the pN/pS they \
                 feed, can be wrong. The bound is a floor, so the real number can only be higher. \
                 Supply per-sample genotypes to have the codon each sample carries scored instead. \
                 Run with -vv to list the codons.",
                cooccurring_bound
            ),
        ));
    }
    // Exact negatives: carriers were available and no sample carries two of them, so the
    // independent per-SNP classification is not merely defensible here, it is correct.
    let exact_apart = d.codons_with_carriers.saturating_sub(d.cooccurring_exact);
    if exact_apart > 0 {
        out.push((
            log::Level::Info,
            format!(
                "{} codon(s) carry more than one SNP but NO sample carries two of them (checked \
                 against the per-sample genotypes), so classifying each against the reference \
                 codon is exact for those codons.",
                exact_apart
            ),
        ));
    }
    // The genuine unknowns: no carriers, and frequencies that force nothing.
    let phase_unknown = d
        .multi_snp_codons
        .saturating_sub(d.codons_with_carriers)
        .saturating_sub(cooccurring_bound);
    if phase_unknown > 0 {
        out.push((
            log::Level::Info,
            format!(
                "{} further codon(s) carry more than one SNP with neither per-sample genotypes nor \
                 a forcing allele-frequency bound, so whether one genome carries two of them is \
                 unknown and each was classified against the reference codon (run with -vv to list \
                 them).",
                phase_unknown
            ),
        ));
    }
    out
}

/// What is known about the SNPs sharing one codon, and on what evidence.
///
/// The distinction is the point of this stage: `Observed` answers the question, in both
/// directions, from named samples; the other two are what an input without per-sample
/// genotypes leaves you with.
enum CodonShare<'a> {
    /// Every position in the codon had a carrier set. Each entry of `positions` is a
    /// distinct genomic position paired with the samples carrying a variant base there,
    /// so any single allele can be asked whether it shares a sample with another
    /// position; `shared` is that question answered over the codon as a whole.
    Observed { positions: Vec<(usize, &'a CarrierSet)>, shared: bool },
    /// No carriers, but the allele frequencies leave no room: some sample must carry two
    /// of them (see [`forced_cooccurrence`]).
    ForcedByFrequency,
    /// No carriers and no forcing bound: whether one genome carries two is unknown.
    Unknown,
}

impl CodonShare<'_> {
    /// The `-vv` sentence explaining what this codon's classification is worth.
    fn explanation(&self) -> &'static str {
        match self {
            CodonShare::Observed { shared: true, .. } => {
                "At least one sample carries two of them, so this residue was scored as the codon \
                 those samples carry, by Nei-Gojobori pathway averaging, rather than one SNP at a \
                 time against the reference codon."
            }
            CodonShare::Observed { shared: false, .. } => {
                "No sample carries more than one of them, so scoring each against the reference \
                 codon is exact here."
            }
            CodonShare::ForcedByFrequency => {
                "There are no per-sample genotypes, but their allele frequencies leave no room to \
                 avoid one another, so some sample carries two. With no genotypes there is no \
                 codon to score jointly, so the per-SNP calls reported for this residue are not \
                 what it encodes."
            }
            CodonShare::Unknown => {
                "There are no per-sample genotypes and the frequencies force nothing, so each SNP \
                 was scored against the reference codon on its own."
            }
        }
    }
}

/// Do the allele frequencies of the SNPs sharing one codon force them onto the same
/// haplotype?
///
/// This is the fallback for an input with **no per-sample genotypes**. Where carriers
/// survived parsing, [`compute_pn_ps`] intersects them instead and gets the exact answer
/// in both directions; this bound is a floor, so it misses real co-occurrences whose
/// frequencies happen not to prove themselves. It is kept, unweakened, because for an
/// AF-only VCF it is the only sound inference available.
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
    /// Codons carrying two or more retained SNPs at distinct positions.
    pub multi_snp_codons: usize,
    /// The subset of `multi_snp_codons` where the SNPs provably occur together in one
    /// genome. With per-sample genotypes these are scored as the codon those samples
    /// carry; without them, they are the codons whose joint change stays unevaluated.
    pub cooccurring_codons: usize,
    /// Genes contributing at least one `cooccurring_codons` codon.
    pub genes_with_cooccurring: usize,
    /// Of `cooccurring_codons`, those established by intersecting per-sample carrier
    /// sets. The remainder come from the allele-frequency bound, which is available for
    /// an AF-only VCF and is a floor rather than a count.
    pub cooccurring_exact: usize,
    /// Multi-SNP codons where every position had a carrier set, so both the positive and
    /// the negative answer are exact. `multi_snp_codons - codons_with_carriers` is how
    /// many were left to the frequency bound.
    pub codons_with_carriers: usize,
    /// ALT alleles whose carriers hold a codon differing from the reference at more than
    /// one base, so the allele was scored as part of a multi-nucleotide change by
    /// Nei-Gojobori pathway averaging rather than on its own.
    pub mnv_alleles: usize,
    /// Of those, how many were classified into pN/pS but kept out of the
    /// McDonald-Kreitman 2x2 and the SFS bins, which need whole alleles.
    pub mnv_excluded_from_mk: usize,
    /// Codons whose joint change is a premature **stop** that none of their SNPs creates
    /// on its own. Scored one SNP at a time these read as ordinary missense or synonymous
    /// changes, which is the most damaging form the old classification took.
    pub mnv_new_stops: usize,
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
    ///
    /// Every caller filters on something unique to its own gene name: the sink is
    /// process-wide, so a needle that another test could also emit would make this an
    /// assertion about the whole test binary rather than about one run.
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

    /// The run-level messages for a diagnostics record, at one level, as plain strings.
    fn messages(d: &ComputeDiagnostics, want: log::Level) -> Vec<String> {
        same_codon_messages(d)
            .into_iter()
            .filter(|(level, _)| *level == want)
            .map(|(_, m)| m)
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

    /// Single-ALT SNP builder with no per-sample identity: the AF-only input, where the
    /// frequency bound is all there is.
    fn snp(pos: usize, ref_allele: u8, alt: u8, af: f64) -> VcfSnp {
        VcfSnp {
            chrom: "chr1".to_string(),
            pos,
            ref_allele,
            alt_alleles: vec![alt],
            alt_freqs: vec![af],
            gt_counts: None,
            carriers: None,
            filter: "PASS".to_string(),
            depth: None,
        }
    }

    /// The same SNP with its carriers retained: `samples` are the indices of the haploid
    /// samples carrying the ALT, out of `n`. `af` stays a free parameter so a test can
    /// make the frequencies and the genotypes tell different stories, which is exactly
    /// the case the exact detector exists for.
    fn snp_carried(
        pos: usize,
        ref_allele: u8,
        alt: u8,
        af: f64,
        samples: &[usize],
        n: usize,
    ) -> VcfSnp {
        let mut set = CarrierSet::new(n);
        for &s in samples {
            set.insert(s);
        }
        VcfSnp { carriers: Some(vec![set]), ..snp(pos, ref_allele, alt, af) }
    }

    #[test]
    fn fixed_snps_in_one_codon_warn_and_are_still_scored_apart() {
        // Reference CDS ATG CTT TAA with chr1:4 C>T and chr1:6 T>A, both at AF 1.0 and no
        // genotype columns. Scored one at a time against the reference codon CTT: C>T
        // gives TTT (Phe, missense L2F) and T>A gives CTA (Leu, synonymous). The haplotype
        // every sample actually carries is TTA, which is Leu: one synonymous
        // multi-nucleotide change. Calling it is out of scope, so the run must SAY so, and
        // must say it as the frequency BOUND it is rather than as an observation.
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

        // Detection: one codon, forced by the frequencies, in one gene, and NOT counted
        // as exact, because nothing about a named sample was known.
        assert_eq!(diag.multi_snp_codons, 1, "codon 2 carries both SNPs");
        assert_eq!(diag.cooccurring_codons, 1, "AF 1.0 + 1.0 forces co-occurrence");
        assert_eq!(diag.genes_with_cooccurring, 1);
        assert_eq!(diag.cooccurring_exact, 0, "no carriers, so nothing was observed");
        assert_eq!(diag.codons_with_carriers, 0);

        // Exactly one warning, and it is stated as the bound it is: an AF-only input must
        // never claim a sample was seen carrying both.
        let warns = messages(&diag, log::Level::Warn);
        assert_eq!(warns.len(), 1, "expected exactly one same-codon warning, got {:?}", warns);
        assert!(warns[0].starts_with("1 codon(s) carry"), "warning text: {}", warns[0]);
        assert!(warns[0].contains("allele frequencies leave no room"), "{}", warns[0]);
        assert!(warns[0].contains("floor"), "the bound must be labelled a floor: {}", warns[0]);
        assert!(
            !warns[0].contains("carries TOGETHER"),
            "an AF-only input must not claim a sample was observed carrying both: {}",
            warns[0]
        );
        assert!(messages(&diag, log::Level::Info).is_empty(), "nothing left to report");

        // -vv detail names the residue and both positions, so the user can act on it.
        let details = captured(from, "DEBUG", "Gene leu_fixed residue");
        assert_eq!(details.len(), 1, "expected one debug detail line, got {:?}", details);
        assert!(
            details[0].contains("residue 2")
                && details[0].contains("chr1:4 (AF 1.000)")
                && details[0].contains("6 (AF 1.000)"),
            "debug detail: {}",
            details[0]
        );

        // Every row of the codon is marked, since without carriers the overlap cannot be
        // attributed to one ALT and under-warning is the expensive mistake.
        assert_eq!(res[0].variants.len(), 2);
        assert!(res[0].variants.iter().all(|v| v.codon_shared == Some(true)));
    }

    #[test]
    fn a_shared_sample_is_detected_where_the_frequency_bound_gives_up() {
        // The whole point of retaining carriers. Two SNPs in codon 2 at AF 0.3 and 0.2:
        // the pigeonhole bound proves nothing (0.5 <= 1), and before carriers existed the
        // run said "phase unknown". But samples 0-3 are listed for BOTH, and the genotypes
        // are haploid, so those genomes carry TTA. On the bundled 20-sample toy VCF this
        // case is 5 of the 13 real co-occurrences, the 38% the bound was missing.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_shared");
        let snps = vec![
            snp_carried(4, b'C', b'T', 0.3, &[0, 1, 2, 3, 4, 5], 20),
            snp_carried(6, b'T', b'A', 0.2, &[0, 1, 2, 3], 20),
        ];

        let _serialised = LOG_TESTS.lock().unwrap_or_else(|e| e.into_inner());
        let from = start_capture();
        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        assert!(
            !forced_cooccurrence(&[0.3, 0.2]),
            "the frequency bound must genuinely fail here, or this test proves nothing"
        );
        assert_eq!(diag.multi_snp_codons, 1);
        assert_eq!(diag.cooccurring_codons, 1, "the carriers overlap, so it is detected");
        assert_eq!(diag.cooccurring_exact, 1, "and detected exactly, not by a bound");
        assert_eq!(diag.codons_with_carriers, 1);
        assert_eq!(diag.genes_with_cooccurring, 1);

        // Both alleles are now scored on the codon their carriers hold. Samples 0-3 carry
        // TTA (Leu, k = 2, pathway score (1, 1), so half a unit each way per base);
        // samples 4 and 5 carry C>T alone, which is TTT (Phe, nonsynonymous).
        //   C>T: (4 * 0.5 + 2 * 0) / 6 syn, (4 * 0.5 + 2 * 1) / 6 nonsyn
        //   T>A: carried only alongside C>T, so 0.5 / 0.5
        assert!((res[0].syn_snps - (1.0 / 3.0 + 0.5)).abs() < 1e-12, "syn {}", res[0].syn_snps);
        assert!(
            (res[0].nonsyn_snps - (2.0 / 3.0 + 0.5)).abs() < 1e-12,
            "nonsyn {}",
            res[0].nonsyn_snps
        );
        // Two alleles, so two units, whatever the split: this is the invariant the
        // per-nucleotide site denominator depends on.
        assert!((res[0].total_snps - 2.0).abs() < 1e-12, "total {}", res[0].total_snps);
        assert_eq!(res[0].n_snps, 2);

        // The report is stated as the observation it is, and is INFO rather than a
        // warning: the joint change was evaluated, so there is nothing to caveat.
        assert!(
            messages(&diag, log::Level::Warn).is_empty(),
            "a codon eskaks scored correctly must not warn: {:?}",
            messages(&diag, log::Level::Warn)
        );
        let infos = messages(&diag, log::Level::Info);
        assert_eq!(infos.len(), 1, "expected one observed-co-occurrence line, got {:?}", infos);
        assert!(infos[0].starts_with("1 codon(s) in 1 gene(s)"), "text: {}", infos[0]);
        assert!(infos[0].contains("carries TOGETHER"), "{}", infos[0]);
        assert!(
            infos[0].contains("scored as the codon those samples actually carry"),
            "the message must say the joint change was evaluated: {}",
            infos[0]
        );
        assert!(
            !infos[0].contains("allele frequencies leave no room"),
            "an observation must not be reported as a frequency bound: {}",
            infos[0]
        );
        let details = captured(from, "DEBUG", "Gene leu_shared residue");
        assert!(
            details[0].contains("chr1:4 (AF 0.300, 6 carrier(s))")
                && details[0].contains("6 (AF 0.200, 4 carrier(s))"),
            "debug detail: {}",
            details[0]
        );

        // Both rows are marked: each ALT is carried by a sample that also carries the other.
        assert!(res[0].variants.iter().all(|v| v.codon_shared == Some(true)));
        // Each row reports the codon most of THAT allele's carriers hold. C>T is carried
        // by 6 samples, 4 of them alongside T>A, so its row reports the joint TTA (Leu)
        // rather than the TTT (Phe) that only 2 of them have -- and its count still keeps
        // the TTT minority, which is why its split is not a clean half.
        let by_pos = |pos: usize| {
            res[0].variants.iter().find(|v| v.pos == pos).expect("row exists").clone()
        };
        assert_eq!(&by_pos(4).alt_codon, b"TTA");
        assert_eq!(by_pos(4).alt_aa, b'L');
        assert_eq!(&by_pos(6).alt_codon, b"TTA");
        assert_eq!(by_pos(6).alt_aa, b'L');
        assert!(
            !res[0].variants.iter().any(|v| v.alt_aa == b'F'),
            "no row may report a residue a minority of carriers has"
        );
    }

    #[test]
    fn the_headline_case_reports_no_change_the_genomes_do_not_carry() {
        // The bug report, end to end. Reference CDS ATG CTT TAA with chr1:4 C>T and
        // chr1:6 T>A, both at AF 1.0 and both carried by every one of the 20 samples.
        // Scored apart that is a missense L2F plus a synonymous L2L. The codon every
        // genome holds is TTA, which is Leu: no amino-acid change at all, and the L2F is a
        // change no sample carries.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_joint");
        let all: Vec<usize> = (0..20).collect();
        let snps = vec![
            snp_carried(4, b'C', b'T', 1.0, &all, 20),
            snp_carried(6, b'T', b'A', 1.0, &all, 20),
        ];
        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        // No row claims a phenylalanine.
        assert_eq!(res[0].variants.len(), 2);
        for v in &res[0].variants {
            assert_eq!(&v.alt_codon, b"TTA", "the codon every sample carries");
            assert_eq!(v.alt_aa, b'L', "TTA is Leu, so the residue is unchanged");
            assert_eq!(v.effect, SnpEffect::Synonymous);
            assert_eq!(v.codon_shared, Some(true));
        }
        assert!(
            !res[0].variants.iter().any(|v| v.alt_aa == b'F'),
            "the missense L2F no genome carries must be gone"
        );

        // The counts: Nei-Gojobori pathway averaging over the two orders the changes could
        // have happened in gives (1 syn, 1 nonsyn) for CTT -> TTA. That is 2 units for 2
        // nucleotide differences, which is what keeps the numerator commensurate with the
        // per-nucleotide site denominator. Counting the MNV as one event would give 1.
        assert!((res[0].syn_snps - 1.0).abs() < 1e-12, "syn {}", res[0].syn_snps);
        assert!((res[0].nonsyn_snps - 1.0).abs() < 1e-12, "nonsyn {}", res[0].nonsyn_snps);
        assert!((res[0].total_snps - 2.0).abs() < 1e-12, "Nd + Sd must be k = 2");

        // Both alleles leave the McDonald-Kreitman table and the SFS, and say so.
        assert_eq!(res[0].mk_dn + res[0].mk_ds + res[0].mk_pn + res[0].mk_ps, 0);
        assert_eq!(res[0].mk_mnv_excluded, 2);
        assert_eq!(res[0].sfs_nonsyn.iter().sum::<u32>() + res[0].sfs_syn.iter().sum::<u32>(), 0);
        assert_eq!(diag.mnv_alleles, 2);
        assert_eq!(diag.mnv_excluded_from_mk, 2);
        assert_eq!(diag.mnv_new_stops, 0, "TTA is Leu, not a stop");
    }

    #[test]
    fn a_joint_stop_is_reported_as_nonsense_and_leaves_the_counts() {
        // The most damaging shape the old classification took. Reference CDS ATG CAT TAA:
        // chr1:4 C>T alone gives TAT (Tyr, missense) and chr1:6 T>A alone gives CAA (Gln,
        // missense), so a per-SNP run reports two ordinary missense changes. Every sample
        // carries both, and together they give TAA: a premature stop.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"ATGCATTAA".to_vec());
        let gene = Gene {
            name: "his_stop".to_string(),
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
        let all: Vec<usize> = (0..20).collect();
        let snps = vec![
            snp_carried(4, b'C', b'T', 1.0, &all, 20),
            snp_carried(6, b'T', b'A', 1.0, &all, 20),
        ];
        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        for v in &res[0].variants {
            assert_eq!(&v.alt_codon, b"TAA");
            assert_eq!(v.effect, SnpEffect::Nonsense, "the genomes carry a truncation");
        }
        // A change to a stop has always been excluded from the pN/pS counts, so the gene
        // now has no classified coding SNP at all rather than two missense ones.
        assert_eq!(res[0].nonsyn_snps, 0.0);
        assert_eq!(res[0].syn_snps, 0.0);
        assert_eq!(res[0].n_snps, 0);
        assert_eq!(diag.mnv_new_stops, 1, "one codon gains a stop only jointly");
        assert_eq!(diag.mnv_excluded_from_mk, 0, "they were never in the pN/pS counts");

        // And it is worth a warning: nothing else in the run would tell a reader that a
        // residue annotated as two missense changes is really a truncation.
        let warns = messages(&diag, log::Level::Warn);
        assert_eq!(warns.len(), 1, "expected the joint-stop warning, got {:?}", warns);
        assert!(warns[0].contains("premature STOP"), "{}", warns[0]);
    }

    #[test]
    fn disjoint_carriers_in_one_codon_are_not_a_multi_nucleotide_change() {
        // The other direction, which the frequency bound can never give: two SNPs at AF
        // 0.6 and 0.6 in one codon, whose bound SAYS they must overlap (1.2 > 1), but
        // whose genotypes put them in disjoint samples. The bound is a statement about
        // the whole cohort; here the cohort is 20 and the two sets of 12 and 12 would
        // have to overlap -- so this is the case where the declared AF disagrees with the
        // genotypes, and the genotypes are the ones that decide.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_apart");
        let snps = vec![
            snp_carried(4, b'C', b'T', 0.6, &[0, 1, 2, 3], 20),
            snp_carried(6, b'T', b'A', 0.6, &[10, 11, 12], 20),
        ];

        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        assert!(forced_cooccurrence(&[0.6, 0.6]), "the bound would have flagged this codon");
        assert_eq!(diag.multi_snp_codons, 1);
        assert_eq!(diag.cooccurring_codons, 0, "no sample carries both, so nothing is wrong");
        assert_eq!(diag.cooccurring_exact, 0);
        assert_eq!(diag.codons_with_carriers, 1, "the negative is exact, not 'unknown'");
        assert!(
            messages(&diag, log::Level::Warn).is_empty(),
            "an observed non-overlap must not warn"
        );
        let infos = messages(&diag, log::Level::Info);
        assert_eq!(infos.len(), 1);
        assert!(
            infos[0].contains("NO sample carries two of them"),
            "it is reported as the exact negative it is: {}",
            infos[0]
        );
        assert!(
            res[0].variants.iter().all(|v| v.codon_shared == Some(false)),
            "no row of this codon may be marked"
        );
    }

    #[test]
    fn only_the_alleles_that_actually_share_a_sample_are_marked() {
        // Three SNP positions in one codon is the case where a per-codon mark would be
        // too coarse: positions 4 and 5 share sample 0, position 6 is carried by nobody
        // else. Its row is scored against the reference codon correctly, so marking it
        // would send a reader chasing a variant that is fine.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_partial");
        let snps = vec![
            snp_carried(4, b'C', b'T', 0.1, &[0, 1], 20),
            snp_carried(5, b'T', b'C', 0.1, &[0, 2], 20),
            snp_carried(6, b'T', b'A', 0.1, &[7], 20),
        ];
        let (res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        assert_eq!(diag.cooccurring_codons, 1, "the codon as a whole is affected");
        let marks: Vec<(usize, Option<bool>)> =
            res[0].variants.iter().map(|v| (v.pos, v.codon_shared)).collect();
        assert_eq!(
            marks,
            vec![(4, Some(true)), (5, Some(true)), (6, Some(false))],
            "only the two alleles sharing sample 0 are marked"
        );
    }

    #[test]
    fn one_position_split_across_records_keeps_every_carrier() {
        // A VCF may list one site on more than one record. The duplicates collapse to a
        // single position, and their carriers must be UNIONED into the survivor: sample 9
        // appears only on the second chr1:4 record, and it is the sample that also carries
        // chr1:6. Keeping the highest-AF record and discarding the other one's carriers
        // would silently lose the only overlap there is.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = leu_gene("leu_split_record");
        let snps = vec![
            snp_carried(4, b'C', b'T', 0.2, &[0, 1], 20),
            snp_carried(4, b'C', b'T', 0.1, &[9], 20),
            snp_carried(6, b'T', b'A', 0.1, &[9], 20),
        ];
        let (_res, diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);

        assert_eq!(diag.multi_snp_codons, 1, "chr1:4 twice is one position, not two");
        assert_eq!(diag.codons_with_carriers, 1);
        assert_eq!(
            diag.cooccurring_exact, 1,
            "sample 9 carries both, and it is listed only on the lower-AF chr1:4 record"
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
        assert_eq!(diag.codons_with_carriers, 0, "there was nothing to check exactly");
        assert!(
            messages(&diag, log::Level::Warn).is_empty(),
            "a pair with no carriers and no forcing bound must not warn"
        );
        let infos = messages(&diag, log::Level::Info);
        assert_eq!(infos.len(), 1);
        assert!(
            infos[0].contains("neither per-sample genotypes nor a forcing"),
            "it is reported as the genuine unknown it is: {}",
            infos[0]
        );
        // The -vv detail says the same thing, scoped to this gene.
        let details = captured(from, "DEBUG", "Gene leu_rare residue");
        assert_eq!(details.len(), 1, "expected one debug detail line, got {:?}", details);
        assert!(details[0].contains("frequencies force nothing"), "detail: {}", details[0]);
    }

    #[test]
    fn a_mixed_run_reports_both_bases_separately_and_adds_up() {
        // The four groups are disjoint and must partition the multi-SNP codons, or a
        // reader adding the lines up gets a number that is not the one at the top. Built
        // from a diagnostics record rather than a run, so every combination is reachable:
        // 10 multi-SNP codons, 6 with carriers (4 of them shared), 2 more forced by the
        // frequency bound, leaving 2 genuinely unknown.
        let d = ComputeDiagnostics {
            multi_snp_codons: 10,
            cooccurring_codons: 6,
            cooccurring_exact: 4,
            codons_with_carriers: 6,
            genes_with_cooccurring: 3,
            mnv_alleles: 9,
            mnv_excluded_from_mk: 7,
            ..ComputeDiagnostics::default()
        };
        let msgs = same_codon_messages(&d);
        assert_eq!(msgs.len(), 4, "all four groups are non-empty here: {:?}", msgs);
        assert!(msgs[0].1.starts_with("4 codon(s) in 3 gene(s)"), "{}", msgs[0].1);
        assert!(msgs[1].1.starts_with("2 codon(s) carry"), "{}", msgs[1].1);
        assert!(msgs[2].1.starts_with("2 codon(s) carry more than one SNP but NO"), "{}", msgs[2].1);
        assert!(msgs[3].1.starts_with("2 further codon(s)"), "{}", msgs[3].1);
        // 4 observed + 2 bounded + 2 exact-apart + 2 unknown = the 10 at the top.
        assert_eq!((4 + 2) + 2 + 2, d.multi_snp_codons);
        // Only the frequency-bound group warns: it is the only one whose joint change was
        // not evaluated, so it is the only one whose numbers can be wrong.
        assert_eq!(msgs[0].0, log::Level::Info);
        assert_eq!(msgs[1].0, log::Level::Warn);
        assert_eq!(msgs[2].0, log::Level::Info);
        assert_eq!(msgs[3].0, log::Level::Info);
        // The observed line carries the allele counts, so a reader can reconcile the
        // McDonald-Kreitman totals with the pN/pS ones without running anything again.
        assert!(msgs[0].1.contains("9 allele(s)") && msgs[0].1.contains("7 are excluded"),
            "{}", msgs[0].1);

        // The joint-stop warning is its own group, so it cannot be crowded out.
        let with_stops = ComputeDiagnostics { mnv_new_stops: 3, ..d };
        let msgs = same_codon_messages(&with_stops);
        assert_eq!(msgs.len(), 5);
        assert_eq!(msgs[1].0, log::Level::Warn);
        assert!(msgs[1].1.starts_with("3 codon(s) encode a premature STOP"), "{}", msgs[1].1);
    }

    #[test]
    fn a_clean_run_says_nothing_about_same_codon_snps() {
        assert!(same_codon_messages(&ComputeDiagnostics::default()).is_empty());
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
            carriers: None,
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
