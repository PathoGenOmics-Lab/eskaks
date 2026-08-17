//! Scoring the codon each sample actually carries.
//!
//! eskaks used to build every alternate codon by copying the **reference** codon and
//! mutating one position, so two SNPs in one codon were each scored against the
//! reference and their joint change was never evaluated. This module removes that
//! assumption for any input that carries per-sample genotypes.
//!
//! The whole thing rests on one fact: eskaks assumes **haploid** genotypes (the
//! *M. tuberculosis* setting it was written for). Two alleles listed for one sample are
//! therefore two alleles in one molecule, by definition, so "which codon does this
//! sample carry?" is answered by a set lookup over sample indices. No phasing, no reads,
//! no BAM. Where the input has no genotypes there is nothing to look up, and every
//! allele is scored alone exactly as before ([`solo_call`]).
//!
//! # What is counted
//!
//! For a codon carrying variants at several positions, the codons actually **observed**
//! across the cohort are enumerated ([`observed_codons`]): a sample carrying variants at
//! positions 1 and 3 has that doubly-mutated codon; a sample carrying only one of them
//! has the singly-mutated one. Each observed codon is then scored against the reference
//! with Nei-Gojobori pathway averaging ([`pathway_diffs`]), which is exactly how
//! [`crate::models::nei`] handles a multi-difference codon pair in the FASTA path:
//! average over the orders in which the changes could have happened, skipping any order
//! whose intermediate codon is a stop.
//!
//! The invariant that keeps the counts commensurate with the per-nucleotide site
//! denominator in [`crate::vcf_analysis::count_sites`] is `Nd + Sd = k`, where `k` is the
//! number of nucleotide differences: a two-nucleotide change contributes **two**, not
//! one. Pathway averaging gives that for free (each of the `k` steps of a pathway is
//! classified, so every pathway scores `k`), and the per-codon total is split evenly over
//! the `k` alleles that make it up, so one allele still contributes exactly one unit and
//! a gene's total is still its number of classified alleles.
//!
//! Where one allele is carried on more than one background (some samples carry it
//! alone, others carry it with a neighbour), its share is averaged over those observed
//! codons, weighted by how many samples carry each. That is what "the codon each sample
//! actually carries" means when the samples disagree.

use super::*;
use std::collections::BTreeMap;

/// One ALT allele of a codon, as the per-codon scorer needs it.
pub(crate) struct CodonAllele<'a> {
    /// Which base of the codon this ALT changes (0, 1 or 2), on the **coding** strand.
    pub pos_in_codon: usize,
    /// The alternate base on the coding strand (already complemented for a minus gene).
    pub alt_in_cds: u8,
    /// Samples carrying this ALT, when the input had genotypes. `None` anywhere in the
    /// codon drops the whole codon back to the reference-codon scoring, which is the
    /// only thing an AF-only VCF supports.
    pub carriers: Option<&'a CarrierSet>,
}

/// How one ALT allele was classified, once the codon its carriers actually hold is known.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AlleleCall {
    /// The alternate codon this allele's carriers hold, on the coding strand. Equal to
    /// the single-mutant codon unless a sample carries a neighbouring SNP too; when the
    /// carriers disagree, this is the codon the most of them hold.
    pub alt_codon: [u8; 3],
    /// Amino acid encoded by `alt_codon`, and the effect against the reference codon.
    pub alt_aa: u8,
    pub effect: SnpEffect,
    /// Nucleotide differences between `alt_codon` and the reference codon: 1 for an
    /// ordinary SNP, 2 or 3 when most of this allele's carriers hold a multi-nucleotide
    /// change.
    pub diffs: u8,
    /// Does **any** sample carrying this ALT hold it inside a multi-nucleotide codon
    /// change? Wider than `diffs > 1`, which speaks only for the majority: an allele
    /// carried by 300 samples of whom 50 also carry a neighbour has `diffs == 1` and
    /// `joint == true`, and its `syn`/`nonsyn` split is a blend of the two backgrounds
    /// rather than a whole 1 in one class. That is the distinction the McDonald-Kreitman
    /// 2x2 and the SFS bins turn on, since both need whole alleles: they take `!joint`.
    pub joint: bool,
    /// This allele's share of its codon's Nei-Gojobori pathway score. `syn + nonsyn == 1`
    /// when the allele is counted, `0` when it is not (a change to or from a stop, which
    /// the pN/pS counts have always excluded).
    pub syn: f64,
    pub nonsyn: f64,
    /// Does this allele enter the pN/pS numerator at all?
    pub counted: bool,
}

impl AlleleCall {
    /// A call carrying no weight: the allele is shown in the variants table but excluded
    /// from every count, which is what a nonsense or stop-loss change has always been.
    fn excluded(alt_codon: [u8; 3], alt_aa: u8, effect: SnpEffect, diffs: u8, joint: bool) -> Self {
        AlleleCall {
            alt_codon,
            alt_aa,
            effect,
            diffs,
            joint,
            syn: 0.0,
            nonsyn: 0.0,
            counted: false,
        }
    }
}

/// Classify the coding effect of `alt_aa` against `ref_aa`, the one place the four
/// effect labels are decided.
fn effect_of(ref_aa: u8, alt_aa: u8) -> SnpEffect {
    if ref_aa == alt_aa {
        SnpEffect::Synonymous
    } else if alt_aa == b'*' {
        SnpEffect::Nonsense
    } else if ref_aa == b'*' {
        SnpEffect::StopLoss
    } else {
        SnpEffect::Missense
    }
}

/// Nucleotide differences between two codons.
fn hamming(a: &[u8; 3], b: &[u8; 3]) -> u8 {
    (0..3).filter(|&i| a[i] != b[i]).count() as u8
}

/// Score one ALT allele on its own against the reference codon: the classification
/// eskaks has always applied, and the only one an input without per-sample genotypes
/// supports. Returns `None` for an ambiguous alternate codon, which the caller drops.
fn solo_call(ref_codon: &[u8; 3], ref_aa: u8, a: &CodonAllele, gc: &GeneticCode) -> Option<AlleleCall> {
    let mut alt_codon = *ref_codon;
    alt_codon[a.pos_in_codon] = a.alt_in_cds;
    let alt_aa = codon_to_aa(&alt_codon, gc)?;
    let effect = effect_of(ref_aa, alt_aa);
    let diffs = hamming(ref_codon, &alt_codon);
    // A change to or from a stop is excluded from the pN/pS counts (the Nei-Gojobori
    // exclude-nonsense convention this crate uses throughout), but still reported.
    if ref_aa == b'*' || alt_aa == b'*' {
        return Some(AlleleCall::excluded(alt_codon, alt_aa, effect, diffs, false));
    }
    let syn = if ref_aa == alt_aa { 1.0 } else { 0.0 };
    Some(AlleleCall {
        alt_codon,
        alt_aa,
        effect,
        diffs,
        joint: false,
        syn,
        nonsyn: 1.0 - syn,
        counted: true,
    })
}

/// The distinct alternate codons the cohort actually carries at one codon, each with the
/// number of samples carrying it, ascending by codon so the result is deterministic.
///
/// A sample's codon is the reference codon with every position it carries a variant at
/// replaced. Samples carrying nothing here hold the reference and are not listed. Where
/// two records describe the same position (a VCF may split one site), the first listed
/// allele that lists the sample wins, so a duplicate record cannot silently overwrite
/// the base a sample was already given.
///
/// Iteration is over the union of the carrier sets, so this costs the number of samples
/// carrying *something* at this codon, not the cohort size.
pub(crate) fn observed_codons(ref_codon: &[u8; 3], alleles: &[CodonAllele]) -> Vec<([u8; 3], usize)> {
    let mut union = CarrierSet::default();
    for a in alleles {
        if let Some(cs) = a.carriers {
            union.union_with(cs);
        }
    }
    let mut counts: BTreeMap<[u8; 3], usize> = BTreeMap::new();
    for sample in union.samples() {
        let mut codon = *ref_codon;
        let mut filled = [false; 3];
        for a in alleles {
            if filled[a.pos_in_codon] {
                continue;
            }
            if a.carriers.is_some_and(|cs| cs.contains(sample)) {
                codon[a.pos_in_codon] = a.alt_in_cds;
                filled[a.pos_in_codon] = true;
            }
        }
        if codon != *ref_codon {
            *counts.entry(codon).or_insert(0) += 1;
        }
    }
    counts.into_iter().collect()
}

/// Classify every ALT allele of one codon by the codon its carriers actually hold.
///
/// Returns one entry per input allele, in the same order; `None` where the alternate
/// codon is ambiguous and the allele is dropped, exactly as the per-SNP path dropped it.
/// `ref_codon` must translate (the caller checks), and its amino acid is `ref_aa`.
///
/// Two fast paths return the historic per-SNP classification bit for bit, and between
/// them they cover every input that is not a genuine multi-nucleotide change:
///
/// * **no carriers** anywhere in the codon: an AF-only VCF has nothing to intersect, so
///   the reference-codon scoring is all there is (and stays the documented behaviour);
/// * **a single variant position**: alternative alleles of one site exclude one another
///   in a haploid genome, so no sample can hold two of them and no MNV is possible. This
///   is also what keeps a genome-scale run affordable: the carrier walk below runs only
///   at the rare codon carrying variants at two or three of its bases.
pub(crate) fn call_codon(
    ref_codon: &[u8; 3],
    ref_aa: u8,
    alleles: &[CodonAllele],
    gc: &GeneticCode,
) -> Vec<Option<AlleleCall>> {
    let solo = |a: &CodonAllele| solo_call(ref_codon, ref_aa, a, gc);
    let single_position = alleles.iter().all(|a| a.pos_in_codon == alleles[0].pos_in_codon);
    let all_carried = alleles.iter().all(|a| a.carriers.is_some());
    if alleles.is_empty() || single_position || !all_carried {
        return alleles.iter().map(solo).collect();
    }

    let observed = observed_codons(ref_codon, alleles);
    alleles
        .iter()
        .map(|a| {
            // The observed codons that contain THIS allele, with their carrier counts.
            let contexts: Vec<&([u8; 3], usize)> =
                observed.iter().filter(|(z, _)| z[a.pos_in_codon] == a.alt_in_cds).collect();
            if contexts.is_empty() {
                // Nobody carries it (an ALT listed with an all-reference genotype column,
                // or one shadowed by a duplicate record at the same position), so there is
                // no observed codon to score and the reference codon is the only context.
                return solo(a);
            }
            // The codon most of its carriers hold, ties broken toward the lower codon so
            // the reported row never depends on iteration order.
            let (dom_codon, _) = contexts
                .iter()
                .copied()
                .reduce(|best, c| if c.1 > best.1 { c } else { best })
                .expect("contexts is non-empty");
            let alt_aa = codon_to_aa(dom_codon, gc)?;
            let effect = effect_of(ref_aa, alt_aa);
            let diffs = hamming(ref_codon, dom_codon);
            // Any background at all, not just the majority one: an allele carried on both
            // a solo and a joint background contributes a blended split, which the
            // McDonald-Kreitman 2x2 and the SFS cannot take.
            let joint = contexts.iter().any(|(z, _)| hamming(ref_codon, z) > 1);
            // That dominant codon decides whether the allele is counted at all: a change to
            // or from a stop is excluded from the pN/pS counts, as it always has been, and
            // "this allele produces a stop" has to be one answer rather than a fraction,
            // because the row's `Effect` is one label. So an allele whose carriers mostly
            // hold a truncated protein leaves the counts even if a minority of them hold a
            // coding codon, and the reverse case keeps the allele and drops only the stop
            // background from the average below.
            if ref_aa == b'*' || alt_aa == b'*' {
                return Some(AlleleCall::excluded(*dom_codon, alt_aa, effect, diffs, joint));
            }
            // Average this allele's share over the backgrounds it is carried on, weighted
            // by how many samples carry each. Each observed codon's pathway score sums to
            // its own k, so dividing by k makes every allele contribute exactly one unit
            // and keeps Nd + Sd equal to the number of nucleotide differences.
            let (mut syn, mut nonsyn, mut weight) = (0.0f64, 0.0f64, 0.0f64);
            for (z, n) in contexts {
                let Some(z_aa) = codon_to_aa(z, gc) else { continue };
                // A background whose joint codon is a stop is excluded from the counts,
                // like any other change to a stop, and its weight goes with it.
                if z_aa == b'*' {
                    continue;
                }
                let k = hamming(ref_codon, z) as f64;
                if k == 0.0 {
                    continue;
                }
                let (sd, nd) = pathway_diffs(ref_codon, z, gc);
                syn += *n as f64 * sd / k;
                nonsyn += *n as f64 * nd / k;
                weight += *n as f64;
            }
            if weight <= 0.0 {
                return Some(AlleleCall::excluded(*dom_codon, alt_aa, effect, diffs, joint));
            }
            Some(AlleleCall {
                alt_codon: *dom_codon,
                alt_aa,
                effect,
                diffs,
                joint,
                syn: syn / weight,
                nonsyn: nonsyn / weight,
                counted: true,
            })
        })
        .collect()
}

/// Nei-Gojobori (1986) pathway analysis for one pair of codons: `(synonymous,
/// nonsynonymous)` differences, summing to the number of differing positions.
///
/// This is the same convention as [`crate::models::nei`]'s precomputed table, deliberately
/// so, because one crate must not hold two answers for "how many synonymous differences separate
/// these codons". For codons differing at one position the classification is direct; at
/// two positions it averages the 2 orders, at three the 6, and any order whose
/// intermediate codon is a **stop** is skipped. When every order is skipped the same
/// fallbacks are used as there: `(0.5, 1.5)` for a two-difference pair and `(1.0, 2.0)`
/// for a three-difference one. `cov_pathway_diffs_match_the_fasta_path` checks the two
/// agree over all 4,096 codon pairs of several genetic codes, so they cannot drift.
pub(crate) fn pathway_diffs(a: &[u8; 3], b: &[u8; 3], gc: &GeneticCode) -> (f64, f64) {
    let diffs: Vec<usize> = (0..3).filter(|&i| a[i] != b[i]).collect();
    let (Some(aa_a), Some(aa_b)) = (codon_to_aa(a, gc), codon_to_aa(b, gc)) else {
        return (0.0, 0.0);
    };
    match diffs.len() {
        0 => (0.0, 0.0),
        1 => {
            if aa_a == aa_b {
                (1.0, 0.0)
            } else {
                (0.0, 1.0)
            }
        }
        2 => {
            let orders = [[diffs[0], diffs[1]], [diffs[1], diffs[0]]];
            average_paths(a, b, aa_a, aa_b, &orders, gc).unwrap_or((0.5, 1.5))
        }
        _ => {
            let orders = [
                [0usize, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            average_paths(a, b, aa_a, aa_b, &orders, gc).unwrap_or((1.0, 2.0))
        }
    }
}

/// Average `(syn, nonsyn)` over the mutational orders in `orders`, skipping any order
/// that passes through a stop codon. `None` when every order does, which is the caller's
/// cue to use the fallback.
///
/// Each order lists the positions to change in sequence; the final step needs no
/// intermediate check, so only the first `len - 1` are examined, matching
/// [`crate::models::nei`]'s `pathway_2diff` / `pathway_3diff`.
fn average_paths<const N: usize>(
    a: &[u8; 3],
    b: &[u8; 3],
    aa_a: u8,
    aa_b: u8,
    orders: &[[usize; N]],
    gc: &GeneticCode,
) -> Option<(f64, f64)> {
    let (mut total_syn, mut total_nonsyn, mut valid) = (0.0f64, 0.0f64, 0u32);
    'orders: for order in orders {
        let mut cur = *a;
        let mut prev_aa = aa_a;
        let (mut syn, mut nonsyn) = (0.0f64, 0.0f64);
        for (step, &pos) in order.iter().enumerate() {
            cur[pos] = b[pos];
            let aa = if step + 1 == order.len() {
                aa_b
            } else {
                match codon_to_aa(&cur, gc) {
                    Some(aa) if aa != b'*' => aa,
                    // An intermediate stop (or an untranslatable intermediate) makes this
                    // order impossible; drop the whole order, not just the step.
                    _ => continue 'orders,
                }
            };
            if aa == prev_aa {
                syn += 1.0;
            } else {
                nonsyn += 1.0;
            }
            prev_aa = aa;
        }
        total_syn += syn;
        total_nonsyn += nonsyn;
        valid += 1;
    }
    (valid > 0).then(|| (total_syn / valid as f64, total_nonsyn / valid as f64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::nei::NeiTables;

    fn gc() -> &'static GeneticCode {
        crate::genetic_code::get_table(1).expect("standard code")
    }

    /// A carrier set over `n` haploid samples.
    fn carriers(samples: &[usize], n: usize) -> CarrierSet {
        let mut set = CarrierSet::new(n);
        for &s in samples {
            set.insert(s);
        }
        set
    }

    fn allele<'a>(pos_in_codon: usize, alt: u8, cs: Option<&'a CarrierSet>) -> CodonAllele<'a> {
        CodonAllele { pos_in_codon, alt_in_cds: alt, carriers: cs }
    }

    #[test]
    fn the_headline_case_is_scored_as_the_codon_samples_carry() {
        // Reference codon CTT (Leu) with C>T at base 1 and T>A at base 3, both carried by
        // every sample. Scored apart, that is a missense L->F plus a synonymous L->L. The
        // codon every sample really holds is TTA, which is Leu: no amino acid changes.
        let all = carriers(&(0..20).collect::<Vec<_>>(), 20);
        let alleles = [allele(0, b'T', Some(&all)), allele(2, b'A', Some(&all))];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());
        for c in &calls {
            let c = c.as_ref().expect("both alleles classify");
            assert_eq!(&c.alt_codon, b"TTA", "the codon the cohort carries");
            assert_eq!(c.alt_aa, b'L');
            assert_eq!(c.effect, SnpEffect::Synonymous, "Leu -> Leu encodes no change");
            assert_eq!(c.diffs, 2, "a two-nucleotide change");
            assert!(c.counted);
        }
        // Pathway averaging over the two orders: CTT->TTT->TTA is (nonsyn, nonsyn) and
        // CTT->CTA->TTA is (syn, syn), so the pair averages to (1, 1) over 2 differences,
        // and each allele carries half of it. Nd + Sd is still 2, not 1.
        let total_syn: f64 = calls.iter().flatten().map(|c| c.syn).sum();
        let total_nonsyn: f64 = calls.iter().flatten().map(|c| c.nonsyn).sum();
        assert!((total_syn - 1.0).abs() < 1e-12, "syn {}", total_syn);
        assert!((total_nonsyn - 1.0).abs() < 1e-12, "nonsyn {}", total_nonsyn);
        assert!((total_syn + total_nonsyn - 2.0).abs() < 1e-12, "Nd + Sd must equal k = 2");
    }

    #[test]
    fn a_joint_stop_is_nonsense_even_though_neither_snp_is() {
        // Reference CAT (His). Alone, C>T at base 1 gives TAT (Tyr) and T>A at base 3
        // gives CAA (Gln): two ordinary missense changes, both counted. Carried together
        // they give TAA, a premature STOP that neither change creates on its own. This is
        // the case the bug report singles out: reported as a missense plus a missense
        // where the genome carries a truncation.
        let all = carriers(&(0..10).collect::<Vec<_>>(), 10);
        let alleles = [allele(0, b'T', Some(&all)), allele(2, b'A', Some(&all))];
        let singles = call_codon(b"CAT", b'H', &[allele(0, b'T', None), allele(2, b'A', None)], gc());
        for c in singles.iter().flatten() {
            assert_eq!(c.effect, SnpEffect::Missense, "each SNP alone is a missense");
            assert!(c.counted, "and each alone would be counted");
        }
        let calls = call_codon(b"CAT", b'H', &alleles, gc());
        for c in calls.iter().flatten() {
            assert_eq!(&c.alt_codon, b"TAA");
            assert_eq!(c.alt_aa, b'*');
            assert_eq!(c.effect, SnpEffect::Nonsense, "the codon carried is a premature stop");
            assert!(!c.counted, "a change to a stop is excluded from the pN/pS counts");
            assert_eq!((c.syn, c.nonsyn), (0.0, 0.0));
        }
    }

    #[test]
    fn a_stop_background_decides_by_who_carries_it() {
        // The live shape of the toy genome's one joint stop: CAT (His) with C>T at base 1
        // carried by 10 samples, 6 of whom also carry T>A at base 3. Those 6 hold TAA, a
        // stop; the other 4 hold TAT (Tyr). The rule has to give ONE answer, because the
        // row carries one `Effect` label, and it follows the codon most of the carriers
        // hold: here the stop wins, so the allele is reported as nonsense and leaves the
        // counts entirely, exactly as a plain nonsense SNP does.
        let first = carriers(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 20);
        let second = carriers(&[0, 1, 2, 3, 4, 5], 20);
        let alleles = [allele(0, b'T', Some(&first)), allele(2, b'A', Some(&second))];
        let calls = call_codon(b"CAT", b'H', &alleles, gc());
        let a0 = calls[0].as_ref().expect("classified");
        assert_eq!(&a0.alt_codon, b"TAA");
        assert_eq!(a0.effect, SnpEffect::Nonsense);
        assert!(!a0.counted, "6 of its 10 carriers hold a truncation");
        assert_eq!((a0.syn, a0.nonsyn), (0.0, 0.0));

        // Tip the majority the other way and the allele comes back: 4 carriers on the
        // stop, 6 on TAT, so it is a missense and only the stop background is dropped from
        // the average (the surviving background is a plain 1-base change, so 1 nonsyn).
        let second = carriers(&[0, 1, 2, 3], 20);
        let alleles = [allele(0, b'T', Some(&first)), allele(2, b'A', Some(&second))];
        let calls = call_codon(b"CAT", b'H', &alleles, gc());
        let a0 = calls[0].as_ref().expect("classified");
        assert_eq!(&a0.alt_codon, b"TAT");
        assert_eq!(a0.effect, SnpEffect::Missense);
        assert!(a0.counted);
        assert!((a0.nonsyn - 1.0).abs() < 1e-12 && a0.syn == 0.0, "{:?}", a0);
        // The T>A allele is on the stop background for all of its carriers, so it stays out
        // either way: no codon it produces is countable.
        let a1 = calls[1].as_ref().expect("classified");
        assert_eq!(a1.effect, SnpEffect::Nonsense);
        assert!(!a1.counted);
    }

    #[test]
    fn disjoint_carriers_in_one_codon_are_two_ordinary_snps() {
        // Two SNPs in one codon that NO sample carries together are two single-nucleotide
        // changes, and must be scored exactly as they always were.
        let first = carriers(&[0, 1, 2], 20);
        let second = carriers(&[10, 11], 20);
        let alleles = [allele(0, b'T', Some(&first)), allele(2, b'A', Some(&second))];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());
        let solo = call_codon(b"CTT", b'L', &[allele(0, b'T', None), allele(2, b'A', None)], gc());
        assert_eq!(calls, solo, "no shared sample means nothing changes");
        assert_eq!(calls[0].as_ref().unwrap().effect, SnpEffect::Missense, "L2F stands");
        assert_eq!(calls[1].as_ref().unwrap().effect, SnpEffect::Synonymous);
    }

    #[test]
    fn an_allele_on_two_backgrounds_is_averaged_over_its_carriers() {
        // Three samples carry C>T alone (codon TTT, Phe: missense) and one carries it with
        // T>A (codon TTA, Leu: no change). The row reports the codon most of them hold,
        // and the count is the carrier-weighted average of the two backgrounds.
        let first = carriers(&[0, 1, 2, 3], 20);
        let second = carriers(&[3], 20);
        let alleles = [allele(0, b'T', Some(&first)), allele(2, b'A', Some(&second))];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());

        let a0 = calls[0].as_ref().expect("classified");
        assert_eq!(&a0.alt_codon, b"TTT", "3 of its 4 carriers hold TTT");
        assert_eq!(a0.effect, SnpEffect::Missense);
        assert_eq!(a0.diffs, 1);
        // The distinction the McDonald-Kreitman table turns on: this allele is mostly a
        // plain SNP, so `diffs` is 1, but one of its carriers holds it jointly and its
        // count below is a blend, so it is not a whole allele in either class.
        assert!(a0.joint, "diffs speaks for the majority; joint speaks for anyone");
        // 3 carriers on TTT (nonsynonymous, k = 1) and 1 on TTA (half of the (1,1)
        // pathway score, k = 2): nonsyn = (3*1 + 1*0.5)/4 = 0.875.
        assert!((a0.nonsyn - 0.875).abs() < 1e-12, "nonsyn {}", a0.nonsyn);
        assert!((a0.syn - 0.125).abs() < 1e-12, "syn {}", a0.syn);
        assert!((a0.syn + a0.nonsyn - 1.0).abs() < 1e-12, "an allele is always one unit");

        // The T>A allele has a single background, TTA, so it takes half of (1, 1).
        let a1 = calls[1].as_ref().expect("classified");
        assert_eq!(&a1.alt_codon, b"TTA");
        assert_eq!(a1.diffs, 2);
        assert!(a1.joint);
        assert!((a1.syn - 0.5).abs() < 1e-12 && (a1.nonsyn - 0.5).abs() < 1e-12);

        // Neither allele is a whole allele in one class, so neither can go into a
        // McDonald-Kreitman cell or an SFS bin.
        assert!(calls.iter().flatten().all(|c| c.joint));
    }

    #[test]
    fn three_snps_in_one_codon_use_the_six_pathway_average() {
        // Every base of GGG (Gly) mutated to CCC (Pro) in the same sample. All 6 orders
        // are stop-free and each has exactly one synonymous step, so the pair scores
        // (1, 2) over 3 differences, a third of it per allele.
        let all = carriers(&[0, 1], 4);
        let alleles = [
            allele(0, b'C', Some(&all)),
            allele(1, b'C', Some(&all)),
            allele(2, b'C', Some(&all)),
        ];
        let calls = call_codon(b"GGG", b'G', &alleles, gc());
        let syn: f64 = calls.iter().flatten().map(|c| c.syn).sum();
        let nonsyn: f64 = calls.iter().flatten().map(|c| c.nonsyn).sum();
        assert!((syn - 1.0).abs() < 1e-12, "syn {syn}");
        assert!((nonsyn - 2.0).abs() < 1e-12, "nonsyn {nonsyn}");
        assert!((syn + nonsyn - 3.0).abs() < 1e-12, "Nd + Sd must equal k = 3");
        for c in calls.iter().flatten() {
            assert_eq!(c.diffs, 3);
            assert_eq!(&c.alt_codon, b"CCC");
            assert!((c.syn + c.nonsyn - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn an_input_without_genotypes_keeps_the_reference_codon_scoring() {
        // No carriers anywhere: there is nothing to intersect, so every allele is scored
        // against the reference codon. This is the AF-only path and it must not move.
        let alleles = [allele(0, b'T', None), allele(2, b'A', None)];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());
        assert_eq!(calls[0].as_ref().unwrap().effect, SnpEffect::Missense);
        assert_eq!(calls[1].as_ref().unwrap().effect, SnpEffect::Synonymous);
        assert!(calls.iter().flatten().all(|c| c.diffs == 1));

        // One position missing its carriers is enough: a codon cannot be half-checked.
        let one = carriers(&[0], 4);
        let mixed = [allele(0, b'T', Some(&one)), allele(2, b'A', None)];
        assert_eq!(call_codon(b"CTT", b'L', &mixed, gc()), calls);
    }

    #[test]
    fn two_alts_at_one_position_never_form_a_haplotype() {
        // A multiallelic site is ONE position: a haploid sample carries at most one of its
        // ALTs, so however the carriers fall the codon can only be singly mutated.
        let a = carriers(&[0, 1], 4);
        let b = carriers(&[2, 3], 4);
        let alleles = [allele(0, b'T', Some(&a)), allele(0, b'A', Some(&b))];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());
        assert!(calls.iter().flatten().all(|c| c.diffs == 1));
        assert_eq!(&calls[0].as_ref().unwrap().alt_codon, b"TTT");
        assert_eq!(&calls[1].as_ref().unwrap().alt_codon, b"ATT");
    }

    #[test]
    fn an_allele_no_sample_carries_falls_back_to_the_reference_codon() {
        // An ALT listed with an all-reference genotype column has no observed codon at
        // all. Dropping it would lose the row; scoring it alone is the honest fallback.
        let none = carriers(&[], 20);
        let other = carriers(&[5], 20);
        let alleles = [allele(0, b'T', Some(&none)), allele(2, b'A', Some(&other))];
        let calls = call_codon(b"CTT", b'L', &alleles, gc());
        assert_eq!(&calls[0].as_ref().unwrap().alt_codon, b"TTT");
        assert_eq!(calls[0].as_ref().unwrap().effect, SnpEffect::Missense);
        assert_eq!(&calls[1].as_ref().unwrap().alt_codon, b"CTA");
    }

    #[test]
    fn an_ambiguous_alternate_codon_is_dropped_as_it_always_was() {
        let all = carriers(&[0], 4);
        let alleles = [allele(0, b'N', Some(&all))];
        assert!(call_codon(b"CTT", b'L', &alleles, gc())[0].is_none());
    }

    #[test]
    fn observed_codons_are_the_codons_samples_hold_and_nothing_else() {
        let first = carriers(&[0, 1, 2], 8);
        let second = carriers(&[2, 3], 8);
        let alleles = [allele(0, b'T', Some(&first)), allele(2, b'A', Some(&second))];
        // 0 and 1 hold TTT, 2 holds TTA, 3 holds CTA; 4..8 hold the reference and are
        // not listed at all.
        assert_eq!(
            observed_codons(b"CTT", &alleles),
            vec![(*b"CTA", 1), (*b"TTA", 1), (*b"TTT", 2)]
        );
    }

    #[test]
    fn one_position_split_across_records_does_not_overwrite_a_sample() {
        // Two records describe chr:4, listing different ALTs for the same sample (a
        // malformed or contradictory input). The first listed allele wins, so the codon a
        // sample gets never depends on record order beyond that stated rule.
        let a = carriers(&[0], 4);
        let b = carriers(&[0], 4);
        let alleles = [allele(0, b'T', Some(&a)), allele(0, b'A', Some(&b))];
        assert_eq!(observed_codons(b"CTT", &alleles), vec![(*b"TTT", 1)]);
    }

    #[test]
    fn pathway_diffs_matches_the_hand_derived_nei_gojobori_values() {
        // The values nei.rs's own tests derive by hand, recomputed by this scorer.
        assert_eq!(pathway_diffs(b"GCT", b"GCC", gc()), (1.0, 0.0)); // Ala -> Ala
        assert_eq!(pathway_diffs(b"GCT", b"ATT", gc()), (0.0, 2.0)); // Ala -> Ile, both paths nonsyn
        assert_eq!(pathway_diffs(b"CCC", b"GGG", gc()), (1.0, 2.0)); // Pro -> Gly, 6 paths
        // TGG (Trp) vs TAA: both intermediates (TAG, TGA) are stops, so every order is
        // skipped and the documented fallback applies.
        assert_eq!(pathway_diffs(b"TGG", b"TAA", gc()), (0.5, 1.5));
        assert_eq!(pathway_diffs(b"GCT", b"GCT", gc()), (0.0, 0.0));
    }

    #[test]
    fn cov_pathway_diffs_match_the_fasta_path() {
        // One crate must not hold two answers for "how many synonymous differences
        // separate these codons". Compared over all 4,096 pairs, for the standard code and
        // for two whose extra stops make the all-stop pathway fallbacks reachable.
        //
        // Nei index = 16*remap(b2) + 4*remap(b1) + remap(b3), remap A->2 C->1 G->3 T->0,
        // so its inverse maps slot values back to bases through "TCAG".
        let inv = b"TCAG";
        for id in [1u8, 2, 11] {
            let code = crate::genetic_code::get_table(id).expect("table exists");
            let tables = NeiTables::with_genetic_code(code);
            for ia in 0..64usize {
                for ib in 0..64usize {
                    let un = |idx: usize| {
                        let (s0, s1, s2) = (idx / 16, (idx / 4) % 4, idx % 4);
                        [inv[s1], inv[s0], inv[s2]]
                    };
                    let (a, b) = (un(ia), un(ib));
                    let (syn, nonsyn) = pathway_diffs(&a, &b, code);
                    let (want_syn, want_nonsyn) = tables.diff_entry(ia, ib);
                    assert!(
                        (syn - want_syn).abs() < 1e-6 && (nonsyn - want_nonsyn).abs() < 1e-6,
                        "code {id}, {} vs {}: got ({syn}, {nonsyn}), nei table has \
                         ({want_syn}, {want_nonsyn})",
                        std::str::from_utf8(&a).expect("ascii"),
                        std::str::from_utf8(&b).expect("ascii"),
                    );
                }
            }
        }
    }

    #[test]
    fn every_codon_pair_splits_into_exactly_its_nucleotide_differences() {
        // The invariant the per-nucleotide site denominator depends on: Nd + Sd = k for
        // every pair of translatable, non-stop codons, under every supported code.
        let bases = b"ACGT";
        for (id, _name) in crate::genetic_code::list_tables() {
            let code = crate::genetic_code::get_table(id).expect("listed table exists");
            for i in 0..64usize {
                for j in 0..64usize {
                    let un = |x: usize| [bases[x / 16], bases[(x / 4) % 4], bases[x % 4]];
                    let (a, b) = (un(i), un(j));
                    if codon_to_aa(&a, code) == Some(b'*') || codon_to_aa(&b, code) == Some(b'*') {
                        continue;
                    }
                    let k = hamming(&a, &b) as f64;
                    let (syn, nonsyn) = pathway_diffs(&a, &b, code);
                    assert!(
                        (syn + nonsyn - k).abs() < 1e-9,
                        "code {id}: {} vs {} split into {} + {} != {k}",
                        std::str::from_utf8(&a).expect("ascii"),
                        std::str::from_utf8(&b).expect("ascii"),
                        syn,
                        nonsyn
                    );
                }
            }
        }
    }
}
