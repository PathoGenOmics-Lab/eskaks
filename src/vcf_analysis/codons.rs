//! Per-codon recurrence scan: which residues collect an unusual number of
//! *distinct* nonsynonymous alleles.
//!
//! # What is tested, and what is deliberately not
//!
//! The statistic is **allelic multiplicity**: `X_c`, the number of distinct missense
//! ALT alleles observed at codon `c`. Each distinct allele is at minimum one
//! independent mutational event, which is the only recurrence claim a tool with no
//! phylogeny can make. Under the null, a codon's `X_c` is Binomial(`A_c`, `theta`),
//! where `A_c` is the number of single-nucleotide changes that codon could undergo to
//! a different amino acid (an integer from 0 to 9, stop-creating changes excluded, see
//! [`codon_change_split`]) and `theta = Σ X_c / Σ A_c` is the plug-in genome-wide rate
//! per possible change. The reported p-value is the one-sided upper tail `P(X >= x_c)`.
//!
//! Two things this is **not**:
//!
//! * **Not a per-codon pN/pS test.** A codon has at most 9 possible SNP alleles ever,
//!   and its null nonsynonymous fraction `N_c/(N_c+S_c)` has median 0.667 and minimum
//!   0.500 over the 61 sense codons of code 11. Enumerating every sense codon and
//!   every achievable (k, n), the smallest two-sided mid-p reachable in the
//!   nonsynonymous-excess direction is **0.0529** (CTA, 5 alleles of 5): **no codon can
//!   ever reach even a nominal 0.05**, let alone survive any multiple-testing
//!   correction. ATG and TGG have `S_c == 0`, so they have no null at all (embB M306
//!   is ATG). Such a column would be significant nowhere, so none is emitted.
//! * **Not a test on carrier counts.** Carrier columns are reported and never tested:
//!   without a phylogeny one old mutation carried by an expanded lineage cannot be
//!   told from many independent origins. `X_c` is the quantity clonal expansion does
//!   not inflate, which is exactly why it is the tested one.
//!
//! `theta` is small (1e-3 to 4e-2 at realistic cohort sizes), unlike the codon's own
//! nonsynonymous fraction, so the tail is steep and the test has real power: at a
//! 1000-isolate scale it recovers gyrA 94 (p ~ 4e-13) and embB 306 (p ~ 8e-8)
//! comfortably, katG 315 marginally, and does NOT recover rpoB S450L or rpsL K43R,
//! whose evidence is carrier recurrence rather than allelic multiplicity.

use super::*;
use std::collections::BTreeMap;

/// One codon carrying at least one coding SNP: its observed distinct alleles, its
/// combinatorial opportunity, and the recurrence test built from the two.
#[derive(Debug, Clone)]
pub struct CodonRow {
    /// Parent gene, its contig and strand (copied from the parent [`GenePnPs`]).
    pub gene: String,
    pub chrom: String,
    pub strand: char,
    /// 1-based residue number within the protein, and the reference residue.
    pub aa_pos: usize,
    pub ref_aa: u8,
    /// The reference codon, on the coding strand.
    pub ref_codon: [u8; 3],
    /// `A_c` / `S_c`: possible nonsynonymous and synonymous single-nucleotide changes
    /// at this codon (stop-creating and ambiguous changes excluded).
    pub poss_nonsyn: u32,
    pub poss_syn: u32,
    /// Distinct ALT alleles observed here, by effect. `nonsyn_alleles` is the tested
    /// statistic `X_c`; the other two are reported only.
    pub nonsyn_alleles: u32,
    pub syn_alleles: u32,
    /// Distinct alleles that create a stop codon (nonsense) or replace a reference
    /// stop (stop-loss). Excluded from `X_c` and from `A_c`, exactly as `count_sites`
    /// excludes changes to a stop, but shown because a loss of function matters.
    pub nonsense_alleles: u32,
    /// Distinct alternate amino acids over all amino-acid-changing alleles, and the
    /// changes themselves (`S315I;S315N;S315T`, sorted). `distinct_aa` can be below
    /// `nonsyn_alleles` when two different nucleotide alleles encode the same residue.
    pub distinct_aa: u32,
    pub aa_changes: String,
    /// Carrier bounds over this codon's alleles: the largest single-allele carrier
    /// count (a lower bound on codon carriers) and their sum (an upper bound, since a
    /// sample carrying two of them is counted twice). `None` when the input has
    /// neither genotypes nor a sample size, in which case the test still runs.
    pub carriers_max: Option<usize>,
    pub carriers_sum: Option<usize>,
    /// Highest allele frequency among this codon's alleles.
    pub max_af: f64,
    /// `A_c · theta`: how many distinct nonsynonymous alleles this codon was expected
    /// to collect by chance. The effect-size anchor for the p-value.
    pub exp_nonsyn_alleles: f64,
    /// One-sided upper-tail p-value `P(X >= nonsyn_alleles)` and its BH q-value over
    /// the whole coding genome. Both NaN where the test is not meaningful.
    pub p_recurrence: f64,
    pub q_recurrence: f64,
    /// True when some sample carries two of this codon's SNPs: they are then a single
    /// multi-nucleotide event, not independent origins, so the codon is not tested and is
    /// dropped from the family. Read from the per-sample genotypes where the input has
    /// them, and from the allele-frequency pigeonhole bound where it does not.
    pub cooccurring: bool,
    /// Parent-gene context, so a globally elevated gene is not mistaken for a
    /// codon-specific hit. `gene_q_bh` is the gene's neutrality-test BH q-value.
    pub repetitive: bool,
    pub gene_pn_ps: f64,
    pub gene_q_bh: f64,
    /// The parent gene's 1-based start. Not written to the table; it disambiguates the
    /// join back to the per-gene results when two genes share a name.
    pub gene_start: usize,
}

/// A whole-run per-codon scan: the rows, plus the null it was tested against.
#[derive(Debug, Clone)]
pub struct CodonScan {
    /// One row per codon carrying at least one coding SNP, already sorted.
    pub rows: Vec<CodonRow>,
    /// Plug-in rate per possible nonsynonymous change, `Σ X_c / Σ A_c` over the family.
    pub theta: f64,
    /// Multiple-testing family size: every codon of every analysed gene with `A_c > 0`
    /// (minus repetitive genes under `--exclude-repetitive`, minus phase-suppressed
    /// codons), NOT the number of rows written.
    pub family_m: usize,
    /// Denominator and numerator of `theta`.
    pub family_poss_nonsyn: u64,
    pub observed_nonsyn_alleles: u64,
    /// Codons suppressed because their SNPs are provably one multi-nucleotide event.
    pub suppressed_cooccurring: usize,
    /// Effective sample size used for the carrier columns (0 = unknown).
    pub n_samples: usize,
}

impl CodonScan {
    /// Rows significant at the given FDR threshold.
    pub fn significant(&self, fdr: f64) -> usize {
        self.rows.iter().filter(|r| r.q_recurrence.is_finite() && r.q_recurrence < fdr).count()
    }

    /// Join each row's parent-gene BH q-value from the (possibly `--min-snps`-filtered)
    /// per-gene results, once the gene-level correction has run. A gene dropped from
    /// the per-gene table keeps `NA` here: it was never in the gene test's family.
    pub fn attach_gene_qvalues(&mut self, results: &[GenePnPs]) {
        let mut q: HashMap<(&str, &str, usize), f64> = HashMap::new();
        for g in results {
            q.insert((g.name.as_str(), g.chrom.as_str(), g.genome_start), g.q_value);
        }
        for r in &mut self.rows {
            r.gene_q_bh = q
                .get(&(r.gene.as_str(), r.chrom.as_str(), r.gene_start))
                .copied()
                .unwrap_or(f64::NAN);
        }
    }
}

/// Carrier count behind one ALT allele: the exact genotype-derived count when the VCF
/// has genotypes, else `round(af · n)` over the effective sample size, exactly as the
/// diversity path does. `None` when the input offers neither.
fn carriers_of(v: &Variant, n_samples: usize) -> Option<usize> {
    match v.gt_derived {
        Some(k) => Some(k),
        None if n_samples > 0 && v.af.is_finite() => {
            Some((v.af * n_samples as f64).round().max(0.0) as usize)
        }
        None => None,
    }
}

/// Build the per-codon recurrence scan from the per-gene results.
///
/// Call this on the **unfiltered** results, before `--min-snps` drops low-count genes:
/// the family is the whole coding genome, exactly as the genome-wide pooled pN/pS is
/// pooled over all genes. `n_samples` is the effective sample size for the carrier
/// columns (0 when the input has no genotypes and no sample count), and
/// `exclude_repetitive` removes PE/PPE/PGRS/IS genes from the family and from `theta`
/// (their SNP calls are frequently mapping artefacts, and a mutational or mapping
/// hotspot is exactly what this test cannot distinguish from selection).
pub fn compute_codon_scan(
    results: &[GenePnPs],
    gc: &GeneticCode,
    n_samples: usize,
    exclude_repetitive: bool,
) -> CodonScan {
    let in_family = |g: &GenePnPs| !(exclude_repetitive && g.repetitive);

    // The family before phase suppression: every codon of every analysed gene that
    // could show a nonsynonymous allele, whether or not a SNP landed on it.
    let mut family_m: i64 = 0;
    let mut family_poss: i64 = 0;
    for g in results {
        if in_family(g) {
            family_m += g.scan_codons as i64;
            family_poss += g.scan_poss_nonsyn as i64;
        }
    }

    let mut rows: Vec<CodonRow> = Vec::new();
    for g in results {
        if g.variants.is_empty() {
            continue;
        }
        // BTreeMap, so codons come out in residue order whatever the strand.
        let mut by_codon: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        for v in &g.variants {
            by_codon.entry(v.aa_pos).or_default().push(v);
        }

        for (aa_pos, mut vars) in by_codon {
            // Collapse duplicate (position, ALT) records to one allele. The merge path
            // already keys on (CHROM, POS, ALT), but a directly parsed VCF may split one
            // allele across two records, and counting it twice would manufacture a hit.
            // The descending-AF tiebreak keeps the record with the highest frequency,
            // matching how the same-codon phase check collapses repeated positions.
            vars.sort_by(|a, b| {
                a.pos
                    .cmp(&b.pos)
                    .then(a.alt_allele.cmp(&b.alt_allele))
                    .then(b.af.total_cmp(&a.af))
            });
            vars.dedup_by(|a, b| a.pos == b.pos && a.alt_allele == b.alt_allele);

            let ref_codon = vars[0].ref_codon;
            let ref_aa = vars[0].ref_aa;
            let (poss_nonsyn, poss_syn) = codon_change_split(&ref_codon, gc);

            let (mut nonsyn_alleles, mut syn_alleles, mut nonsense_alleles) = (0u32, 0u32, 0u32);
            let mut changes: Vec<String> = Vec::new();
            let mut alt_aas: Vec<u8> = Vec::new();
            let mut max_af = f64::NEG_INFINITY;
            let (mut carriers_known, mut c_max, mut c_sum) = (true, 0usize, 0usize);

            for v in &vars {
                match v.effect {
                    SnpEffect::Synonymous => syn_alleles += 1,
                    SnpEffect::Missense => nonsyn_alleles += 1,
                    SnpEffect::Nonsense | SnpEffect::StopLoss => nonsense_alleles += 1,
                }
                if v.effect != SnpEffect::Synonymous {
                    changes.push(format!("{}{}{}", ref_aa as char, aa_pos, v.alt_aa as char));
                    alt_aas.push(v.alt_aa);
                }
                if v.af.is_finite() {
                    max_af = max_af.max(v.af);
                }
                match carriers_of(v, n_samples) {
                    Some(c) => {
                        c_max = c_max.max(c);
                        c_sum += c;
                    }
                    None => carriers_known = false,
                }
            }
            changes.sort();
            changes.dedup();
            alt_aas.sort_unstable();
            alt_aas.dedup();

            // A codon whose SNPs one sample carries together is ONE multi-nucleotide event,
            // not several independent origins, so testing it would count a single mutation
            // several times. `Variant::codon_shared` is that verdict, already decided per
            // allele: read from the per-sample genotypes where they exist (so a codon whose
            // SNPs sit in disjoint samples is correctly NOT suppressed, which the
            // allele-frequency floor could never tell), and from the frequency bound alone
            // for an AF-only VCF, which is what this used to be everywhere.
            let cooccurring = vars.iter().any(|v| v.codon_shared == Some(true));

            rows.push(CodonRow {
                gene: g.name.clone(),
                chrom: g.chrom.clone(),
                strand: g.strand,
                aa_pos,
                ref_aa,
                ref_codon,
                poss_nonsyn,
                poss_syn,
                nonsyn_alleles,
                syn_alleles,
                nonsense_alleles,
                distinct_aa: alt_aas.len() as u32,
                aa_changes: if changes.is_empty() { String::new() } else { changes.join(";") },
                carriers_max: carriers_known.then_some(c_max),
                carriers_sum: carriers_known.then_some(c_sum),
                max_af: if max_af.is_finite() { max_af } else { f64::NAN },
                exp_nonsyn_alleles: f64::NAN,
                p_recurrence: f64::NAN,
                q_recurrence: f64::NAN,
                cooccurring,
                repetitive: g.repetitive,
                gene_pn_ps: g.pn_ps,
                gene_q_bh: f64::NAN,
                gene_start: g.genome_start,
            });
        }
    }

    // Phase-suppressed codons leave the family entirely: their alleles are one
    // multi-nucleotide event, so they are neither a test nor evidence about theta.
    let mut suppressed = 0usize;
    let mut observed: u64 = 0;
    for r in &rows {
        // A stop or ambiguous reference codon has A_c == 0 and was never counted into
        // the family by `codon_space`, so it needs no adjustment here.
        if r.poss_nonsyn == 0 || (exclude_repetitive && r.repetitive) {
            continue;
        }
        if r.cooccurring {
            suppressed += 1;
            family_m -= 1;
            family_poss -= r.poss_nonsyn as i64;
            continue;
        }
        observed += r.nonsyn_alleles as u64;
    }
    let family_m = family_m.max(0) as usize;
    let family_poss = family_poss.max(0) as u64;

    let theta = if family_poss > 0 { observed as f64 / family_poss as f64 } else { f64::NAN };

    for r in rows.iter_mut() {
        if theta.is_finite() && r.poss_nonsyn > 0 {
            r.exp_nonsyn_alleles = r.poss_nonsyn as f64 * theta;
        }
        // No null where the codon can undergo no nonsynonymous change, where the
        // alleles are provably one event, or where theta is unestimable.
        if r.cooccurring || r.poss_nonsyn == 0 || !theta.is_finite() {
            continue;
        }
        r.p_recurrence =
            crate::stats::binomial_upper_tail_p(r.nonsyn_alleles as u64, r.poss_nonsyn as u64, theta);
    }

    // BH over the whole coding genome. Codons outside the family contribute NaN, the
    // same out-of-family convention `apply_multiple_testing` uses for repetitive genes:
    // their p-value is still shown, their q-value is not.
    let q_in: Vec<f64> = rows
        .iter()
        .map(|r| if exclude_repetitive && r.repetitive { f64::NAN } else { r.p_recurrence })
        .collect();
    let qvals = crate::stats::benjamini_hochberg_with_m(&q_in, family_m);
    for (r, q) in rows.iter_mut().zip(qvals) {
        r.q_recurrence = q;
    }

    // Deterministic ranking: p ascending with NA last, then the raw statistic, then
    // the (untested) carrier bound, then gene and residue.
    rows.sort_by(|a, b| {
        let key = |p: f64| if p.is_nan() { f64::INFINITY } else { p };
        key(a.p_recurrence)
            .total_cmp(&key(b.p_recurrence))
            .then(b.nonsyn_alleles.cmp(&a.nonsyn_alleles))
            .then(b.carriers_max.unwrap_or(0).cmp(&a.carriers_max.unwrap_or(0)))
            .then(a.gene.cmp(&b.gene))
            .then(a.aa_pos.cmp(&b.aa_pos))
            .then(a.chrom.cmp(&b.chrom))
    });

    CodonScan {
        rows,
        theta,
        family_m,
        family_poss_nonsyn: family_poss,
        observed_nonsyn_alleles: observed,
        suppressed_cooccurring: suppressed,
        n_samples,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff::{CdsExon, Strand};
    use crate::vcf::VcfSnp;

    /// Reference CDS `TCG GCT TTT` = Ser Ala Phe, one plus-strand gene on chr1.
    ///
    /// Its codon-level opportunity, by hand from the standard table:
    ///   TCG (Ser): pos1 -> ACG/CCG/GCG (3 nonsyn), pos2 -> TAG stop (excluded) +
    ///              TGG/TTG (2 nonsyn), pos3 -> TCA/TCC/TCT (3 syn). A = 5, S = 3.
    ///   GCT (Ala): third position is 4-fold, so A = 6, S = 3.
    ///   TTT (Phe): only TTC is synonymous, so A = 8, S = 1.
    /// Family: 3 codons, 5 + 6 + 8 = 19 possible nonsynonymous changes.
    fn ser_gene() -> (HashMap<String, Vec<u8>>, Gene) {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"TCGGCTTTT".to_vec());
        let gene = Gene {
            name: "g1".to_string(),
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

    /// Run the real pipeline, so the scan is exercised through the same `Variant`
    /// records the CLI builds rather than hand-made ones.
    fn scan_of(snps: &[VcfSnp], n_samples: usize) -> CodonScan {
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = ser_gene();
        let (results, _diag) = compute_pn_ps(&reference, &[gene], snps, gc, false, 1.0, 0.95);
        assert_eq!(results.len(), 1, "the toy gene must be analysed");
        assert_eq!(results[0].scan_codons, 3, "3 codons with A_c > 0");
        assert_eq!(results[0].scan_poss_nonsyn, 19, "5 + 6 + 8 possible nonsyn changes");
        compute_codon_scan(&results, gc, n_samples, false)
    }

    #[test]
    fn a_codon_mixing_synonymous_and_nonsynonymous_snps_is_scored_by_hand() {
        // Codon 1 (TCG) collects three distinct alleles:
        //   chr1:1 T>G -> GCG  Ser->Ala  missense   (S1A)
        //   chr1:2 C>G -> TGG  Ser->Trp  missense   (S1W)
        //   chr1:3 G>A -> TCA  Ser->Ser  synonymous (not counted in X_c)
        // Codon 2 (GCT) collects one:
        //   chr1:4 G>A -> ACT  Ala->Thr  missense   (A2T)
        // Frequencies stay low, so the phase bound never fires (0.3+0.2+0.2 <= 2).
        let snps = vec![
            snp(1, b'T', b'G', 0.3),
            snp(2, b'C', b'G', 0.2),
            snp(3, b'G', b'A', 0.2),
            snp(4, b'G', b'A', 0.4),
        ];
        let scan = scan_of(&snps, 10);

        // theta counts DISTINCT NONSYNONYMOUS alleles only: 2 at codon 1 + 1 at codon 2,
        // over all 19 possible nonsynonymous changes in the gene. The synonymous allele
        // is reported and never enters the null.
        assert_eq!(scan.observed_nonsyn_alleles, 3);
        assert_eq!(scan.family_poss_nonsyn, 19);
        assert_eq!(scan.family_m, 3, "m is every codon, not the 2 rows written");
        assert!((scan.theta - 3.0 / 19.0).abs() < 1e-15, "theta {}", scan.theta);

        assert_eq!(scan.rows.len(), 2, "codon 3 carries no SNP, so it has no row");
        let c1 = &scan.rows[0];
        assert_eq!((c1.aa_pos, c1.ref_aa, &c1.ref_codon), (1, b'S', b"TCG"));
        assert_eq!((c1.poss_nonsyn, c1.poss_syn), (5, 3));
        assert_eq!((c1.nonsyn_alleles, c1.syn_alleles, c1.nonsense_alleles), (2, 1, 0));
        assert_eq!(c1.distinct_aa, 2);
        assert_eq!(c1.aa_changes, "S1A;S1W", "sorted, synonymous alleles excluded");
        // Carriers span every allele at the codon, synonymous included: the max is a
        // lower bound on the codon's carriers and the sum an upper bound.
        assert_eq!((c1.carriers_max, c1.carriers_sum), (Some(3), Some(7)));
        assert!((c1.max_af - 0.3).abs() < 1e-12);
        assert!(!c1.cooccurring);

        // Statistic: P(X >= 2) for X ~ Binomial(5, 3/19), whole atom, one-sided.
        //   q = 16/19; P(X=0) = q^5 = 0.4234787...,  P(X=1) = 5(3/19)q^4 = 0.3970113...
        //   P(X >= 2) = 1 - 0.8204900... = 0.1795099...
        let (t, q) = (3.0f64 / 19.0, 16.0f64 / 19.0);
        let want_c1 = 1.0 - q.powi(5) - 5.0 * t * q.powi(4);
        assert!((c1.p_recurrence - want_c1).abs() < 1e-12, "p {}", c1.p_recurrence);
        assert!((c1.p_recurrence - 0.179_509_9).abs() < 1e-6, "p {}", c1.p_recurrence);
        assert!((c1.exp_nonsyn_alleles - 5.0 * t).abs() < 1e-12, "A_c * theta");

        let c2 = &scan.rows[1];
        assert_eq!((c2.aa_pos, c2.nonsyn_alleles, c2.syn_alleles), (2, 1, 0));
        assert_eq!(c2.aa_changes, "A2T");
        let want_c2 = 1.0 - q.powi(6);
        assert!((c2.p_recurrence - want_c2).abs() < 1e-12, "p {}", c2.p_recurrence);

        // Multiple testing over m = 3, NOT over the 2 printed rows: the untested codon
        // sits at p = 1.0 and takes the last rank, so the step-up minima are
        //   rank 2: min(1, p2 * 3/2),  rank 1: min(that, p1 * 3/1).
        let q2 = (want_c2 * 3.0 / 2.0).min(1.0);
        let q1 = q2.min(want_c1 * 3.0);
        assert!((c1.q_recurrence - q1).abs() < 1e-12, "q1 {} want {}", c1.q_recurrence, q1);
        assert!((c2.q_recurrence - q2).abs() < 1e-12, "q2 {} want {}", c2.q_recurrence, q2);
        assert!((q1 - 0.538_528_1).abs() < 1e-6, "q1 drifted: {q1}");
        // Correcting over the 2 rows alone would have been visibly laxer.
        assert!(q1 > (want_c1 * 2.0).min(1.0), "m must be the genome, not the table");
        assert_eq!(scan.significant(0.05), 0, "nothing survives a 3-codon family");
    }

    #[test]
    fn a_single_snp_in_a_codon_is_a_real_but_unremarkable_row() {
        // The edge case the whole design turns on: one allele is one mutational event,
        // which is exactly what the null expects, so it must produce a row with a
        // finite, distinctly unimpressive p-value rather than a hit or an NA.
        let scan = scan_of(&[snp(4, b'G', b'A', 0.4)], 10);
        assert_eq!(scan.rows.len(), 1);
        assert_eq!(scan.observed_nonsyn_alleles, 1);
        assert!((scan.theta - 1.0 / 19.0).abs() < 1e-15);

        let r = &scan.rows[0];
        assert_eq!((r.aa_pos, r.nonsyn_alleles, r.syn_alleles), (2, 1, 0));
        assert_eq!((r.poss_nonsyn, r.poss_syn), (6, 3));
        assert_eq!(r.distinct_aa, 1);
        // P(X >= 1) = 1 - (18/19)^6 = 0.2770...; over m = 3 that is 0.831..., nowhere
        // near significance, which is the point.
        let want = 1.0 - (18.0f64 / 19.0).powi(6);
        assert!((r.p_recurrence - want).abs() < 1e-12, "p {}", r.p_recurrence);
        assert!((r.p_recurrence - 0.277_041_4).abs() < 1e-6, "p {}", r.p_recurrence);
        assert!((r.q_recurrence - (want * 3.0).min(1.0)).abs() < 1e-12);
        assert_eq!(scan.significant(0.05), 0);
    }

    #[test]
    fn a_synonymous_only_codon_gets_p_equal_to_one_not_na() {
        // X_c = 0 is a real observation, so P(X >= 0) = 1 exactly. NA here would say
        // "not tested", which would be a lie.
        let scan = scan_of(&[snp(3, b'G', b'A', 0.2), snp(4, b'G', b'A', 0.4)], 10);
        let syn_row = scan.rows.iter().find(|r| r.aa_pos == 1).expect("codon 1 row");
        assert_eq!((syn_row.nonsyn_alleles, syn_row.syn_alleles), (0, 1));
        assert_eq!(syn_row.p_recurrence, 1.0);
        assert_eq!(syn_row.aa_changes, "", "no amino-acid change to list");
        assert!(syn_row.q_recurrence.is_finite(), "still inside the family");
    }

    #[test]
    fn a_phase_forced_codon_is_suppressed_and_leaves_the_family() {
        // Two SNPs in codon 1 at AF 1.0: every sample carries both, so they are one
        // multi-nucleotide event, not two origins. Testing them would count a single
        // mutation twice, which is the false-positive hole this closes.
        let snps = vec![snp(1, b'T', b'G', 1.0), snp(2, b'C', b'G', 1.0)];
        let scan = scan_of(&snps, 10);
        assert_eq!(scan.suppressed_cooccurring, 1);

        let r = &scan.rows[0];
        assert!(r.cooccurring);
        assert_eq!(r.nonsyn_alleles, 2, "the alleles are still counted and shown");
        assert!(r.p_recurrence.is_nan() && r.q_recurrence.is_nan(), "not tested");
        // The codon leaves the family on both counts: m drops from 3 to 2, its 5
        // possible changes leave theta's denominator, and its 2 alleles leave the
        // numerator (which is then 0, since no other codon carries one).
        assert_eq!(scan.family_m, 2);
        assert_eq!(scan.family_poss_nonsyn, 19 - 5);
        assert_eq!(scan.observed_nonsyn_alleles, 0);
        assert_eq!(scan.theta, 0.0);
    }

    #[test]
    fn two_alts_at_one_position_are_two_independent_events() {
        // A multiallelic site is ONE position: its ALTs exclude one another, so they can
        // never be a haplotype, and each is a separate mutational origin. Both count.
        let multi = VcfSnp {
            chrom: "chr1".to_string(),
            pos: 1,
            ref_allele: b'T',
            alt_alleles: vec![b'G', b'A'],
            alt_freqs: vec![0.6, 0.4],
            gt_counts: None,
            carriers: None,
            filter: "PASS".to_string(),
            depth: None,
        };
        // T>G gives GCG (Ala), T>A gives ACG (Thr): two distinct missense alleles.
        let scan = scan_of(&[multi], 10);
        assert_eq!(scan.rows.len(), 1);
        let r = &scan.rows[0];
        assert!(!r.cooccurring, "one position can never be forced onto a haplotype");
        assert_eq!(r.nonsyn_alleles, 2);
        assert_eq!(r.aa_changes, "S1A;S1T");
        assert!(r.p_recurrence.is_finite());
    }

    #[test]
    fn one_allele_split_across_two_records_is_not_two_events() {
        // A directly parsed VCF can list the same (position, ALT) twice. Counting it
        // twice would manufacture recurrence out of a duplicated line, so the scan
        // collapses on (position, ALT) before counting.
        let dup = vec![snp(1, b'T', b'G', 0.3), snp(1, b'T', b'G', 0.3)];
        let scan = scan_of(&dup, 10);
        assert_eq!(scan.rows.len(), 1);
        assert_eq!(scan.rows[0].nonsyn_alleles, 1, "one allele, however often listed");
        assert_eq!(scan.observed_nonsyn_alleles, 1);
    }

    #[test]
    fn carrier_columns_are_na_without_genotypes_but_the_test_still_runs() {
        // The scan needs allele identity, never counts, so an AF-only VCF with no
        // sample size loses only the descriptive columns.
        let scan = scan_of(&[snp(4, b'G', b'A', 0.4)], 0);
        let r = &scan.rows[0];
        assert_eq!((r.carriers_max, r.carriers_sum), (None, None));
        assert!(r.p_recurrence.is_finite(), "the test is unaffected");
        assert_eq!(r.nonsyn_alleles, 1);
    }

    #[test]
    fn repetitive_genes_leave_the_family_under_exclude_repetitive() {
        // Same convention as the per-gene test: the p-value is still shown, the
        // q-value is not, and the codon contributes to neither m nor theta.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (mut reference, mut gene) = ser_gene();
        gene.name = "PPE_toy".to_string();
        reference.insert("chr1".to_string(), b"TCGGCTTTT".to_vec());
        let snps = vec![snp(4, b'G', b'A', 0.4)];
        let (results, _d) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);
        assert!(results[0].repetitive, "PPE_toy must be flagged repetitive");

        let kept = compute_codon_scan(&results, gc, 10, false);
        assert_eq!(kept.family_m, 3);
        assert!(kept.rows[0].q_recurrence.is_finite());

        let dropped = compute_codon_scan(&results, gc, 10, true);
        assert_eq!(dropped.family_m, 0, "the only gene left the family");
        assert_eq!(dropped.family_poss_nonsyn, 0);
        assert!(dropped.theta.is_nan(), "no family, no plug-in rate");
        assert!(dropped.rows[0].repetitive);
        assert!(dropped.rows[0].q_recurrence.is_nan(), "out of the correction family");
    }

    #[test]
    fn the_scan_is_indifferent_to_af_weighting_and_to_kappa() {
        // X_c is an integer count of alleles, so neither the AF weighting that turns
        // pN/pS into piN/piS nor the ts/tv site model can move it. That is the clearest
        // statement that this measures something the pN/pS columns do not.
        let gc = crate::genetic_code::get_table(1).expect("standard code");
        let (reference, gene) = ser_gene();
        let snps = vec![snp(1, b'T', b'G', 0.3), snp(4, b'G', b'A', 0.4)];
        let base = {
            let (r, _) = compute_pn_ps(&reference, std::slice::from_ref(&gene), &snps, gc, false, 1.0, 0.95);
            compute_codon_scan(&r, gc, 10, false)
        };
        for (af_weighted, kappa) in [(true, 1.0), (false, 4.0), (true, 4.0)] {
            let (r, _) =
                compute_pn_ps(&reference, std::slice::from_ref(&gene), &snps, gc, af_weighted, kappa, 0.95);
            let other = compute_codon_scan(&r, gc, 10, false);
            assert_eq!(other.family_poss_nonsyn, base.family_poss_nonsyn);
            assert_eq!(other.observed_nonsyn_alleles, base.observed_nonsyn_alleles);
            assert_eq!(other.rows.len(), base.rows.len());
            for (a, b) in other.rows.iter().zip(&base.rows) {
                assert_eq!((a.aa_pos, a.nonsyn_alleles, a.poss_nonsyn), (b.aa_pos, b.nonsyn_alleles, b.poss_nonsyn));
                assert!((a.p_recurrence - b.p_recurrence).abs() < 1e-15);
            }
        }
    }
}
