use proptest::prelude::*;

/// Compare two f64 values, treating NaN==NaN as equal.
fn approx_eq(a: f64, b: f64) -> bool {
    if a.is_nan() && b.is_nan() { return true; }
    (a - b).abs() < 1e-10
}

/// Generate a random codon index sequence (valid indices 0-63, or 255 for gap).
fn arb_codon_seq(len: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop_oneof![
            9 => 0u8..64,   // 90% valid codons
            1 => Just(255u8), // 10% gaps
        ],
        len,
    )
}

// ─── Nei model properties ────────────────────────────────────────────────────

mod nei {
    use super::*;

    fn tables() -> Box<eskaks::models::nei::NeiTables> {
        eskaks::models::nei::NeiTables::new()
    }

    proptest! {
        #[test]
        fn identical_seqs_give_zero(seq in arb_codon_seq(100)) {
            let t = tables();
            let (dn, ds) = t.compute_pair(&seq, &seq);
            prop_assert!((dn - 0.0).abs() < 1e-10, "dn should be 0 for identical seqs, got {}", dn);
            prop_assert!((ds - 0.0).abs() < 1e-10, "ds should be 0 for identical seqs, got {}", ds);
        }

        #[test]
        fn symmetry(seq1 in arb_codon_seq(50), seq2 in arb_codon_seq(50)) {
            let t = tables();
            let (dn_ab, ds_ab) = t.compute_pair(&seq1, &seq2);
            let (dn_ba, ds_ba) = t.compute_pair(&seq2, &seq1);
            prop_assert!(approx_eq(dn_ab, dn_ba), "dn not symmetric: {} vs {}", dn_ab, dn_ba);
            prop_assert!(approx_eq(ds_ab, ds_ba), "ds not symmetric: {} vs {}", ds_ab, ds_ba);
        }

        #[test]
        fn dn_ds_non_negative(seq1 in arb_codon_seq(50), seq2 in arb_codon_seq(50)) {
            let t = tables();
            let (dn, ds) = t.compute_pair(&seq1, &seq2);
            // NaN (from saturation) is acceptable; negative values are not
            prop_assert!(dn.is_nan() || dn >= 0.0, "dn should be >= 0, got {}", dn);
            prop_assert!(ds.is_nan() || ds >= 0.0, "ds should be >= 0, got {}", ds);
        }

        #[test]
        fn all_gaps_give_nan(len in 1usize..100) {
            let t = tables();
            let gaps = vec![255u8; len];
            let (dn, ds) = t.compute_pair(&gaps, &gaps);
            prop_assert!(dn.is_nan(), "dn should be NaN for all-gap sequences, got {}", dn);
            prop_assert!(ds.is_nan(), "ds should be NaN for all-gap sequences, got {}", ds);
        }
    }
}

// ─── Li model properties ─────────────────────────────────────────────────────

mod li {
    use super::*;

    fn tables() -> Box<eskaks::models::li::LiTables> {
        eskaks::models::li::LiTables::new()
    }

    proptest! {
        #[test]
        fn identical_seqs_give_zero(seq in arb_codon_seq(100)) {
            let t = tables();
            let (dn, ds) = t.compute_pair(&seq, &seq);
            prop_assert!((dn - 0.0).abs() < 1e-10, "dn should be 0 for identical seqs, got {}", dn);
            prop_assert!((ds - 0.0).abs() < 1e-10, "ds should be 0 for identical seqs, got {}", ds);
        }

        #[test]
        fn symmetry(seq1 in arb_codon_seq(50), seq2 in arb_codon_seq(50)) {
            let t = tables();
            let (dn_ab, ds_ab) = t.compute_pair(&seq1, &seq2);
            let (dn_ba, ds_ba) = t.compute_pair(&seq2, &seq1);
            prop_assert!(approx_eq(dn_ab, dn_ba), "dn not symmetric: {} vs {}", dn_ab, dn_ba);
            prop_assert!(approx_eq(ds_ab, ds_ba), "ds not symmetric: {} vs {}", ds_ab, ds_ba);
        }

        #[test]
        fn dn_ds_non_negative(seq1 in arb_codon_seq(50), seq2 in arb_codon_seq(50)) {
            let t = tables();
            let (dn, ds) = t.compute_pair(&seq1, &seq2);
            prop_assert!(dn.is_nan() || dn >= 0.0, "dn should be >= 0, got {}", dn);
            prop_assert!(ds.is_nan() || ds >= 0.0, "ds should be >= 0, got {}", ds);
        }

        #[test]
        fn all_gaps_give_nan(len in 1usize..100) {
            let t = tables();
            let gaps = vec![255u8; len];
            let (dn, ds) = t.compute_pair(&gaps, &gaps);
            prop_assert!(dn.is_nan(), "dn should be NaN for all-gap sequences, got {}", dn);
            prop_assert!(ds.is_nan(), "ds should be NaN for all-gap sequences, got {}", ds);
        }
    }
}

// ─── Cross-model consistency ─────────────────────────────────────────────────

mod cross_model {
    use super::*;

    proptest! {
        /// Both models should agree on identical sequences (both zero).
        #[test]
        fn both_models_zero_for_identical(seq in arb_codon_seq(100)) {
            let nei = eskaks::models::nei::NeiTables::new();
            let li = eskaks::models::li::LiTables::new();
            let (dn_nei, ds_nei) = nei.compute_pair(&seq, &seq);
            let (dn_li, ds_li) = li.compute_pair(&seq, &seq);
            prop_assert!((dn_nei - 0.0).abs() < 1e-10);
            prop_assert!((ds_nei - 0.0).abs() < 1e-10);
            prop_assert!((dn_li - 0.0).abs() < 1e-10);
            prop_assert!((ds_li - 0.0).abs() < 1e-10);
        }

        /// Both models should agree that identical sequences have zero dN/dS,
        /// and both should produce NaN for all-gap sequences.
        /// (Nei and Li CAN diverge on nonzero values — different counting methods.)
        #[test]
        fn both_nan_for_gaps(len in 1usize..50) {
            let nei = eskaks::models::nei::NeiTables::new();
            let li = eskaks::models::li::LiTables::new();
            let gaps = vec![255u8; len];
            let (dn_nei, ds_nei) = nei.compute_pair(&gaps, &gaps);
            let (dn_li, ds_li) = li.compute_pair(&gaps, &gaps);
            prop_assert!(dn_nei.is_nan() && ds_nei.is_nan(), "nei should be NaN for gaps");
            prop_assert!(dn_li.is_nan() && ds_li.is_nan(), "li should be NaN for gaps");
        }
    }
}

// ─── Population-genetics diversity properties ─────────────────────────────────

mod diversity {
    use super::*;
    use eskaks::stats::{tajimas_d, theta_pi_varn, theta_watterson};

    // Uniform-n sites for the fixed-sample-size property tests: pair every derived
    // count with the same n, the (k, n_i) shape theta_pi_varn consumes.
    fn uniform_sites(n: usize, counts: &[usize]) -> Vec<(usize, usize)> {
        counts.iter().map(|&k| (k, n)).collect()
    }

    // An explicit haploid genotype matrix: `n` samples over `sites` sites, flattened
    // row-major into a bit vector (0 = ancestral, 1 = derived).
    fn arb_genotype_matrix() -> impl Strategy<Value = (usize, usize, Vec<u8>)> {
        (2usize..8, 1usize..12).prop_flat_map(|(n, sites)| {
            prop::collection::vec(0u8..2, n * sites).prop_map(move |bits| (n, sites, bits))
        })
    }

    proptest! {
        // Differential oracle: nucleotide diversity from the tool's site-frequency
        // formula must equal the brute-force average number of pairwise differences
        // over the same explicit genotype matrix (two independent computations).
        #[test]
        fn pi_matches_brute_force_pairwise((n, sites, bits) in arb_genotype_matrix()) {
            let gt = |sample: usize, site: usize| bits[sample * sites + site];

            // Tool: per-site derived count -> theta_pi over the segregating sites.
            let mut counts = Vec::new();
            for site in 0..sites {
                let k = (0..n).filter(|&s| gt(s, site) == 1).count();
                if k > 0 && k < n {
                    counts.push(k);
                }
            }
            let tool = theta_pi_varn(&uniform_sites(n, &counts));

            // Oracle: total pairwise differences / number of pairs.
            let mut diffs = 0usize;
            for i in 0..n {
                for j in (i + 1)..n {
                    diffs += (0..sites).filter(|&site| gt(i, site) != gt(j, site)).count();
                }
            }
            let oracle = diffs as f64 / (n * (n - 1) / 2) as f64;
            prop_assert!((tool - oracle).abs() < 1e-9, "pi: tool={tool} oracle={oracle}");
        }

        // Nucleotide diversity is invariant to which allele is called "derived":
        // folding every site's count k -> n-k must give the identical value.
        #[test]
        fn pi_is_polarization_invariant(
            n in 2usize..40,
            raw in prop::collection::vec(0usize..40, 0..25),
        ) {
            let counts: Vec<usize> = raw.iter().map(|&k| k % (n + 1)).collect();
            let folded: Vec<usize> = counts.iter().map(|&k| n - k).collect();
            prop_assert!(approx_eq(
                theta_pi_varn(&uniform_sites(n, &counts)),
                theta_pi_varn(&uniform_sites(n, &folded))
            ));
        }

        // pi is a finite, non-negative quantity for any segregating-count vector.
        #[test]
        fn pi_is_finite_non_negative(
            n in 2usize..50,
            raw in prop::collection::vec(0usize..50, 0..30),
        ) {
            let counts: Vec<usize> = raw.iter().map(|&k| k % (n + 1)).collect();
            let pi = theta_pi_varn(&uniform_sites(n, &counts));
            prop_assert!(pi.is_finite() && pi >= 0.0, "pi = {pi}");
        }

        // Watterson theta is linear in the number of segregating sites and >= 0.
        #[test]
        fn watterson_is_linear_in_s(n in 2usize..40, s in 0usize..200) {
            let base = theta_watterson(n, s);
            let double = theta_watterson(n, 2 * s);
            prop_assert!(base >= 0.0);
            prop_assert!((double - 2.0 * base).abs() <= 1e-9 * (1.0 + double.abs()));
        }

        // Tajima's D is never an infinity: it is a finite number, or NaN when
        // undefined (n < 2, or no segregating sites), never garbage.
        #[test]
        fn tajimas_d_is_finite_or_nan(n in 0usize..60, s in 0usize..200, pi in 0.0f64..50.0) {
            prop_assert!(!tajimas_d(n, s, pi).is_infinite());
        }

        #[test]
        fn tajimas_d_undefined_regimes_are_nan(pi in 0.0f64..10.0) {
            prop_assert!(tajimas_d(1, 5, pi).is_nan(), "n < 2 must be NaN");
            prop_assert!(tajimas_d(10, 0, pi).is_nan(), "s = 0 must be NaN");
        }
    }
}

// ─── Statistical-distribution properties ──────────────────────────────────────

mod distributions {
    use super::*;
    use eskaks::stats::{benjamini_hochberg, binomial_two_sided_p, wilson_interval};

    proptest! {
        // A Wilson interval stays inside [0, 1] and brackets the point estimate.
        #[test]
        fn wilson_brackets_point_estimate(k in 0u64..1000, extra in 0u64..1000) {
            let n = k + extra + 1; // n >= 1 and k <= n
            let (lo, hi) = wilson_interval(k, n, 0.95);
            let phat = k as f64 / n as f64;
            prop_assert!(0.0 <= lo && lo <= hi && hi <= 1.0, "interval [{lo}, {hi}]");
            prop_assert!(lo - 1e-9 <= phat && phat <= hi + 1e-9, "phat {phat} outside [{lo}, {hi}]");
        }

        // A two-sided binomial p-value is a probability (or NaN at a degenerate p0).
        #[test]
        fn binomial_p_is_a_probability(k in 0u64..400, extra in 0u64..400, p0 in 0.01f64..0.99) {
            let p = binomial_two_sided_p(k, k + extra, p0);
            prop_assert!(p.is_nan() || (0.0..=1.0).contains(&p), "p = {p}");
        }

        // Benjamini-Hochberg adjustment preserves length, stays in [0, 1], and never
        // reports a q-value below its raw p-value (the adjustment only inflates).
        #[test]
        fn bh_q_values_are_valid(pvals in prop::collection::vec(0.0f64..=1.0, 1..60)) {
            let q = benjamini_hochberg(&pvals);
            prop_assert_eq!(q.len(), pvals.len());
            for (qi, pi) in q.iter().zip(&pvals) {
                prop_assert!((0.0..=1.0).contains(qi), "q out of range: {qi}");
                prop_assert!(*qi >= *pi - 1e-9, "BH q {qi} below raw p {pi}");
            }
        }
    }
}

// ─── Parser robustness: no input may ever panic ───────────────────────────────

mod parser_robustness {
    use super::*;
    use std::io::Write;

    fn write_tmp(bytes: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(bytes).unwrap();
        f.flush().unwrap();
        f
    }

    // A VCF data line with adversarial-but-plausible fields (huge/empty/negative POS,
    // odd REF/ALT, malformed AF) to probe the field parsing deeply.
    fn arb_vcf_line() -> impl Strategy<Value = String> {
        (
            "[A-Za-z0-9_.]{0,10}",
            "-?[0-9]{0,12}",
            "[ACGTN.,<>*]{0,6}",
            "[ACGTN.,*]{0,6}",
            "[A-Za-z0-9.,=;eE+-]{0,20}",
        )
            .prop_map(|(c, p, r, a, info)| format!("{c}\t{p}\t.\t{r}\t{a}\t60\tPASS\t{info}"))
    }

    proptest! {
        // Arbitrary bytes must never panic the VCF parser: a clean Ok or Err only.
        #[test]
        fn parse_vcf_never_panics(bytes in prop::collection::vec(any::<u8>(), 0..4096)) {
            let f = write_tmp(&bytes);
            let _ = eskaks::vcf::parse_vcf(f.path());
        }

        #[test]
        fn parse_gff3_never_panics(bytes in prop::collection::vec(any::<u8>(), 0..4096)) {
            let f = write_tmp(&bytes);
            let _ = eskaks::gff::parse_gff3(f.path());
        }

        // Structured-but-hostile VCF records exercise POS/AF/allele parsing paths.
        #[test]
        fn parse_vcf_structured_never_panics(lines in prop::collection::vec(arb_vcf_line(), 0..40)) {
            let mut body =
                String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
            for l in &lines {
                body.push_str(l);
                body.push('\n');
            }
            let f = write_tmp(body.as_bytes());
            let _ = eskaks::vcf::parse_vcf(f.path());
        }
    }
}
