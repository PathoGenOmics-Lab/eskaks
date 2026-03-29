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
