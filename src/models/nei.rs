use crate::codon::INVALID_CODON;

/// Jukes-Cantor saturation threshold: p-values >= this indicate saturation.
const JC_SATURATION_THRESHOLD: f64 = 0.749;

const AA_ARRAY: [char; 64] = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', '*', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G',
];

const SYN_SITE_ARRAY: [usize; 64] = [
    1, 1, 2, 2, 3, 3, 4, 4, 2, 2, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    3, 3, 4, 4, 1, 1, 2, 2, 3, 3, 3, 3,
];

/// Calculates dN and dS using the Nei-Gojobori model with Jukes-Cantor correction.
#[inline]
pub fn calculate_syn_nonsyn_from_indices(codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
    let mut count_valid_codons = 0u32;
    let mut syn_diffs = 0.0;
    let mut nonsyn_diffs = 0.0;
    let mut sum_syn_sites_seq1 = 0.0;
    let mut sum_syn_sites_seq2 = 0.0;

    // Process 4 codons at a time for branch-free fast path
    let chunks1 = codon_indices1.chunks_exact(4);
    let chunks2 = codon_indices2.chunks_exact(4);
    let rem1 = chunks1.remainder();
    let rem2 = chunks2.remainder();

    for (ch1, ch2) in chunks1.zip(chunks2) {
        let a1 = ch1[0] as usize; let a2 = ch1[1] as usize;
        let a3 = ch1[2] as usize; let a4 = ch1[3] as usize;
        let b1 = ch2[0] as usize; let b2 = ch2[1] as usize;
        let b3 = ch2[2] as usize; let b4 = ch2[3] as usize;

        if (a1 | a2 | a3 | a4 | b1 | b2 | b3 | b4) < INVALID_CODON as usize {
            // SAFETY: all indices verified < 64, which is within AA_ARRAY and SYN_SITE_ARRAY bounds
            unsafe {
                count_valid_codons += 4;
                sum_syn_sites_seq1 += *SYN_SITE_ARRAY.get_unchecked(a1) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(a2) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(a3) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(a4) as f64;
                sum_syn_sites_seq2 += *SYN_SITE_ARRAY.get_unchecked(b1) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(b2) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(b3) as f64
                    + *SYN_SITE_ARRAY.get_unchecked(b4) as f64;
                if a1 != b1 { if *AA_ARRAY.get_unchecked(a1) == *AA_ARRAY.get_unchecked(b1) { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; } }
                if a2 != b2 { if *AA_ARRAY.get_unchecked(a2) == *AA_ARRAY.get_unchecked(b2) { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; } }
                if a3 != b3 { if *AA_ARRAY.get_unchecked(a3) == *AA_ARRAY.get_unchecked(b3) { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; } }
                if a4 != b4 { if *AA_ARRAY.get_unchecked(a4) == *AA_ARRAY.get_unchecked(b4) { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; } }
            }
        } else {
            for (&c1, &c2) in ch1.iter().zip(ch2.iter()) {
                let i1 = c1 as usize;
                let i2 = c2 as usize;
                if i1 >= INVALID_CODON as usize || i2 >= INVALID_CODON as usize { continue; }
                count_valid_codons += 1;
                sum_syn_sites_seq1 += SYN_SITE_ARRAY[i1] as f64;
                sum_syn_sites_seq2 += SYN_SITE_ARRAY[i2] as f64;
                if i1 != i2 {
                    if AA_ARRAY[i1] == AA_ARRAY[i2] { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; }
                }
            }
        }
    }

    // Handle remainder codons (0-3)
    for (&c1, &c2) in rem1.iter().zip(rem2.iter()) {
        let idx1 = c1 as usize;
        let idx2 = c2 as usize;
        if idx1 >= INVALID_CODON as usize || idx2 >= INVALID_CODON as usize { continue; }
        count_valid_codons += 1;
        sum_syn_sites_seq1 += SYN_SITE_ARRAY[idx1] as f64;
        sum_syn_sites_seq2 += SYN_SITE_ARRAY[idx2] as f64;
        if idx1 != idx2 {
            if AA_ARRAY[idx1] == AA_ARRAY[idx2] { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; }
        }
    }

    if count_valid_codons == 0 {
        return (f64::NAN, f64::NAN);
    }

    let potential_syn_sites = (sum_syn_sites_seq1 / 3.0 + sum_syn_sites_seq2 / 3.0) / 2.0;
    let potential_nonsyn_sites = (count_valid_codons as f64) * 3.0 - potential_syn_sites;
    let ps = if potential_syn_sites > 0.0 { syn_diffs / potential_syn_sites } else { 0.0 };
    let pn = if potential_nonsyn_sites > 0.0 { nonsyn_diffs / potential_nonsyn_sites } else { 0.0 };

    let mut ds = if ps < JC_SATURATION_THRESHOLD {
        -0.75 * (1.0 - (4.0 / 3.0) * ps).ln()
    } else {
        f64::NAN
    };
    let mut dn = if pn < JC_SATURATION_THRESHOLD {
        -0.75 * (1.0 - (4.0 / 3.0) * pn).ln()
    } else {
        f64::NAN
    };

    if ds.is_finite() && ds < 0.0 { ds = 0.0; }
    if dn.is_finite() && dn < 0.0 { dn = 0.0; }

    (dn, ds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::fasta_to_codon_indices;
    use crate::models::Model;

    /// Helper: run the Nei model on two ASCII DNA sequences.
    fn nei(seq1: &[u8], seq2: &[u8]) -> (f64, f64) {
        let idx1 = fasta_to_codon_indices(seq1, Model::Nei);
        let idx2 = fasta_to_codon_indices(seq2, Model::Nei);
        calculate_syn_nonsyn_from_indices(&idx1, &idx2)
    }

    const EPSILON: f64 = 1e-5;

    // --- Trivial cases ---

    #[test]
    fn nei_identical_sequences_give_zero() {
        // dN = dS = 0 when sequences are identical
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGGCTGCT");
        assert_eq!(dn, 0.0);
        assert_eq!(ds, 0.0);
    }

    #[test]
    fn nei_all_ambiguous_gives_nan() {
        // All-N sequences → no valid codons → NaN
        let (dn, ds) = nei(b"NNNNNNNNN", b"NNNNNNNNN");
        assert!(dn.is_nan(), "expected dN=NaN, got {}", dn);
        assert!(ds.is_nan(), "expected dS=NaN, got {}", ds);
    }

    // --- One synonymous change ---
    // seq1: ATG GCT GCT  (Met, Ala, Ala)
    // seq2: ATG GCC GCT  (Met, Ala, Ala)  — GCT→GCC at codon 2 (Ala→Ala)
    // Hand-calculated:
    //   S = (6/3 + 6/3) / 2 = 2.0,  N = 7.0
    //   pS = 1/2 = 0.5,  pN = 0/7 = 0
    //   dS = -0.75·ln(1 - 4/3·0.5) = -0.75·ln(1/3) ≈ 0.82396
    //   dN = 0.0
    #[test]
    fn nei_one_synonymous_change() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        assert!((dn - 0.0).abs() < EPSILON, "dN should be 0, got {}", dn);
        assert!((ds - 0.82396).abs() < EPSILON, "dS should be ~0.82396, got {}", ds);
    }

    // --- One nonsynonymous change ---
    // seq1: ATG GCT GCT  (Met, Ala, Ala)
    // seq3: ATG ATT GCT  (Met, Ile, Ala)  — GCT→ATT at codon 2 (Ala→Ile)
    // Hand-calculated:
    //   S = (6/3 + 5/3) / 2 = 11/6 ≈ 1.8333,  N = 43/6 ≈ 7.1667
    //   pS = 0,  pN = 6/43 ≈ 0.13953
    //   dS = 0.0
    //   dN = -0.75·ln(1 - 4/3·6/43) ≈ 0.15439
    #[test]
    fn nei_one_nonsynonymous_change() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGATTGCT");
        assert!((ds - 0.0).abs() < EPSILON, "dS should be 0, got {}", ds);
        assert!((dn - 0.15439).abs() < EPSILON, "dN should be ~0.15439, got {}", dn);
    }

    // --- Two nonsynonymous changes ---
    // seq1: ATG GCT GCT  (Met, Ala, Ala)
    // seq4: ATG ATT ATG  (Met, Ile, Met)  — 2 nonsyn changes
    // Hand-calculated:
    //   S = (6/3 + 2/3) / 2 = 4/3 ≈ 1.3333,  N = 23/3 ≈ 7.6667
    //   pN = 2/(23/3) = 6/23 ≈ 0.26087
    //   dN ≈ 0.32058,  dS = 0.0
    #[test]
    fn nei_two_nonsynonymous_changes() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGATTATG");
        assert!((ds - 0.0).abs() < EPSILON, "dS should be 0, got {}", ds);
        assert!((dn - 0.32058).abs() < EPSILON, "dN should be ~0.32058, got {}", dn);
    }

    // --- dN/dS > 0 for mixed changes ---
    // Both synonymous and nonsynonymous differences present.
    // seq1: ATG GCT GCT GGT  (Met, Ala, Ala, Gly) — 4 codons, 12 bp
    // seq2: ATG ATT GCC GGT  (Met, Ile, Ala, Gly)
    //   Differences: GCT→ATT (Ala→Ile, nonsyn) + GCT→GCC (Ala→Ala, syn)
    // S=(9/3+8/3)/2=2.833, pS=1/2.833=0.353, pN=1/9.167=0.109  → both dN,dS > 0
    #[test]
    fn nei_mixed_changes_both_dn_ds_positive() {
        let (dn, ds) = nei(b"ATGGCTGCTGGT", b"ATGATTGCCGGT");
        assert!(dn > 0.0, "dN should be > 0, got {}", dn);
        assert!(ds > 0.0, "dS should be > 0, got {}", ds);
        assert!((dn - 0.11789).abs() < EPSILON, "dN should be ~0.11789, got {}", dn);
        assert!((ds - 0.47699).abs() < EPSILON, "dS should be ~0.47699, got {}", ds);
    }

    // --- Symmetry: swap seq1 and seq2 should give the same result ---
    #[test]
    fn nei_result_is_symmetric() {
        let (dn_ab, ds_ab) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        let (dn_ba, ds_ba) = nei(b"ATGGCCGCT", b"ATGGCTGCT");
        assert!((dn_ab - dn_ba).abs() < EPSILON, "dN should be symmetric");
        assert!((ds_ab - ds_ba).abs() < EPSILON, "dS should be symmetric");
    }

    // --- RNA input (U instead of T) ---
    #[test]
    fn nei_rna_input_gives_same_as_dna() {
        let (dn_dna, ds_dna) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        let (dn_rna, ds_rna) = nei(b"AUGGCUGCU", b"AUGGCCGCU");
        assert!((dn_dna - dn_rna).abs() < EPSILON);
        assert!((ds_dna - ds_rna).abs() < EPSILON);
    }

    // --- Jukes-Cantor saturation: pS >= 0.749 → dS = NaN ---
    // GCTGCT vs GCAGCA: 2 synonymous changes (GCT→GCA, Ala→Ala at each codon)
    // S_avg = 2.0, pS = 2/2.0 = 1.0 >= JC_SATURATION_THRESHOLD (0.749) → dS = NaN
    #[test]
    fn nei_saturated_ps_gives_nan_ds() {
        let (dn, ds) = nei(b"GCTGCT", b"GCAGCA");
        assert!(ds.is_nan(), "saturated pS should give NaN dS, got {}", ds);
        assert!((dn).abs() < EPSILON, "dN should be 0 for purely synonymous pair, got {}", dn);
    }
}
