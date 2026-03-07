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
