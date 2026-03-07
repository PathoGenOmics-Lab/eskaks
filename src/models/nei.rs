use crate::codon::INVALID_CODON;

/// Umbral de saturacion de Jukes-Cantor: valores p >= este indican saturacion.
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

/// Calcula dN y dS usando el modelo Nei-Gojobori con correccion Jukes-Cantor.
pub fn calculate_syn_nonsyn_from_indices(codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
    let mut count_valid_codons = 0;
    let mut syn_diffs = 0.0;
    let mut nonsyn_diffs = 0.0;
    let mut sum_syn_sites_seq1 = 0.0;
    let mut sum_syn_sites_seq2 = 0.0;

    for k in 0..codon_indices1.len() {
        let idx1 = codon_indices1[k] as usize;
        let idx2 = codon_indices2[k] as usize;
        if idx1 >= INVALID_CODON as usize || idx2 >= INVALID_CODON as usize {
            continue;
        }
        count_valid_codons += 1;
        sum_syn_sites_seq1 += SYN_SITE_ARRAY[idx1] as f64;
        sum_syn_sites_seq2 += SYN_SITE_ARRAY[idx2] as f64;
        if idx1 != idx2 {
            let aa1 = AA_ARRAY[idx1];
            let aa2 = AA_ARRAY[idx2];
            if aa1 == aa2 { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; }
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
