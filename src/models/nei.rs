use crate::codon::INVALID_CODON;

/// Jukes-Cantor saturation threshold: p-values >= this indicate saturation.
const JC_SATURATION_THRESHOLD: f64 = 0.749;

/// Standard genetic code amino acids, indexed by Nei codon index.
/// Index: 16*remap(b2) + 4*remap(b1) + remap(b3), where remap: A→2, C→1, G→3, T→0.
const AA_ARRAY: [char; 64] = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', '*', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G',
];

/// Synonymous sites per codon (x3). Divide by 3 to get actual synonymous sites.
const SYN_SITE_ARRAY: [usize; 64] = [
    1, 1, 2, 2, 3, 3, 4, 4, 2, 2, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    3, 3, 4, 4, 1, 1, 2, 2, 3, 3, 3, 3,
];

/// Precomputed table of (synonymous_diffs, nonsynonymous_diffs) for each pair of
/// Nei codon indices, using Nei-Gojobori (1986) pathway analysis.
///
/// For codons differing at 1 position: direct syn/nonsyn classification.
/// For codons differing at 2 positions: average over 2 pathways.
/// For codons differing at 3 positions: average over 6 pathways.
/// Pathways through stop codon intermediates are excluded; if all pathways
/// are invalid, a fallback is used (0.5/1.5 for 2-diff, 1.0/2.0 for 3-diff).
///
/// Table is 4096 entries * 8 bytes = 32 KB, fits in L1 cache.
pub struct NeiTables {
    diff_table: [(f32, f32); 4096],
}

/// Reconstruct a Nei codon index from 3 slot values.
#[inline(always)]
fn slots_to_index(s: &[u8; 3]) -> usize {
    s[0] as usize * 16 + s[1] as usize * 4 + s[2] as usize
}

impl NeiTables {
    /// Build the pathway analysis diff table for all 4096 codon pairs.
    pub fn new() -> Box<NeiTables> {
        let mut tables = Box::new(NeiTables {
            diff_table: [(0.0f32, 0.0f32); 4096],
        });

        for idx_a in 0u8..64 {
            for idx_b in (idx_a + 1)..64 {
                // Decompose into 3 slots
                let slots_a = [idx_a / 16, (idx_a / 4) % 4, idx_a % 4];
                let slots_b = [idx_b / 16, (idx_b / 4) % 4, idx_b % 4];

                // Find differing positions
                let mut diffs: [usize; 3] = [0; 3];
                let mut ndiffs = 0;
                for pos in 0..3 {
                    if slots_a[pos] != slots_b[pos] {
                        diffs[ndiffs] = pos;
                        ndiffs += 1;
                    }
                }

                if ndiffs == 0 { continue; }

                let (sd, nd) = match ndiffs {
                    1 => {
                        // Single nucleotide change, no intermediates
                        if AA_ARRAY[idx_a as usize] == AA_ARRAY[idx_b as usize] {
                            (1.0f64, 0.0f64)
                        } else {
                            (0.0f64, 1.0f64)
                        }
                    }
                    2 => Self::pathway_2diff(&slots_a, &slots_b, &diffs, idx_a, idx_b),
                    3 => Self::pathway_3diff(&slots_a, &slots_b, idx_a, idx_b),
                    _ => unreachable!(),
                };

                let entry = (sd as f32, nd as f32);
                // Table is symmetric
                tables.diff_table[idx_a as usize * 64 + idx_b as usize] = entry;
                tables.diff_table[idx_b as usize * 64 + idx_a as usize] = entry;
            }
        }

        tables
    }

    /// Pathway analysis for codons differing at exactly 2 positions.
    /// 2 pathways (permutations of 2 diff positions).
    fn pathway_2diff(
        slots_a: &[u8; 3], slots_b: &[u8; 3], diffs: &[usize; 3],
        idx_a: u8, idx_b: u8,
    ) -> (f64, f64) {
        let perms = [(diffs[0], diffs[1]), (diffs[1], diffs[0])];
        let mut total_sd = 0.0f64;
        let mut total_nd = 0.0f64;
        let mut valid_paths = 0u32;

        for &(first, _second) in &perms {
            // Intermediate codon: change `first` position
            let mut inter = *slots_a;
            inter[first] = slots_b[first];
            let inter_idx = slots_to_index(&inter);

            // Skip pathway if intermediate is a stop codon
            if AA_ARRAY[inter_idx] == '*' { continue; }

            // Step 1: codon_a -> intermediate
            if AA_ARRAY[idx_a as usize] == AA_ARRAY[inter_idx] {
                total_sd += 1.0;
            } else {
                total_nd += 1.0;
            }
            // Step 2: intermediate -> codon_b
            if AA_ARRAY[inter_idx] == AA_ARRAY[idx_b as usize] {
                total_sd += 1.0;
            } else {
                total_nd += 1.0;
            }
            valid_paths += 1;
        }

        if valid_paths > 0 {
            (total_sd / valid_paths as f64, total_nd / valid_paths as f64)
        } else {
            // Fallback: all pathways through stop codons
            (0.5, 1.5)
        }
    }

    /// Pathway analysis for codons differing at all 3 positions.
    /// 6 pathways (3! = 6 permutations).
    fn pathway_3diff(
        slots_a: &[u8; 3], slots_b: &[u8; 3],
        idx_a: u8, idx_b: u8,
    ) -> (f64, f64) {
        let perms: [(usize, usize, usize); 6] = [
            (0, 1, 2), (0, 2, 1), (1, 0, 2),
            (1, 2, 0), (2, 0, 1), (2, 1, 0),
        ];
        let mut total_sd = 0.0f64;
        let mut total_nd = 0.0f64;
        let mut valid_paths = 0u32;

        for &(first, second, _third) in &perms {
            let mut cur = *slots_a;

            // Intermediate 1: change `first` position
            cur[first] = slots_b[first];
            let inter1_idx = slots_to_index(&cur);
            if AA_ARRAY[inter1_idx] == '*' { continue; }

            // Intermediate 2: change `second` position
            cur[second] = slots_b[second];
            let inter2_idx = slots_to_index(&cur);
            if AA_ARRAY[inter2_idx] == '*' { continue; }

            // All intermediates valid, classify 3 steps
            let aa_a = AA_ARRAY[idx_a as usize];
            let aa_1 = AA_ARRAY[inter1_idx];
            let aa_2 = AA_ARRAY[inter2_idx];
            let aa_b = AA_ARRAY[idx_b as usize];

            let mut path_sd = 0.0;
            let mut path_nd = 0.0;
            if aa_a == aa_1 { path_sd += 1.0; } else { path_nd += 1.0; }
            if aa_1 == aa_2 { path_sd += 1.0; } else { path_nd += 1.0; }
            if aa_2 == aa_b { path_sd += 1.0; } else { path_nd += 1.0; }

            total_sd += path_sd;
            total_nd += path_nd;
            valid_paths += 1;
        }

        if valid_paths > 0 {
            (total_sd / valid_paths as f64, total_nd / valid_paths as f64)
        } else {
            // Fallback: all pathways through stop codons
            (1.0, 2.0)
        }
    }

    /// Calculates dN and dS using the Nei-Gojobori model with Jukes-Cantor correction.
    /// Uses the precomputed pathway analysis diff table for proper handling of
    /// codons differing at multiple nucleotide positions.
    #[inline]
    pub fn compute_pair(&self, codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
        let mut count_valid_codons = 0u32;
        let mut syn_diffs = 0.0f64;
        let mut nonsyn_diffs = 0.0f64;
        let mut sum_syn_sites_seq1 = 0.0f64;
        let mut sum_syn_sites_seq2 = 0.0f64;

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
                // SAFETY: all indices verified < 64, which is within bounds for
                // SYN_SITE_ARRAY (64 entries) and diff_table (4096 entries, max index = 63*64+63 = 4095)
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
                    if a1 != b1 { let e = *self.diff_table.get_unchecked(a1 * 64 + b1); syn_diffs += e.0 as f64; nonsyn_diffs += e.1 as f64; }
                    if a2 != b2 { let e = *self.diff_table.get_unchecked(a2 * 64 + b2); syn_diffs += e.0 as f64; nonsyn_diffs += e.1 as f64; }
                    if a3 != b3 { let e = *self.diff_table.get_unchecked(a3 * 64 + b3); syn_diffs += e.0 as f64; nonsyn_diffs += e.1 as f64; }
                    if a4 != b4 { let e = *self.diff_table.get_unchecked(a4 * 64 + b4); syn_diffs += e.0 as f64; nonsyn_diffs += e.1 as f64; }
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
                        let e = self.diff_table[i1 * 64 + i2];
                        syn_diffs += e.0 as f64;
                        nonsyn_diffs += e.1 as f64;
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
                let e = self.diff_table[idx1 * 64 + idx2];
                syn_diffs += e.0 as f64;
                nonsyn_diffs += e.1 as f64;
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::fasta_to_codon_indices;
    use crate::models::Model;

    /// Helper: run the Nei model on two ASCII DNA sequences.
    fn nei(seq1: &[u8], seq2: &[u8]) -> (f64, f64) {
        let tables = NeiTables::new();
        let idx1 = fasta_to_codon_indices(seq1, Model::Nei);
        let idx2 = fasta_to_codon_indices(seq2, Model::Nei);
        tables.compute_pair(&idx1, &idx2)
    }

    const EPSILON: f64 = 1e-4;

    // --- Trivial cases ---

    #[test]
    fn nei_identical_sequences_give_zero() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGGCTGCT");
        assert_eq!(dn, 0.0);
        assert_eq!(ds, 0.0);
    }

    #[test]
    fn nei_all_ambiguous_gives_nan() {
        let (dn, ds) = nei(b"NNNNNNNNN", b"NNNNNNNNN");
        assert!(dn.is_nan(), "expected dN=NaN, got {}", dn);
        assert!(ds.is_nan(), "expected dS=NaN, got {}", ds);
    }

    // --- One synonymous change (1 nucleotide position differs) ---
    // seq1: ATG GCT GCT  (Met, Ala, Ala)
    // seq2: ATG GCC GCT  (Met, Ala, Ala)  -- GCT->GCC at codon 2, 1 position differs
    // Pathway analysis: 1-diff case, same AA (Ala->Ala) -> sd=1, nd=0
    // S = (6/3 + 6/3) / 2 = 2.0,  N = 7.0
    // pS = 1/2 = 0.5,  pN = 0/7 = 0
    // dS = -0.75*ln(1 - 4/3*0.5) = -0.75*ln(1/3) ~ 0.82396
    // dN = 0.0
    #[test]
    fn nei_one_synonymous_change() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        assert!((dn - 0.0).abs() < EPSILON, "dN should be 0, got {}", dn);
        assert!((ds - 0.82396).abs() < EPSILON, "dS should be ~0.82396, got {}", ds);
    }

    // --- One amino acid change but 2 nucleotide positions differ ---
    // seq1: ATG GCT GCT  (Met, Ala, Ala)
    // seq3: ATG ATT GCT  (Met, Ile, Ala)  -- GCT->ATT, 2 positions differ (pos 0: G->A, pos 1: C->T)
    // Pathway analysis for GCT->ATT (2 diffs, 2 pathways):
    //   Path 1: GCT->ACT->ATT  (Ala->Thr->Ile: nd, nd)  sd=0, nd=2
    //   Path 2: GCT->GTT->ATT  (Ala->Val->Ile: nd, nd)   sd=0, nd=2
    //   Average: sd=0, nd=2
    // S = (6/3 + 5/3) / 2 = 11/6,  N = 43/6
    // pS = 0,  pN = 2/(43/6) = 12/43
    // dS = 0.0
    // dN = -0.75*ln(1 - 4/3*12/43) ~ 0.34888
    #[test]
    fn nei_two_position_nonsynonymous_change() {
        let (dn, ds) = nei(b"ATGGCTGCT", b"ATGATTGCT");
        assert!((ds - 0.0).abs() < EPSILON, "dS should be 0, got {}", ds);
        assert!((dn - 0.34902).abs() < EPSILON, "dN should be ~0.34902, got {}", dn);
    }

    // --- Symmetry: swap seq1 and seq2 should give the same result ---
    #[test]
    fn nei_result_is_symmetric() {
        let (dn_ab, ds_ab) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        let (dn_ba, ds_ba) = nei(b"ATGGCCGCT", b"ATGGCTGCT");
        assert!((dn_ab - dn_ba).abs() < EPSILON, "dN should be symmetric");
        assert!((ds_ab - ds_ba).abs() < EPSILON, "dS should be symmetric");
    }

    // --- Symmetry for multi-position diffs ---
    #[test]
    fn nei_result_is_symmetric_multi_diff() {
        let (dn_ab, ds_ab) = nei(b"ATGGCTGCT", b"ATGATTGCT");
        let (dn_ba, ds_ba) = nei(b"ATGATTGCT", b"ATGGCTGCT");
        assert!((dn_ab - dn_ba).abs() < EPSILON, "dN should be symmetric for multi-diff");
        assert!((ds_ab - ds_ba).abs() < EPSILON, "dS should be symmetric for multi-diff");
    }

    // --- RNA input (U instead of T) ---
    #[test]
    fn nei_rna_input_gives_same_as_dna() {
        let (dn_dna, ds_dna) = nei(b"ATGGCTGCT", b"ATGGCCGCT");
        let (dn_rna, ds_rna) = nei(b"AUGGCUGCU", b"AUGGCCGCU");
        assert!((dn_dna - dn_rna).abs() < EPSILON);
        assert!((ds_dna - ds_rna).abs() < EPSILON);
    }

    // --- Jukes-Cantor saturation: pS >= 0.749 -> dS = NaN ---
    // GCTGCT vs GCAGCA: 2 synonymous changes (GCT->GCA, Ala->Ala at each codon)
    // Each is a 1-diff case. sd=2.
    // S_avg = 2.0, pS = 2/2.0 = 1.0 >= JC_SATURATION_THRESHOLD (0.749) -> dS = NaN
    #[test]
    fn nei_saturated_ps_gives_nan_ds() {
        let (dn, ds) = nei(b"GCTGCT", b"GCAGCA");
        assert!(ds.is_nan(), "saturated pS should give NaN dS, got {}", ds);
        assert!((dn).abs() < EPSILON, "dN should be 0 for purely synonymous pair, got {}", dn);
    }

    // --- Pathway analysis: verify diff table entries directly ---
    #[test]
    fn nei_diff_table_one_position_syn() {
        // GCT (Ala) -> GCC (Ala): 1 position differs, synonymous
        let tables = NeiTables::new();
        let idx_gct = fasta_to_codon_indices(b"GCT", Model::Nei)[0] as usize;
        let idx_gcc = fasta_to_codon_indices(b"GCC", Model::Nei)[0] as usize;
        let e = tables.diff_table[idx_gct * 64 + idx_gcc];
        assert!((e.0 - 1.0).abs() < 1e-6, "sd should be 1.0, got {}", e.0);
        assert!((e.1 - 0.0).abs() < 1e-6, "nd should be 0.0, got {}", e.1);
    }

    #[test]
    fn nei_diff_table_two_position_all_nonsyn() {
        // GCT (Ala) -> ATT (Ile): 2 positions differ, all pathways nonsynonymous
        let tables = NeiTables::new();
        let idx_gct = fasta_to_codon_indices(b"GCT", Model::Nei)[0] as usize;
        let idx_att = fasta_to_codon_indices(b"ATT", Model::Nei)[0] as usize;
        let e = tables.diff_table[idx_gct * 64 + idx_att];
        assert!((e.0 - 0.0).abs() < 1e-6, "sd should be 0.0, got {}", e.0);
        assert!((e.1 - 2.0).abs() < 1e-6, "nd should be 2.0, got {}", e.1);
    }

    #[test]
    fn nei_diff_table_is_symmetric() {
        let tables = NeiTables::new();
        for a in 0..64usize {
            for b in 0..64usize {
                let ab = tables.diff_table[a * 64 + b];
                let ba = tables.diff_table[b * 64 + a];
                assert!((ab.0 - ba.0).abs() < 1e-6 && (ab.1 - ba.1).abs() < 1e-6,
                    "diff_table not symmetric for ({}, {}): ({}, {}) vs ({}, {})", a, b, ab.0, ab.1, ba.0, ba.1);
            }
        }
    }

    // --- dN/dS > 0 for mixed changes (longer sequence) ---
    // seq1: ATG GCT GCT GGT  (Met, Ala, Ala, Gly)
    // seq2: ATG ATT GCC GGT  (Met, Ile, Ala, Gly)
    //   GCT->ATT: 2-diff pathway -> sd=0, nd=2
    //   GCT->GCC: 1-diff -> sd=1, nd=0
    //   Total: sd=1, nd=2
    // S=(9/3+8/3)/2=17/6, N=55/6
    // pS=6/17, pN=12/55
    // Both dN and dS > 0
    #[test]
    fn nei_mixed_changes_both_dn_ds_positive() {
        let (dn, ds) = nei(b"ATGGCTGCTGGT", b"ATGATTGCCGGT");
        assert!(dn > 0.0, "dN should be > 0, got {}", dn);
        assert!(ds > 0.0, "dS should be > 0, got {}", ds);
    }
}
