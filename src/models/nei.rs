use crate::codon::INVALID_CODON;
use crate::genetic_code::{self, GeneticCode};

/// Jukes-Cantor saturation threshold: p-values >= this indicate saturation.
const JC_SATURATION_THRESHOLD: f64 = 0.749;

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
    /// Amino acid array (Nei indexing) for the active genetic code.
    aa_array: [char; 64],
    /// Synonymous SITES per codon (Nei indexing), Nei-Gojobori exclude-changes-to-stop
    /// convention (S = 3·syn/(syn+nonsyn); N = 3 − S). Matches `count_sites` in the vcf path.
    syn_site_array: [f64; 64],
    /// Whether each codon index is a stop. Stop codons are excluded from both site and
    /// diff counting (as the vcf `count_sites` does), so a terminal stop does not add
    /// phantom nonsynonymous sites and deflate dN and the dN/dS ratio.
    is_stop: [bool; 64],
}

/// Reconstruct a Nei codon index from 3 slot values.
#[inline(always)]
fn slots_to_index(s: &[u8; 3]) -> usize {
    s[0] as usize * 16 + s[1] as usize * 4 + s[2] as usize
}

impl NeiTables {
    /// Build the pathway analysis diff table using the standard genetic code (table 1).
    #[allow(dead_code)]
    pub fn new() -> Box<NeiTables> {
        let gc = genetic_code::get_table(1).unwrap();
        Self::with_genetic_code(gc)
    }

    /// Build the pathway analysis diff table for any NCBI genetic code.
    pub fn with_genetic_code(gc: &GeneticCode) -> Box<NeiTables> {
        let aa_array = genetic_code::li_to_nei_aa(&gc.aa_table);
        let syn_site_array = genetic_code::compute_syn_sites(&aa_array);
        let is_stop: [bool; 64] = std::array::from_fn(|i| aa_array[i] == '*');

        let mut tables = Box::new(NeiTables {
            diff_table: [(0.0f32, 0.0f32); 4096],
            aa_array,
            syn_site_array,
            is_stop,
        });

        let aa = &tables.aa_array;

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
                        if aa[idx_a as usize] == aa[idx_b as usize] {
                            (1.0f64, 0.0f64)
                        } else {
                            (0.0f64, 1.0f64)
                        }
                    }
                    2 => Self::pathway_2diff(aa, &slots_a, &slots_b, &diffs, idx_a, idx_b),
                    3 => Self::pathway_3diff(aa, &slots_a, &slots_b, idx_a, idx_b),
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
        aa: &[char; 64],
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
            if aa[inter_idx] == '*' { continue; }

            // Step 1: codon_a -> intermediate
            if aa[idx_a as usize] == aa[inter_idx] {
                total_sd += 1.0;
            } else {
                total_nd += 1.0;
            }
            // Step 2: intermediate -> codon_b
            if aa[inter_idx] == aa[idx_b as usize] {
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
        aa: &[char; 64],
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
            if aa[inter1_idx] == '*' { continue; }

            // Intermediate 2: change `second` position
            cur[second] = slots_b[second];
            let inter2_idx = slots_to_index(&cur);
            if aa[inter2_idx] == '*' { continue; }

            // All intermediates valid, classify 3 steps
            let aa_a = aa[idx_a as usize];
            let aa_1 = aa[inter1_idx];
            let aa_2 = aa[inter2_idx];
            let aa_b = aa[idx_b as usize];

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
        let (dn, ds, _, _) = self.compute_pair_stats(codon_indices1, codon_indices2);
        (dn, ds)
    }

    /// Like [`compute_pair`] but also returns the Nei-Gojobori analytic
    /// sampling variances `(var_dn, var_ds)` from the Jukes-Cantor delta method:
    /// `V(d) = p(1-p) / [L · (1 - 4/3·p)²]`, where `p` is the difference
    /// proportion and `L` the site count. NaN where the distance is NaN
    /// (saturation) or the sites are zero.
    #[inline]
    pub fn compute_pair_stats(&self, codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64, f64, f64) {
        // Compare strictly position-by-position over the common length. If the inputs
        // differ in length (an imperfect / differently-trimmed alignment — the loader
        // only warns, never bails), the chunks_exact(4) fast path below would otherwise
        // desync: chunks1.zip(chunks2) drops the longer sequence's extra full chunks and
        // the two `.remainder()` slices start at different codon offsets, mis-registering
        // codons and reporting an arbitrary wrong dN/dS. Clamping to min() keeps both the
        // chunk stream and the remainder aligned to the same absolute codon index.
        let n = codon_indices1.len().min(codon_indices2.len());
        let codon_indices1 = &codon_indices1[..n];
        let codon_indices2 = &codon_indices2[..n];

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

            let all_valid = (a1 | a2 | a3 | a4 | b1 | b2 | b3 | b4) < INVALID_CODON as usize;
            // Fast path only when all 8 codons are valid AND none is a stop; the `&&`
            // short-circuits so the is_stop reads happen only once every index is < 64.
            // A chunk containing a stop falls to the scalar loop below, which skips it.
            let fast = all_valid
                && !(self.is_stop[a1] || self.is_stop[a2] || self.is_stop[a3] || self.is_stop[a4]
                    || self.is_stop[b1] || self.is_stop[b2] || self.is_stop[b3] || self.is_stop[b4]);
            if fast {
                // SAFETY: All indices verified < 64 (INVALID_CODON) by the bitwise OR check above.
                // syn_site_array has 64 entries, diff_table has 4096 entries (max index = 63*64+63 = 4095).
                unsafe {
                    count_valid_codons += 4;
                    sum_syn_sites_seq1 += *self.syn_site_array.get_unchecked(a1)
                        + *self.syn_site_array.get_unchecked(a2)
                        + *self.syn_site_array.get_unchecked(a3)
                        + *self.syn_site_array.get_unchecked(a4);
                    sum_syn_sites_seq2 += *self.syn_site_array.get_unchecked(b1)
                        + *self.syn_site_array.get_unchecked(b2)
                        + *self.syn_site_array.get_unchecked(b3)
                        + *self.syn_site_array.get_unchecked(b4);
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
                    if self.is_stop[i1] || self.is_stop[i2] { continue; } // exclude stop codons
                    count_valid_codons += 1;
                    sum_syn_sites_seq1 += self.syn_site_array[i1];
                    sum_syn_sites_seq2 += self.syn_site_array[i2];
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
            if self.is_stop[idx1] || self.is_stop[idx2] { continue; } // exclude stop codons
            count_valid_codons += 1;
            sum_syn_sites_seq1 += self.syn_site_array[idx1];
            sum_syn_sites_seq2 += self.syn_site_array[idx2];
            if idx1 != idx2 {
                let e = self.diff_table[idx1 * 64 + idx2];
                syn_diffs += e.0 as f64;
                nonsyn_diffs += e.1 as f64;
            }
        }

        if count_valid_codons == 0 {
            return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
        }

        // syn_site_array already holds per-codon synonymous SITES (not change counts),
        // so average the two sequences directly; N is the per-codon complement (3 − S).
        let potential_syn_sites = (sum_syn_sites_seq1 + sum_syn_sites_seq2) / 2.0;
        let potential_nonsyn_sites = (count_valid_codons as f64) * 3.0 - potential_syn_sites;
        // A zero-site denominator means the quantity cannot be estimated (e.g. a pair of
        // all-Met/Trp codons has no synonymous sites), so it is undefined, not 0.0. A
        // spurious 0.0 would flow through JC into dS = 0 and read as strong purifying
        // selection; NaN matches the Li model on the same input. (Nonsyn sites are always
        // > 0 for real codons, but keep the guard symmetric.)
        let ps = if potential_syn_sites > 0.0 { syn_diffs / potential_syn_sites } else { f64::NAN };
        let pn = if potential_nonsyn_sites > 0.0 { nonsyn_diffs / potential_nonsyn_sites } else { f64::NAN };

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

        let var_ds = jc_variance(ps, potential_syn_sites);
        let var_dn = jc_variance(pn, potential_nonsyn_sites);

        (dn, ds, var_dn, var_ds)
    }
}

/// The pathway-analysis entry for one pair of Nei codon indices, so the VCF path's
/// per-codon scorer can be checked against this table rather than trusted to agree with
/// it (see `vcf_analysis::mnv::pathway_diffs`). Test-only: nothing outside the tests
/// needs to reach into the table, and exposing it in a normal build would read as dead
/// code in the binary target.
#[cfg(test)]
impl NeiTables {
    pub(crate) fn diff_entry(&self, idx_a: usize, idx_b: usize) -> (f64, f64) {
        let e = self.diff_table[idx_a * 64 + idx_b];
        (e.0 as f64, e.1 as f64)
    }
}

/// Jukes-Cantor delta-method variance of a corrected distance:
/// `V(d) = p(1-p) / [L · (1 - 4/3·p)²]`. NaN at/above saturation or L == 0.
#[inline]
fn jc_variance(p: f64, sites: f64) -> f64 {
    // NaN p falls through to the arithmetic below, which yields NaN.
    if sites <= 0.0 || p >= JC_SATURATION_THRESHOLD {
        return f64::NAN;
    }
    let denom = 1.0 - (4.0 / 3.0) * p;
    p * (1.0 - p) / (sites * denom * denom)
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


    // === cov_ hardening tests (appended) ===

    // Degenerate pathway: a 2-difference codon pair whose BOTH single-step
    // intermediates are stop codons uses the documented (0.5, 1.5) fallback
    // (pathway_2diff, nei.rs lines 133-135). TGG(Trp) vs TAA(stop) route only
    // through TAG and TGA, both stop codons => valid_paths == 0 => (0.5, 1.5).
    #[test]
    fn cov_nei_two_diff_all_stop_pathway_fallback() {
        let tables = NeiTables::new();
        let idx_tgg = fasta_to_codon_indices(b"TGG", Model::Nei)[0] as usize;
        let idx_taa = fasta_to_codon_indices(b"TAA", Model::Nei)[0] as usize;
        let e = tables.diff_table[idx_tgg * 64 + idx_taa];
        assert!((e.0 - 0.5).abs() < 1e-6, "sd fallback should be 0.5, got {}", e.0);
        assert!((e.1 - 1.5).abs() < 1e-6, "nd fallback should be 1.5, got {}", e.1);
    }

    // 3-difference pathway analysis (pathway_3diff), value derived from first
    // principles. CCC(Pro) vs GGG(Gly): all 6 stop-free pathways yield
    // (sd, nd) = (1, 2) -- each path has exactly one synonymous step (Pro->Pro or
    // Gly->Gly) and two non-synonymous steps -- so the averaged entry is (1, 2).
    #[test]
    fn cov_nei_three_diff_pathway_value() {
        let tables = NeiTables::new();
        let idx_ccc = fasta_to_codon_indices(b"CCC", Model::Nei)[0] as usize;
        let idx_ggg = fasta_to_codon_indices(b"GGG", Model::Nei)[0] as usize;
        let e = tables.diff_table[idx_ccc * 64 + idx_ggg];
        assert!((e.0 - 1.0).abs() < 1e-6, "sd should be 1.0, got {}", e.0);
        assert!((e.1 - 2.0).abs() < 1e-6, "nd should be 2.0, got {}", e.1);
    }

    // Saturation of BOTH pN and pS => both dN and dS are NaN. CCC..CCC vs
    // GGG..GGG: per-codon (sd, nd) = (1, 2) and S = 1 (one 4-fold third position).
    // For k codons pS = k/k = 1 and pN = 2k/2k = 1, both >= JC_SATURATION_THRESHOLD
    // (0.749) => NaN.
    #[test]
    fn cov_nei_saturated_pn_and_ps_give_nan() {
        let (dn, ds) = nei(b"CCCCCC", b"GGGGGG");
        assert!(dn.is_nan(), "saturated pN should give NaN dN, got {}", dn);
        assert!(ds.is_nan(), "saturated pS should give NaN dS, got {}", ds);
    }

    #[test]
    fn cov_nei_unequal_length_compares_overlap_not_misregistered() {
        // Unequal-length inputs (an imperfect / differently-trimmed alignment) must be
        // compared strictly position-by-position over the common prefix, not mis-registered
        // by the chunks_exact fast path. Result must equal truncating the longer sequence.
        let long = b"ATGATGATGATGAAACCCATGATGGGGTTT"; // 10 codons
        let short = b"ATGATGATGATGGGGTTT"; //             6 codons (differs at codons 5,6)
        let trunc = b"ATGATGATGATGAAACCC"; //             long's first 6 codons
        let (dn_u, ds_u) = nei(long, short);
        let (dn_t, ds_t) = nei(trunc, short);
        assert!((dn_u - dn_t).abs() < 1e-12 || (dn_u.is_nan() && dn_t.is_nan()), "dN {} vs {}", dn_u, dn_t);
        assert!((ds_u - ds_t).abs() < 1e-12 || (ds_u.is_nan() && ds_t.is_nan()), "dS {} vs {}", ds_u, ds_t);
        // The real overlap has nonsynonymous divergence — the old code silently reported 0.
        assert!(dn_u > 0.0, "overlap divergence must not be hidden, got dN={}", dn_u);
    }

    #[test]
    fn cov_nei_zero_synonymous_sites_gives_nan_ds() {
        // A pair of only Met (ATG) and Trp (TGG) codons has zero synonymous sites, so dS
        // cannot be estimated: it must be NaN, not a spurious 0.0 (which flows through JC
        // to dS = 0 and reads as strong purifying selection). Matches the Li model, which
        // returns Ks = NaN for the same input. dN stays defined (nonsynonymous sites exist).
        let (dn, ds) = nei(b"ATGATGATG", b"ATGTGGATG");
        assert!(ds.is_nan(), "zero synonymous sites should give NaN dS, got {}", ds);
        assert!(dn.is_finite(), "dN is still defined, got {}", dn);
    }

    // Jukes-Cantor delta-method variance invariants (compute_pair_stats +
    // jc_variance). V(d) = p(1-p) / [L*(1 - 4/3 p)^2] is finite and >= 0 for
    // 0 <= p < 0.749, and NaN at/above saturation.
    #[test]
    fn cov_nei_variance_nonneg_and_nan_at_saturation() {
        let tables = NeiTables::new();
        // Moderate: one synonymous change -> pS = 0.5 (var_ds > 0), pN = 0 (var_dn = 0).
        let a1 = fasta_to_codon_indices(b"ATGGCTGCT", Model::Nei);
        let a2 = fasta_to_codon_indices(b"ATGGCCGCT", Model::Nei);
        let (dn, ds, var_dn, var_ds) = tables.compute_pair_stats(&a1, &a2);
        assert!(dn.abs() < EPSILON, "dN should be 0, got {}", dn);
        assert!(ds > 0.0, "dS should be > 0, got {}", ds);
        assert!(var_ds.is_finite() && var_ds > 0.0, "var_ds must be finite & > 0, got {}", var_ds);
        assert!(var_dn.is_finite() && var_dn >= 0.0, "var_dn must be finite & >= 0, got {}", var_dn);

        // Saturated: pS = pN = 1 -> distances and variances all NaN.
        let s1 = fasta_to_codon_indices(b"CCCCCC", Model::Nei);
        let s2 = fasta_to_codon_indices(b"GGGGGG", Model::Nei);
        let (dn2, ds2, var_dn2, var_ds2) = tables.compute_pair_stats(&s1, &s2);
        assert!(dn2.is_nan() && ds2.is_nan(), "saturated distances must be NaN");
        assert!(var_dn2.is_nan(), "var_dn must be NaN at saturation, got {}", var_dn2);
        assert!(var_ds2.is_nan(), "var_ds must be NaN at saturation, got {}", var_ds2);
    }
}
