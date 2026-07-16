use crate::codon::INVALID_CODON;

// --- Li model constants ---

/// Minimum denominator for Kimura correction stability.
/// Below this threshold the denominators (1-2p-q, 1-2q) are near-zero, indicating
/// biological saturation; the correction is unreliable and NaN is returned instead.
const LI_EPSILON: f64 = 1e-6;

// --- Genetic code for stop codon detection ---
// The genetic code is now stored per-instance in LiTables.gc,
// allowing alternative translation tables (mitochondrial, etc.).

// --- Internal biological functions ---

fn codon_to_index(cod: &[char; 3]) -> usize {
    let base_to_num = |c: char| match c { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'U' => 3, _ => 4 };
    let n1 = base_to_num(cod[0]);
    let n2 = base_to_num(cod[1]);
    let n3 = base_to_num(cod[2]);
    if n1 == 4 || n2 == 4 || n3 == 4 { INVALID_CODON as usize } else { 16 * n1 + 4 * n2 + n3 }
}

fn decode_codon(n: usize) -> [char; 3] {
    let b1 = n / 16;
    let r1 = n % 16;
    let b2 = r1 / 4;
    let b3 = r1 % 4;
    let map = |x: usize| match x { 0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'X' };
    [map(b1), map(b2), map(b3)]
}

/// Returns the amino acid for a codon index (0-63). Stop codons return b'!'.
fn get_amino_acid(gc: &[u8; 64], idx: usize) -> u8 {
    gc[idx]
}

/// Returns true if the codon is a stop codon.
fn is_stop_codon(gc: &[u8; 64], cod: &[char; 3]) -> bool {
    let idx = codon_to_index(cod);
    idx >= 64 || get_amino_acid(gc, idx) == b'!'
}

/// Classify degeneracy of a codon position (0, 2, or 4-fold).
/// Matches KaKs_Calculator's getCodonClass: counts how many of the 3 alternative
/// bases at this position produce the same amino acid.
/// Returns 0 (0-fold), 1 (2-fold), or 2 (4-fold) as array index.
fn get_codon_class(gc: &[u8; 64], cod: &[char; 3], pos: usize) -> usize {
    let idx = codon_to_index(cod);
    if idx >= 64 { return 0; }
    let aa = get_amino_acid(gc, idx);
    if aa == b'!' { return 0; }

    let bases = ['A', 'C', 'G', 'T'];
    let mut count = 0;
    for &b in &bases {
        if b != cod[pos] {
            let mut alt = *cod;
            alt[pos] = b;
            let alt_idx = codon_to_index(&alt);
            if alt_idx < 64 {
                let alt_aa = get_amino_acid(gc, alt_idx);
                if alt_aa != b'!' && alt_aa == aa {
                    count += 1;
                }
            }
        }
    }
    match count {
        0 => 0,     // 0-fold
        3 => 2,     // 4-fold
        _ => 1,     // 2-fold (1 or 2 synonymous alternatives)
    }
}

/// Classify a nucleotide change as transition ('i') or transversion ('v').
fn is_transition(nt1: char, nt2: char) -> bool {
    matches!((nt1, nt2), ('A','G')|('G','A')|('C','T')|('T','C'))
}

/// Count transitions/transversions for a single nucleotide change between
/// two codons that differ at exactly one position. Classifies the change
/// by the degeneracy class of the differing position and adds 0.5 to
/// the appropriate ti/tv bucket for each codon's classification.
///
/// Special cases from KaKs_Calculator (following Prof. Li's kakstools.py):
/// - CGA<->AGA and CGG<->AGG at position 1: these are Arg->Arg changes
///   that look like transversions (C->A) but are classified as synonymous
///   transitions. In KaKs_Calculator LWL85, CGA<->AGA and CGG<->AGG at pos 0
///   get special handling: Si_temp[class1] += 0.5 and Vi_temp[class2] += 0.5.
fn count_1diff_titv(gc: &[u8; 64], cod1: &[char; 3], cod2: &[char; 3], weight: f64,
                     ti: &mut [f64; 3], tv: &mut [f64; 3]) {
    for pos in 0..3 {
        if cod1[pos] != cod2[pos] {
            let class1 = get_codon_class(gc, cod1, pos);
            let class2 = get_codon_class(gc, cod2, pos);

            // Special case: CGA<->AGA and CGG<->AGG at position 0 (Arg synonymous)
            // In KaKs_Calculator LWL85, this gets Si_temp[class1]+=0.5, Vi_temp[class2]+=0.5.
            // Only apply it when the two codons are actually SYNONYMOUS (and non-stop)
            // under the active genetic code — in some mitochondrial codes AGA is Ser or
            // a stop, so this transversion is nonsynonymous and must NOT get the split.
            let is_arg_special = pos == 0 && (
                (cod1 == &['C','G','A'] && cod2 == &['A','G','A']) ||
                (cod1 == &['A','G','A'] && cod2 == &['C','G','A']) ||
                (cod1 == &['C','G','G'] && cod2 == &['A','G','G']) ||
                (cod1 == &['A','G','G'] && cod2 == &['C','G','G'])
            ) && {
                let a1 = get_amino_acid(gc, codon_to_index(cod1));
                a1 != b'!' && a1 == get_amino_acid(gc, codon_to_index(cod2))
            };

            if is_arg_special {
                ti[class1] += 0.5 * weight;
                tv[class2] += 0.5 * weight;
            } else if is_transition(cod1[pos], cod2[pos]) {
                ti[class1] += 0.5 * weight;
                ti[class2] += 0.5 * weight;
            } else {
                tv[class1] += 0.5 * weight;
                tv[class2] += 0.5 * weight;
            }
            break; // only one difference
        }
    }
}

// --- LiTables: precomputed AoS lookup tables for cache-friendly access ---

/// Per-codon-pair precomputed data, stored contiguously for cache locality.
/// Each entry is 72 bytes (9 × f64), fitting in 1-2 cache lines.
#[repr(C)]
#[derive(Clone, Copy)]
struct CodonPairData {
    l: [f64; 3],
    ti: [f64; 3],
    tv: [f64; 3],
}

impl Default for CodonPairData {
    fn default() -> Self {
        CodonPairData { l: [0.0; 3], ti: [0.0; 3], tv: [0.0; 3] }
    }
}

/// Precomputed lookup tables for the Li (1993) / Pamilo-Bianchi (1993) model.
/// Matches KaKs_Calculator 2.0 LPB method exactly:
/// - Sites (L) counted from original codons only (averaged)
/// - Substitutions counted via equal-weight pathway analysis excluding stop codons
/// - Kimura two-parameter correction: A_k = -0.5*ln(1-2P-Q) + 0.25*ln(1-2Q),
///   B_k = -0.5*ln(1-2Q)
///
/// Uses AoS layout: one contiguous [CodonPairData; 4096] for cache-friendly access.
/// `same_l` is a compact 64-entry table (1.5KB) for identical-codon fast path.
pub struct LiTables {
    data: [CodonPairData; 4096],
    same_l: [[f64; 3]; 64],
    /// The genetic code used for this instance.
    gc: [u8; 64],
}

#[inline(always)]
fn flat_idx(i: usize, j: usize) -> usize {
    i * 64 + j
}

/// Accumulate a single codon pair's data into running sums (different codons).
#[inline(always)]
fn accumulate(d: &CodonPairData, l_sum: &mut [f64; 3], ti_sum: &mut [f64; 3], tv_sum: &mut [f64; 3]) {
    l_sum[0] += d.l[0]; l_sum[1] += d.l[1]; l_sum[2] += d.l[2];
    ti_sum[0] += d.ti[0]; ti_sum[1] += d.ti[1]; ti_sum[2] += d.ti[2];
    tv_sum[0] += d.tv[0]; tv_sum[1] += d.tv[1]; tv_sum[2] += d.tv[2];
}

/// Fast path for identical codons: only l[] values are non-zero (ti=tv=0).
/// The same_l table (1.5KB) fits in L1 cache, avoiding the 288KB main table.
#[inline(always)]
fn accumulate_same_l(l_vals: &[f64; 3], l_sum: &mut [f64; 3]) {
    l_sum[0] += l_vals[0]; l_sum[1] += l_vals[1]; l_sum[2] += l_vals[2];
}

impl LiTables {
    /// Builds the lookup tables matching KaKs_Calculator 2.0 LPB methodology.
    ///
    /// For each codon pair (i, j):
    /// 1. L[class] = average of both codons' site classifications (all 3 positions)
    /// 2. Si/Vi = transitions/transversions via pathway analysis:
    ///    - 0 diff: nothing
    ///    - 1 diff: direct classification
    ///    - 2 diff: enumerate 2 pathways, exclude stop intermediates, equal weight
    ///    - 3 diff: enumerate 6 pathways, exclude stop intermediates, equal weight
    ///
    /// Build lookup tables using the standard genetic code (table 1).
    #[allow(dead_code)]
    pub fn new() -> Box<Self> {
        use crate::genetic_code;
        let gc = genetic_code::get_table(1).unwrap();
        Self::with_genetic_code(&gc.aa_table)
    }

    /// Build lookup tables for any NCBI genetic code.
    /// The input table uses `b'*'` for stop codons (NCBI convention);
    /// internally we convert to `b'!'` for Li model compatibility.
    pub fn with_genetic_code(gc_table: &[u8; 64]) -> Box<Self> {
        // SAFETY: LiTables is 256 KB+, too large for stack. We allocate on the heap
        // and zero-initialize. All fields (data: [CodonPairData; 4096], same_l, gc)
        // are valid when zero-filled (f64 zero, [u8] zero). Layout matches LiTables exactly.
        let mut tables: Box<LiTables> = unsafe {
            let layout = std::alloc::Layout::new::<LiTables>();
            let ptr = std::alloc::alloc_zeroed(layout) as *mut LiTables;
            if ptr.is_null() { std::alloc::handle_alloc_error(layout); }
            Box::from_raw(ptr)
        };
        // Convert * → ! for stop codons (Li convention)
        tables.gc = *gc_table;
        for aa in tables.gc.iter_mut() {
            if *aa == b'*' { *aa = b'!'; }
        }
        let gc = &tables.gc as *const [u8; 64];
        // SAFETY: gc points into tables which lives for the duration of this function
        // and the returned Box. All reads happen within this scope.
        let gc = unsafe { &*gc };

        for i in 0..64usize {
            for j in i..64usize {
                let cod1 = decode_codon(i);
                let cod2 = decode_codon(j);

                // Step 1: Count sites from original codons (KaKs_Calculator approach)
                // Each codon contributes its classification; average of both (each × 0.5)
                let mut l_accum = [0.0; 3];
                for pos in 0..3 {
                    let c1 = get_codon_class(gc, &cod1, pos);
                    let c2 = get_codon_class(gc, &cod2, pos);
                    l_accum[c1] += 0.5;
                    l_accum[c2] += 0.5;
                }

                // Step 2: Count substitutions via pathway analysis
                let mut ti_accum = [0.0; 3];
                let mut tv_accum = [0.0; 3];

                let mut diff_positions = [0usize; 3];
                let mut ndiff = 0;
                for pos in 0..3 {
                    if cod1[pos] != cod2[pos] {
                        diff_positions[ndiff] = pos;
                        ndiff += 1;
                    }
                }

                match ndiff {
                    0 => {} // identical codons: no substitutions
                    1 => {
                        count_1diff_titv(gc, &cod1, &cod2, 1.0, &mut ti_accum, &mut tv_accum);
                    }
                    2 => {
                        // 2 pathways: change pos[0] first or pos[1] first
                        let d0 = diff_positions[0];
                        let d1 = diff_positions[1];

                        // Pathway 1: cod1 -> int1 -> cod2 (change d0 first)
                        let mut int1 = cod1;
                        int1[d0] = cod2[d0];

                        // Pathway 2: cod1 -> int2 -> cod2 (change d1 first)
                        let mut int2 = cod1;
                        int2[d1] = cod2[d1];

                        let p1_valid = !is_stop_codon(gc, &int1);
                        let p2_valid = !is_stop_codon(gc, &int2);

                        let valid_count = p1_valid as u32 + p2_valid as u32;
                        if valid_count > 0 {
                            let w = 1.0 / valid_count as f64;
                            if p1_valid {
                                count_1diff_titv(gc, &cod1, &int1, w, &mut ti_accum, &mut tv_accum);
                                count_1diff_titv(gc, &int1, &cod2, w, &mut ti_accum, &mut tv_accum);
                            }
                            if p2_valid {
                                count_1diff_titv(gc, &cod1, &int2, w, &mut ti_accum, &mut tv_accum);
                                count_1diff_titv(gc, &int2, &cod2, w, &mut ti_accum, &mut tv_accum);
                            }
                        }
                    }
                    3 => {
                        // 6 pathways: all permutations of (0,1,2)
                        let perms: [[usize; 3]; 6] = [
                            [0, 1, 2], [0, 2, 1], [1, 0, 2],
                            [1, 2, 0], [2, 0, 1], [2, 1, 0],
                        ];

                        let mut valid_count = 0u32;
                        let mut path_valid = [false; 6];

                        for (pi, perm) in perms.iter().enumerate() {
                            // cod1 -> int1 -> int2 -> cod2
                            let mut int1 = cod1;
                            int1[perm[0]] = cod2[perm[0]];
                            let mut int2 = int1;
                            int2[perm[1]] = cod2[perm[1]];

                            if !is_stop_codon(gc, &int1) && !is_stop_codon(gc, &int2) {
                                path_valid[pi] = true;
                                valid_count += 1;
                            }
                        }

                        if valid_count > 0 {
                            let w = 1.0 / valid_count as f64;
                            for (pi, perm) in perms.iter().enumerate() {
                                if !path_valid[pi] { continue; }
                                let mut int1 = cod1;
                                int1[perm[0]] = cod2[perm[0]];
                                let mut int2 = int1;
                                int2[perm[1]] = cod2[perm[1]];

                                count_1diff_titv(gc, &cod1, &int1, w, &mut ti_accum, &mut tv_accum);
                                count_1diff_titv(gc, &int1, &int2, w, &mut ti_accum, &mut tv_accum);
                                count_1diff_titv(gc, &int2, &cod2, w, &mut ti_accum, &mut tv_accum);
                            }
                        }
                    }
                    _ => unreachable!(),
                }

                let entry = CodonPairData { l: l_accum, ti: ti_accum, tv: tv_accum };
                tables.data[flat_idx(i, j)] = entry;
                tables.data[flat_idx(j, i)] = entry;
            }
        }

        // Build L1-cache-resident same-codon table from diagonal of main table
        for i in 0..64 {
            tables.same_l[i] = tables.data[flat_idx(i, i)].l;
        }

        tables
    }

    /// Computes Ka and Ks for a pair of sequences encoded as codon indices.
    /// Uses AoS layout + chunk-of-4 processing + unsafe unchecked indexing
    /// with same-codon fast path: ~95% of codons are identical in typical
    /// alignments and only need the L1-cache-resident same_l table (1.5KB)
    /// instead of the full 288KB data table.
    #[inline]
    pub fn compute_pair(&self, codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
        // Compare strictly position-by-position over the common length; see the note in
        // nei.rs::compute_pair_stats. Without this clamp the chunks_exact(4) fast path
        // mis-registers codons for unequal-length inputs, yielding a wrong Ka/Ks.
        let n = codon_indices1.len().min(codon_indices2.len());
        let codon_indices1 = &codon_indices1[..n];
        let codon_indices2 = &codon_indices2[..n];

        let mut l_sum = [0.0; 3];
        let mut ti_sum = [0.0; 3];
        let mut tv_sum = [0.0; 3];

        let data = &self.data;
        let same_l = &self.same_l;

        // Process 4 codons at a time
        let chunks1 = codon_indices1.chunks_exact(4);
        let chunks2 = codon_indices2.chunks_exact(4);
        let rem1 = chunks1.remainder();
        let rem2 = chunks2.remainder();

        for (ch1, ch2) in chunks1.zip(chunks2) {
            let a1 = ch1[0] as usize; let a2 = ch1[1] as usize;
            let a3 = ch1[2] as usize; let a4 = ch1[3] as usize;
            let b1 = ch2[0] as usize; let b2 = ch2[1] as usize;
            let b3 = ch2[2] as usize; let b4 = ch2[3] as usize;

            // Chunk-level all-identical check: with ~95% identical codons,
            // ~81% of chunks have all 4 pairs identical (0.95^4).
            // These use only the L1-resident same_l table (1.5KB),
            // completely skipping the 288KB data table.
            if ((a1 ^ b1) | (a2 ^ b2) | (a3 ^ b3) | (a4 ^ b4)) == 0
                && (a1 | a2 | a3 | a4) < INVALID_CODON as usize {
                // SAFETY: all identical AND all < 64 (INVALID_CODON), so indices are in bounds
                unsafe {
                    accumulate_same_l(same_l.get_unchecked(a1), &mut l_sum);
                    accumulate_same_l(same_l.get_unchecked(a2), &mut l_sum);
                    accumulate_same_l(same_l.get_unchecked(a3), &mut l_sum);
                    accumulate_same_l(same_l.get_unchecked(a4), &mut l_sum);
                }
            } else if (a1 | a2 | a3 | a4 | b1 | b2 | b3 | b4) < INVALID_CODON as usize {
                // Some differ, all valid: use full 288KB table
                // SAFETY: all indices verified < 64, so flat index < 4096
                unsafe {
                    accumulate(data.get_unchecked(a1 * 64 + b1), &mut l_sum, &mut ti_sum, &mut tv_sum);
                    accumulate(data.get_unchecked(a2 * 64 + b2), &mut l_sum, &mut ti_sum, &mut tv_sum);
                    accumulate(data.get_unchecked(a3 * 64 + b3), &mut l_sum, &mut ti_sum, &mut tv_sum);
                    accumulate(data.get_unchecked(a4 * 64 + b4), &mut l_sum, &mut ti_sum, &mut tv_sum);
                }
            } else {
                // Some invalid codons: check individually
                for (&c1, &c2) in ch1.iter().zip(ch2.iter()) {
                    let n1 = c1 as usize;
                    let n2 = c2 as usize;
                    if n1 >= INVALID_CODON as usize || n2 >= INVALID_CODON as usize { continue; }
                    // SAFETY: n1 < 64 && n2 < 64 guarantees flat < 4096
                    unsafe { accumulate(data.get_unchecked(n1 * 64 + n2), &mut l_sum, &mut ti_sum, &mut tv_sum); }
                }
            }
        }

        // Handle remainder codons (0-3)
        for (&c1, &c2) in rem1.iter().zip(rem2.iter()) {
            let n1 = c1 as usize;
            let n2 = c2 as usize;
            if n1 >= INVALID_CODON as usize || n2 >= INVALID_CODON as usize { continue; }
            // SAFETY: n1 < 64 && n2 < 64 guarantees flat < 4096
            unsafe { accumulate(data.get_unchecked(n1 * 64 + n2), &mut l_sum, &mut ti_sum, &mut tv_sum); }
        }

        // Kimura two-parameter correction per degeneracy class
        let mut p_k = [0.0; 3];
        let mut q_k = [0.0; 3];
        for k in 0..3 {
            if l_sum[k] > 0.0 {
                p_k[k] = ti_sum[k] / l_sum[k];
                q_k[k] = tv_sum[k] / l_sum[k];
            }
        }

        let mut a_k = [0.0f64; 3];
        let mut b_k = [0.0f64; 3];
        for k in 0..3 {
            if l_sum[k] == 0.0 { continue; }
            let denom_a = 1.0 - 2.0 * p_k[k] - q_k[k];
            let denom_b = 1.0 - 2.0 * q_k[k];
            // Kimura two-parameter correction (Li 1993):
            //   A_k = -0.5*ln(1-2P-Q) + 0.25*ln(1-2Q)   (needs BOTH terms finite)
            //   B_k = -0.5*ln(1-2Q)                      (depends only on Q)
            // Gate each independently: when transitions saturate (denom_a ≤ ε) but
            // transversions do not, B_k is still recoverable and must not be zeroed.
            if denom_b > LI_EPSILON {
                let b_val = -0.5 * denom_b.ln();
                b_k[k] = if b_val >= 0.0 { b_val } else { 0.0 };
                if denom_a > LI_EPSILON {
                    let a_val = -0.5 * denom_a.ln() + 0.25 * denom_b.ln();
                    a_k[k] = if a_val >= 0.0 { a_val } else { 0.0 };
                } else {
                    // Transitions saturated (1-2P-Q ≤ 0): the A_k Kimura correction is
                    // undefined. Propagate NaN so Ks does not silently collapse to 0 and
                    // report a spurious dN/dS = ∞ (the module doc promises NaN here).
                    a_k[k] = f64::NAN;
                }
            } else {
                // Transversions saturated (1-2Q ≤ 0): both A_k and B_k are undefined.
                a_k[k] = f64::NAN;
                b_k[k] = f64::NAN;
            }
        }

        // LPB93 Ka/Ks formulas (Li 1993, Pamilo & Bianchi 1993)
        let l0 = l_sum[0]; let l2 = l_sum[1]; let l4 = l_sum[2];

        // Ka = A₀ + (L₀·B₀ + L₂·B₂) / (L₀ + L₂)
        let mut ka_val = f64::NAN;
        if (l0 + l2) > 0.0 {
            ka_val = a_k[0] + (l0 * b_k[0] + l2 * b_k[1]) / (l0 + l2);
        }

        // Ks = B₄ + (L₂·A₂ + L₄·A₄) / (L₂ + L₄)
        let mut ks_val = f64::NAN;
        if (l2 + l4) > 0.0 {
            ks_val = b_k[2] + (l2 * a_k[1] + l4 * a_k[2]) / (l2 + l4);
        }

        if ka_val.is_finite() && ka_val < 0.0 { ka_val = 0.0; }
        if ks_val.is_finite() && ks_val < 0.0 { ks_val = 0.0; }

        (ka_val, ks_val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codon::fasta_to_codon_indices;
    use crate::models::Model;

    const EPSILON: f64 = 1e-4;

    fn li(seq1: &[u8], seq2: &[u8]) -> (f64, f64) {
        let idx1 = fasta_to_codon_indices(seq1, Model::Li);
        let idx2 = fasta_to_codon_indices(seq2, Model::Li);
        LiTables::new().compute_pair(&idx1, &idx2)
    }

    #[test]
    fn li_identical_seqs_give_zero() {
        let (dn, ds) = li(b"ATGGCTGCT", b"ATGGCTGCT");
        assert!((dn).abs() < EPSILON, "dN for identical seqs should be 0, got {}", dn);
        assert!((ds).abs() < EPSILON, "dS for identical seqs should be 0, got {}", ds);
    }

    #[test]
    fn li_unequal_length_compares_overlap_not_misregistered() {
        // See the Nei twin: unequal-length inputs must be compared over the common prefix,
        // not mis-registered by the chunked fast path.
        let long = b"ATGATGATGATGAAACCCATGATGGGGTTT"; // 10 codons
        let short = b"ATGATGATGATGGGGTTT"; //             6 codons
        let trunc = b"ATGATGATGATGAAACCC"; //             long's first 6 codons
        let (dn_u, ds_u) = li(long, short);
        let (dn_t, ds_t) = li(trunc, short);
        assert!((dn_u - dn_t).abs() < 1e-12 || (dn_u.is_nan() && dn_t.is_nan()), "dN {} vs {}", dn_u, dn_t);
        assert!((ds_u - ds_t).abs() < 1e-12 || (ds_u.is_nan() && ds_t.is_nan()), "dS {} vs {}", ds_u, ds_t);
        assert!(dn_u > 0.0, "overlap divergence must not be hidden, got dN={}", dn_u);
    }

    #[test]
    fn li_one_synonymous_change_dn_is_zero() {
        // 4 codons (12bp): ATGGCTGCTGCT vs ATGGCCGCTGCT — GCT→GCC (Ala→Ala) at codon 2
        let (dn, ds) = li(b"ATGGCTGCTGCT", b"ATGGCCGCTGCT");
        assert!((dn).abs() < EPSILON, "dN should be 0 for purely synonymous change, got {}", dn);
        assert!(ds > 0.0, "dS should be > 0 for synonymous change, got {}", ds);
    }

    #[test]
    fn li_one_nonsynonymous_change_ds_is_zero() {
        // ATGGCTGCT vs ATGATTGCT: GCT→ATT (Ala→Ile), 2 nucleotide differences
        // but both are nonsynonymous via pathway analysis
        let (dn, ds) = li(b"ATGGCTGCTGCTGCT", b"ATGATTGCTGCTGCT");
        assert!(ds < 0.01, "dS should be ~0 for nonsynonymous change, got {}", ds);
        assert!(dn > 0.0, "dN should be > 0 for nonsynonymous change, got {}", dn);
    }

    #[test]
    fn li_symmetry() {
        // li(A,B) == li(B,A) — use 4-codon sequences to avoid Kimura saturation
        let (dn_ab, ds_ab) = li(b"ATGGCTGCTGCT", b"ATGATTGCTGCT");
        let (dn_ba, ds_ba) = li(b"ATGATTGCTGCT", b"ATGGCTGCTGCT");
        assert!((dn_ab - dn_ba).abs() < EPSILON, "dN not symmetric: {} vs {}", dn_ab, dn_ba);
        assert!((ds_ab - ds_ba).abs() < EPSILON, "dS not symmetric: {} vs {}", ds_ab, ds_ba);
    }

    #[test]
    fn li_rna_input_same_as_dna() {
        // U→T substitution should give same results — use 4-codon sequences
        let (dn_dna, ds_dna) = li(b"ATGGCTGCTGCT", b"ATGATTGCTGCT");
        let (dn_rna, ds_rna) = li(b"AUGGCUGCUGCU", b"AUGAUUGCUGCU");
        assert!((dn_dna - dn_rna).abs() < EPSILON, "RNA dN differs from DNA: {} vs {}", dn_dna, dn_rna);
        assert!((ds_dna - ds_rna).abs() < EPSILON, "RNA dS differs from DNA: {} vs {}", ds_dna, ds_rna);
    }


    // === cov_ hardening tests (appended) ===

    // Coverage of CodonPairData::default() (li.rs lines 141-143): the Default
    // impl body sets l/ti/tv each to [0.0; 3]. A zero entry must also be a no-op
    // when accumulated (conservation).
    #[test]
    fn cov_default_codon_pair_data_is_zero() {
        let d = CodonPairData::default();
        assert_eq!(d.l, [0.0; 3]);
        assert_eq!(d.ti, [0.0; 3]);
        assert_eq!(d.tv, [0.0; 3]);

        let mut l_sum = [0.0f64; 3];
        let mut ti_sum = [0.0f64; 3];
        let mut tv_sum = [0.0f64; 3];
        accumulate(&d, &mut l_sum, &mut ti_sum, &mut tv_sum);
        assert_eq!(l_sum, [0.0; 3]);
        assert_eq!(ti_sum, [0.0; 3]);
        assert_eq!(tv_sum, [0.0; 3]);
    }

    // Invariant: an alignment whose only differences are 4-fold synonymous
    // substitutions gives dN == 0 and dS > 0. Such substitutions land only in the
    // 4-fold (class-2) counters, which feed Ks. LPB reduction:
    //   Ka = A0 + (L0*B0 + L2*B2)/(L0+L2); with no 0-fold/2-fold substitutions
    //   A0 = B0 = B2 = 0  =>  Ka = 0.
    #[test]
    fn cov_li_all_synonymous_fourfold_dn_is_zero() {
        // Met Ala Gly Ala Gly Ala; codon 2 GCT->GCC and codon 5 GGT->GGC
        // (both synonymous, 4-fold, third position). No non-synonymous difference.
        let (dn, ds) = li(b"ATGGCTGGTGCTGGTGCT", b"ATGGCCGGTGCTGGCGCT");
        assert!(dn.abs() < EPSILON, "dN must be 0 for purely 4-fold synonymous diffs, got {}", dn);
        assert!(ds > 0.0, "dS must be > 0 for synonymous diffs, got {}", ds);
    }

    // Invariant: a single 0-fold non-synonymous substitution gives dS == 0 and
    // dN > 0. No substitution in the 2-fold/4-fold classes => A1 = A2 = B2 = 0 =>
    //   Ks = B4 + (L2*A2 + L4*A4)/(L2+L4) = 0.
    #[test]
    fn cov_li_zerofold_nonsynonymous_ds_is_zero() {
        // codon 2 GCT(Ala)->GAT(Asp): single change at the middle base (C->A), 0-fold.
        let (dn, ds) = li(b"ATGGCTGCTGCT", b"ATGGATGCTGCT");
        assert!(ds.abs() < EPSILON, "dS must be 0 for a 0-fold non-synonymous change, got {}", ds);
        assert!(dn > 0.0, "dN must be > 0 for a non-synonymous change, got {}", dn);
    }

    // Invariant: symmetry d(A,B) == d(B,A) and non-negativity/finiteness for a
    // mixed (synonymous + non-synonymous) alignment. The estimator is symmetric
    // in its two sequences and clamps negative corrected distances to 0.
    #[test]
    fn cov_li_symmetry_and_nonnegative_mixed() {
        // codon 2 GCT->GCC (synonymous), codon 3 GGT(Gly)->GAT(Asp) (non-synonymous).
        let (dn_ab, ds_ab) = li(b"ATGGCTGGTGCT", b"ATGGCCGATGCT");
        let (dn_ba, ds_ba) = li(b"ATGGCCGATGCT", b"ATGGCTGGTGCT");
        assert!((dn_ab - dn_ba).abs() < EPSILON, "dN not symmetric: {} vs {}", dn_ab, dn_ba);
        assert!((ds_ab - ds_ba).abs() < EPSILON, "dS not symmetric: {} vs {}", ds_ab, ds_ba);
        assert!(dn_ab.is_finite() && dn_ab >= 0.0, "dN must be finite & >= 0, got {}", dn_ab);
        assert!(ds_ab.is_finite() && ds_ab >= 0.0, "dS must be finite & >= 0, got {}", ds_ab);
        assert!(dn_ab > 0.0, "dN must be > 0 (non-synonymous change present), got {}", dn_ab);
        assert!(ds_ab > 0.0, "dS must be > 0 (synonymous change present), got {}", ds_ab);
    }

    // Branch coverage + invariant: when transversions saturate the Kimura
    // denominator (1 - 2Q <= LI_EPSILON), the correction is suppressed (A_k/B_k
    // stay 0) rather than producing NaN or a negative value. Three codons all
    // GCT(Ala)->GCA(Ala): a synonymous transversion at every 4-fold position =>
    // Q4 = 3/3 = 1 => denom_b = 1 - 2Q = -1 <= LI_EPSILON => B4 = A4 = 0 => dS = 0.
    // No 0-fold/2-fold substitution => dN = 0.
    #[test]
    fn cov_li_saturation_returns_nan_not_zero() {
        // All-transversion synonymous divergence saturates the 4-fold class
        // (1-2Q ≤ 0), so the Kimura correction for dS is undefined and must be NaN,
        // NOT silently collapse to 0 (which would report a spurious dN/dS = ∞).
        let (dn, ds) = li(b"GCTGCTGCT", b"GCAGCAGCA");
        assert!(ds.is_nan(), "dS must be NaN under transversion saturation, got {}", ds);
        // No non-synonymous change at the (unsaturated) 0-/2-fold positions -> dN = 0.
        assert!(dn.is_finite() && dn.abs() < EPSILON, "dN must be finite 0, got {}", dn);
    }

    #[test]
    fn cov_li_synonymous_transition_saturation_does_not_give_infinite_ratio() {
        // Regression: heavily transition-skewed synonymous divergence used to zero the
        // A_k term and collapse Ks to exactly 0, giving dN/dS = ∞ (false positive
        // selection). Ks must now be NaN (saturated), so the ratio is not a finite > 1.
        let a: Vec<u8> = b"GCT".iter().cycle().take(60).copied().collect::<Vec<_>>()
            .into_iter().chain(b"TTT".iter().cycle().take(60).copied()).collect();
        let b: Vec<u8> = b"GCC".iter().cycle().take(60).copied().collect::<Vec<_>>()
            .into_iter()
            .chain(b"CTT".iter().cycle().take(12).copied())
            .chain(b"TTT".iter().cycle().take(48).copied())
            .collect();
        let (dn, ds) = li(&a, &b);
        assert!(dn.is_finite() && dn >= 0.0, "dN finite, got {}", dn);
        assert!(!(ds.is_finite() && ds == 0.0), "dS must not silently be 0 under saturation, got {}", ds);
    }
}
