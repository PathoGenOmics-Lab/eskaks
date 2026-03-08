use crate::models::Model;

/// Sentinel value indicating an invalid/ambiguous codon index.
pub const INVALID_CODON: u8 = 64;

/// Converts a base byte (ASCII) to a compact number (0-4).
/// A=0, C=1, G=2, T/U=3, Other=4 (Ambiguous/N)
#[inline(always)]
pub fn base_to_dna5(b: u8) -> u8 {
    const LUT: [u8; 256] = {
        let mut t = [4; 256];
        t[b'A' as usize] = 0; t[b'a' as usize] = 0;
        t[b'C' as usize] = 1; t[b'c' as usize] = 1;
        t[b'G' as usize] = 2; t[b'g' as usize] = 2;
        t[b'T' as usize] = 3; t[b't' as usize] = 3;
        t[b'U' as usize] = 3; t[b'u' as usize] = 3;
        t
    };
    LUT[b as usize]
}

/// Converts an ASCII byte sequence to a compact DNA5 sequence.
pub fn seq_to_dna5(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().map(|&b| base_to_dna5(b)).collect()
}

/// Converts a compact DNA5 sequence to a list of codon indices.
pub fn dna5_to_codon_indices(dna5_bases: &[u8], model: Model) -> Vec<u8> {
    let mut v = Vec::with_capacity(dna5_bases.len() / 3);
    for chunk in dna5_bases.chunks_exact(3) {
        let (b1, b2, b3) = (chunk[0], chunk[1], chunk[2]);

        if b1 > 3 || b2 > 3 || b3 > 3 {
            v.push(INVALID_CODON);
        } else {
            let index = match model {
                Model::Li => 16 * b1 + 4 * b2 + b3,
                Model::Nei => {
                    let remap = |b: u8| match b { 0 => 2, 1 => 1, 2 => 3, 3 => 0, _ => unreachable!() };
                    16 * remap(b2) + 4 * remap(b1) + remap(b3)
                }
            };
            v.push(index);
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::Model;

    /// Helper: convert an ASCII codon string to a Nei model codon index.
    fn nei_index(codon: &[u8; 3]) -> u8 {
        let dna5 = seq_to_dna5(codon);
        dna5_to_codon_indices(&dna5, Model::Nei)[0]
    }

    /// Helper: convert an ASCII codon string to a Li model codon index.
    fn li_index(codon: &[u8; 3]) -> u8 {
        let dna5 = seq_to_dna5(codon);
        dna5_to_codon_indices(&dna5, Model::Li)[0]
    }

    #[test]
    fn nei_atg_is_met() {
        // ATG encodes Met (M) at AA_ARRAY index 11
        assert_eq!(nei_index(b"ATG"), 11);
    }

    #[test]
    fn nei_tgg_is_trp() {
        // TGG encodes Trp (W) at AA_ARRAY index 51
        assert_eq!(nei_index(b"TGG"), 51);
    }

    #[test]
    fn nei_ttt_is_phe() {
        // TTT encodes Phe (F) at AA_ARRAY index 0
        assert_eq!(nei_index(b"TTT"), 0);
    }

    #[test]
    fn nei_taa_is_stop() {
        // TAA is a stop codon at AA_ARRAY index 34
        assert_eq!(nei_index(b"TAA"), 34);
    }

    #[test]
    fn nei_ggg_is_gly() {
        // GGG encodes Gly (G) at AA_ARRAY index 63
        assert_eq!(nei_index(b"GGG"), 63);
    }

    #[test]
    fn li_aaa_index_zero() {
        // Li model: AAA → 16*0 + 4*0 + 0 = 0
        assert_eq!(li_index(b"AAA"), 0);
    }

    #[test]
    fn li_ttt_index_63() {
        // Li model: TTT → 16*3 + 4*3 + 3 = 63
        assert_eq!(li_index(b"TTT"), 63);
    }

    #[test]
    fn ambiguous_base_gives_invalid_codon() {
        let dna5 = seq_to_dna5(b"ATN");
        let indices = dna5_to_codon_indices(&dna5, Model::Nei);
        assert_eq!(indices[0], INVALID_CODON);
    }

    #[test]
    fn multiple_codons() {
        let dna5 = seq_to_dna5(b"ATGTTT");
        let indices = dna5_to_codon_indices(&dna5, Model::Nei);
        assert_eq!(indices.len(), 2);
        assert_eq!(indices[0], 11); // ATG = Met
        assert_eq!(indices[1], 0);  // TTT = Phe
    }

    #[test]
    fn trailing_bases_ignored() {
        // 7 bases → only 2 complete codons
        let dna5 = seq_to_dna5(b"ATGTTTG");
        let indices = dna5_to_codon_indices(&dna5, Model::Nei);
        assert_eq!(indices.len(), 2);
    }
}
