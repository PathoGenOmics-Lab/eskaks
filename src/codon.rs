use crate::models::Model;

/// Sentinel value indicating an invalid/ambiguous codon index.
pub const INVALID_CODON: u8 = 64;

/// Converts a base byte (ASCII) to a compact number (0-4).
/// A=0, C=1, G=2, T/U=3, Other=4 (Ambiguous/N)
#[inline(always)]
fn base_to_dna5(b: u8) -> u8 {
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

/// Converts FASTA bytes directly to codon indices, skipping the DNA5 intermediate.
/// This avoids allocating a full-length DNA5 Vec (saves L bytes per sequence).
pub fn fasta_to_codon_indices(fasta_bytes: &[u8], model: Model) -> Vec<u8> {
    let mut v = Vec::with_capacity(fasta_bytes.len() / 3);
    for chunk in fasta_bytes.chunks_exact(3) {
        let b1 = base_to_dna5(chunk[0]);
        let b2 = base_to_dna5(chunk[1]);
        let b3 = base_to_dna5(chunk[2]);

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

/// Counts gap characters ('-' and '.') in raw FASTA bytes.
pub fn count_gaps(fasta_bytes: &[u8]) -> usize {
    fasta_bytes.iter().filter(|&&b| b == b'-' || b == b'.').count()
}

/// Extracts the grouping key from a sequence ID.
/// If `by_first_letter` is true, returns the first character (UTF-8 aware);
/// otherwise returns the part before the first '_', or the full ID if no '_' is present.
pub fn extract_group_key(id: &str, by_first_letter: bool) -> &str {
    if by_first_letter {
        &id[..id.chars().next().map(|c| c.len_utf8()).unwrap_or(0)]
    } else {
        id.split('_').next().unwrap_or(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::Model;

    /// Helper: convert an ASCII codon string to a Nei model codon index.
    fn nei_index(codon: &[u8; 3]) -> u8 {
        fasta_to_codon_indices(codon, Model::Nei)[0]
    }

    /// Helper: convert an ASCII codon string to a Li model codon index.
    fn li_index(codon: &[u8; 3]) -> u8 {
        fasta_to_codon_indices(codon, Model::Li)[0]
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
        let indices = fasta_to_codon_indices(b"ATN", Model::Nei);
        assert_eq!(indices[0], INVALID_CODON);
    }

    #[test]
    fn multiple_codons() {
        let indices = fasta_to_codon_indices(b"ATGTTT", Model::Nei);
        assert_eq!(indices.len(), 2);
        assert_eq!(indices[0], 11); // ATG = Met
        assert_eq!(indices[1], 0);  // TTT = Phe
    }

    #[test]
    fn trailing_bases_ignored() {
        // 7 bases → only 2 complete codons
        let indices = fasta_to_codon_indices(b"ATGTTTG", Model::Nei);
        assert_eq!(indices.len(), 2);
    }

    #[test]
    fn cov_extract_group_key_first_letter_ascii() {
        // by_first_letter = true returns the first UTF-8 character (ASCII: 1 byte).
        assert_eq!(extract_group_key("A_foo_bar", true), "A");
        assert_eq!(extract_group_key("XYZ", true), "X");
    }

    #[test]
    fn cov_extract_group_key_first_letter_multibyte() {
        // Multi-byte first char: len_utf8() drives the slice end so it lands on a
        // char boundary. U+00E9 (e-acute) = 2 bytes, U+4E2D (CJK) = 3 bytes.
        assert_eq!(extract_group_key("\u{e9}mile", true), "\u{e9}");
        assert_eq!(extract_group_key("\u{4e2d}x", true), "\u{4e2d}");
    }

    #[test]
    fn cov_extract_group_key_first_letter_empty() {
        // Empty id: chars().next() is None => len_utf8 defaults to 0 => &id[..0] == "".
        assert_eq!(extract_group_key("", true), "");
    }

    #[test]
    fn cov_extract_group_key_by_underscore() {
        // by_first_letter = false returns the prefix before the first '_'.
        assert_eq!(extract_group_key("grp_sample1", false), "grp");
        // No '_' present => split yields the whole string as its single element.
        assert_eq!(extract_group_key("wholeid", false), "wholeid");
        assert_eq!(extract_group_key("", false), "");
    }
}
