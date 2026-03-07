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
                    16 * remap(b1) + 4 * remap(b2) + remap(b3)
                }
            };
            v.push(index);
        }
    }
    v
}
