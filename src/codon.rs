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
                    16 * remap(b1) + 4 * remap(b2) + remap(b3)
                }
            };
            v.push(index);
        }
    }
    v
}
