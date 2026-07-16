//! SplitMix64: a small seeded PRNG for reproducible bootstrap resampling.

/// Minimal seeded PRNG (SplitMix64) for reproducible bootstrap resampling —
/// avoids pulling in the `rand` crate. Not cryptographic.
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `[0, bound)` (modulo bias is negligible for the
    /// gene/site counts used in bootstrap resampling).
    #[inline]
    pub fn below(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }
}

