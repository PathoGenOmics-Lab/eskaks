//! Which samples carry each ALT allele, as a fixed-width bitset over sample indices.

/// The samples carrying one ALT allele.
///
/// eskaks assumes **haploid** genotypes (the *M. tuberculosis* setting it was written
/// for). That assumption is what makes this worth keeping: two alleles found in the
/// same sample sit on the same molecule by definition, so "do these two SNPs ever
/// occur together?" is a set intersection over named samples, not an inference from
/// allele frequencies, and it needs neither phasing nor reads.
///
/// One `u64` per 64 samples per allele. A 5,000-sample by 50,000-site cohort costs
/// about 31 MB, against the ~1.8 GB such a run already holds, so the dense form is
/// kept: a sparse or hybrid representation would be complexity bought with nothing.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CarrierSet {
    /// Bit `i % 64` of word `i / 64` is set when sample `i` carries the allele.
    words: Vec<u64>,
    /// How many sample indices this set was sized for. `insert` ignores anything at or
    /// above it, so a malformed genotype column cannot silently grow the set.
    n_samples: usize,
}

impl CarrierSet {
    /// An empty set over `n_samples` sample indices.
    pub fn new(n_samples: usize) -> Self {
        Self { words: vec![0u64; n_samples.div_ceil(64)], n_samples }
    }

    /// Record that sample `i` carries the allele. Out-of-range indices are dropped.
    pub fn insert(&mut self, sample: usize) {
        if sample < self.n_samples {
            self.words[sample / 64] |= 1u64 << (sample % 64);
        }
    }

    /// How many samples carry the allele.
    pub fn len(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Does any single sample carry **both** alleles? Under the haploid assumption that
    /// is exactly co-occurrence on one molecule, which is the whole point of keeping
    /// carriers instead of counts. Sets built from different records may have different
    /// word lengths; the shared prefix is the only part where both can be set, since a
    /// sample index beyond the shorter set is absent from it by construction.
    pub fn intersects(&self, other: &Self) -> bool {
        self.words.iter().zip(other.words.iter()).any(|(a, b)| a & b != 0)
    }

    /// Fold `other`'s carriers into this set, widening it if `other` covers more samples.
    pub fn union_with(&mut self, other: &Self) {
        if other.words.len() > self.words.len() {
            self.words.resize(other.words.len(), 0);
        }
        self.n_samples = self.n_samples.max(other.n_samples);
        for (a, b) in self.words.iter_mut().zip(other.words.iter()) {
            *a |= *b;
        }
    }

    /// Heap bytes held by this set, so a run can report what retaining carriers cost
    /// instead of asking the reader to take the estimate on trust.
    pub fn heap_bytes(&self) -> usize {
        self.words.len() * std::mem::size_of::<u64>()
    }

    /// Does sample `i` carry the allele?
    ///
    /// A set built for fewer samples than the index asked about simply does not carry
    /// it, so this is safe to call with the union's sample indices against a narrower
    /// per-allele set.
    pub fn contains(&self, sample: usize) -> bool {
        sample < self.n_samples && ((self.words[sample / 64] >> (sample % 64)) & 1) == 1
    }

    /// The sample indices carrying the allele, ascending.
    ///
    /// Walks set bits (`x & (x - 1)` clears the lowest one) rather than every index, so
    /// enumerating the codons a cohort actually carries costs the number of carriers,
    /// not the number of samples. That distinction is what keeps the per-codon pass
    /// affordable on a cohort of thousands where almost every allele is rare.
    pub fn samples(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(w, &word)| {
            let base = w * 64;
            let mut rest = word;
            std::iter::from_fn(move || {
                if rest == 0 {
                    return None;
                }
                let bit = rest.trailing_zeros() as usize;
                rest &= rest - 1; // clear the lowest set bit
                Some(base + bit)
            })
        })
    }
}

/// Introspection the bitset's own tests assert on (and `is_empty`, which clippy wants
/// beside `len`). The analysis paths need only `insert`, `len`, `intersects`,
/// `union_with`, `contains`, `samples` and `heap_bytes`, so without this these read as
/// dead code in the binary target, which compiles the module a second time.
#[allow(dead_code)]
impl CarrierSet {
    /// How many sample indices this set covers (not how many carry the allele).
    pub fn capacity(&self) -> usize {
        self.n_samples
    }

    /// True when no sample carries the allele.
    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_bitset_records_and_intersects_exactly() {
        let mut a = CarrierSet::new(200);
        let mut b = CarrierSet::new(200);
        for i in [0usize, 63, 64, 199] {
            a.insert(i);
        }
        for i in [1usize, 65, 128] {
            b.insert(i);
        }
        assert_eq!(a.len(), 4);
        assert_eq!(b.len(), 3);
        assert!(a.contains(63) && a.contains(64) && !a.contains(65));
        // Word boundaries are where an off-by-one would hide, so they are checked from
        // both sides: disjoint sets that straddle 64 must NOT report an intersection.
        assert!(!a.intersects(&b), "disjoint sets across word boundaries");
        b.insert(64);
        assert!(a.intersects(&b) && b.intersects(&a), "a shared sample is symmetric");
    }

    #[test]
    fn out_of_range_samples_are_dropped_not_wrapped() {
        // A malformed genotype line must not set a bit belonging to another sample.
        let mut s = CarrierSet::new(10);
        s.insert(10);
        s.insert(usize::MAX);
        assert!(s.is_empty(), "no index outside 0..10 may set a bit");
        assert!(!s.contains(10));
    }

    #[test]
    fn union_widens_and_keeps_both_sides() {
        let mut wide = CarrierSet::new(70);
        wide.insert(69);
        let mut narrow = CarrierSet::new(5);
        narrow.insert(2);
        narrow.union_with(&wide);
        assert_eq!(narrow.capacity(), 70);
        assert!(narrow.contains(2) && narrow.contains(69));
        assert_eq!(narrow.len(), 2);
    }

    #[test]
    fn samples_walks_set_bits_in_order_across_words() {
        // The per-codon pass iterates carriers, not samples, so this iterator has to be
        // exact at the word boundaries and to cost nothing on an empty word.
        let mut s = CarrierSet::new(200);
        for i in [0usize, 63, 64, 65, 128, 199] {
            s.insert(i);
        }
        assert_eq!(s.samples().collect::<Vec<_>>(), vec![0, 63, 64, 65, 128, 199]);
        assert_eq!(s.samples().count(), s.len());
        // An all-zero word must yield nothing rather than index 0 (the `x & (x - 1)`
        // walk would underflow if it ever stepped past zero).
        let empty = CarrierSet::new(200);
        assert_eq!(empty.samples().next(), None);
        // And a full word yields every index in it.
        let mut full = CarrierSet::new(64);
        for i in 0..64 {
            full.insert(i);
        }
        assert_eq!(full.samples().collect::<Vec<_>>(), (0..64).collect::<Vec<_>>());
    }

    #[test]
    fn memory_is_one_word_per_64_samples() {
        // The costing behind the decision to keep the dense form: 5,000 samples is 79
        // words, so 50,000 alleles come to about 31 MB.
        let s = CarrierSet::new(5_000);
        assert_eq!(s.heap_bytes(), 79 * 8);
        assert!((s.heap_bytes() * 50_000) < 33 * 1024 * 1024);
    }
}
