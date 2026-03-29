# FAQ

## Why is eskaks so much faster?

Three main reasons:

1. **Precomputed lookup tables**: All 4,096 codon pairs have their diffs/sites precomputed at startup. Each pairwise comparison is just table lookups + a final correction formula.

2. **Cache optimization**: The Li model's lookup table (288 KB) fits in L2 cache. Identical codons (~95% in typical alignments) use a separate 1.5 KB table that fits in L1 cache.

3. **Parallel + streaming**: Rayon distributes pairs across CPU cores, and a dedicated writer thread handles I/O without blocking computation.

## Why do I get NaN values?

NaN means the substitution proportion has reached **saturation** (p ≥ 0.749 for Nei-Gojobori, or denominator ≤ 0 for Li). This happens when sequences are too divergent — the correction formula breaks down because there have been so many substitutions that the signal is lost.

**What to do**: This is biologically correct. Very divergent sequences simply can't be reliably compared at the nucleotide level. Consider using protein-level methods instead.

## Which model should I use?

- **Nei-Gojobori**: Simpler, faster, good for quick scans. Use when you want results comparable to most published analyses.
- **Li (1993)**: More accurate because it accounts for the well-known transition/transversion bias. Use when accuracy matters more than speed (though the speed difference is minimal).

## My sequences have internal stop codons

eskaks warns about these automatically. Common causes:
- **Wrong reading frame**: Shift your sequences by 1 or 2 bases
- **Pseudogene**: The gene has been inactivated
- **Frameshift mutation**: An indel has disrupted the reading frame
- **Wrong genetic code**: Try `--genetic-code 2` (mitochondrial) etc.

## Can I use eskaks with non-standard genetic codes?

Yes! Use `--genetic-code <N>` with any of the 20 supported NCBI tables. Run `eskaks --list-codes` to see all options.

## How do I cite eskaks?

See the [CITATION.cff](https://github.com/PathoGenOmics-Lab/eskaks/blob/main/CITATION.cff) file, or use GitHub's "Cite this repository" button.

## Can I use eskaks as a library?

Yes. eskaks exposes a `lib.rs` with public modules. Add it as a dependency in your `Cargo.toml`:

```toml
[dependencies]
eskaks = { git = "https://github.com/PathoGenOmics-Lab/eskaks.git" }
```

Then use `eskaks::models::nei::NeiTables` or `eskaks::models::li::LiTables` directly.
