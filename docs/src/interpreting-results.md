# Interpreting Results

## dN/dS (ω) ratios

The dN/dS ratio (also called ω or Ka/Ks) measures the balance between nonsynonymous (amino acid-changing) and synonymous (silent) substitutions:

| dN/dS | Interpretation |
|-------|----------------|
| **ω < 1** | **Purifying selection**: Most amino acid changes are deleterious and removed by natural selection. This is by far the most common result for functional genes. |
| **ω ≈ 1** | **Neutral evolution**: Amino acid changes are neither beneficial nor deleterious. May indicate a pseudogene or relaxed constraint. |
| **ω > 1** | **Positive selection**: Amino acid changes are favored. Rare in whole-gene comparisons; more common in specific codons or domains. |

## Common values

- **Housekeeping genes**: ω ≈ 0.01–0.1 (strong purifying selection)
- **Immune genes**: ω ≈ 0.5–2.0 (variable, some under positive selection)
- **Pseudogenes**: ω ≈ 1.0 (no functional constraint)
- **Viral surface proteins**: ω > 1 common in specific sites

## Special values

- **NaN**: The Jukes-Cantor or Kimura correction reached saturation. This means the sequences are too divergent for reliable estimation. Typical when `p ≥ 0.75` (more than 75% of sites differ).
- **0.0 / 0.0 = 0.0**: Identical sequences. No substitutions observed.
- **dN > 0, dS = 0**: All observed changes are nonsynonymous. The ratio is technically infinity, reported as `Inf` (TSV/CSV) or `null` (JSON).

## Caveats

1. **Pairwise ω is a genome-wide average**. Positive selection at specific sites can be masked by purifying selection elsewhere. Use site-level methods (PAML, HyPhy) for finer resolution.

2. **Short sequences**: dN/dS estimates become unreliable with few codons. Use `--min-codons` to filter very short sequences.

3. **Saturation**: Very divergent sequences (>75% divergence) yield NaN. This is biologically correct, the signal is saturated and the estimate is unreliable.

4. **Internal stop codons**: eskaks warns about these. They usually indicate frameshifts, pseudogenes, or incorrect reading frames. Consider excluding these sequences.

5. **Recombination**: dN/dS assumes a single phylogenetic history. Recombination can bias estimates. Consider using sliding windows (`--window-size`) to detect mosaic patterns.

## Sliding window interpretation

Window analysis reveals variation along the alignment:

- **Peaks (ω > 1)**: Potential positive selection hotspots (e.g., surface-exposed domains)
- **Valleys (ω ≈ 0)**: Highly conserved regions (e.g., catalytic sites, structural cores)
- **Noisy windows**: Short windows or few differences produce unreliable estimates. Use windows of at least 50–100 codons.
