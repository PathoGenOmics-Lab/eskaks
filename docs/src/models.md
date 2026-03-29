# Models

eskaks implements two classical models for estimating synonymous (dS) and nonsynonymous (dN) substitution rates.

## Nei-Gojobori (1986)

The default model. A straightforward counting approach:

1. **Count sites**: Each codon position is classified as synonymous or nonsynonymous based on how many single-nucleotide changes at that position would change the amino acid.
2. **Count differences**: For each codon pair, differences are classified as synonymous or nonsynonymous via pathway analysis:
   - 1-position difference: direct classification
   - 2-position differences: average over 2 pathways (excluding stop codon intermediates)
   - 3-position differences: average over 6 pathways (excluding stop codon intermediates)
3. **Jukes-Cantor correction**: Corrects for multiple substitutions at the same site.
   - Formula: `d = -3/4 × ln(1 - 4p/3)` where `p` is the proportion of differences
   - Saturates when `p ≥ 0.749` (returns NaN)

**When to use**: Fast exploratory analyses, large datasets, when simplicity is preferred.

## Li (1993) / LPB93

A more sophisticated model that accounts for transition/transversion bias:

1. **Classify sites**: Each codon position is classified as 0-fold, 2-fold, or 4-fold degenerate.
2. **Count substitutions**: Transitions and transversions are counted separately for each degeneracy class, using the same pathway analysis as Nei-Gojobori.
3. **Kimura two-parameter correction**: Separately corrects for transitions and transversions:
   - `A_k = -0.5 × ln(1 - 2P - Q) + 0.25 × ln(1 - 2Q)`
   - `B_k = -0.5 × ln(1 - 2Q)`
4. **LPB93 formulas**: Combines the corrections across degeneracy classes:
   - `Ka = A₀ + (L₀×B₀ + L₂×B₂) / (L₀ + L₂)`
   - `Ks = B₄ + (L₂×A₂ + L₄×A₄) / (L₂ + L₄)`

**When to use**: More accurate estimates, especially when transition/transversion ratios are unequal (which is almost always the case in real data).

## Comparison

| Aspect | Nei-Gojobori | Li (1993) |
|--------|-------------|-----------|
| Correction | Jukes-Cantor (equal rates) | Kimura 2-parameter (ti/tv) |
| Site classification | Syn/nonsyn | 0-fold/2-fold/4-fold |
| Speed | Slightly faster | Slightly slower |
| Accuracy | Good for similar sequences | Better for divergent sequences |
| Reference tool | KaKs_Calculator NG | KaKs_Calculator LPB |

Both models handle all 20 NCBI genetic code tables via `--genetic-code`.
