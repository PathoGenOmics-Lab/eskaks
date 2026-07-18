# Models

eskaks implements two classical models for estimating synonymous (dS) and nonsynonymous (dN) substitution rates.

## Nei-Gojobori (1986)

The default model. A straightforward counting approach:

1. **Count sites**: For each codon, its nine single-nucleotide changes are classified as synonymous or nonsynonymous, *excluding* changes that create a stop codon; the codon's three sites are then split in that ratio (`S = 3 × syn/(syn+nonsyn)`, `N = 3 − S`). This is the same exclude-nonsense convention the VCF pN/pS path uses, so both modes agree codon-for-codon.
2. **Count differences**: For each codon pair, differences are classified as synonymous or nonsynonymous via pathway analysis:
   - 1-position difference: direct classification
   - 2-position differences: average over 2 pathways (excluding stop codon intermediates)
   - 3-position differences: average over 6 pathways (excluding stop codon intermediates)
3. **Jukes-Cantor correction** for multiple hits at the same site, where \(p\) is the proportion of differences (\(p_N\) for dN, \(p_S\) for dS):

    \[ d = -\frac{3}{4}\,\ln\!\left(1 - \frac{4}{3}\,p\right) \]

The ratio is then \(\text{dN/dS} = d_N / d_S\): below 1 signals purifying selection, above 1 diversifying/positive selection.

!!! warning "Saturation"
    The Jukes-Cantor correction diverges as \(p \to 0.75\). eskaks returns `NaN` (not an unstable or infinite value) once \(p \ge 0.749\), so a saturated pair is reported as undefined rather than silently wrong.

!!! tip "When to use"
    Fast exploratory analyses, large datasets, and when you want agreement with `KaKs_Calculator`'s NG implementation. It is the default (`--model nei`).

## Li (1993) / LPB93

A more sophisticated model that accounts for transition/transversion bias:

1. **Classify sites**: Each codon position is classified as 0-fold, 2-fold, or 4-fold degenerate.
2. **Count substitutions**: Transitions and transversions are counted separately for each degeneracy class, using the same pathway analysis as Nei-Gojobori.
3. **Kimura two-parameter correction**, separately for transitions (\(P\)) and transversions (\(Q\)) in each degeneracy class \(k\):

    \[ A_k = -\tfrac{1}{2}\ln(1 - 2P - Q) + \tfrac{1}{4}\ln(1 - 2Q), \qquad B_k = -\tfrac{1}{2}\ln(1 - 2Q) \]

4. **LPB93 combination** across the 0-fold (\(L_0\)), 2-fold (\(L_2\)) and 4-fold (\(L_4\)) site classes:

    \[ K_a = A_0 + \frac{L_0 B_0 + L_2 B_2}{L_0 + L_2}, \qquad K_s = B_4 + \frac{L_2 A_2 + L_4 A_4}{L_2 + L_4} \]

!!! tip "When to use"
    More accurate estimates, especially when the transition/transversion ratio is far from 1 (almost always, in real data). Validated at \(R^2 = 1.000\) against `KaKs_Calculator`'s LPB. Select with `--model li`.

## Comparison

| Aspect | Nei-Gojobori | Li (1993) |
|--------|-------------|-----------|
| Correction | Jukes-Cantor (equal rates) | Kimura 2-parameter (ti/tv) |
| Site classification | Syn/nonsyn | 0-fold/2-fold/4-fold |
| Speed | Slightly faster | Slightly slower |
| Accuracy | Good for similar sequences | Better for divergent sequences |
| Reference tool | KaKs_Calculator NG | KaKs_Calculator LPB |

Both models handle all 20 NCBI genetic code tables via `--genetic-code`.
