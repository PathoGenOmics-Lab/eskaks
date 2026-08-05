---
description: >-
  The Nei-Gojobori (1986) and Li/LPB93 (1993) substitution models eskaks uses for
  dN/dS, plus the Ina (1995) transition-bias-aware site-counting correction.
---

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

Both models handle all 20 NCBI genetic code tables via `--genetic-code` (see
[Genetic codes](genetic-codes.md) for the full list and how to pick one).

## Transition-bias-aware site counting (Ina 1995) { #ina-1995 }

Classic Nei-Gojobori counts synonymous and nonsynonymous **sites** by treating all
nine single-nucleotide changes at a codon as equally likely. Real mutational
spectra are not uniform — most genomes, and *M. tuberculosis* especially, are
strongly **transition-biased** (A↔G, C↔T changes dominate). Because the synonymous
change at a 2-fold degenerate third position is almost always a transition,
equal-rate counting **under-counts synonymous sites** and biases pN/pS downward.

[Ina (1995)](https://doi.org/10.1007/BF00173196) generalised the Nei-Gojobori
site count to weight each candidate change by its relative mutation rate — `κ` for
a transition, `1` for a transversion. eskaks exposes this in the VCF path through
`--kappa`: a 2-fold synonymous site's contribution moves from `1/3` (κ=1) to
`κ/(κ+2)`, 4-fold degenerate sites stay synonymous regardless (κ-invariant), and
each codon still contributes exactly three sites — only the **S/N split** shifts.
`--kappa 1` reproduces the classic equal-rate count bit-for-bit. See
[Mutation-spectrum-aware site counting](vcf-analysis.md#kappa) for the full
derivation and worked direction.
