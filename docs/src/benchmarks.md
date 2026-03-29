# Benchmarks

## Reproducing benchmarks

```bash
make benchmark
```

This runs the full pipeline:
1. Generates synthetic datasets (small/medium/large)
2. Runs eskaks, KaKs_Calculator, BioPython, and PAML yn00
3. Compares accuracy and generates plots

### Requirements

- eskaks (built with `make release`)
- Python 3 with matplotlib, numpy
- [KaKs_Calculator](https://github.com/kullrich/kakscalculator2) (optional)
- [BioPython](https://biopython.org/) (optional)
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) yn00 (optional)

Tools that aren't installed are simply skipped.

## Performance results

Wall-clock time for pairwise dN/dS computation:

| Dataset | eskaks (4t) | KaKs_Calc | PAML yn00 | BioPython |
|---------|-------------|-----------|-----------|-----------|
| 20 seq × 300 bp | 2 ms | 34 ms | 8 ms | 610 ms |
| 100 seq × 3000 bp | 6 ms | 7,703 ms | 697 ms | 111,619 ms |
| 500 seq × 3000 bp | 74 ms | 195,456 ms | — | — |

## Accuracy results

Validated on 20 sequences, 300 codons, 190 pairwise comparisons:

| Comparison | dN R² | dS R² |
|-----------|-------|-------|
| eskaks Li vs KaKs_Calculator LPB | 1.000000 | 1.000000 |
| eskaks Nei vs KaKs_Calculator NG | 0.999397 | 0.995155 |
| eskaks Nei vs BioPython NG86 | 0.998169 | 0.996981 |

The Li model achieves **exact agreement** (R² = 1.0) with KaKs_Calculator's LPB method.
Small Nei differences are due to minor pathway-counting heuristics and are within the inter-tool variation observed between KaKs_Calculator and BioPython.
