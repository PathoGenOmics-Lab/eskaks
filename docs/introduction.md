# eskaks

**eskaks** is a high-performance command-line tool for evolutionary rate analysis. It calculates pairwise dN/dS (Ka/Ks) from codon-aligned sequences and per-gene pN/pS from VCF files. It implements the Nei-Gojobori (1986) and Li (1993) models with precomputed lookup tables, achieving over 1,000× speedup compared to existing tools.

> **New to this?** Start with the [**getting-started tutorial**](tutorial.md) — it runs a full analysis on bundled example data with no background assumed. Every unfamiliar term is in the [glossary](glossary.md).

**In one line:** eskaks measures natural selection on genes by comparing how fast *amino-acid-changing* mutations accumulate versus *silent* ones. A ratio below 1 means selection is removing harmful changes (a healthy, constrained gene); a ratio above 1 flags genes where change is being favoured (drug targets, antigens, immune genes).

## Why eskaks?

| Feature | eskaks | KaKs_Calculator | BioPython | PAML yn00 |
|---------|--------|-----------------|-----------|-----------|
| Speed (100 seqs) | **6 ms** | 7,703 ms | 111,619 ms | 697 ms |
| Li model R² | **1.000** | reference | — | — |
| Nei model R² | **0.999** | reference | 0.996 | — |
| Genetic codes | **20 tables** | 1 | 1 | limited |
| JSON output | ✅ | ❌ | ❌ | ❌ |
| Stdin pipe | ✅ | ❌ | ❌ | ❌ |
| Parallel | ✅ (rayon) | ❌ | ❌ | ❌ |

## Key Features

- **Two models**: Nei-Gojobori (1986) with Jukes-Cantor correction, Li (1993)/LPB93 with Kimura two-parameter correction
- **20 NCBI genetic code tables**: Standard, mitochondrial, plastid, and more
- **Output modes**: Pairwise, lineage summary, group average, sliding window
- **Output formats**: TSV, CSV, JSON
- **SVG plots**: Histograms, window plots, group bar charts
- **Pipeline-friendly**: Stdin support, JSON output, non-zero exit on errors

## Citation

If you use eskaks in your research, please cite:

```
Ruiz-Rodriguez P, Coscolla M (2026). eskaks: Fast pairwise dN/dS calculation.
https://github.com/PathoGenOmics-Lab/eskaks
```
