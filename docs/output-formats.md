---
description: >-
  Every file eskaks writes — the pairwise dN/dS tables, the per-gene pN/pS,
  variants, diversity and MK tables, SVG plots and the interactive HTML report —
  in TSV, CSV or JSON.
---

# Output Formats

`--format tsv|csv|json` applies to **both** subcommands and **every** output mode:
the pairwise dN/dS tables below, and the [`eskaks vcf` tables](#eskaks-vcf-outputs)
further down. The examples in this first section use `eskaks fasta`, but the format
rules (quoting, special-value normalisation) are identical for the VCF path.

## TSV (default)

Tab-separated values. Standard bioinformatics format.

```bash
eskaks fasta genes.fasta -o results
# → results_pairwise_results.tsv
```

```
Seq1	Seq2	dN	dS	dN/dS
gene_A	gene_B	0.052300	0.321400	0.162700
```

## CSV

Comma-separated values. Compatible with Excel and pandas.

```bash
eskaks fasta genes.fasta --format csv -o results
# → results_pairwise_results.csv
```

## JSON

JSON array of objects. Best for programmatic parsing.

```bash
eskaks fasta genes.fasta --format json -o results
# → results_pairwise_results.json
```

```json
[
{"seq1":"gene_A","seq2":"gene_B","dN":0.052300,"dS":0.321400,"dN_dS":0.162700},
{"seq1":"gene_A","seq2":"gene_C","dN":0.089100,"dS":null,"dN_dS":null}
]
```

Special values:
- `NaN` (saturation, or no comparable codons) → `null`
- `Infinity` (dS=0, dN>0) → `null`
- `-0.0` → `0.0`

`--format` applies to **every** output mode, not just the default pairwise table:
`--lineage`, `--group-average`, and `--window-size` also emit a valid JSON array (or
TSV/CSV) with the fields for that mode.

In TSV/CSV, any sequence id or gene name that contains the delimiter, a quote, or a
newline is quoted (RFC 4180), so the columns never shift.

## Output modes

| Mode | Flag | Output file |
|------|------|------------|
| Pairwise (default) | - | `<prefix>_pairwise_results.<ext>` |
| Lineage summary | `--lineage` | `<prefix>_lineage_summary.<ext>` |
| Group average | `--group-average` | `<prefix>_group_avg_dn_ds.<ext>` |
| Sliding window | `--window-size N` | `<prefix>_pairwise_windows.<ext>` |

## SVG Plots (`eskaks fasta`)

Add `--plot` to generate SVG visualizations:

- **Pairwise mode**: dN/dS histogram (`<prefix>_dnds_histogram.svg`)
- **Window mode**: dN/dS along the alignment (`<prefix>_window_plot.svg`)
- **Group mode**: Bar chart with CI (`<prefix>_group_dnds.svg`)
- **Lineage mode**: Bar chart by lineage (`<prefix>_lineage_dnds.svg`)

## `eskaks vcf` outputs

The VCF path writes one always-on table plus opt-in files per flag. All tables
honour `--format tsv|csv|json` and the same quoting / special-value rules as above.

| File | Flag | What it holds |
|------|------|---------------|
| `<prefix>_pnps.<ext>` | (default) | Per-gene pN/pS: counts, fractional N/S sites, the ratio, the neutrality-test p-value and BH/Bonferroni corrections, plus gene coordinates. [Full column list](vcf-analysis.md#output). |
| `<prefix>_variants.<ext>` | `--variants` | One row per coding SNP: position, base change, `S315T`-style amino-acid change, AF, and effect (synonymous / missense / nonsense / stop_loss). [Details](vcf-analysis.md#per-variant-table-variants). |
| `<prefix>_diversity.<ext>` | `--diversity` | Per-gene πN/πS, Watterson θ and Tajima's D (needs a sample size — a multi-sample VCF or `--vcf-list`). [Details](vcf-analysis.md#population-diversity-diversity). |
| `<prefix>_mk.<ext>` | `--mk` | Per-gene McDonald-Kreitman 2×2 table, Neutrality Index, α and a Fisher p/q. [Details](vcf-analysis.md#mcdonald-kreitman-test). |
| `<prefix>_pnps_manhattan.svg`, `<prefix>_pvalue_manhattan.svg` | `--plot` | Manhattan plots of pN/pS and of `−log10(p)` by genome position, significant genes outlined. |
| `<prefix>_report.html` | `--report` | Self-contained [interactive dashboard](vcf-analysis.md#interactive-html-report) — no CDN or internet needed. See the [live example](example-report.md). |

JSON emits a proper array of typed objects (numbers stay numeric, `NA`/undefined
become `null`, `-0.0` normalises to `0.0`), so the tables drop straight into
`pandas.read_json` or `jq`. The `--divergence <FILE>` input is read, not written —
it feeds the report's polymorphism-vs-divergence panel.
