# Output Formats

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

## SVG Plots

Add `--plot` to generate SVG visualizations:

- **Pairwise mode**: dN/dS histogram (`<prefix>_dnds_histogram.svg`)
- **Window mode**: dN/dS along the alignment (`<prefix>_window_plot.svg`)
- **Group mode**: Bar chart with CI (`<prefix>_group_dnds.svg`)
- **Lineage mode**: Bar chart by lineage (`<prefix>_lineage_dnds.svg`)
