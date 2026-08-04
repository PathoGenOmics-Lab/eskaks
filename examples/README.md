# Example data

Small, self-contained datasets so you can run eskaks immediately. They are used in
the [getting-started tutorial](../docs/tutorial.md).

## `genes.fasta` for `eskaks fasta` (dN/dS)

Six codon-aligned versions of one gene (180 bp) from six strains. All 15 pairs come
out well under purifying selection (dN/dS ≈ 0.08–0.34, i.e. always < 1).

```bash
eskaks fasta examples/genes.fasta -o first_run
```

## `lineages.fasta` for grouped `eskaks fasta` (`--lineage` / `--group-average`)

The same 180 bp gene from six isolates, but named by **lineage**: `Lineage2`,
`Lineage4`, and `Bovis`, two isolates each. Isolates within a lineage share a block
of substitutions, so they cluster. Use it to compare groups instead of individual
pairs:

```bash
# mean dN/dS between every pair of lineages (and within each)
eskaks fasta examples/lineages.fasta -o lin --group-average

# group by the first letter of the ID instead of splitting on '_':
# Lineage2 + Lineage4 merge into "L", Bovis stays "B"
eskaks fasta examples/lineages.fasta -o lin_fl --group-average --first-letter-lineage
```

## `toy_genome/` for `eskaks vcf` (pN/pS)

A miniature genome (12 genes) with the three inputs `eskaks vcf` needs, plus a
divergence table for the report's reconciliation panel:

| File | What it is |
|---|---|
| `reference.fasta` | the genome sequence (contig `chr1`) |
| `genes.gff3` | gene annotation (CDS features) |
| `variants.vcf` | the SNPs, with `AF` and `DP` in the INFO field |
| `variants_mixed.vcf` | 12 `gene01` SNPs with a mix of `PASS` and `LowQual` FILTERs, to try `--pass-only` |
| `divergence.tsv` | a per-gene dN/dS table (`gene <TAB> dN/dS`) |

```bash
eskaks vcf \
  --ref examples/toy_genome/reference.fasta \
  --gff examples/toy_genome/genes.gff3 \
  --vcf examples/toy_genome/variants.vcf \
  --genetic-code 11 --report --plot \
  --divergence examples/toy_genome/divergence.tsv \
  -o toy_scan
# then open toy_scan_report.html in a browser
```

> These are **synthetic** datasets built only to demonstrate the tool; the numbers
> are not biologically meaningful. Two genes are named `PPE_toy1` / `PE_PGRS_toy2`
> to show how repetitive genes are flagged.
