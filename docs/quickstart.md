# Quick Start

## Pairwise dN/dS (FASTA)

```bash
# Compute pairwise dN/dS with Nei model (default)
eskaks fasta aligned_genes.fasta -o results

# Use the Li (1993) model with 8 threads
eskaks fasta aligned_genes.fasta --model li --workers 8 -o results

# Read from stdin
cat aligned_genes.fasta | eskaks fasta - -o results
```

## pN/pS per gene (VCF)

```bash
# One VCF per sample: AF computed as fraction of samples
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf sample1.vcf --vcf sample2.vcf --vcf sample3.vcf \
  --af-weighted --genetic-code 11 -o population_pnps

# Or use a list file with VCF paths
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf-list samples.txt --af-weighted --plot -o population_pnps

# Single multi-sample VCF also works
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --pass-only --min-af 0.05 -o filtered
```

## Input requirements

Your input must be a **codon-aligned** FASTA file:

- All sequences should be the same length
- Sequence length should be a multiple of 3 (complete codons)
- Standard DNA/RNA alphabet (A, C, G, T/U)
- Gaps (`-`, `.`) are treated as ambiguous and skipped
- Ambiguous bases (N, etc.) produce invalid codons (also skipped)

> **Tip**: If your sequences are not codon-aligned, use tools like [MAFFT](https://mafft.cbrc.jp/) + [PAL2NAL](http://www.bork.embl.de/pal2nal/) or [MACSE](https://bioweb.supagro.inra.fr/macse/) first.

## Output

By default, eskaks produces a TSV file with pairwise results:

```
Seq1    Seq2    dN      dS      dN/dS
gene_A  gene_B  0.0523  0.3214  0.1627
gene_A  gene_C  0.0891  0.4102  0.2172
...
```

Each run ends with a confirmation of what was done and where the output went:

```
── Done ───────────────────────────────────
  Sequences:  6 (6 unique) from aligned_genes.fasta
  Model:      Nei-Gojobori
  Output:
    results_pairwise_results.tsv
────────────────────────────────────────────
```

Add `--summary` for the full statistics block, or `--quiet` to suppress this.

## Common workflows

### Positive selection scan

```bash
# Compute pairwise dN/dS with summary statistics
eskaks fasta genes.fasta --model nei --summary --plot -o scan

# Filter pairs with dN/dS > 1 (positive selection)
awk -F'\t' 'NR>1 && $5>1' scan_pairwise_results.tsv
```

### Sliding window analysis

```bash
# 50-codon windows, stepping by 10
eskaks fasta genes.fasta --window-size 50 --window-step 10 --plot -o windows
```

### Group comparisons

```bash
# Sequences named like "lineageA_gene1", "lineageB_gene2"
eskaks fasta genes.fasta --group-average -o groups

# Or group by first letter
eskaks fasta genes.fasta --lineage --first-letter-lineage -o lineages
```

### JSON for pipelines

```bash
# Machine-readable output
eskaks fasta genes.fasta --format json -o results
cat results_pairwise_results.json | python3 -c "
import json, sys
data = json.load(sys.stdin)
pos = [r for r in data if r['dN_dS'] and r['dN_dS'] > 1]
print(f'{len(pos)} pairs under positive selection')
"
```
