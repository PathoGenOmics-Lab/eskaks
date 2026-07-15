# VCF Analysis (pN/pS)

eskaks can compute **pN/pS per gene** directly from a VCF file, a reference genome, and a GFF3 annotation.

## What is pN/pS?

While dN/dS measures *fixed* substitutions between species, **pN/pS** measures the ratio of *polymorphic* nonsynonymous to synonymous variants within a population:

| Metric | What it measures | Input |
|---|---|---|
| dN/dS | Divergence (fixed substitutions) | Aligned sequences (FASTA) |
| pN/pS | Polymorphism (segregating variants) | VCF + reference + annotation |

| pN/pS | Interpretation |
|---|---|
| **< 1** | Purifying selection — most amino acid changes are removed |
| **≈ 1** | Neutral evolution or pseudogene |
| **> 1** | Possible positive/diversifying selection |

## Usage

### One VCF per sample (typical workflow)

```bash
# Provide each sample's VCF individually
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf sample1.vcf --vcf sample2.vcf --vcf sample3.vcf \
  --af-weighted --genetic-code 11 -o population_pnps

# Or use a file with one VCF path per line
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf-list samples.txt --af-weighted -o population_pnps
```

When multiple VCFs are provided, eskaks **merges** them and computes the allele frequency as the fraction of samples carrying each variant. For example, if 30 out of 100 samples have a SNP → AF = 0.3.

### Single multi-sample VCF

```bash
eskaks vcf --ref reference.fasta --gff annotation.gff3 --vcf population.vcf -o results
```

When a single VCF is provided, allele frequencies are taken from INFO/AF or calculated from GT fields.

### Required arguments

| Flag | Description |
|---|---|
| `--ref <FASTA>` | Reference genome in FASTA format |
| `--gff <GFF3>` | Gene annotation in GFF3 format (CDS features) |
| `--vcf <VCF>` | VCF file(s) — use multiple times for per-sample VCFs |

### Options

| Flag | Description | Default |
|---|---|---|
| `--vcf-list <FILE>` | File with one VCF path per line (alternative to multiple `--vcf`) | none |
| `--af-weighted` | Weight SNP counts by allele frequency (πN/πS instead of pN/pS) | off |
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table | `1` |
| `--pass-only` | Only include FILTER=PASS (or `.`) variants | off |
| `--min-af <FLOAT>` | Minimum allele frequency (0.0–1.0) | none |
| `--max-af <FLOAT>` | Maximum allele frequency — use 0.99 to exclude fixed variants | none |
| `--min-depth <INT>` | Minimum read depth (INFO/DP) | none |
| `--plot` | Generate Manhattan-style SVG plot | off |

## Output

File: `<prefix>_pnps.tsv` (or `.csv` / `.json`)

| Column | Description |
|---|---|
| Gene | Gene name (from GFF3 Name/gene/locus_tag attribute) |
| Length_bp | Total CDS length in base pairs |
| N_sites | Nonsynonymous sites (fractional) |
| S_sites | Synonymous sites (fractional) |
| pN | Proportion of nonsynonymous polymorphisms (nonsyn_SNPs / N_sites) |
| pS | Proportion of synonymous polymorphisms (syn_SNPs / S_sites) |
| pN/pS | Ratio of pN to pS |
| Nonsyn_SNPs | Count of nonsynonymous SNPs in this gene |
| Syn_SNPs | Count of synonymous SNPs in this gene |
| Total_SNPs | Total SNPs falling in this gene |

## Genome-wide (pooled) pN/pS

After writing the per-gene table, eskaks prints a summary to stderr that ends
with a **genome-wide** estimate pooled across every analyzed gene:

```
── pN/pS Summary ──────────────────────────
  Genes analyzed:      2
  Genes with SNPs:     2
  Total synonymous:    2.00
  Total nonsynonymous: 1.00
  ── Genome-wide (pooled) ──────────────────
  N / S sites:         36.6 / 11.4
  Overall pN / pS:     0.027335 / 0.175182
  Overall pN/pS:       0.156036
  Selection:           purifying selection (pN/pS < 1)
───────────────────────────────────────────
```

The pooled ratio sums counts and sites over all genes **before** dividing:

```
pN = Σ nonsyn_SNPs / Σ N_sites
pS = Σ syn_SNPs    / Σ S_sites
```

This is deliberately *not* the mean of the per-gene pN/pS column. Averaging
ratios gives a single noisy gene (few sites, extreme ratio) the same weight as
a whole chromosome; pooling weights each gene by its number of sites, which is
the standard way to report an overall signal of selection. Under `--af-weighted`
the pooled figure is πN/πS.

The `Selection:` line is a coarse convenience label (`< 0.9` purifying, `0.9–1.1`
near-neutral, `> 1.1` diversifying), not a statistical test — formal inference
needs an explicit null model.

## pN/pS vs πN/πS

| Mode | Flag | How SNPs count | Best for |
|---|---|---|---|
| **pN/pS** | (default) | Each SNP counts as 1 | Presence/absence of variants |
| **πN/πS** | `--af-weighted` | Each SNP weighted by AF | Population diversity analysis |

**Example**: A SNP at AF=0.3 contributes 1.0 to pN/pS but only 0.3 to πN/πS.

## Examples

```bash
# Population πN/πS from per-sample VCFs (M. tuberculosis)
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf-list samples.txt --af-weighted --genetic-code 11 \
  --min-af 0.01 --max-af 0.99 --plot -o mtb_pnps

# Simple pN/pS from a single multi-sample VCF
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --pass-only --min-depth 10 -o filtered

# JSON output
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --format json -o results
```

## How it works

1. **Load reference**: Parse FASTA into a sequence map
2. **Parse GFF3**: Extract CDS features, group by gene (Parent/gene_id), handle multi-exon, strand, phase
3. **Parse VCF(s)**: Extract SNPs (skip indels). Multiple per-sample VCFs are merged; AF = fraction of samples with the variant. Single VCFs use INFO/AF or GT fields.
4. **Apply filters**: PASS-only, min/max AF, minimum depth
5. **For each gene**:
   - Extract the CDS sequence from the reference (handling exon order, reverse complement for minus strand, phase offset)
   - **Count sites**: For each reference codon, enumerate all 9 possible single-nucleotide changes. Classify each as synonymous or nonsynonymous. S_sites = syn_changes/3, N_sites = nonsyn_changes/3.
   - **Classify SNPs**: For each SNP within the gene's CDS, reconstruct the reference and alternate codons. Look up amino acids → synonymous or nonsynonymous.
   - **Compute**: pN = nonsyn_SNPs / N_sites, pS = syn_SNPs / S_sites, pN/pS = pN / pS

## Input format details

### VCF
- Standard VCFv4.x format
- Only single-base SNPs are used (indels and MNPs are skipped)
- Multi-allelic sites are supported (each ALT allele classified independently)
- **One VCF per sample** (recommended): provide via `--vcf sample1.vcf --vcf sample2.vcf` or `--vcf-list samples.txt`. AF is computed as fraction of samples carrying the variant.
- **Single multi-sample VCF**: AF parsed from `INFO/AF`, or calculated from `GT` fields
- Read depth: parsed from `INFO/DP`
- REF allele is verified against the reference genome (mismatches are warned and skipped)

### GFF3
- Standard GFF3 format
- Only `CDS` feature types are used
- Multi-exon genes: grouped by `Parent=` attribute (or `gene_id=` for GTF-style)
- Gene name: extracted from `gene=`, `Name=`, or `locus_tag=` attributes
- Phase (column 8): used to adjust reading frame for first exon

### Reference FASTA
- Standard FASTA format
- Sequence names must match the `CHROM` column in the VCF and `seqid` column in the GFF3

## Caveats

1. **pN/pS ≠ dN/dS**: pN/pS does not correct for multiple substitutions (no Jukes-Cantor or Kimura). It measures raw polymorphism proportions, not evolutionary rates.
2. **Low SNP counts**: Genes with very few SNPs produce unreliable per-gene pN/pS estimates. Consider filtering genes with < 5 total SNPs, or rely on the [genome-wide pooled estimate](#genome-wide-pooled-pnps) for the overall signal.
3. **Overlapping genes**: Each SNP is assigned to all genes whose CDS regions overlap its position. Overlapping genes on opposite strands will classify the same SNP differently.
4. **Fixed variants**: By default, fixed variants (AF=1.0) are included. Use `--max-af 0.99` to exclude them and analyze only segregating polymorphisms.
5. **pN/pS vs πN/πS**: Without `--af-weighted`, all SNPs count equally (pN/pS). With `--af-weighted`, rare variants contribute less than common ones (πN/πS). Choose based on your question.
