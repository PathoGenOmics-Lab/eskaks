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

```bash
eskaks vcf --ref reference.fasta --gff annotation.gff3 --vcf variants.vcf -o results
```

### Required arguments

| Flag | Description |
|---|---|
| `--ref <FASTA>` | Reference genome in FASTA format |
| `--gff <GFF3>` | Gene annotation in GFF3 format (CDS features) |
| `--vcf <VCF>` | Variants in VCF format (SNPs) |

### Options

| Flag | Description | Default |
|---|---|---|
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table | `1` |
| `--pass-only` | Only include FILTER=PASS (or `.`) variants | off |
| `--min-af <FLOAT>` | Minimum allele frequency (0.0–1.0) | none |
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

## Examples

```bash
# Basic pN/pS for M. tuberculosis
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf population.vcf \
  --genetic-code 11 -o mtb_pnps

# Filter: only PASS variants with AF ≥ 0.05 and depth ≥ 10
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --pass-only --min-af 0.05 --min-depth 10 -o filtered

# JSON output with Manhattan plot
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --format json --plot -o results
```

## How it works

1. **Load reference**: Parse FASTA into a sequence map
2. **Parse GFF3**: Extract CDS features, group by gene (Parent/gene_id), handle multi-exon, strand, phase
3. **Parse VCF**: Extract SNPs (skip indels), parse allele frequencies from INFO/AF or GT fields
4. **Apply filters**: PASS-only, minimum AF, minimum depth
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
- Allele frequency: parsed from `INFO/AF`, or calculated from `GT` fields if AF is absent
- Read depth: parsed from `INFO/DP`

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
2. **Low SNP counts**: Genes with very few SNPs produce unreliable pN/pS estimates. Consider filtering genes with < 5 total SNPs.
3. **Overlapping genes**: Each SNP is assigned to all genes whose CDS regions overlap its position. Overlapping genes on opposite strands will classify the same SNP differently.
4. **Allele frequency weighting**: Currently, each SNP is counted once regardless of its allele frequency. Future versions may support AF-weighted counting.
