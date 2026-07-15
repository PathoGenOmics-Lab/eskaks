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
| `--kappa <FLOAT>` | Transition/transversion rate ratio for [spectrum-aware site counting](#mutation-spectrum-aware-site-counting---kappa) | `1.0` (equal rates) |
| `--min-snps <INT>` | Drop genes with fewer SNPs from the per-gene table, plot, and [test](#per-gene-neutrality-test) (the pooled estimate still uses all genes) | `0` |
| `--fdr <FLOAT>` | FDR threshold for calling genes significant in the [neutrality test](#per-gene-neutrality-test) | `0.05` |
| `--mk` | Also run a per-gene [McDonald-Kreitman test](#mcdonald-kreitman-test) (writes `<prefix>_mk.<ext>`) | off |
| `--mk-fixed-af <FLOAT>` | AF at/above which a variant is "fixed" (divergence) vs polymorphic in the MK test | `0.99` |
| `--bootstrap <INT>` | Bootstrap replicates for a 95% CI on the genome-wide pooled pN/pS (0 = off) | `0` |
| `--seed <INT>` | Seed for reproducible bootstrap resampling | `42` |
| `--workers <INT>` | Parallel threads for the per-gene computation | `4` |
| `--plot` | Generate Manhattan-style SVG plot (significant genes outlined) | off |

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
| Chrom / Start / End / Strand | Gene location (1-based), so the table is self-contained and joinable |
| Exp_N_frac | Expected nonsynonymous fraction under neutrality, `N/(N+S)` (the test's null) |
| P_value | Two-sided exact-binomial p-value for H0: pN/pS = 1 (`NA` if untested) |
| Q_value_BH | Benjamini-Hochberg FDR q-value across all tested genes |
| P_Bonferroni | Bonferroni-corrected p-value across all tested genes |

## Per-gene neutrality test

Under strict neutrality the nonsynonymous fraction of a gene's SNPs equals its
**mutational opportunity** `N/(N+S)` (the `Exp_N_frac` column). eskaks tests each
gene against this null with a **two-sided exact binomial test** (observed
nonsynonymous SNPs vs `N/(N+S)`), giving a `P_value` per gene. Because a bacterial
genome has thousands of genes, the p-values are corrected for multiple testing two
ways: **Benjamini-Hochberg** FDR q-values (`Q_value_BH`) and the more conservative
**Bonferroni** (`P_Bonferroni`). The summary reports how many genes are significant
at `--fdr` (default 0.05), and `--plot` outlines those genes on the Manhattan plot.

```bash
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --kappa 2 --min-snps 5 --fdr 0.05 --plot -o mtb_scan
```

- Use `--min-snps` to drop low-count genes whose p-values are unreliable; the
  correction is then applied only over the genes actually tested.
- With `--kappa`, the null `N/(N+S)` is computed under the same ts/tv model as the
  observed counts, so the test and the site model stay consistent.
- The test is **skipped under `--af-weighted`** (πN/πS uses fractional counts, which
  a binomial test can't take).

> **Caveat:** the binomial null assumes SNPs are independent draws over sites. In
> highly clonal or strongly linked bacterial populations this is violated and the
> p-values are anti-conservative — treat them as a ranking aid, and lean on the FDR
> column rather than raw p-values.

## McDonald-Kreitman test

`--mk` writes `<prefix>_mk.<ext>` with a per-gene McDonald-Kreitman test. eskaks
already classifies every SNP as synonymous/nonsynonymous and knows each variant's
allele frequency, so the fixed-vs-polymorphic split is computed with no extra
input: ALTs with `AF >= --mk-fixed-af` (default 0.99) are treated as **fixed**
(divergence), the rest as **polymorphic**. Each gene gets the 2×2 table
`[Dn, Ds; Pn, Ps]` plus:

- **NI** (Neutrality Index) = `(Pn/Ps)/(Dn/Ds)` — > 1 suggests purifying selection, < 1 adaptive.
- **alpha** = `1 − (Ds·Pn)/(Dn·Ps)` — the estimated proportion of adaptive substitutions.
- **Fisher_p** — a two-sided Fisher exact test on the table, with a **Fisher_q_BH** FDR q-value across genes.

```bash
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --mk --mk-fixed-af 0.95 -o mtb_mk
```

> **Caveat:** this is a **reference-polarized** MK test — "fixed" means high AF
> *within your sample*, which conflates high-frequency derived alleles with true
> between-species divergence. It is a fast screen for adaptation, not a substitute
> for an outgroup-based MK test.

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

Add `--bootstrap N` (with `--seed`) for a reproducible **95% confidence interval**
on the pooled ratio, obtained by resampling genes with replacement. A wide interval
warns that a few gene-rich loci dominate the pooled estimate.

## Mutation-spectrum-aware site counting (`--kappa`)

By default eskaks counts synonymous (S) and nonsynonymous (N) sites the classic
Nei-Gojobori way: every possible single-nucleotide change at a codon is treated
as equally likely. Real mutational spectra are **not** uniform — most genomes,
and *M. tuberculosis* especially, are strongly **transition-biased** (A↔G, C↔T
mutations are far more frequent than transversions).

This matters because the synonymous change at a 2-fold degenerate third position
is almost always a transition. Counting it with equal weights **under-counts
synonymous sites**, which inflates N, deflates S, and biases pN/pS **downward** —
potentially making genuinely conserved genes look neutral and masking selection.

`--kappa <ratio>` corrects this by weighting each candidate change by its
relative mutation rate — `kappa` for a transition, `1` for a transversion —
when counting sites (the [modified Nei-Gojobori / Ina 1995 correction](models.md)):

```bash
# Transition/transversion ratio of ~2 (a typical bacterial value)
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --kappa 2 -o mtb_pnps
```

- A 2-fold synonymous site's contribution moves from `1/3` (κ=1) to `κ/(κ+2)`.
- 4-fold degenerate sites are synonymous regardless, so they are **κ-invariant**.
- The total number of sites per codon is unchanged; only the **S/N split** moves.
- Weighting keeps the same codon-level normalisation as the default, so `κ` is the
  *only* thing that changes — the correction is not confounded with a scheme switch.
- **Across the coding genome**, transition bias (`κ>1`) generally raises total S,
  lowers total N, and so raises pN/pS relative to the equal-rates estimate.
  Individual codons can move either way — a codon whose transition-reachable
  changes are mostly *nonsynonymous* (e.g. some 4-fold codons) shifts the other
  way — so read the direction genome-wide, not gene-by-gene.

Only the site **denominators** change. The observed synonymous/nonsynonymous SNP
counts (the numerators) are read directly from the VCF and are never rate-weighted,
so `--kappa` is orthogonal to `--af-weighted`. `--kappa 1` (the default) reproduces
the classic equal-rates counting exactly (bit-for-bit).

> A single `κ` captures the transition/transversion contrast but not base-specific
> asymmetries (e.g. C→T ≠ A→G) or context effects. Supply a value from the
> literature or estimate it from your own data (e.g. from 4-fold degenerate sites).

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
   - **Count sites**: For each reference codon, enumerate all 9 possible single-nucleotide changes, *excluding* changes to stop codons, and classify each as synonymous or nonsynonymous. Each codon contributes exactly 3 sites, split proportionally: S_sites = 3 × syn/(syn+nonsyn). With [`--kappa`](#mutation-spectrum-aware-site-counting---kappa) each change is weighted by its transition/transversion rate before this split.
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
