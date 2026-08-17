---
description: >-
  Plain-language definitions of every term used across eskaks and its output —
  dN/dS, pN/pS, synonymous vs nonsynonymous, FDR, Tajima's D and more.
---

# Glossary

Plain-language definitions of the terms used across eskaks and its output. If a
concept in the report or docs is unfamiliar, look here first.

## The core idea

- **Codon**: a triplet of DNA bases that codes for one amino acid. Genes are read
  three bases at a time.
- **Synonymous change (S)**: a DNA change that leaves the amino acid unchanged
  ("silent"). Selection mostly ignores these, so they accumulate freely.
- **Nonsynonymous change (N)**: a DNA change that alters the amino acid. Selection
  acts on these.
- **dN**: the rate of nonsynonymous change, per nonsynonymous **site**.
- **dS**: the rate of synonymous change, per synonymous **site**.
- **dN/dS (ω, Ka/Ks)**: the ratio of the two, comparing **species** (fixed
  differences). `< 1` purifying, `≈ 1` neutral, `> 1` positive selection.
- **pN / pS / pN/pS**: the same idea but from **polymorphism within a population**
  (variants still segregating), computed from a VCF instead of an alignment.
- **N sites / S sites**: the *number of opportunities* for nonsynonymous vs
  synonymous change in a gene, given its codons. Used to normalize the raw counts
  (so a long gene isn't automatically "more selected" than a short one).
- **πN/πS**: pN/pS weighted by each variant's allele frequency (turned on with
  `--af-weighted`).

## Kinds of selection

- **Purifying (negative) selection**: amino-acid changes are harmful and get
  removed. `pN/pS < 1`. The normal state of a working gene.
- **Neutral evolution**: changes are neither good nor bad. `pN/pS ≈ 1`.
- **Positive / diversifying selection**: amino-acid changes are favoured. `pN/pS > 1`.
  Look here for drug targets, antigens, immune-evasion genes.
- **Relaxed constraint**: a gene that used to be constrained but no longer is
  (e.g. a pseudogene).

## Inputs

- **FASTA**: a text file of sequences. For `eskaks fasta` it must be
  **codon-aligned**: every sequence the same length, in frame, gaps in multiples of 3.
- **VCF**: Variant Call Format: the list of positions where your samples differ
  from the reference (the SNPs).
- **GFF3**: an annotation file saying where the genes (CDS features) are on the genome.
- **Reference FASTA**: the genome sequence the VCF positions refer to.
- **Contig name**: the chromosome/sequence identifier (e.g. `chr1`, `NC_000962.3`).
  It must be **identical** across the VCF, reference, and GFF3, or eskaks can't line
  them up.
- **Allele frequency (AF)**: how common a variant is (0–1). Read from the VCF
  `INFO/AF` field.

## The statistics

- **Neutrality test**: for each gene, a two-sided **mid-p binomial test** asking
  "is the nonsynonymous fraction of SNPs different from what neutral evolution
  predicts (`N/(N+S)`)?" Gives a **p-value** per gene.
- **mid-p**: the correction that makes that discrete test usable. A binomial count
  is a whole number, so the textbook version (which counts the observed count fully
  on both sides) is left with far more than the nominal 5% in hand and rarely calls
  anything on a small gene. Mid-p counts the observed count only **half** on each
  side, which centres the p-values on 0.5 under a genuine null instead of ~0.8.
  The trade: the test is **no longer exact**, and on very small genes its real error
  rate can run slightly above the nominal one (see
  [per-gene neutrality test](vcf-analysis.md#per-gene-neutrality-test)).
- **p-value**: the probability of seeing a result this extreme if the gene were
  neutral. Small = surprising = evidence of selection.
- **Multiple testing**: testing thousands of genes produces false positives by
  chance. Two corrections:
  - **FDR / Benjamini-Hochberg (q-value)**: controls the *proportion* of false
    discoveries. More permissive; good for screening. **Use this one by default.**
  - **Bonferroni**: controls the chance of *any* false positive. Stricter; good for
    a confident shortlist.
- **Genome-wide (pooled) pN/pS**: one overall ratio, pooling SNP and site counts
  across all genes, as a summary of the whole coding genome.
- **Bootstrap 95% CI**: a confidence interval on the genome-wide estimate, obtained
  by resampling genes many times (`--bootstrap`).
- **Wilson confidence interval**: a per-gene confidence interval on pN/pS; if it
  excludes 1, the gene departs from neutrality.

## Advanced terms (in the report / advanced flags)

- **κ (kappa)**: the transition/transversion rate ratio. Some organisms (like *M.
  tuberculosis*) mutate transitions much more often; `--kappa` corrects the site
  counts for this so pN/pS isn't biased.
- **McDonald-Kreitman (MK) test**: contrasts **fixed** vs **polymorphic** changes to
  detect adaptation; reports the Neutrality Index and **α** (fraction of adaptive
  substitutions).
- **Genomic inflation factor (λ)**: how far the whole set of p-values departs from
  chance. `λ ≫ 1` means many strong signals, real widespread selection, *or*
  systematic bias / clonal linkage. `λ ≈ 1` is the normal reading, but it is **not**
  proof of calibration: the test is discrete, so a genuine null already lands a
  little below 1 (about **0.90** for genes with a handful of SNPs, about **0.97**
  once they carry tens). Read λ upwards, as an inflation flag, never downwards as a
  clean bill of health.
- **Genomic control**: an optional correction (`--genomic-control`) that divides
  each test statistic by λ, for when the inflation is artefactual.
- **z(N)**: a standardized measure of how far a gene's observed nonsynonymous count
  is from expectation; a power-aware effect size.
- **Site frequency spectrum (SFS)**: the distribution of variants by allele
  frequency. Purifying selection keeps harmful variants rare, so pN/pS **falling** as
  frequency rises is a signature of constraint.
- **PE/PPE, PGRS, IS elements**: repetitive, hard-to-map gene families (especially in
  *M. tuberculosis*) whose SNP calls are often artefacts. `--exclude-repetitive`
  drops them from the pooled estimate and the test.
- **Same-codon SNPs (MNV)**: two or more SNPs inside one codon. eskaks classifies each
  one against the *reference* codon, never against its neighbour, so the joint
  amino-acid change a genome actually carries is not the one reported. Such codons are
  counted in every run, and warned about when the allele frequencies put their SNPs on
  one haplotype. See [Same-codon SNPs](vcf-analysis.md#same-codon-snps).

## Special output values

- **NaN**: "not a number": the estimate is undefined (e.g. sequences too divergent
  and the correction saturated, or a gene with no variation). Reported as `null` in JSON.
- **inf**: infinity: e.g. `dN > 0` but `dS = 0`, so the ratio is unbounded. Reported as `null` in JSON.
- **NA**: the value was not computed for this row (e.g. an untested gene's p-value).
