---
description: >-
  What eskaks' dN/dS and pN/pS numbers mean, the traps (saturation, low power,
  PE/PPE mapping artefacts, clonality) and a panel-by-panel guide to the
  interactive report.
---

# Interpreting Results

## dN/dS (ω) ratios

The dN/dS ratio (also called ω or Ka/Ks) measures the balance between nonsynonymous (amino acid-changing) and synonymous (silent) substitutions:

| dN/dS | Interpretation |
|-------|----------------|
| **ω < 1** | **Purifying selection**: Most amino acid changes are deleterious and removed by natural selection. This is by far the most common result for functional genes. |
| **ω ≈ 1** | **Neutral evolution**: Amino acid changes are neither beneficial nor deleterious. May indicate a pseudogene or relaxed constraint. |
| **ω > 1** | **Positive selection**: Amino acid changes are favored. Rare in whole-gene comparisons; more common in specific codons or domains. |

## Common values

- **Housekeeping genes**: ω ≈ 0.01–0.1 (strong purifying selection)
- **Immune genes**: ω ≈ 0.5–2.0 (variable, some under positive selection)
- **Pseudogenes**: ω ≈ 1.0 (no functional constraint)
- **Viral surface proteins**: ω > 1 common in specific sites

## Special values

- **NaN**: undefined. Either the Jukes-Cantor / Kimura correction reached saturation (the sequences are too divergent to estimate reliably, typical when `p ≥ 0.75` — more than 75% of sites differ), or the pair has no comparable codons at all (e.g. an all-`N` / all-gap sequence). Reported as `null` in JSON.
- **0.0 / 0.0 = 0.0**: Identical sequences. No substitutions observed.
- **dN > 0, dS = 0**: All observed changes are nonsynonymous. The ratio is technically infinity, reported as `inf` (TSV/CSV) or `null` (JSON).

## Caveats

1. **Pairwise ω is a genome-wide average**. Positive selection at specific sites can be masked by purifying selection elsewhere. Use site-level methods (PAML, HyPhy) for finer resolution.

2. **Short sequences**: dN/dS estimates become unreliable with few codons. Use `--min-codons` to filter very short sequences.

3. **Saturation**: Very divergent sequences (>75% divergence) yield NaN. This is biologically correct, the signal is saturated and the estimate is unreliable.

4. **Internal stop codons**: eskaks warns about these. They usually indicate frameshifts, pseudogenes, or incorrect reading frames. Consider excluding these sequences.

5. **Recombination**: dN/dS assumes a single phylogenetic history. Recombination can bias estimates. Consider using sliding windows (`--window-size`) to detect mosaic patterns.

## Sliding window interpretation

Window analysis reveals variation along the alignment:

- **Peaks (ω > 1)**: Potential positive selection hotspots (e.g., surface-exposed domains)
- **Valleys (ω ≈ 0)**: Highly conserved regions (e.g., catalytic sites, structural cores)
- **Noisy windows**: Short windows or few differences produce unreliable estimates. Use windows of at least 50–100 codons.

## pN/pS: reading within-species selection

`eskaks vcf` reports **pN/pS**, the polymorphism analogue of dN/dS. It reads the
same way — below 1 is purifying selection, ≈ 1 is neutral/underpowered, above 1
hints at diversifying selection — but polymorphism data has its own traps that
routinely produce false hits. Read a per-gene pN/pS through these three filters
**before** you believe it:

1. **Is it significant, and does the CI exclude 1?** A ratio is a point estimate.
   Most bacterial genes carry only a handful of SNPs, so a pN/pS of 1.8 with two
   nonsynonymous SNPs is noise. Lean on the `Q_value_BH` column and the Wilson CI
   whiskers in the report, not the raw ratio — and remember a **non-significant
   gene is underpowered, not proven neutral**.

2. **Is it a repetitive / hard-to-map gene?** In *M. tuberculosis* the **PE/PPE/PGRS**
   families, IS elements and maturases are the single most common source of spurious
   high pN/pS: paralogous reads misalign, manufacturing high-AF "nonsynonymous"
   calls that are mapping artefacts, not selection. A pN/pS > 1 in a `PE`/`PPE`
   gene is a **red flag to distrust**, not a discovery. Use `--exclude-repetitive`
   to drop them from the pooled estimate and the test family, and check the report's
   core-vs-repetitive comparison.

3. **Is the population clonal?** Genome-wide linkage in clonal organisms breaks the
   per-gene test's independence assumption and inflates significance (watch the
   **λ** card). See the clonality note under [Limitations](#limitations).

!!! warning "pN/pS is not dN/dS"
    pN/pS measures *raw polymorphism proportions* with **no** multiple-hit
    correction (no Jukes-Cantor/Kimura), so it is not on the same scale as a
    divergence dN/dS and the two should never be compared numerically. To separate
    selection from demography you need diversity statistics (`--diversity`: πN/πS,
    Watterson θ, Tajima's D) or an outgroup-based test; a single pN/pS cannot.

## Reading the interactive report

`eskaks vcf --report` (and `eskaks fasta --report`) writes a self-contained HTML
dashboard — see the [**live example**](example-report.md). Every panel also carries
an **"i"** button repeating this guidance in-app. Here is what each one answers:

- **Verdict banner & summary cards** — the headline: the genome-wide **pooled**
  pN/pS (sites-weighted, not the mean of the column), its bootstrap CI if you ran
  `--bootstrap`, the significant-gene count at your FDR, and the genomic-inflation
  **λ**. Start here, then drill down.
- **Manhattan** — every gene laid out by genome position, height = significance
  (`−log10 p`), with the FDR/Bonferroni line drawn in. Points above the line are
  your candidate hits. Toggle the metric to `pN/pS` or `z(N)` to see effect size
  rather than significance.
- **Volcano** — effect (pN/pS) on x against significance on y. The top corners are
  what matter: strong departure *and* significant. A gene high on the y-axis but
  near pN/pS = 1 is significant but small-effect.
- **p-value QQ plot (λ)** — observed vs expected p-values. On the diagonal = well
  calibrated; a curve lifting off the diagonal (λ ≫ 1) means far more low p-values
  than chance — either pervasive selection or, in clonal data, inflation. This is
  the panel that tells you whether to trust the raw p-values at all.
- **Power funnel** — pN/pS against SNP count, with per-gene Wilson CI whiskers.
  Low-count genes fan out into wide CIs at the left; only genes whose CI **excludes
  1** are departing from neutrality. Guards you against over-reading small genes.
- **Observed vs expected** — each gene's observed nonsynonymous SNP count against
  its neutral expectation `N/(N+S)`. Points below the line are conserved (fewer
  amino-acid changes than expected), above the line diversifying.
- **Allele-frequency spectrum (SFS)** — pN/pS binned by allele frequency. A profile
  that **falls** as frequency rises is the classic signature of purifying selection
  keeping deleterious nonsynonymous variants rare. Needs a multi-sample cohort; a
  single sample shows an empty-state note.
- **pN/pS distribution** — the histogram of per-gene ratios, centred well below 1
  for a genome under broad constraint. The long right tail is where candidate
  targets sit.
- **McDonald-Kreitman** (with `--mk`) — per-gene `[Dn Ds; Pn Ps]` table with the
  Neutrality Index and α. NI > 1 / α < 0 suggests purifying selection, NI < 1 /
  α > 0 adaptive. Remember this is a **reference-polarized** proxy, not an
  outgroup MK test.
- **Polymorphism vs divergence** (with `--divergence`) — reconciles within-species
  pN/pS against between-species dN/dS per gene. Genes far off the diagonal are
  where polymorphism and divergence disagree — the interesting selection stories.
- **Per-gene table** — sortable, filterable, exportable (CSV/JSON). Click any row
  and the gene lights up across every panel above; this is how you go from a dot
  on the Manhattan to its actual `S315T`-style variants.

## Limitations

- `eskaks fasta` expects **codon-aligned** input (in-frame, gap lengths multiples
  of 3); align with MAFFT + PAL2NAL or MACSE first.
- `eskaks vcf` uses **SNPs only**; indels and multi-nucleotide variants are not
  codon-annotated (see [get_MNV](https://github.com/PathoGenOmics-Lab/get_MNV) for MNVs).
- pN/pS is estimated from **within-sample polymorphism**, so many genes are
  **underpowered** in low-diversity organisms; a non-significant result is not
  evidence of neutrality.
- The per-gene neutrality test assumes independent SNPs. In **clonal** organisms
  (e.g. *M. tuberculosis*) genome-wide linkage inflates significance;
  `--genomic-control` is a pragmatic correction, but a high λ can also reflect
  **real** pervasive selection, so apply it only when you suspect systematic bias.
- Contig names must match across the VCF, reference FASTA, and GFF3.
