---
description: >-
  A zero-to-results, no-background-assumed walkthrough of eskaks on its bundled
  example data — pairwise dN/dS, group averages, per-gene pN/pS from a VCF and the
  interactive report.
---

# Getting started: a hands-on tutorial

New to eskaks? This walkthrough takes you from zero to a finished analysis using
**example data that ships with the tool**: no biology background assumed. Copy and
paste each command; your numbers should match the ones shown here (rows may appear
in a different order).

!!! tip "In a hurry? Run the demo"
    `eskaks --demo` runs **both** analyses end to end on bundled data - no input
    files, no flags: pairwise dN/dS *and* the full per-gene pN/pS scan (including
    `--variants` and `--diversity`), each writing an interactive report. It's the
    fastest way to confirm your build works before following the steps below.

!!! tip "Prefer a notebook?"
    The same analysis is available as a runnable [**Jupyter notebook**](tutorial-notebook.ipynb)
    that loads the outputs with `pandas` and plots them with `matplotlib` — read it
    online or download the `.ipynb` and run it yourself.

## What does eskaks actually measure?

Genes are written in **codons** (triplets of DNA). A mutation in a codon is either:

- **synonymous**: it changes the DNA but **not** the amino acid (a "silent" change), or
- **nonsynonymous**: it **does** change the amino acid.

Natural selection mostly cares about the protein, so it acts on nonsynonymous
changes and largely ignores synonymous ones. eskaks counts both and takes their
ratio, correcting for how many sites of each kind exist:

- **dN/dS** (from aligned sequences) compares **species**: fixed differences.
- **pN/pS** (from a VCF) compares individuals **within a population**: variants
  still segregating.

Read the ratio like a thermostat for selection:

| Ratio | Meaning |
|---|---|
| **≈ 0** (well below 1) | **Purifying selection**: amino-acid changes are being removed. The normal state of most functional genes. |
| **≈ 1** | **Neutral**: no net selection (or too little data to tell). |
| **> 1** | **Positive / diversifying selection**: amino-acid changes are favoured. Interesting: drug targets, antigens, immune genes. |

That's the whole idea. Now let's run it.

## 0. Install

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks
make release
cp target/release/eskaks ~/.local/bin/   # or run ./target/release/eskaks
eskaks --version
```

(Needs [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70. See
[installation](installation.md) for details.)

The example files below live in the `examples/` folder of the repository.

## 1. Your first run: dN/dS from sequences

`examples/genes.fasta` holds six versions of the same gene from six strains,
already **codon-aligned** (in frame, all the same length). Compute dN/dS for every
pair:

```bash
eskaks fasta examples/genes.fasta -o first_run
```

eskaks confirms what it did and where the output went (add `--quiet` to hide this):

```text
── Done ───────────────────────────────────
  Sequences:  6 (6 unique) from examples/genes.fasta
  Model:      Nei-Gojobori
  Output:
    first_run_pairwise_results.tsv
────────────────────────────────────────────
```

Open `first_run_pairwise_results.tsv`: one row per pair of strains (order may vary,
the values don't):

```text
Seq1      Seq2      dN        dS        dN/dS
strain_B  strain_C  0.038372  0.423770  0.090548
strain_B  strain_F  0.038421  0.349819  0.109830
strain_E  strain_F  0.054217  0.319140  0.169885
...
```

**How to read it:** every `dN/dS` is well below 1 (≈ 0.09–0.22). These genes are
under **purifying selection**: exactly what you expect for a functioning gene:
silent changes accumulate freely (`dS` is high), but amino-acid changes are held
back (`dN` is low). A value above 1 would have flagged positive selection.

> 💡 Your own sequences must be **codon-aligned** first (in frame, gaps in
> multiples of 3). Align them with [MAFFT](https://mafft.cbrc.jp/) +
> [PAL2NAL](http://www.bork.embl.de/pal2nal/) or [MACSE](https://bioweb.supagro.inra.fr/macse/).

## 2. Compare groups, not just pairs

When your sequences fall into groups (lineages, clades, treatment arms), you usually
want the mean dN/dS *between* groups rather than one row per pair.
`examples/lineages.fasta` holds the same gene from six isolates named by lineage
(`Lineage2`, `Lineage4`, `Bovis`, two isolates each):

```bash
eskaks fasta examples/lineages.fasta -o lineages --group-average
```

`lineages_group_avg_dn_ds.tsv` has one row per pair of groups (plus each group
against itself). The between-lineage rows are the interesting ones:

```text
Group1    Group2    NumComparisons  Mean_dN/dS
Lineage2  Lineage4  4               0.608113
Lineage2  Bovis     4               0.700671
Lineage4  Bovis     4               1.207500
```

eskaks reads the group from the part of each ID before the first `_`. To group by the
first letter of the ID instead (so `Lineage2` and `Lineage4` merge into a single `L`
group and `Bovis` becomes `B`), add `--first-letter-lineage`:

```bash
eskaks fasta examples/lineages.fasta -o lineages_fl --group-average --first-letter-lineage
```

For a per-genome view (each isolate's mean dN/dS against every lineage), swap
`--group-average` for `--lineage`.

## 3. Per-gene pN/pS from a VCF

`examples/toy_genome/` is a miniature genome with three files, the three inputs
`eskaks vcf` always needs:

| File | What it is |
|---|---|
| `reference.fasta` | the genome sequence |
| `genes.gff3` | where the genes are (annotation) |
| `variants.vcf` | the SNPs found in your samples |

Run the full analysis and write an interactive report:

```bash
eskaks vcf \
  --ref examples/toy_genome/reference.fasta \
  --gff examples/toy_genome/genes.gff3 \
  --vcf examples/toy_genome/variants.vcf \
  --genetic-code 11 \
  --report --plot \
  --divergence examples/toy_genome/divergence.tsv \
  -o toy_scan
```

`--genetic-code 11` is the bacterial code (run `eskaks --list-codes` for the list).
eskaks prints a summary to the screen, then the list of files written:

```text
── pN/pS Summary ──────────────────────────
  Genes analyzed:      12
  Genes with SNPs:     12
  SNPs used (in CDS):  264 of 264 parsed
  Total synonymous:    115.00
  Total nonsynonymous: 149.00
  ── Genome-wide (pooled) ──────────────────
  N / S sites:         3200.3 / 1101.7
  Overall pN / pS:     0.046558 / 0.104384
  Overall pN/pS:       0.446028
  Selection:           purifying selection (pN/pS < 1)
  ── Neutrality test (pN/pS = 1) ───────────
  Genes tested:        12
  Significant genes:   7  (BH-FDR < 0.05)
```

The **genome-wide** pN/pS is 0.45, the genome as a whole is under purifying
selection, and 7 of 12 genes reject the "no selection" null. The `264 of 264`
accounting confirms every SNP landed inside a CDS: if many were dropped, that line
would flag a contig-name or coordinate mismatch between your VCF, GFF and reference.

## 4. Read the per-gene table

Open `toy_scan_pnps.tsv` (a plain table, Excel, R, pandas, or `column -t` all work):

```text
Gene       pN/pS     Nonsyn  Syn  SNPs  P_value   Q_value_BH
gene01     0.2259    10      15   25    4.97e-4   ...
gene03     0.1235    4       11   15    2.65e-4   ...
PPE_toy1   1.6171    18      4    22    0.54      ...
gene05     0.8884    19      7    26    0.94      ...
```

- `gene01`, `gene03`: `pN/pS` well below 1 with a **small q-value** → significant
  **purifying** selection.
- `PPE_toy1`: `pN/pS` above 1, but its **q-value is not significant** (0.54), and
  the name (`PPE`) marks it as a **repetitive** gene, where SNP calls are often
  mapping artefacts. Treat it with caution, not excitement.
- `gene05`: `pN/pS ≈ 1` and not significant → no evidence of selection here.

> ⚠️ A **non-significant** gene is usually just **underpowered** (too few SNPs);
> it is *not* proof of neutrality.

## 5. Explore the interactive report

Open `toy_scan_report.html` in any browser (it is fully self-contained, no
internet needed). Use the **table of contents on the left** to jump between panels,
and click the small **"i"** on any panel for a plain-language explanation of it.

Things to try:

- The **Manhattan** and **Volcano** plots, significant genes stand out; click one to
  highlight it everywhere.
- The **"How to read this report"** box at the top, a glossary of every metric.
- The **👁 CVD** button (colour-blind-safe palette + shapes) and the **◑ Theme**
  (light/dark) button in the top bar.
- Export the filtered table with **⤓ CSV** / **⤓ JSON**, or **🖶 Print** to PDF.

## 6. Now use your own data

You need the same three inputs, with **matching contig names** across all of them:

- a **reference FASTA** of your organism,
- a **GFF3** annotation (CDS features),
- one or more **VCF** files of your variants (one per sample, or a multi-sample VCF).

A realistic *M. tuberculosis* run:

```bash
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --kappa 2 \
  --min-snps 5 --mk --bootstrap 1000 --report -o mtb_scan
```

- `--kappa 2` corrects for the transition-heavy mutation spectrum of TB
  ([why](vcf-analysis.md#kappa)).
- `--min-snps 5` drops genes with too few SNPs to test.
- `--mk` adds a McDonald-Kreitman test; `--bootstrap 1000` adds a confidence
  interval on the genome-wide estimate.

## Where to go next

- [Interpreting results](interpreting-results.md), what the numbers mean, and the pitfalls
- [VCF analysis](vcf-analysis.md), every pN/pS option, the neutrality test, the report
- [Output formats](output-formats.md), every file eskaks writes and how to parse it
- [Glossary](glossary.md), every term in one place
- [CLI reference](cli-reference.md), all flags
- [FAQ](faq.md), common questions and errors
