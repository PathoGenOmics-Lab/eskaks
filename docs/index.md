---
title: eskaks
description: >-
  Fast pairwise dN/dS and per-gene pN/pS for molecular evolution and
  Mycobacterium tuberculosis selection analysis. Nei-Gojobori and Li (1993),
  1000× faster than existing tools.
# The sidebar is NOT hidden here, matching get_MNV. Hiding it gives the landing
# page a wider hero, at the cost of the one thing a first-time reader needs from
# a landing page: the shape of the documentation. Arriving with no visible map
# and having to click before seeing what exists is a worse trade than a narrower
# hero. Only the table of contents is hidden, since the hero and the card grid
# below it are the page's own navigation.
hide:
  - toc
---

<div class="eskaks-hero" markdown>

![eskaks](assets/logo.svg#only-light){ .eskaks-wordmark }
![eskaks](assets/logo-dark.svg#only-dark){ .eskaks-wordmark }

# eskaks

<p class="eskaks-hero__lead">
Fast pairwise <strong>dN/dS</strong> and per-gene <strong>pN/pS</strong> for
molecular evolution and <em>Mycobacterium tuberculosis</em> selection analysis.
Nei-Gojobori and Li (1993) with precomputed lookup tables, <strong>1000×</strong>
faster than existing tools.
</p>

<div class="eskaks-hero__actions" markdown>
[Get started :material-rocket-launch:](tutorial.md){ .md-button .md-button--primary }
[Install :material-download:](installation.md){ .md-button }
[View on GitHub :fontawesome-brands-github:](https://github.com/PathoGenOmics-Lab/eskaks){ .md-button }
</div>

</div>

!!! tip "New to selection analysis?"
    Start with the [**getting-started tutorial**](tutorial.md) — it runs a full analysis on the bundled example data with no background assumed. Every unfamiliar term is defined in the [glossary](glossary.md), and hovering a dotted abbreviation like dN/dS shows its meaning.

**In one line:** eskaks measures natural selection on genes by comparing how fast *amino-acid-changing* mutations accumulate versus *silent* ones. A ratio below 1 means selection is removing harmful changes (a conserved gene); a ratio above 1 flags genes where change is favoured (drug targets, antigens).

## What can it do?

<div class="grid cards" markdown>

-   :material-dna:{ .lg .middle } **Pairwise dN/dS**

    ---

    Codon-aligned FASTA in, dN/dS out — Nei-Gojobori or Li/LPB93, sliding windows, per-lineage and per-group summaries.

    [:octicons-arrow-right-24: Models](models.md)

-   :material-file-chart:{ .lg .middle } **Per-gene pN/pS**

    ---

    VCF + reference + GFF3 → per-gene pN/pS, an exact neutrality test, FDR correction, an MK screen, and π / Watterson θ / Tajima's D.

    [:octicons-arrow-right-24: VCF analysis](vcf-analysis.md)

-   :material-chart-scatter-plot:{ .lg .middle } **Interactive report**

    ---

    A single self-contained HTML dashboard — Manhattan, volcano, QQ, MK — with a colour-blind mode and CSV/JSON/Print export. No internet needed.

    [:octicons-arrow-right-24: See the live report](example-report.md)

-   :material-lightning-bolt:{ .lg .middle } **Fast & pipeline-ready**

    ---

    Precomputed lookup tables + rayon parallelism, stdin support, JSON output, and a non-zero exit on error.

    [:octicons-arrow-right-24: Performance](performance.md)

</div>

## Try it in one command

=== "Pairwise dN/dS (FASTA)"

    ```bash
    eskaks fasta alignment.fasta --report -o results
    ```

    Reads codon-aligned CDS, writes `results_pairwise_results.tsv` and an interactive `results_report.html`.

=== "Per-gene pN/pS (VCF)"

    ```bash
    eskaks vcf --ref genome.fasta --gff genes.gff3 --vcf variants.vcf \
      --genetic-code 11 --report --variants --diversity -o scan
    ```

    Writes per-gene pN/pS, a per-variant table (`S315T`-style keys), diversity statistics, and the report.

=== "No data yet?"

    ```bash
    eskaks --demo
    ```

    Runs both the dN/dS and per-gene pN/pS analyses on bundled example data - no input files needed.

## How does it compare?

| Feature | eskaks | KaKs_Calculator | BioPython | PAML yn00 |
|---|---|---|---|---|
| Speed (100 seqs) | **6 ms** | 7,703 ms | 111,619 ms | 697 ms |
| Li model R² vs KaKs | **1.000** | reference | — | — |
| Nei model R² vs KaKs | **0.999** | reference | 0.996 | — |
| Genetic codes | **20 tables** | 1 | 1 | limited |
| pN/pS from VCF | :material-check: | :material-close: | :material-close: | :material-close: |
| JSON output | :material-check: | :material-close: | :material-close: | :material-close: |
| Stdin pipe | :material-check: | :material-close: | :material-close: | :material-close: |
| Parallel | :material-check: (rayon) | :material-close: | :material-close: | :material-close: |

!!! quote "Citation"
    Ruiz-Rodriguez P, Coscolla M (2026). *eskaks: fast pairwise dN/dS and per-gene pN/pS.*
    <https://github.com/PathoGenOmics-Lab/eskaks>
