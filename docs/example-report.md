---
title: Example report
description: >-
  A live, fully-interactive eskaks pN/pS report embedded in the docs — Manhattan,
  volcano, QQ, McDonald-Kreitman and a per-gene table, built from the bundled toy
  genome.
hide:
  - toc
---

# The interactive report

`eskaks vcf --report` (and `eskaks fasta --report`) writes a **single,
self-contained HTML file** — every style and script inlined, no internet needed —
with the full interactive dashboard: a genome-wide verdict, stat cards, Manhattan,
volcano, QQ, McDonald-Kreitman, a power funnel, the site-frequency spectrum, and
a per-gene table. Below is a **live example** built from the bundled toy genome.

!!! tip "It's fully interactive — try it"
    Click a gene in any panel to highlight it everywhere · toggle **FDR ↔ Bonferroni**
    stringency · switch the **colour-blind–safe** palette · flip **light/dark** ·
    export **CSV / JSON** or print to PDF. Hover the **ⓘ** buttons for how to read
    each panel. [Open it full-screen :material-open-in-new:](assets/example-report.html){ target="_blank" rel="noopener" }

<iframe
  class="example-report-frame"
  src="../assets/example-report.html"
  title="eskaks interactive pN/pS report — toy genome example"
  loading="lazy"
  style="width:100%; height:84vh; min-height:520px; border:1px solid var(--md-default-fg-color--lightest); border-radius:.5rem;">
</iframe>

This example was generated with `--report --divergence` (the divergence file adds
the polymorphism-vs-divergence panel). See [interpreting results](interpreting-results.md)
for how to read each panel, and [VCF analysis](vcf-analysis.md) for the full command.
