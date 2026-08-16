---
title: Example reports
description: >-
  Two live, fully-interactive eskaks dashboards embedded in the docs: the
  per-gene pN/pS report built from a VCF, and the pairwise dN/dS report built
  from a codon alignment.
hide:
  - toc
---

# The interactive reports

`--report` writes a **single, self-contained HTML file**: every style and script
inlined, no internet needed, so it opens by double-clicking and can be emailed or
archived alongside the tables. There are **two such reports**, one per analysis,
and they are different dashboards rather than two skins of the same one. They
take different inputs, count different things, and answer different questions.

|  | pN/pS report | dN/dS report |
|---|---|---|
| **Written by** | `eskaks vcf --report` | `eskaks fasta --report` |
| **Input** | VCF + reference FASTA + GFF | one codon-aligned FASTA |
| **Question** | which genes depart from neutrality *within* a population | how strong is selection *between* sequences, pair by pair |
| **Unit of analysis** | one row per gene | one point per pair of sequences |
| **Panels in the example** | verdict, stat cards, selection regimes, Manhattan, volcano, QQ, polymorphism vs divergence, power funnel, allele-frequency spectrum, per-gene table | verdict, stat cards, sliding window, per-lineage scatter, dN vs dS, pairwise distribution |
| **Live example** | [go to the pN/pS dashboard](#pnps-report) | [go to the dN/dS dashboard](#dnds-report) |

Both dashboards below are live and both are the real output of the command shown
above them, run on the datasets bundled in `examples/`.

## Per-gene pN/pS, from a VCF { #pnps-report }

```bash
eskaks vcf \
  --ref examples/toy_genome/reference.fasta \
  --gff examples/toy_genome/genes.gff3 \
  --vcf examples/toy_genome/variants.vcf \
  --divergence examples/toy_genome/divergence.tsv \
  --genetic-code 11 --report --bootstrap 500 -o toy_scan
# writes toy_scan_report.html
```

Twelve genes of the toy genome, scored one by one. `--divergence` is what adds
the polymorphism-versus-divergence panel, which sets each gene's within-sample
pN/pS against its long-term dN/dS; `--bootstrap 500` puts a confidence interval
on the genome-wide ratio. Adding `--mk` would put a McDonald-Kreitman panel and
its columns alongside them.

!!! tip "It is fully interactive, try it"
    Click a gene in any panel to highlight it everywhere · search for a gene by
    name · toggle **FDR ↔ Bonferroni** stringency · switch the
    **colour-blind–safe** palette · flip **light/dark** · export **CSV / JSON**
    or print to PDF. Hover the **ⓘ** buttons for how to read each panel.
    [Open it full-screen :material-open-in-new:](assets/example-report.html){ target="_blank" rel="noopener" }

<iframe
  class="example-report-frame"
  src="../assets/example-report.html"
  title="eskaks interactive pN/pS report, toy genome example"
  loading="lazy"
  style="width:100%; height:84vh; min-height:520px; border:1px solid var(--md-default-fg-color--lightest); border-radius:.5rem;">
</iframe>

## Pairwise dN/dS, from an alignment { #dnds-report }

```bash
eskaks fasta examples/lineages.fasta --lineage --report -o lineage_demo
# writes lineage_demo_report.html
```

The same 180 bp gene from six isolates, two each from `Lineage2`, `Lineage4` and
`Bovis`, so all 15 pairs are compared and then read back by lineage. The sliding
window, the dN-versus-dS scatter and the distribution of per-pair ratios come
with every `--report`; `--lineage` is what fills the per-lineage panel, with one
point per genome and a bar at each lineage mean. Its alternative,
`--group-average`, cannot be combined with it: that mode replaces the same panel
with one group mean and a 95% confidence interval per group pair.

!!! tip "Also fully interactive, with a smaller control set"
    Hover any point, bar or window to read its numbers · flip **light/dark** ·
    export **CSV / JSON** or print to PDF. Hover the **ⓘ** buttons for how to
    read each panel. There is no gene search or per-gene table here, because the
    unit is a pair of sequences rather than a gene.
    [Open it full-screen :material-open-in-new:](assets/example-report-fasta.html){ target="_blank" rel="noopener" }

<iframe
  class="example-report-frame"
  src="../assets/example-report-fasta.html"
  title="eskaks interactive dN/dS report, lineage example"
  loading="lazy"
  style="width:100%; height:84vh; min-height:520px; border:1px solid var(--md-default-fg-color--lightest); border-radius:.5rem;">
</iframe>

See [interpreting results](interpreting-results.md) for how to read each panel,
[VCF analysis](vcf-analysis.md) for the full pN/pS command and its options, and
the [quick start](quickstart.md) for the dN/dS modes these panels are built from.
