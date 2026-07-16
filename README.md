<p align="center">
  <img src="img/esKaKs.svg" height="200" alt="eskaks logo" />
</p>

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-%23dea584?style=flat-square)](https://www.rust-lang.org/)
[![Version](https://img.shields.io/badge/version-1.4.0-%23149389?style=flat-square)](https://github.com/PathoGenOmics-Lab/eskaks/releases)
[![CI](https://github.com/PathoGenOmics-Lab/eskaks/actions/workflows/ci.yml/badge.svg)](https://github.com/PathoGenOmics-Lab/eskaks/actions/workflows/ci.yml)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square)](https://github.com/PathoGenOmics-Lab)

**Pairwise dN/dS (Ka/Ks) and per-gene pN/pS from codon-aligned sequences or VCF files.**
**Pure Rust · no C dependencies · self-contained interactive HTML reports.**

[Tutorial](docs/tutorial.md) · [Quick Start](#quick-start) · [Report](#interactive-html-report) · [Docs](docs/) · [Citation](#citation)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## What is eskaks?

> [!TIP]
> **New here?** The [**getting-started tutorial**](docs/tutorial.md) walks you from
> install to a finished analysis using example data that ships with the tool, no
> background assumed. Unfamiliar with a term? See the [glossary](docs/glossary.md).

**In plain terms:** genes are read in **codons** (DNA triplets). A mutation is
either *synonymous* (silent: the protein is unchanged) or *nonsynonymous* (it
changes the protein). Natural selection acts on the protein, so comparing the two
rates tells you what selection is doing: a ratio **below 1** means harmful changes
are being removed (*purifying selection*, the normal state of a working gene), and a
ratio **above 1** means changes are being favoured (*positive selection*: drug
targets, antigens, immune genes).

eskaks measures the strength and direction of that selection on protein-coding
genes, in two modes:

- **`eskaks fasta`**: **pairwise dN/dS** (Ka/Ks) from codon-aligned sequences,
  using the Nei-Gojobori (1986) or Li (1993)/LPB93 model.
- **`eskaks vcf`**: **per-gene pN/pS** from population variants (VCF + reference
  FASTA + GFF3), with a genome-wide selection scan and a self-contained interactive
  HTML report.

It implements the classical substitution models with precomputed lookup tables,
achieving a **1,280× speedup** over KaKs_Calculator while staying numerically
accurate (R² = 1.0 for the Li model).

> [!NOTE]
> `eskaks vcf` estimates pN/pS from **within-population polymorphism** (variants
> segregating in your samples), which answers a different question than
> between-species dN/dS divergence. The report can reconcile the two side by side
> when you supply a divergence table with `--divergence`.

## Features

| Feature | Description |
|---|---|
| 🧬 Two dN/dS models | [Nei-Gojobori (1986) and Li (1993)/LPB93](docs/models.md), precomputed lookup tables |
| 🧪 Per-gene pN/pS | From VCF + reference + GFF3, plus a pooled **genome-wide** estimate |
| 🔬 Neutrality test | Exact binomial vs `pN/pS = 1`, with Benjamini-Hochberg FDR and Bonferroni |
| 📐 Spectrum-aware sites | `--kappa` ts/tv weighting (modified Nei-Gojobori) for biased genomes like *M. tuberculosis* |
| 🧮 McDonald-Kreitman | `--mk` per-gene fixed/polymorphic table, Neutrality Index, α, Fisher exact p |
| 📊 Confidence intervals | Bootstrap 95% CI on genome-wide pN/pS, and per-gene Wilson CIs |
| 🧭 Genomic control | `--genomic-control` inflation-factor (λ) correction for clonal linkage |
| 🧹 Core-genome mode | `--exclude-repetitive` drops PE/PPE/PGRS/IS from the pooled estimate and the test |
| 🖥️ Interactive report | `--report` writes a self-contained dashboard (see [below](#interactive-html-report)) |
| 🎨 Genetic codes | [20 NCBI translation tables](docs/genetic-codes.md) (standard, mitochondrial, plastid, …) |
| 📁 Flexible outputs | TSV, CSV, JSON (`null` for NaN/Infinity), SVG plots, HTML report |
| ⚡ Parallel & fast | Rayon multi-threading; deterministic output regardless of thread count |
| 🔗 Pipeline-friendly | Stdin support, JSON output, non-zero exit on errors |

## Installation

Requires [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70.0.

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks
make release
cp target/release/eskaks ~/.local/bin/
```

`make release` enables native-CPU optimizations; `cargo install --path .` also
works. See [docs/installation.md](docs/installation.md) for details.

## Quick Start

**Try it right now** with the datasets in [`examples/`](examples/) (no data of your
own needed):

```bash
# dN/dS from aligned sequences
eskaks fasta examples/genes.fasta -o first_run

# per-gene pN/pS from a VCF + an interactive report
eskaks vcf --ref examples/toy_genome/reference.fasta \
  --gff examples/toy_genome/genes.gff3 \
  --vcf examples/toy_genome/variants.vcf \
  --genetic-code 11 --report -o toy_scan
# → open toy_scan_report.html in a browser
```

The [tutorial](docs/tutorial.md) explains every step and how to read the output; the
[quick-start guide](docs/quickstart.md) and [CLI reference](docs/cli-reference.md)
list more workflows and every flag.

> [!TIP]
> If your sequences are not codon-aligned, align them with [MAFFT](https://mafft.cbrc.jp/) + [PAL2NAL](http://www.bork.embl.de/pal2nal/) or [MACSE](https://bioweb.supagro.inra.fr/macse/) first.

## Interactive HTML report

`--report` writes a single, self-contained `.html` file (no network, no CDN) that
turns the per-gene table into a linked, explorable dashboard with a sticky table of
contents: a genome-wide verdict banner, summary cards, an interactive Manhattan /
volcano / p-value QQ, a McDonald-Kreitman panel, a polymorphism-vs-divergence
reconciliation, a power funnel, and the pN/pS distribution, each with a small
**"i"** button explaining how to read it. A colour-blind mode (Okabe-Ito palette +
direction shapes), light/dark themes, CSV/JSON export, and canvas rendering for
whole genomes are built in.

See the [VCF analysis guide](docs/vcf-analysis.md#interactive-html-report) for the
full panel-by-panel tour.

## Example Output

```text
Gene     Chrom  Start   Strand  pN     pS     pN/pS  N    S   SNPs  P_value   Q_value_BH
katG     chr1   2153889 -       0.021  0.010  2.05   12   3   15    3.1e-04   1.2e-02
rpoB     chr1   759807  +       0.004  0.031  0.14   4    22  26    2.6e-05   1.5e-03
Rv0001   chr1   1       +       0.009  0.012  0.75   6    9   15    0.34      0.61
```

**Reading it:** `pN/pS > 1` with a significant q-value flags **diversifying /
positive** selection (e.g. drug targets, antigens); `pN/pS < 1` flags **purifying**
selection (conserved, essential genes); a non-significant gene is usually just
**underpowered** (too few SNPs), not neutral. More in
[interpreting results](docs/interpreting-results.md).

## Performance

<p align="center">
  <img src="benchmarks/plots/performance_bars.png" width="700" alt="Performance comparison">
</p>

Up to **2,641× faster** than KaKs_Calculator, with **R² = 1.0** accuracy on the Li
model. Full benchmarks and a feature-by-feature comparison are in
[docs/performance.md](docs/performance.md).

## Documentation

| Document | Description |
|---|---|
| [**Tutorial**](docs/tutorial.md) | **Start here**: a hands-on walkthrough with example data |
| [Quick Start](docs/quickstart.md) | Copy-paste commands for common workflows |
| [Glossary](docs/glossary.md) | Plain-language definitions of every term |
| [VCF Analysis (pN/pS)](docs/vcf-analysis.md) | pN/pS per gene, neutrality test, MK, genomic control, SFS, the report |
| [Models](docs/models.md) | Nei-Gojobori vs Li: formulas, differences, when to use each |
| [Genetic Codes](docs/genetic-codes.md) | 20 NCBI translation tables with examples |
| [Interpreting Results](docs/interpreting-results.md) | What the numbers mean, pitfalls, and limitations |
| [Output Formats](docs/output-formats.md) | TSV, CSV, JSON and SVG plots |
| [CLI Reference](docs/cli-reference.md) | All flags, exit codes, examples |
| [Performance & Accuracy](docs/performance.md) | Benchmarks and feature comparison |
| [Development](docs/development.md) | Build, test, and source layout |
| [FAQ](docs/faq.md) | Speed, NaN, stop codons, library usage |
| [Changelog](CHANGELOG.md) | Version history |

## Citation

If you use eskaks in your research, please cite:

> Ruiz-Rodriguez P, Coscollá M. **eskaks: fast pairwise dN/dS and per-gene pN/pS from sequences or VCFs.** https://github.com/PathoGenOmics-Lab/eskaks

```bibtex
@software{ruiz-rodriguez_eskaks_2026,
  title   = {eskaks: fast pairwise dN/dS and per-gene pN/pS from sequences or VCFs},
  author  = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year    = {2026},
  url     = {https://github.com/PathoGenOmics-Lab/eskaks},
  version = {1.4.0},
  license = {GPL-3.0}
}
```

## License

[GNU General Public License v3.0](LICENSE)

---

<h2 id="contributors" align="center">

✨ Contributors
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
eskaks is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
