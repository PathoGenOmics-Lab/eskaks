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
genes. It runs in two modes:

- **`eskaks fasta`**: **pairwise dN/dS** (Ka/Ks) from codon-aligned sequences,
  using the Nei-Gojobori (1986) or Li (1993)/LPB93 model.
- **`eskaks vcf`**: **per-gene pN/pS** from population variants (VCF + reference
  FASTA + GFF3), with a genome-wide scan for selection and a full interactive
  report.

It implements the classical substitution models with precomputed lookup tables,
achieving a **1,280× speedup** over KaKs_Calculator while staying numerically
accurate (R² = 1.0 for the Li model).

The tool takes:

- **Sequences**: a codon-aligned FASTA (`fasta` mode), **or**
- **Variants**: one or more VCFs + a reference FASTA + a GFF3 annotation (`vcf` mode)

It writes tables (TSV/CSV/JSON), SVG plots, and a self-contained interactive HTML
report: no internet, no CDN, works offline on an HPC login node.

> [!NOTE]
> `eskaks vcf` estimates pN/pS from **within-population polymorphism** (variants
> segregating in your samples), which answers a different question than
> between-species dN/dS divergence. The report can reconcile the two side by side
> when you supply a divergence table with `--divergence`.

**Main features:**

- Pairwise dN/dS with two classical models and 20 NCBI genetic codes
- Per-gene pN/pS with a mutation-spectrum-aware (`--kappa`) site count
- A per-gene neutrality test (exact binomial) with FDR/Bonferroni correction
- McDonald-Kreitman test, bootstrap CIs, and an optional genomic-control correction
- A self-contained, interactive HTML report with a colour-blind mode
- Scales from a handful of genes to whole genomes and large cohorts

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

`make release` enables native-CPU optimizations automatically. You can also build
with plain Cargo:

```bash
cargo install --path .
```

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

The [tutorial](docs/tutorial.md) explains every step and how to read the output.

### Pairwise dN/dS (FASTA)

```bash
# Basic pairwise dN/dS (Nei model, 4 threads)
eskaks fasta input.fasta -o results

# Li model, 8 threads
eskaks fasta input.fasta --model li --workers 8 -o results

# Sliding window with an SVG plot
eskaks fasta input.fasta --window-size 100 --window-step 10 --plot -o windows

# From stdin
cat input.fasta | eskaks fasta - -o results
```

### Per-gene pN/pS (VCF)

```bash
# One VCF per sample, AF-weighted (πN/πS)
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf sample1.vcf --vcf sample2.vcf --vcf sample3.vcf \
  --af-weighted --genetic-code 11 -o population_pnps

# A single multi-sample VCF also works
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf population.vcf -o results
```

### Full selection scan with the interactive report

```bash
eskaks vcf \
  --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --kappa 2 \
  --mk --bootstrap 1000 --seed 42 \
  --report --plot -o mtb_scan
# → mtb_scan_pnps.tsv, mtb_scan_mk.tsv, mtb_scan_report.html, plots …
```

Run `eskaks fasta --help` or `eskaks vcf --help` for the complete list of options.

> [!TIP]
> If your sequences are not codon-aligned, align them with [MAFFT](https://mafft.cbrc.jp/) + [PAL2NAL](http://www.bork.embl.de/pal2nal/) or [MACSE](https://bioweb.supagro.inra.fr/macse/) first.

## Common Arguments

**`eskaks vcf`** (per-gene pN/pS):

| Argument | What it does | Default |
|---|---|---|
| `--ref <FILE>` | Reference FASTA. Contig names must match the VCF and GFF. | - |
| `--gff <FILE>` | Gene annotation (GFF3). | - |
| `--vcf <FILE>` | Variant file; repeat for one VCF per sample. | - |
| `--vcf-list <FILE>` | A text file listing VCF paths, one per line. | - |
| `--genetic-code <N>` | NCBI translation table (`11` for bacteria). | `1` |
| `--kappa <F>` | ts/tv rate ratio for spectrum-aware site counting. | `1.0` |
| `--af-weighted` | Weight SNPs by allele frequency (reports πN/πS). | off |
| `--min-af / --max-af / --min-depth` | Variant filters. | - |
| `--pass-only` | Keep only `FILTER=PASS` records. | off |
| `--min-snps <N>` | Drop genes with fewer than N SNPs from the table, plot, and test. | `0` |
| `--fdr <F>` | Benjamini-Hochberg threshold for calling genes significant. | `0.05` |
| `--mk` / `--mk-fixed-af <F>` | Run the McDonald-Kreitman test; "fixed" AF cutoff. | off / `0.99` |
| `--bootstrap <N>` / `--seed <N>` | Replicates for the genome-wide 95% CI; RNG seed. | `0` / `42` |
| `--genomic-control` | Divide each χ² by the inflation factor λ and re-test. | off |
| `--exclude-repetitive` | Drop PE/PPE/PGRS/IS genes from the pooled estimate and test. | off |
| `--divergence <FILE>` | Per-gene dN/dS table for the report's polymorphism-vs-divergence panel. | - |
| `--report` | Write a self-contained interactive HTML report. | off |
| `--plot` | Write SVG Manhattan / p-value plots. | off |
| `--format <tsv\|csv\|json>` | Table output format. | `tsv` |
| `--workers <N>` | Threads. Output is deterministic regardless. | `4` |

**`eskaks fasta`** (pairwise dN/dS): `--model <nei\|li>`, `--genetic-code <N>`,
`--window-size` / `--window-step`, `--lineage`, `--group-average`, `--neutrality`
(NG variance Z-test), `--bootstrap` / `--seed` (per-pair CIs), `--report`,
`--plot`, `--format`, `--workers`. See [docs/cli-reference.md](docs/cli-reference.md).

## Outputs

By default, `eskaks vcf` writes the per-gene table:

```text
<prefix>_pnps.tsv
```

Optional outputs: `<prefix>_mk.<ext>` (`--mk`), `<prefix>_report.html`
(`--report`), and `<prefix>_pnps_plot.svg` / `<prefix>_pvalue_manhattan.svg`
(`--plot`). A **genome-wide pooled pN/pS** is always printed to stderr.

The most important per-gene columns:

| Column | Meaning |
|---|---|
| `Gene`, `Chrom`, `Start`, `End`, `Strand` | Gene identity and location |
| `N_sites`, `S_sites` | Nonsynonymous / synonymous site counts (κ-weighted) |
| `pN`, `pS`, `pN/pS` | Per-site nonsyn/syn polymorphism and their ratio |
| `Nonsyn`, `Syn`, `SNPs` | Observed SNP counts |
| `Exp_N_frac` | Expected nonsynonymous fraction `N/(N+S)` under neutrality |
| `P_value`, `Q_value_BH`, `P_Bonferroni` | Neutrality-test p-value and corrections |

With `--mk`, extra columns report `Dn`, `Ds`, `Pn`, `Ps`, the Neutrality Index,
α, and the Fisher exact p-value. See [docs/output-formats.md](docs/output-formats.md).

## Interactive HTML report

`--report` writes a single, self-contained `.html` file (no network, no CDN) that
turns the per-gene table into a linked, explorable dashboard with a **sticky table
of contents** down the left side. Every panel and summary card carries a small
**“i”** button that explains how to read it.

- **Genome-wide verdict** banner, summary cards (pooled pN/pS + bootstrap CI,
  significant-gene count, inflation λ), and a **“How to read this report”** glossary
- **Selection-regime census** and a **significant-hits shortlist** (click to filter)
- **Manhattan** (`−log10(p)` / `pN/pS` / `z(N)` toggle), **volcano**, and a **p-value
  QQ** plot with the genomic-inflation factor λ
- **McDonald-Kreitman** (α vs significance, with `--mk`), a **polymorphism-vs-
  divergence** reconciliation (with `--divergence`), a **power funnel** with per-gene
  CI whiskers, an **observed-vs-expected** diagnostic, a **top-genes lollipop**, an
  **allele-frequency spectrum**, and the **pN/pS distribution**
- Click any point or row to highlight that gene across **every** panel; **↑/↓** step
  through genes; a global **FDR ↔ Bonferroni** toggle repaints everything
- **🎨 Colour-blind mode**: a toggle swaps to a validated Okabe-Ito palette and adds
  direction shapes (▲ diversifying / ▼ purifying / ● not significant) so meaning
  never depends on colour alone; light/dark theme toggle; CSV/JSON export; Print/PDF
- **Scales to whole genomes**: scatter panels switch to canvas rendering and the
  table is virtualized above ~1200 genes, so a full genome stays responsive

The report uses the PathoGenOmics-Lab **mycolorsTB** palette. `eskaks fasta --report`
produces a matching dashboard for the dN/dS workflow.

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
**underpowered** (too few SNPs), not neutral.

## Performance

<p align="center">
  <img src="benchmarks/plots/performance_bars.png" width="700" alt="Performance comparison">
</p>

| Dataset | eskaks (4t) | KaKs_Calculator | PAML yn00 | BioPython | Speedup |
|---|---|---|---|---|---|
| 20 seq × 300 bp | 2 ms | 34 ms | 8 ms | 610 ms | 17× |
| 100 seq × 3 kb | 6 ms | 7,703 ms | 697 ms | 111,619 ms | **1,280×** |
| 500 seq × 3 kb | 74 ms | 195,456 ms | - | - | **2,641×** |

Li model achieves **R² = 1.0** vs KaKs_Calculator LPB. Full accuracy data and
methodology in [benchmarks/](benchmarks/).

## Comparison

| | eskaks | KaKs_Calculator | BioPython | PAML yn00 |
|---|---|---|---|---|
| Nei-Gojobori model | ✅ | ✅ | ✅ | ✅ |
| Li/LPB93 model | ✅ | ✅ | ❌ | ❌ |
| Per-gene pN/pS from VCF | ✅ | ❌ | ❌ | ❌ |
| Neutrality test + FDR | ✅ | ❌ | ❌ | ❌ |
| Interactive HTML report | ✅ | ❌ | ❌ | ❌ |
| Custom genetic codes | ✅ (20 tables) | ❌ | ❌ | Limited |
| JSON output / stdin pipe | ✅ | ❌ | ❌ | ❌ |
| Parallel | ✅ | ❌ | ❌ | ❌ |
| Speed (100 seq) | **6 ms** | 7,703 ms | 111,619 ms | 697 ms |

## Documentation

| Document | Description |
|---|---|
| [**Tutorial**](docs/tutorial.md) | **Start here**: a hands-on walkthrough with example data |
| [Glossary](docs/glossary.md) | Plain-language definitions of every term |
| [VCF Analysis (pN/pS)](docs/vcf-analysis.md) | pN/pS per gene, neutrality test, MK, genomic control, SFS, the report |
| [Models](docs/models.md) | Nei-Gojobori vs Li: formulas, differences, when to use each |
| [Genetic Codes](docs/genetic-codes.md) | 20 NCBI translation tables with examples |
| [Interpreting Results](docs/interpreting-results.md) | What dN/dS and pN/pS mean, common pitfalls |
| [Output Formats](docs/output-formats.md) | TSV, CSV, JSON and SVG plots |
| [CLI Reference](docs/cli-reference.md) | All flags, exit codes, examples |
| [FAQ](docs/faq.md) | Speed, NaN, stop codons, library usage |
| [Changelog](CHANGELOG.md) | Version history |
| [Benchmarks](benchmarks/) | Accuracy validation and performance |

## For Developers

The CLI and library live in `src/`; detailed docs are in `docs/`.

```bash
cargo test          # run the full test suite (268 tests)
cargo clippy --all-targets -- -D warnings
make release        # optimized build
make docs           # build the docs with mdbook
```

```text
src/
├── main.rs           # orchestration + subcommand dispatch
├── cli.rs            # CLI definitions (clap subcommands)
├── input.rs          # FASTA reading, validation, stdin
├── compute.rs        # ComputeEngine (Nei | Li)
├── genetic_code.rs   # 20 NCBI tables
├── stats.rs          # binomial/Wilson/Fisher, FDR, bootstrap, probit
├── vcf.rs / gff.rs   # VCF and GFF3 parsers
├── vcf_analysis.rs   # per-gene pN/pS, neutrality test, MK, genomic control
├── report.rs         # self-contained interactive HTML report
├── plot.rs           # SVG generation
└── models/           # nei.rs, li.rs
```

## Limitations

- `eskaks fasta` expects **codon-aligned** input (in-frame, gap lengths multiples
  of 3); align with MAFFT + PAL2NAL or MACSE first.
- `eskaks vcf` uses **SNPs only**; indels and multi-nucleotide variants are not
  codon-annotated (see [get_MNV](https://github.com/PathoGenOmics-Lab/get_MNV) for MNVs).
- pN/pS is estimated from **within-sample polymorphism**, so many genes are
  **underpowered** in low-diversity organisms; a non-significant result is not
  evidence of neutrality.
- The per-gene neutrality test assumes independent SNPs. In **clonal** organisms
  (e.g. *M. tuberculosis*) genome-wide linkage inflates significance; `--genomic-
  control` is a pragmatic correction, but a high λ can also reflect **real** pervasive
  selection, so apply it only when you suspect systematic bias.
- Contig names must match across the VCF, reference FASTA, and GFF3.

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
