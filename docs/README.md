# Documentation

Detailed documentation for eskaks. Browse on GitHub or build locally with [mdbook](https://rust-lang.github.io/mdBook/).

## Contents

| Document | Description |
|---|---|
| [**Getting started (tutorial)**](tutorial.md) | **Start here**: hands-on walkthrough with bundled example data |
| [Installation](installation.md) | Build from source and install |
| [Quick Start](quickstart.md) | Copy-paste commands for common workflows |
| [Glossary](glossary.md) | Plain-language definitions of every term |
| [Models](models.md) | Nei-Gojobori (1986) vs Li (1993), formulas, differences, when to use each |
| [Genetic Codes](genetic-codes.md) | 20 supported NCBI translation tables with usage examples |
| [Output Formats](output-formats.md) | TSV, CSV, JSON, output modes and SVG plot generation |
| [Interpreting Results](interpreting-results.md) | What dN/dS values mean, common pitfalls, limitations |
| [CLI Reference](cli-reference.md) | Complete flag list, exit codes, environment variables, examples |
| [VCF Analysis (pN/pS)](vcf-analysis.md) | Compute pN/pS per gene from VCF + reference + GFF3 |
| [Performance & Accuracy](performance.md) | Benchmarks and a feature-by-feature comparison |
| [FAQ](faq.md) | Why so fast? NaN values? Internal stop codons? Library usage? |
| [Development](development.md) | Build, test, and source layout |
| [Benchmarks](../benchmarks/) | Accuracy validation and performance comparison |

## Building locally

```bash
make docs          # Build HTML → docs/book/
make docs-serve    # Build + serve with live reload
```

Requires [mdbook](https://rust-lang.github.io/mdBook/):
```bash
cargo install mdbook
```
