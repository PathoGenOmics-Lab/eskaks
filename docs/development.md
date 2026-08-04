---
description: >-
  Contributing to eskaks — repository layout, building and testing, the docs
  workflow, and how the CLI and library fit together.
---

# Development

The CLI and library live in `src/`; the user-facing documentation is in `docs/`.

## Build & test

```bash
cargo test                                  # full test suite (incl. the golden snapshot)
cargo clippy --all-targets -- -D warnings   # lint (CI denies warnings)
make release                                # optimized (native-CPU) build
make docs                                   # build the docs with MkDocs Material
```

Enable the repo git hooks once per clone. A pre-commit check then blocks common
slips (em-dashes, a leftover `dbg!`, a conflict marker) before they land:

```bash
make hooks   # sets core.hooksPath to .githooks
```

`tests/golden.rs` freezes the exact `eskaks vcf` output on the bundled toy genome
(every table, TSV and JSON), so any silent drift fails CI with an exact diff. After
an *intended* output change, regenerate the snapshot and review the diff before
committing:

```bash
BLESS=1 cargo test --test golden
```

## Source layout

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

eskaks can also be used as a Rust library crate rather than through the CLI; see the
[FAQ](faq.md) for a minimal example.
