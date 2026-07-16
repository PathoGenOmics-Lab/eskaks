# Development

The CLI and library live in `src/`; the user-facing documentation is in `docs/`.

## Build & test

```bash
cargo test                                  # run the full test suite
cargo clippy --all-targets -- -D warnings   # lint
make release                                # optimized (native-CPU) build
make docs                                   # build the docs with mdbook
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
