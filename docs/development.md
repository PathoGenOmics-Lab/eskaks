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

## Quality tooling

`cargo test` already runs, beyond the unit and integration suites:

- **Property tests** (`tests/property_tests.rs`, proptest) that assert invariants over
  generated inputs: the dN/dS models, the diversity statistics (pi is
  polarization-invariant, Tajima's D is finite-or-NaN), the distribution helpers
  (Wilson, binomial, Benjamini-Hochberg), and that the VCF and GFF3 parsers never
  panic on arbitrary bytes.
- **A docs-code contract** (`tests/docs_contract.rs`): every pN/pS output column and
  every `eskaks vcf`/`eskaks fasta --help` flag must be documented, so the docs
  cannot silently drift from the binary.

Two deeper tools run out of band (their CI workflows are scheduled weekly, never a
PR gate):

```bash
# Coverage-guided fuzzing of the parsers (needs a nightly toolchain + cargo-fuzz).
cargo +nightly fuzz run parse_vcf     # or parse_gff3
```

```bash
# Mutation testing: inject bugs and confirm the suite kills them (cargo install cargo-mutants).
cargo mutants --file src/stats/diversity.rs   # one module, fast
cargo mutants                                 # whole scientific core (slow)
```

Mutation testing surfaces *test gaps* (a surviving mutant is a bug the suite did not
catch); review each survivor, but expect some to be genuinely *equivalent* mutants
(e.g. flipping `||` to `&&` between two logically equivalent guards) that no test can
kill.

## Source layout

```text
src/
├── main.rs           # subcommand dispatch + --demo
├── cli.rs            # CLI definitions (clap subcommands)
├── run_fasta.rs      # `fasta` orchestration (pairwise / lineage / group / window)
├── run_vcf.rs        # `vcf` orchestration (pN/pS, diversity, MK, report)
├── input.rs          # FASTA reading, validation, stdin
├── compute.rs        # ComputeEngine (Nei | Li)
├── genetic_code.rs   # 20 NCBI tables
├── stats/            # dist.rs (binomial/Wilson/Fisher, FDR, probit), diversity.rs, accum.rs
├── vcf/, gff.rs      # VCF (parse / merge / filter) and GFF3 parsers
├── vcf_analysis/     # per-gene pN/pS, neutrality test, MK, genomic control, diversity
├── report.rs         # self-contained interactive HTML report
├── plot/             # SVG generation (bars, histogram, window)
└── models/           # nei.rs, li.rs
```

eskaks can also be used as a Rust library crate rather than through the CLI; see the
[FAQ](faq.md) for a minimal example.
