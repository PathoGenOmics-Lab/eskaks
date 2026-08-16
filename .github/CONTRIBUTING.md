# Contributing to eskaks

Thank you for your interest in contributing to **eskaks**! We welcome contributions that help improve the tool, fix bugs, or enhance the documentation. This document outlines the process for contributing and the guidelines to follow.

The full developer guide lives on the documentation site, at
[Development](https://pathogenomics-lab.github.io/eskaks/development/). This file is the
short, GitHub-facing version: what to run before you open a pull request, and the few house
rules that are easy to break by accident. Where the two overlap, the site is the fuller
explanation; where they differ, fix the site.

## Table of Contents
- [Project Description](#project-description)
- [How to Contribute](#how-to-contribute)
  - [Reporting Issues](#reporting-issues)
  - [Feature Requests](#feature-requests)
  - [Contributing Code](#contributing-code)
- [Coding Guidelines](#coding-guidelines)
- [Setting Up the Development Environment](#setting-up-the-development-environment)
- [Testing and Quality Gates](#testing-and-quality-gates)
- [Building the Documentation](#building-the-documentation)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Code of Conduct](#code-of-conduct)

## Project Description

`eskaks` measures **natural selection on protein-coding genes**. It has two subcommands:

- `eskaks fasta`: pairwise **dN/dS (Ka/Ks)** from codon-aligned FASTA, with the
  Nei-Gojobori (1986) and Li (1993)/LPB93 models.
- `eskaks vcf`: per-gene **pN/pS** from a VCF plus a reference FASTA and a GFF3
  annotation, with an exact-binomial neutrality test, FDR and Bonferroni correction,
  diversity statistics, a McDonald-Kreitman test, and a self-contained interactive
  HTML report.

It is a **single Rust crate** (package `eskaks`), pure Rust with no C dependencies and no
external tools to install. That is deliberate: the whole point of the project is that a
reviewer can clone it, run `cargo test`, and reproduce every number without a bioinformatics
stack around it. Please keep it that way, and think twice before a change adds a dependency
that pulls in a C toolchain.

### Scope and Limitations

`eskaks` consumes **existing** alignments and **existing** variant calls. It is not a
variant caller, not an aligner, and it does not re-estimate genotypes, ploidy, or phase.
`eskaks fasta` expects sequences that are already codon-aligned (length a multiple of
three, frame preserved); it will not fix a misaligned input for you, it will only report
the selection signal implied by the alignment you gave it.

The statistical limits matter as much as the input limits. pN/pS is a ratio of raw
polymorphism proportions: it does **not** correct for multiple substitutions (no
Jukes-Cantor, no Kimura), so it is not interchangeable with a rate-based dN/dS. Read
[Interpreting results](https://pathogenomics-lab.github.io/eskaks/interpreting-results/)
before proposing a change to what the numbers mean, and
[Models](https://pathogenomics-lab.github.io/eskaks/models/) before touching the estimators
themselves.

## How to Contribute

We appreciate all contributions, whether it is fixing bugs, proposing new features,
improving the documentation, or suggesting a new direction for the tool.

### Reporting Issues

If you encounter a bug, have a question, or want to request a feature, please
[open an issue](https://github.com/PathoGenOmics-Lab/eskaks/issues/new/choose) on our
GitHub repository. When reporting an issue, please include:

- A detailed description of the problem.
- The exact command you ran and the steps to reproduce it.
- Any relevant logs or screenshots. Re-running with `-vv` turns on debug logging and
  usually says far more about *why* a record was skipped than the default output does.
- Version information: paste the output of `eskaks --version` (it embeds the git commit,
  which tells us whether you are on a release or on a local build) and your operating
  system.
- A minimal input that reproduces it, if you can share one. Selection bugs are almost
  always specific to one gene, one codon, or one odd VCF record; a two-line VCF that fails
  is worth more than a whole genome that fails somewhere.

### Feature Requests

We welcome suggestions for new features and improvements! Please open an issue labeled
**Feature Request** and provide as much detail as possible regarding your suggestion and
its potential use cases. For anything statistical, say which paper or method you expect us
to follow: matching a published estimator is a much easier request to act on than a general
wish for "better statistics".

### Contributing Code

If you would like to contribute code, follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix (`git checkout -b feature/new-feature`).
3. Make your changes.
4. Run `make check` (clippy with warnings denied, then the full test suite).
5. Submit a pull request following the guidelines in the **Submitting a Pull Request**
   section.

## Coding Guidelines

- **Do not run `cargo fmt` on this repository.** This is the single easiest way to make a
  pull request unreviewable here. eskaks is not rustfmt-formatted: a repo-wide `cargo fmt`
  rewrites almost every file, and your three-line fix disappears into a thousand lines of
  reflowed whitespace. There is deliberately no `cargo fmt --check` step in CI for the same
  reason. Match the style of the code around your change instead.
- **No em-dashes, anywhere.** Not in code, not in comments, not in docs, not in commit
  messages. Use a comma, a colon, a semicolon, or a full stop. The pre-commit hook
  (see below) rejects a commit that adds one, so this is enforced rather than merely
  requested.
- **English only** for code, comments, documentation, and commit messages.
- **Lint**: `cargo clippy --all-targets -- -D warnings` must pass. CI denies warnings on
  both Linux and macOS, so a clippy lint that only fires on one of them still fails the
  build.
- **Testing**: every behavioural change needs a test. `cargo test` runs everything (see
  the next section); if your change is not observable in any test, that is usually a sign
  the change is either dead code or under-specified.
- **Documentation**: document public functions, structs, and modules with doc comments
  (`///`). If your change adds an output column or a CLI flag, you must also document it in
  `docs/`, because `tests/docs_contract.rs` will fail the build otherwise.
- **MSRV**: the minimum supported Rust version is **1.85** (`rust-version` in
  `Cargo.toml`). Nothing in CI enforces it yet, so please check by hand that you are not
  reaching for a newer language or standard-library feature.

## Setting Up the Development Environment

1. Clone the repository:

   ```bash
   git clone https://github.com/PathoGenOmics-Lab/eskaks.git
   cd eskaks
   ```

2. Enable the repository git hooks. Do this **once per clone**, before your first commit:

   ```bash
   make hooks   # sets core.hooksPath to .githooks
   ```

   This points git at `.githooks/`, whose pre-commit hook runs `scripts/check-hygiene.sh`
   over the **added** lines of your staged diff and refuses the commit if it introduces an
   em-dash, a leftover debugging macro, or a merge-conflict marker. It only looks at added
   lines, so it never trips on pre-existing content, and it excludes the generated golden
   fixtures under `tests/golden/`, which are data rather than prose. If you ever need to get
   past it on purpose, `git commit --no-verify` is the escape hatch, but say why in the pull
   request.

   Note this page names that macro in prose rather than writing it out, because the hook
   reads added lines without knowing which are documentation, so spelling it out here would
   make this very file impossible to commit. Read `scripts/check-hygiene.sh` for the exact
   patterns it rejects.

3. Build:

   ```bash
   cargo build          # debug build
   make release         # optimized build (RUSTFLAGS="-C target-cpu=native")
   ```

   `make release` is what the benchmark numbers are produced with. Note that
   `-C target-cpu=native` produces a binary tuned for *your* machine, so do not ship it to
   a cluster with different hardware and expect it to run.

4. Run the tests to confirm the clone is healthy:

   ```bash
   cargo test
   ```

5. You are ready to start contributing!

## Testing and Quality Gates

`cargo test` is the whole suite: **700 tests** at the time of writing, no feature flags, no
external tools, no test data to download. It covers the unit tests in `src/`, the
integration tests in `tests/`, a golden snapshot, proptest property tests, and a
documentation contract. `make check` runs clippy first and then the same suite, which is
exactly what CI does.

Three parts of the suite are worth understanding before you trip over them:

**The golden snapshot (`tests/golden.rs`).** It runs `eskaks vcf` on the bundled toy genome
and freezes the exact output: every table (pN/pS, variants, diversity,
McDonald-Kreitman), in TSV and in JSON. The non-numeric skeleton (columns, keys, ordering,
formatting) is compared byte for byte, so a shifted column, a renamed key, a changed count,
or a stray token fails with a readable diff. Each numeric token, however, is compared with a
tiny relative tolerance (1e-9). That is not sloppiness: several columns (`p_value`,
`p_bonferroni`, `p_gc`, `q_gc_bh`) are computed through transcendental libm calls whose last
bit differs between macOS and glibc, so an exact byte compare would fail on one half of the
CI matrix even when nothing had actually drifted. The tolerance sits far below any real
change and well above that cross-platform noise.

If the golden test fails, first decide whether the output change was **intended**. If it
was, regenerate the snapshot and read the diff before committing it:

```bash
BLESS=1 cargo test --test golden
```

A blessed diff is part of your pull request and reviewers will read it. Never bless a diff
you cannot explain.

**The docs-code contract (`tests/docs_contract.rs`).** It requires that every pN/pS output
column and every `eskaks vcf` / `eskaks fasta --help` flag appears in the documentation.
The test exists because the docs had already drifted once, leaving four pN/pS columns and
several flags unmentioned. Adding a flag or a column is therefore a two-file change: the
code, and the page in `docs/` that describes it.

**Property tests (`tests/property_tests.rs`, proptest).** These assert invariants over
generated inputs rather than fixed cases: the dN/dS models, the diversity statistics, the
distribution helpers (Wilson, binomial, Benjamini-Hochberg), and the fact that the VCF and
GFF3 parsers never panic on arbitrary bytes. A proptest failure is written to
`tests/property_tests.proptest-regressions`; commit that file, since it pins the shrunk
counterexample so the bug can never come back unnoticed.

Two deeper tools run out of band. Neither is a pull request gate, and neither runs on a
schedule: their workflows (`fuzz.yml`, `mutants.yml`) are `workflow_dispatch` only, started
by hand from the Actions tab, so they cost nothing until someone asks for them. Locally:

```bash
cargo +nightly fuzz run parse_vcf   # or parse_gff3; needs a nightly toolchain + cargo-fuzz
cargo mutants                       # needs cargo install cargo-mutants; slow
```

A surviving mutant is a gap in the tests, not automatically a bug, and some survivors are
genuinely equivalent mutants that no test can kill. Review them, do not chase a score.

## Building the Documentation

The documentation site is built with
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/). The configuration is
`mkdocs.yml` at the repository root and the pages are the Markdown files in `docs/`. There
is **no `docs/requirements.txt`** in this repository: the dependencies are named explicitly
in the install command, both here and in `.github/workflows/docs.yml`, and the two must stay
in step. If you add a plugin, add it to `mkdocs.yml` *and* to the workflow's `pip install`
line, or the build will pass locally and fail in CI.

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install "mkdocs-material[imaging]" \
  mkdocs-git-revision-date-localized-plugin mkdocs-jupyter
.venv-docs/bin/mkdocs serve          # live preview at http://127.0.0.1:8000/eskaks/
.venv-docs/bin/mkdocs build --strict # production build into ./site
```

With the tools already on your `PATH`, `make docs` and `make docs-serve` are the shortcuts.
Be aware that `make docs` runs a plain `mkdocs build`, while the pull request gate runs
`mkdocs build --strict`, which turns warnings into errors: a broken link, a page missing
from the nav, or an unrecognised config key fails there and not in `make docs`. Run the
strict build yourself before pushing a documentation change.

### Publishing

The site is served at `https://pathogenomics-lab.github.io/eskaks/`. The **Docs** workflow
(`.github/workflows/docs.yml`) validates every pull request that touches `docs/`,
`includes/`, `mkdocs.yml`, or the workflow itself with a strict build, and deploys on push
to `main` (it can also be run manually).

Do not run `mkdocs gh-deploy` from your machine. The workflow is what deploys, and a local
deploy would publish whatever happens to be in your working tree, unreviewed and possibly
built with a different plugin set than CI uses.

## Submitting a Pull Request

1. Make sure `make check` passes locally: `cargo clippy --all-targets -- -D warnings` and
   then `cargo test`. CI runs both on Linux **and** macOS, so a change that only works on
   one of them is not done yet.
2. Write a clear commit message: a short subject line, then two or three sentences on what
   the change does and why. No em-dashes.
3. Open the pull request and fill in the template. Include a summary of the changes, why
   they are necessary, and any related issue numbers. If you blessed a new golden snapshot,
   say so and explain the diff.
4. A project maintainer will review your pull request and provide feedback if needed.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By
participating, you agree to abide by its terms. Please be respectful and professional in
all interactions.

We look forward to your contributions!
