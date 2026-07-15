# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **Genome-wide (pooled) pN/pS** in the `eskaks vcf` summary: an overall estimate
  computed by pooling SNP counts and site counts across all genes
  (`pN = Σ nonsyn / Σ N_sites`, `pS = Σ syn / Σ S_sites`) rather than averaging
  per-gene ratios. Respects `--af-weighted` (reported as πN/πS) and includes a
  coarse selection label (purifying / near-neutral / diversifying).
- **`--kappa` mutation-spectrum-aware site counting** for `eskaks vcf`: a
  transition/transversion rate ratio that weights synonymous vs nonsynonymous
  site counting (modified Nei-Gojobori). Corrects the equal-rates bias that
  understates pN/pS under transition-skewed spectra (e.g. *M. tuberculosis*):
  `kappa > 1` raises S and lowers N, moving a 2-fold synonymous site from `1/3`
  to `kappa/(kappa+2)` while leaving 4-fold sites and the observed-SNP numerators
  unchanged. `--kappa 1` (default) reproduces the classic equal-rates counting.

## [1.4.0] - 2026-03-29

### Added
- **VCF subcommand** (`eskaks vcf`): Compute pN/pS per gene from VCF + reference FASTA + GFF3.
  - Classifies SNPs as synonymous or nonsynonymous using the reference codon context
  - Counts N and S sites per gene (fractional, based on all possible single-nucleotide changes)
  - Supports all 20 NCBI genetic codes via `--genetic-code`
  - Filters: `--pass-only`, `--min-af`, `--min-depth`
  - Output: TSV/CSV/JSON with gene-level pN, pS, pN/pS, SNP counts
  - Manhattan-style SVG plot via `--plot`
  - Multi-exon genes, minus strand, phase handling
- New modules: `vcf.rs` (VCF parser), `gff.rs` (GFF3 parser), `vcf_analysis.rs` (pN/pS computation)
- 16 integration tests for VCF subcommand, 12 unit tests in new modules

### Changed
- CLI restructured to subcommands: `eskaks fasta <input>` (original) and `eskaks vcf` (new)
- `eskaks --list-codes` still works at top level
- All existing tests updated for subcommand syntax

## [1.3.0] - 2026-03-29

### Added
- **Stdin support**: Read FASTA from stdin via `-` or `/dev/stdin` (`cat seqs.fasta | eskaks -`).
- **JSON output format**: `--format json` produces a JSON array of objects with `seq1`, `seq2`, `dN`, `dS`, `dN_dS` keys. NaN/Infinity values are serialized as `null`.
- **Internal stop codon warnings**: Detects premature stop codons (excluding the terminal codon) and warns about potential frameshifts or pseudogenes.
- **Makefile**: `make benchmark` runs the full benchmark pipeline (generate → run → plot). Also: `make test`, `make clippy`, `make check`, `make release`.

## [1.2.0] - 2026-03-29

### Changed
- **Architecture refactor**: Split monolithic `main.rs` (391 lines) into 5 focused modules:
  - `cli.rs` — CLI argument definitions (clap derive), cleanly separated from logic
  - `input.rs` — FASTA reading, validation, filtering, deduplication (`SequenceData` struct)
  - `compute.rs` — `ComputeEngine` enum for model-agnostic dispatch (eliminates duplicated closures)
  - `main.rs` — Thin orchestration (~100 lines): parse → load → compute → dispatch
  - `dispatch_output()` and `dispatch_window()` replace inline match blocks
- Total: 11 source modules (was 8), each with a single clear responsibility

## [1.1.0] - 2026-03-29

### Added
- **20 NCBI genetic code tables** (`--genetic-code <N>`): Standard, Vertebrate/Invertebrate
  Mitochondrial, Yeast Mito, Bacterial/Plastid, Ciliate, Echinoderm, Euplotid, and 12 more.
- `--list-codes` flag to list all available translation tables.
- New `genetic_code` module with Li↔Nei index conversion and dynamic synonymous site computation.
- Property-based tests (proptest): symmetry, non-negativity, identity, NaN handling for both models.
- `OutputConfig` struct to reduce function argument count in output module.
- `LineagePlotResult` type alias for cleaner return types.
- SAFETY comments on all 7 `unsafe` blocks.
- `lib.rs` for library-level access to models and utilities.

### Fixed
- Version mismatch: Cargo.toml and CLI now both report `1.1.0`.
- 19 clippy warnings resolved (0 remaining).
- Documented that Nei pathway fallback values (2-diff: 0.5/1.5, 3-diff: 1.0/2.0)
  are unreachable with standard genetic code but triggered by alternative codes.

### Changed
- `NeiTables` and `LiTables` constructors accept any genetic code via `with_genetic_code()`.
- Output functions now take `&OutputConfig` instead of 5+ separate parameters.
- Removed unused `_n_u` parameter from `write_group_average`.

## [1.0.0] - 2026-03-28

### Added
- Initial release with Nei-Gojobori (1986) and Li (1993)/LPB93 models.
- Pairwise, lineage, group average, and sliding window output modes.
- L1-cache-optimized lookup tables (32 KB Nei, 288 KB Li).
- Fast-path for identical codons (~95% of comparisons).
- Multi-threaded computation via rayon.
- TSV/CSV output with optional SVG plots.
- 21 integration tests.
- Benchmarks: 1,280× faster than KaKs_Calculator, 18,600× faster than BioPython.
- Accuracy: R²=1.0 vs KaKs_Calculator (Li), R²≈0.995-0.999 (Nei).
