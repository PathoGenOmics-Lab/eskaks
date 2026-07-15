# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **Genomic-control clonality correction** (`eskaks vcf --genomic-control`): the
  per-gene neutrality χ² is divided by the median-based inflation factor λ
  (never below 1) and re-tested, adding `p_gc`/`q_gc` columns. λ is always
  reported (in the summary and the report's QQ panel/inflation card) as a
  diagnostic — in clonal organisms like *M. tuberculosis* genome-wide linkage
  inflates significance, so the raw per-gene test is anti-conservative.
- **`--exclude-repetitive` core-genome mode**: drops PE/PPE/PGRS, IS elements and
  maturases from the genome-wide pooled estimate and the neutrality-test family
  (their SNP calls are frequently mapping artefacts). The report always shows a
  **core-vs-repetitive** pooled comparison so the gap is visible either way.
- **Per-gene 95% confidence interval on pN/pS** (Wilson score on the
  nonsynonymous SNP fraction, mapped to pN/pS): new `pn_ps_lo`/`pn_ps_hi` output
  and CI whiskers on the report's power funnel.
- **Site-frequency-spectrum panel**: pN/pS split by allele-frequency bin — a
  falling profile is the signature of purifying selection keeping deleterious
  nonsynonymous variants rare (informative for multi-sample cohorts).
- **Log-space −log10(p)**: the neutrality test now also reports a `−log10(p)`
  computed in log space, so genes whose exact p underflows to 0 keep a finite,
  meaningful value.
- **Report scales to whole genomes**: scatter panels render to a `<canvas>` with
  nearest-point hover/click above ~1200 genes, and the per-gene table is
  virtualized — thousands of genes stay responsive.
- **Colour-blind (CVD) mode** in the report: a toggle swaps to a validated
  Okabe-Ito palette **and** adds direction shapes (▲ diversifying / ▼ purifying /
  • not significant) so selection direction never depends on colour alone.
- **Provenance block** in the report: eskaks version, the invoking command line,
  and the input file paths, for reproducibility.
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
- **Per-gene neutrality test** for `eskaks vcf`: a two-sided exact binomial test
  of H0 pN/pS = 1 (observed nonsynonymous SNPs vs the mutational opportunity
  `N/(N+S)`), with **Benjamini-Hochberg FDR** q-values and **Bonferroni**-corrected
  p-values across genes. New `--fdr` threshold; the summary reports significant
  genes and the Manhattan plot outlines them. Skipped under `--af-weighted`.
- **Enriched per-gene output**: added `Chrom`, `Start`, `End`, `Strand`,
  `Exp_N_frac`, `P_value`, `Q_value_BH`, and `P_Bonferroni` columns (appended, so
  existing column positions are unchanged).
- **`--min-snps`** drops low-count genes from the per-gene table, plot, and test
  (the genome-wide pooled estimate still uses all genes).
- **Parallel `eskaks vcf`**: the per-gene computation now runs on `--workers`
  threads (default 4); output is deterministic regardless of thread count.
- **McDonald-Kreitman test** (`--mk`, `--mk-fixed-af`): per-gene 2×2
  fixed/polymorphic table (Dn, Ds, Pn, Ps), Neutrality Index, alpha, and a
  two-sided Fisher exact p-value with BH-FDR q-value (`<prefix>_mk.<ext>`).
  Reference-polarized (fixed = high AF within the sample).
- **Nei-Gojobori analytic variance + Z-test** (`eskaks fasta --neutrality`):
  writes `<output>_pairwise_tests` with dN, dS, their standard errors, the NG
  neutrality Z statistic, and a two-sided p-value (Nei model; NaN for Li).
- **Bootstrap 95% CI** for the genome-wide pooled pN/pS (`eskaks vcf --bootstrap`,
  `--seed`), by resampling genes with replacement.
- **Per-pair bootstrap CIs** for the fasta path (`eskaks fasta --bootstrap`):
  95% CIs on dN, dS, and dN/dS by resampling codon columns — model-agnostic, so
  it covers the Li model too (`<output>_pairwise_bootstrap`).
- **−log10(p) Manhattan plot** (`eskaks vcf --plot` also writes
  `<prefix>_pvalue_manhattan.svg`): genome-position significance scan with a
  Benjamini-Hochberg line at `--fdr` and significant genes highlighted.
- **Interactive HTML report** — a self-contained page (no CDN/network) with
  dynamic, hover-interactive visualizations:
  - `eskaks vcf --report`: a **linked multi-panel dashboard** — a genome-wide
    **verdict banner**, a selection-regime **census**, a **significant-hits
    shortlist**, Manhattan (with a `−log10(p)` / `pN/pS` / **`z(N)`** metric
    toggle), volcano (log2 pN/pS vs −log10 p), a **p-value QQ plot with genomic
    inflation λ**, a McDonald-Kreitman α-vs-significance volcano (with `--mk`), a
    **top-genes lollipop**, a power funnel (pN/pS vs SNP count), an observed-vs-
    expected N-fraction diagnostic, the pN/pS distribution, and an enriched
    per-gene table (now including a power-aware **`z(N)`** standardized
    nonsynonymous-excess column and, with `--mk`, a **`DoS`** Direction-of-
    Selection column). Every panel and summary card carries an **“i”
    interpretation-help popover** (what it shows, how to read the axes, what to
    watch out for), plus a **“How to read this report”** glossary. Clicking any
    point or row highlights that gene across every panel; **↑/↓** step through
    genes and **Esc** clears; a global **FDR(BH) ↔ Bonferroni** stringency toggle
    repaints all panels; light/dark theme toggle; CSV/JSON export (RFC 4180) of
    the filtered view; a **Print / Save-PDF** button; a **sticky table of contents**
    down the left side with scroll-spy; repetitive-gene (PE/PPE) badges; a
    Methods/parameters block. Uses the PathoGenOmics-Lab **mycolorsTB** palette
    (pathogenomics brand + *M. tuberculosis* lineage colors).
  - `eskaks vcf --divergence <FILE>`: adds a **polymorphism-vs-divergence**
    reconciliation panel to the report — per-gene pN/pS (within-sample) vs a
    supplied per-gene dN/dS (divergence), matched by name, with the concordance
    diagonal and quadrant interpretation (past-positive / diversifying /
    purifying / relaxed).
  - `eskaks fasta --report`: a **multi-panel dashboard** — a positional
    **sliding-window dN/dS "Manhattan"** (for aligned input), a **dN-vs-dS
    scatter** (one point per pair), and the **pairwise dN/dS distribution** are
    always shown, plus a **lineage strip-scatter** (one point per genome,
    per-lineage mean as a bar) under `--lineage` and a **group mean ± 95% CI
    scatter** under `--group-average`.

### Performance
- `eskaks vcf` scales to large genomes and cohorts: the per-gene SNP scan is
  now a binary-searched genomic window (O(S + G·log S) instead of O(G·S)), and
  `merge_vcfs` parses per-sample VCFs in parallel — both byte-identical output.

### Fixed
- Data-quality warnings (REF mismatches, skipped genes, saturation) are shown by
  default; added `-v`/`-vv`/`-q`. Validate AF-filter ranges and reconcile
  contig names up front (both previously produced a silent all-NaN output).
  Empty VCFs in a merge no longer abort the run. Aggregate REF-mismatch
  diagnostics instead of one line per SNP; warn on non-multiple-of-3 CDS and on
  genes dropped from the plot. Report the pooled `mean(dN)/mean(dS)` in the
  dN/dS summary instead of only the biased mean of per-pair ratios.

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
