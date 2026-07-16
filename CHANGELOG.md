# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - 2026-07-16

### Added

- **User-experience polish across the CLI:**
  - **`eskaks fasta` now prints a "Done" confirmation** (sequence count, model, and the
    list of files written). Previously a plain `eskaks fasta in.fasta` finished with no
    terminal output at all, leaving the user unsure it had worked or where the results
    went. Suppressed by `--quiet`.
  - **Clean, cargo-style log lines** (`warning: …`, `error: …`, coloured only on a
    terminal) instead of env_logger's default `[<ISO timestamp> LEVEL module::path] …`.
    `-v`/`-vv` and `RUST_LOG` still work.
  - **`--help` for each subcommand now groups flags into sections** (Output / Analysis /
    Statistics / Filtering / Input) and ends with concrete usage **examples**.
  - **Shell completions**: `eskaks --completions bash|zsh|fish|…` prints a completion
    script (via `clap_complete`).
- **Much clearer failure diagnostics**, so an empty or garbage result is never
  mistaken for a clean run:
  - `eskaks vcf` summary now accounts for the SNPs (`SNPs used (in CDS): X of Y`),
    lists the output files written, and warns when SNPs are read but *none* land in a
    CDS (a coordinate/build mismatch), when some VCF contigs match no gene (their SNPs
    are dropped), when many reference genes have an internal stop (likely the wrong
    `--genetic-code`), when `--min-snps` drops every gene, and prints `n/a` (not
    `0.000000`) when there is no pooled estimate.
  - `eskaks fasta` aggregates the internal-stop warnings into one line that names the
    likely wrong `--genetic-code`/frame, and `--lineage`/`--group-average` warn when
    only a single group is detected.
  - Gzip-compressed inputs fail fast with "decompress first" instead of a cryptic
    parse error; a FASTA passed as a VCF/GFF3 is reported as the wrong format (with a
    "is this a FASTA?" hint) instead of an all-NaN "success"; malformed-line warnings
    are capped so a wrong-format file no longer scrolls the terminal.
- The eskaks logo is now shown in the interactive HTML report header, embedded as a
  base64 `data:` URI so the report stays fully self-contained and offline-capable.

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
- **Stdin support**: Read FASTA from stdin via `-` or `/dev/stdin` (`cat seqs.fasta | eskaks -`).
- **JSON output format**: `--format json` produces a JSON array of objects with `seq1`, `seq2`, `dN`, `dS`, `dN_dS` keys. NaN/Infinity values are serialized as `null`.
- **Internal stop codon warnings**: Detects premature stop codons (excluding the terminal codon) and warns about potential frameshifts or pseudogenes.
- **Makefile**: `make benchmark` runs the full benchmark pipeline (generate → run → plot). Also: `make test`, `make clippy`, `make check`, `make release`.
- **20 NCBI genetic code tables** (`--genetic-code <N>`): Standard, Vertebrate/Invertebrate
  Mitochondrial, Yeast Mito, Bacterial/Plastid, Ciliate, Echinoderm, Euplotid, and 12 more.
- `--list-codes` flag to list all available translation tables.
- New `genetic_code` module with Li↔Nei index conversion and dynamic synonymous site computation.
- Property-based tests (proptest): symmetry, non-negativity, identity, NaN handling for both models.
- `OutputConfig` struct to reduce function argument count in output module.
- `LineagePlotResult` type alias for cleaner return types.
- SAFETY comments on all 7 `unsafe` blocks.
- `lib.rs` for library-level access to models and utilities.
- Initial release with Nei-Gojobori (1986) and Li (1993)/LPB93 models.
- Pairwise, lineage, group average, and sliding window output modes.
- L1-cache-optimized lookup tables (32 KB Nei, 288 KB Li).
- Fast-path for identical codons (~95% of comparisons).
- Multi-threaded computation via rayon.
- TSV/CSV output with optional SVG plots.
- 21 integration tests.
- Benchmarks: 1,280× faster than KaKs_Calculator, 18,600× faster than BioPython.
- Accuracy: R²=1.0 vs KaKs_Calculator (Li), R²≈0.995-0.999 (Nei).

### Changed

- CLI restructured to subcommands: `eskaks fasta <input>` (original) and `eskaks vcf` (new)
- `eskaks --list-codes` still works at top level
- All existing tests updated for subcommand syntax
- **Architecture refactor**: Split monolithic `main.rs` (391 lines) into 5 focused modules:
  - `cli.rs` — CLI argument definitions (clap derive), cleanly separated from logic
  - `input.rs` — FASTA reading, validation, filtering, deduplication (`SequenceData` struct)
  - `compute.rs` — `ComputeEngine` enum for model-agnostic dispatch (eliminates duplicated closures)
  - `main.rs` — Thin orchestration (~100 lines): parse → load → compute → dispatch
  - `dispatch_output()` and `dispatch_window()` replace inline match blocks
- Total: 11 source modules (was 8), each with a single clear responsibility
- `NeiTables` and `LiTables` constructors accept any genetic code via `with_genetic_code()`.
- Output functions now take `&OutputConfig` instead of 5+ separate parameters.
- Removed unused `_n_u` parameter from `write_group_average`.

### Fixed

- **UX consistency & safety (from a hands-on UX evaluation):**
  - **`eskaks fasta` now errors on unequal-length (unaligned) input** instead of warning
    and returning a number computed over the truncated overlap — pairwise dN/dS requires
    a codon alignment, and a trustable-but-meaningless result was the most dangerous
    failure mode. The error names the offending sequence and points at aligners. (Window
    mode already required equal lengths; every mode is now consistent.)
  - **`eskaks vcf` now honors `--quiet`** (it previously printed its whole pN/pS summary
    regardless), and **accepts `--summary`** (which used to be a hard "unexpected
    argument" error) — `--summary` forces the summary even under `--quiet`, symmetric
    with `eskaks fasta`.
  - More guardrails and feedback (from the same evaluation): the top-level `--help` now
    ends with usage **examples**; a bare `eskaks alignment.fasta` (subcommand forgotten)
    now **suggests the right form** instead of a bare "unrecognized subcommand"; input
    sequences with **non-standard characters** (not A/C/G/T/U/N or a gap) now **warn**;
    `--window-step` without `--window-size` **warns** it did nothing; a **too-small
    `--bootstrap`** (< 100) **warns** the CI is unreliable; `eskaks vcf` shows a **per-gene
    progress bar** on genome-scale runs; and its end-of-run **output list matches the
    `fasta` "Done" block** (one file per line).
  - Small polish: `--list-codes` now says how to apply a code; `--help` links the docs;
    `--min-af -0.5` gives the tool's clear range error instead of a cryptic clap one; a
    content-error is reported as "could not read as FASTA" rather than "failed to open";
    **duplicate sequence ids** warn; passing a **directory** as input errors clearly; the
    help documents `--workers 0` (= all cores) and that `--seed` only applies with
    `--bootstrap`; and `vcf` help uses `<output>` consistently (was `<prefix>`).
- **Follow-up sweep (adversarial round 4, over the code added earlier in this cycle):**
  - **The HTML report's dN/dS histogram mis-binned identical (zero-divergence) pairs.**
    After the report started aggregating over all id-pairs, `collect_report_pairwise`
    still used a ratio convention where a `dN=0, dS=0` pair became `NaN`, landing every
    duplicate pair in the `[1.0, inf)` (positive-selection) bin instead of `[0.0, 0.2)` —
    so on clonal data the report's distribution panel contradicted `--summary` and the
    table. It now uses the writer's exact ratio rule (`0/0 → 0.0`), so the panels agree.
  - The per-pair **neutrality and bootstrap** tables no longer emit `-0.000000` (they
    now normalise negative zero like the main results file, in TSV/CSV and JSON).
  - **`delim_field` no longer over-quotes TSV fields** that contain a `"` but no tab
    (which corrupted the value for TSV readers); quote-escaping is now CSV-only, while
    tab/newline still force quoting in both formats.
- **Correctness & output-integrity sweep (adversarial round 3, found by running the
  binary on crafted inputs):**
  - **Unequal-length sequences silently produced a wrong dN/dS.** The Nei and Li
    fast paths split each sequence into chunks of 4 codons and handled the remainder
    separately, which desynced when the two inputs differed in length: extra codons of
    the longer sequence were dropped and the remainders were compared out of register,
    so real divergence was hidden (reported as 0) or fabricated. Both models now clamp
    to the common length and compare strictly position-by-position over the overlap.
  - **Reference site-counting ignored RNA `U`.** A `U` in the reference never matched
    the `ACGT` self-skip in `count_sites`, so a `U→T` non-mutation was scored as a
    spurious synonymous change (and `--kappa` weighting collapsed), giving a U reference
    different N/S sites from its DNA equivalent. The CDS is now normalised `U→T`.
  - **A gene id reused on two contigs merged into one gene** (its CDS assembled from a
    single contig, the other gene silently lost). GFF3 grouping now keys on
    `(seqid, gene_id)`.
  - **A duplicate contig id in the reference FASTA silently overwrote** the earlier
    record (last wins), scoring genes against the wrong sequence; it is now a hard error.
  - **Sequence ids / gene names containing the CSV/TSV delimiter, a quote, or a newline
    corrupted the columns** (an extra field or a split row) across every delimited
    writer (FASTA pairwise/lineage/window/tests/bootstrap and VCF pN/pS + MK). Fields
    are now RFC-4180 quoted. **Gene/chrom names with a quote/backslash/control char
    broke the VCF `_pnps.json` / `_mk.json`**; they are now JSON-escaped. Escaping and
    quoting are centralised in one module so no writer drifts out of sync again.
  - **The FASTA loader used the whole header line as the id** (leaking descriptions into
    the `Seq1/Seq2` columns and splitting group keys); it now takes the first
    whitespace token, matching the reference parser and standard FASTA semantics.
  - **`--mk` combined with `--min-af`/`--max-af` silently gutted the McDonald-Kreitman
    table** (the AF filter runs before the fixed/polymorphic split); it now warns.
  - Delimited writers no longer emit a nonsensical `-0.000000` (normalised to `0`,
    matching the JSON formatter).
  - **`--format json` was silently ignored in `--lineage`, `--group-average`, and
    `--window` modes** (they wrote CSV into a `.json` file). All three now emit a real
    JSON array of objects (NaN/Infinity → `null`), consistent with pairwise mode.
  - **The HTML report's headline cards (Pairs, Valid pairs, Pooled, Mean) and the
    dN/dS distribution were computed over de-duplicated unique pairs**, so for a clonal
    dataset (many identical sequences) they contradicted the terminal `--summary` and
    the `_pairwise_results` table. The report now aggregates over all id-pairs
    (weighted by sequence multiplicity), matching them.
- **All-N / all-gap pairs reported dN/dS = 0.0 instead of undefined.** A sequence
  with no comparable codons (entirely `N`, ambiguous, or gap) was rendered as a
  perfect `0.0` self-comparison, which reads as extreme purifying selection rather
  than "no data". The pairwise diagonal now runs through the compute engine, and
  both `compute_pair` and `compute_pair_stats` return `NaN` for an all-invalid
  self-comparison while an identical pair with at least one valid codon still
  reports `0.0`. Sliding-window and group-average paths inherit the fix (all-N
  windows are `NaN`; an all-N group contributes zero comparisons, not a spurious 0).
- **No-data / degenerate inputs → `NaN` sweep (siblings of the all-N pair bug,
  found by an adversarial hunt that ran the binary on crafted inputs):**
  - **Sliding-window mode** hard-coded `dN/dS = 0.0` for identical (de-duplicated)
    all-N/all-gap pairs, diverging from the fixed pairwise path. Each window is now
    `NaN` when it has no comparable codons, and the saturation count includes them.
  - **`--lineage` summary** hard-coded the self-diagonal to `{0,0}` and reused it for
    de-duplicated identical genomes, so an all-N genome contributed a spurious `0.0`
    to its lineage mean. It now runs through the compute engine (`NaN` → excluded),
    matching the pairwise and `--group-average` paths.
  - **dN/dS ratio in the pairwise and window writers** returned `+inf` when `dN` was
    `NaN` (nonsynonymous saturation) over a finite zero `dS`. An undefined numerator
    now yields `NaN`, not a spurious extreme-positive-selection reading.
  - **Nei-Gojobori `dS`** collapsed to `0.0` (not `NaN`) when the pair had zero
    synonymous *sites* (all Met/Trp codons); it now matches the Li model's `NaN`.
  - **Per-gene `pN`/`pS`** collapsed a zero-site denominator to `0.0`; they are now
    `NaN` (undefined density), and a fully untranslatable all-ambiguous/all-N CDS is
    skipped with a warning instead of emitted as an all-zero row.
  - **Genome-wide pooled `pN`/`pS`** (and each bootstrap replicate) reported `0.0`
    when the pooled site denominator was `0`; now `NaN`.
  - **`--exclude-repetitive` summary** printed `Total synonymous/nonsynonymous: 0.00`
    when every gene was filtered out; it now reads `n/a`, consistent with the pooled
    ratio line.
  - **Genomic-control λ** rendered `λ = 0.00` when the median tested-gene χ² was `0`
    (a degenerate, mostly-`p=1` family); it now reports `NaN`/`NA`. The correction is
    unchanged (λ is floored at 1 either way).
- **Adversarial bug-hunt round, second batch:**
  - **`--min-snps` filtered on the AF-weighted fractional total under `--af-weighted`**,
    so it dropped genes that had far more than the threshold in real SNPs. It now
    filters on the raw SNP count (new unweighted `n_snps`).
  - **`--min-af`/`--max-af` filtered whole records**, so at a multi-allelic site a
    sub-threshold allele leaked through and a co-located in-range allele was dropped.
    Filtering is now per-allele; the record is dropped only if no ALT survives.
  - **Bootstrap genome-wide CI** discarded replicates with no synonymous variation
    (`pS = 0`, `pN > 0`) as undefined, truncating the upper tail and biasing the CI
    low; such replicates are now kept as `+∞` toward the upper percentile, and the
    "excluded replicates" warning fires even when *every* replicate is undefined.
  - **A CDS exon extending past the reference end** desynced the SNP→codon offset
    mapping and silently dropped SNPs; such genes are now skipped with a warning.
  - Documented that the `--exclude-repetitive` pooled totals are core-only in the
    summary; corrected the `count_sites` docstring (renormalise-to-3 convention) and
    the genomic-control λ docstring (a discrete-test null yields λ < 1, not ≈ 1).
- **Adversarial bug-hunt round (found by running the binary on crafted inputs):**
  - **Li/LPB93 dS collapsed to 0 (not NaN) at synonymous transition saturation**,
    reporting a spurious dN/dS = ∞ (apparent extreme positive selection) where Nei
    correctly reports strong purifying selection. The undefined Kimura correction now
    propagates NaN, as the module already documented.
  - **A single malformed `INFO/AF` token** (`nan`, `inf`, or a value outside `[0,1]`)
    was parsed into `alt_freqs`, where a `NaN` slipped past `--min-af`/`--max-af`
    (every NaN comparison is false) and poisoned the genome-wide πN/πS to `NaN`.
    Non-finite / out-of-range AF tokens are now rejected (treated as missing).
  - **Out-of-range GT allele indices** (a genotype referencing an undeclared ALT)
    still incremented the allele total, deflating every real allele frequency at the
    site; they are now ignored.
  - **JSON output did not escape sequence IDs**, so an id containing `"` produced
    invalid JSON. IDs are now JSON-escaped in all pairwise writers.
  - **CSV group-average `95%CI` field** embedded an unquoted comma, splitting the
    column and misaligning the row; it is now quoted in CSV mode.
  - **Merged read-depth summed in `u32`** could overflow (panic in debug) on very
    high `DP` across many samples; it now accumulates in `u64`.
  - **Unbounded GFF3 CDS coordinates** could overflow `length_bp` / over-allocate the
    CDS buffer (panic/OOM); an implausibly wide CDS span is now skipped with a warning.
- **Non-deterministic multi-VCF merge (found while hardening test coverage):**
  when merging single-sample VCFs, a position carrying more than one ALT allele
  (multi-allelic sites, or samples disagreeing on the variant base) emitted its
  ALT alleles — and their aligned frequencies — in `HashMap` iteration order,
  which Rust randomizes per run. The merged output therefore varied run-to-run at
  such positions. ALTs are now sorted by base within each position, restoring
  reproducible output.
- **Consistency & robustness (third audit):**
  - A stray internal space/tab in a FASTA sequence line used to occupy a codon slot
    and frameshift every codon after it; internal whitespace is now stripped (gaps
    `-`/`.` are kept) before framing.
  - Under `--exclude-repetitive`, the genome-wide point estimate is core-only but the
    bootstrap 95% CI resampled **all** genes, so the CI need not bound the estimate;
    the CI now resamples the same core-only set.
  - The report's "Genes analyzed"/"With SNPs" cards were computed from the
    post-`--min-snps` slice while the pooled estimate used all genes, so the report
    disagreed with the CLI summary; the report now uses the pre-filter counts.
  - "Genes tested" now counts the actual multiple-testing family (finite q), and the
    CLI "Significant" count uses the GC-corrected q under `--genomic-control`, so the
    CLI and the HTML report agree.
  - GFF3 attribute URL-decoding reassembles multi-byte UTF-8 escapes (`%C3%A9` → `é`)
    instead of emitting Latin-1 bytes; the GT index is read per-record (VCF allows the
    FORMAT order to vary); `merge_vcfs` drops an ALT equal to the merged REF; and
    `--divergence` warns when `--report` is not given.
- **Robustness / escaping (second audit):**
  - A GFF3 CDS with `end < start` underflowed `end - start` (panic in debug, a
    capacity-overflow panic in release, aborting the whole run); such lines are now
    skipped with a warning.
  - SVG plot labels (genome/lineage/group names) and the FASTA report's genome/
    lineage/group names are now XML/HTML-escaped, so a name containing `&`, `<`, `>`
    or `"` no longer produces a blank `.svg` or injects markup into the report.
  - In the pN/pS report, a gene with `pN/pS = ∞` (synonymous count 0, nonsynonymous
    > 0) — the strongest diversifying signal — was serialized as a null ratio and
    then mislabelled **purifying** (blue, ▼) in the census, Manhattan, hits and
    table; it is now correctly classified as diversifying everywhere.
- **Correctness (dN/dS & pN/pS):**
  - Genetic codes **24** and **33** mistranslated `AGA` as Lys; it is **Ser**
    (`AGG` stays Lys), so every codon comparison involving `AGA` and the resulting
    dN/dS was wrong under those two tables.
  - **Multi-exon minus-strand** genes reconstructed their CDS with the exons in
    swapped order, so every SNP was mapped to the wrong codon and most were dropped
    as "REF mismatch"; pN/pS, MK counts, and site totals were all wrong for such
    genes. (Single-exon and plus-strand genes were unaffected.)
  - VCF `INFO/AF`: a missing/malformed token (e.g. `AF=.,0.2,0.3`) was silently
    dropped, shifting every later allele frequency onto the **wrong ALT**; parsing
    now preserves positional alignment.
  - The Li/LPB93 Arg synonymous-transversion special case is now gated on the
    active genetic code (it wrongly applied to non-synonymous `AGA`/`AGG` under
    codes where those are not Arg), and the Kimura `B` term is recovered when only
    transitions saturate.
  - `merge_vcfs` counted duplicate records within one sample more than once when
    computing "fraction of samples"; it now dedupes per sample.
- **Determinism:** the pairwise, sliding-window, lineage, and group-average TSV/JSON
  outputs (and the lineage/group plot data) were written in thread-completion order;
  they are now emitted in a stable, reproducible order regardless of thread count.
- **Robustness:** output files are created in the calling thread (a bad output path
  now returns a clean error instead of panicking a worker) and buffered writers are
  flushed explicitly; `--min-codons` that leaves fewer than two sequences now errors
  clearly instead of proceeding.
- **Report:** gene/chromosome names are HTML/JSON-escaped, so a name containing
  `</script>`, `<`, `&`, or `"` can no longer break the page or inject markup;
  scatter axis bounds are computed without a spread that could overflow the stack on
  whole-genome runs.
- Data-quality warnings (REF mismatches, skipped genes, saturation) are shown by
  default; added `-v`/`-vv`/`-q`. Validate AF-filter ranges and reconcile
  contig names up front (both previously produced a silent all-NaN output).
  Empty VCFs in a merge no longer abort the run. Aggregate REF-mismatch
  diagnostics instead of one line per SNP; warn on non-multiple-of-3 CDS and on
  genes dropped from the plot. Report the pooled `mean(dN)/mean(dS)` in the
  dN/dS summary instead of only the biased mean of per-pair ratios.
- 19 clippy warnings resolved (0 remaining).
- Documented that Nei pathway fallback values (2-diff: 0.5/1.5, 3-diff: 1.0/2.0)
  are unreachable with standard genetic code but triggered by alternative codes.

### Performance

- `eskaks vcf` scales to large genomes and cohorts: the per-gene SNP scan is
  now a binary-searched genomic window (O(S + G·log S) instead of O(G·S)), and
  `merge_vcfs` parses per-sample VCFs in parallel — both byte-identical output.
