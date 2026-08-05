---
description: >-
  The complete eskaks flag reference — global options, `eskaks fasta` (pairwise
  dN/dS) and `eskaks vcf` (per-gene pN/pS), exit codes and environment variables.
---

# CLI Reference

eskaks uses subcommands for its two modes of operation:

```
eskaks <COMMAND> [OPTIONS]
eskaks --list-codes
eskaks --version
```

## Global flags

| Flag | Description |
|---|---|
| `--demo` | Run a quick demo on the bundled example data (no input files needed) and exit |
| `--list-codes` | List available NCBI genetic code tables and exit |
| `--completions <SHELL>` | Print a shell completion script and exit (`bash`, `zsh`, `fish`, `elvish`, `powershell`) |
| `-v, --verbose` | Increase log verbosity: `-v` shows info, `-vv` shows debug. Data-quality warnings show by default |
| `-q, --quiet` | Silence all logs except errors (also hides the run-confirmation block) |
| `-h, --help` | Show help. Per subcommand, flags are grouped into sections (Output / Analysis / Statistics / Filtering / Input) and followed by usage examples |
| `-V, --version` | Show version |

### Run confirmation & logging

Every `eskaks fasta` / `eskaks vcf` run ends with a short confirmation block on
stderr (sequence/gene counts, model, and the list of files written), so a run is
never silent about where its output went. Add `--quiet` to suppress it.

Diagnostics are printed as clean `warning:` / `error:` lines (coloured only when
stderr is a terminal). Set `RUST_LOG` to override the level entirely.

---

## `eskaks fasta`: Pairwise dN/dS

Compute pairwise dN/dS from codon-aligned FASTA sequences.

```
eskaks fasta <INPUT_FILE> [OPTIONS]
```

| Argument | Description |
|---|---|
| `<INPUT_FILE>` | Aligned coding sequences in FASTA format. Use `-` for stdin. |

| Flag | Description | Default |
|---|---|---|
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `-w, --workers <N>` | Number of parallel threads | `4` |
| `--model <nei\|li>` | Substitution model | `nei` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table number | `1` |
| `--lineage` | Lineage summary mode | off |
| `--group-average` | Group average mode | off |
| `--first-letter-lineage` | Group by first character of ID | off |
| `--window-size <N>` | Sliding window size (codons) | off |
| `--window-step <N>` | Window step size | `1` |
| `--min-codons <N>` | Filter sequences with < N valid codons | `0` |
| `--summary` | Print summary statistics to stderr | off |
| `--neutrality` | Write a per-pair Nei-Gojobori neutrality test (`<output>_pairwise_tests`): dN, dS, SEs, Z, p-value (Nei only) | off |
| `--bootstrap <N>` | Per-pair 95% bootstrap CIs on dN, dS, dN/dS (`<output>_pairwise_bootstrap`); resamples codons, works for both models | `0` |
| `--seed <N>` | Seed for reproducible bootstrap resampling | `42` |
| `--plot` | Generate SVG plot(s) | off |
| `--report` | Write an interactive HTML report (`<output>_report.html`): lineage/group scatter with per-group means, or the pairwise dN/dS distribution | off |

### Examples

```bash
# Basic Nei model
eskaks fasta input.fasta

# Li model, 16 threads, CSV output
eskaks fasta input.fasta --model li --workers 16 --format csv -o results

# Vertebrate mitochondrial code, JSON output
eskaks fasta mito.fasta --genetic-code 2 --format json -o mito

# Sliding window with plot
eskaks fasta input.fasta --window-size 100 --window-step 10 --plot -o windows

# From stdin
cat input.fasta | eskaks fasta - --summary -o piped
```

---

## `eskaks vcf`: pN/pS per gene

Compute pN/pS per gene from a VCF file, reference FASTA, and GFF3 annotation.

```
eskaks vcf --ref <FASTA> --gff <GFF3> --vcf <VCF> [OPTIONS]
```

**Inputs & output**

| Flag | Description | Default |
|---|---|---|
| `--ref <FASTA>` | Reference genome in FASTA format | required |
| `--gff <GFF3>` | Gene annotation in GFF3 format | required |
| `--vcf <VCF>` | VCF file(s), use multiple times for per-sample VCFs | one of `--vcf` / `--vcf-list` required |
| `--vcf-list <FILE>` | File with one VCF path per line (alternative to repeated `--vcf`) | one of `--vcf` / `--vcf-list` required |
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table number | `1` |
| `--workers <N>` | Parallel threads (output is deterministic) | `4` |

**Variant filters**

| Flag | Description | Default |
|---|---|---|
| `--pass-only` | Only include FILTER=PASS variants | off |
| `--min-af <FLOAT>` | Minimum allele frequency (0.0–1.0) | none |
| `--max-af <FLOAT>` | Maximum allele frequency (exclude fixed variants) | none |
| `--min-depth <INT>` | Minimum read depth (INFO/DP) | none |
| `--af-weighted` | Weight counts by AF (πN/πS instead of pN/pS) | off |

**Model, test & correction**

| Flag | Description | Default |
|---|---|---|
| `--kappa <FLOAT>` | ts/tv rate ratio for spectrum-aware site counting | `1.0` |
| `--min-snps <N>` | Drop genes with fewer SNPs from the table, plot, and test | `0` |
| `--fdr <FLOAT>` | Benjamini-Hochberg threshold for significant genes | `0.05` |
| `--mk` | Run the McDonald-Kreitman test (`<prefix>_mk.<ext>`) | off |
| `--mk-fixed-af <FLOAT>` | AF at/above which a variant is "fixed" in the MK test | `0.99` |
| `--bootstrap <N>` | Replicates for a 95% CI on the genome-wide pooled pN/pS | `0` |
| `--seed <N>` | Seed for reproducible bootstrap resampling | `42` |
| `--genomic-control` | Divide each χ² by the inflation factor λ and re-test | off |
| `--exclude-repetitive` | Drop PE/PPE/PGRS/IS genes from the pooled estimate and test | off |

**Output & reporting**

| Flag | Description | Default |
|---|---|---|
| `--variants` | Write a per-coding-SNP table (`<prefix>_variants.<ext>`): position, base and amino-acid change (e.g. `S315T`), AF, and effect (synonymous/missense/nonsense/stop_loss) | off |
| `--diversity` | Write per-gene πN/πS, Watterson θ and Tajima's D (`<prefix>_diversity.<ext>`); needs the sample size, so a multi-sample VCF or `--vcf-list` | off |
| `--summary` | Print the pN/pS summary block to stderr | off |
| `--plot` | Generate Manhattan / p-value SVG plots | off |
| `--report` | Write a self-contained interactive HTML report | off |
| `--divergence <FILE>` | Per-gene dN/dS TSV for the report's polymorphism-vs-divergence panel | none |

### Examples

```bash
# Population πN/πS from per-sample VCFs
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf sample1.vcf --vcf sample2.vcf --vcf sample3.vcf \
  --af-weighted --genetic-code 11 -o mtb_pnps

# Full selection scan with the interactive report
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf-list samples.txt \
  --genetic-code 11 --kappa 2 --min-snps 5 \
  --mk --bootstrap 1000 --seed 42 --report --plot -o mtb_scan

# Single multi-sample VCF with filters
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --pass-only --min-af 0.05 --min-depth 10 -o filtered
```

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Error (invalid input, missing file, bad arguments) |
| 2 | Argument parsing error (clap) |

## Environment variables

| Variable | Description |
|---|---|
| `RUST_LOG` | Set log level: `info`, `warn`, `debug` |
