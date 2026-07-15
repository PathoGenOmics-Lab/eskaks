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
| `--list-codes` | List available NCBI genetic code tables and exit |
| `-h, --help` | Show help |
| `-V, --version` | Show version |

---

## `eskaks fasta` — Pairwise dN/dS

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

## `eskaks vcf` — pN/pS per gene

Compute pN/pS per gene from a VCF file, reference FASTA, and GFF3 annotation.

```
eskaks vcf --ref <FASTA> --gff <GFF3> --vcf <VCF> [OPTIONS]
```

| Flag | Description | Default |
|---|---|---|
| `--ref <FASTA>` | Reference genome in FASTA format | required |
| `--gff <GFF3>` | Gene annotation in GFF3 format | required |
| `--vcf <VCF>` | VCF file(s) — use multiple times for per-sample VCFs | required |
| `--vcf-list <FILE>` | File with one VCF path per line | none |
| `--af-weighted` | Weight counts by AF (πN/πS instead of pN/pS) | off |
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table number | `1` |
| `--pass-only` | Only include FILTER=PASS variants | off |
| `--min-af <FLOAT>` | Minimum allele frequency (0.0–1.0) | none |
| `--max-af <FLOAT>` | Maximum allele frequency (exclude fixed variants) | none |
| `--min-depth <INT>` | Minimum read depth (INFO/DP) | none |
| `--plot` | Generate Manhattan-style SVG plot | off |

### Examples

```bash
# Population πN/πS from per-sample VCFs
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf sample1.vcf --vcf sample2.vcf --vcf sample3.vcf \
  --af-weighted --genetic-code 11 -o mtb_pnps

# Or via a list file
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 \
  --vcf-list samples.txt --af-weighted \
  --min-af 0.01 --max-af 0.99 --plot -o mtb_pnps

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
