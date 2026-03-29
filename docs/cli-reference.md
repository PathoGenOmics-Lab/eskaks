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
| `--plot` | Generate SVG plot(s) | off |

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
| `--vcf <VCF>` | Variants in VCF format | required |
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `--format <tsv\|csv\|json>` | Output format | `tsv` |
| `--genetic-code <N>` | NCBI translation table number | `1` |
| `--pass-only` | Only include FILTER=PASS variants | off |
| `--min-af <FLOAT>` | Minimum allele frequency (0.0–1.0) | none |
| `--min-depth <INT>` | Minimum read depth (INFO/DP) | none |
| `--plot` | Generate Manhattan-style SVG plot | off |

### Examples

```bash
# Basic pN/pS for M. tuberculosis
eskaks vcf --ref H37Rv.fasta --gff H37Rv.gff3 --vcf population.vcf \
  --genetic-code 11 -o mtb_pnps

# Filter: only PASS, AF ≥ 0.05, depth ≥ 10
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --pass-only --min-af 0.05 --min-depth 10 -o filtered

# JSON output with plot
eskaks vcf --ref ref.fasta --gff ref.gff3 --vcf calls.vcf \
  --format json --plot -o results
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
