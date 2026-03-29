# CLI Reference

```
eskaks [OPTIONS] <INPUT_FILE>
```

## Arguments

| Argument | Description |
|----------|-------------|
| `<INPUT_FILE>` | Aligned coding sequences in FASTA format. Use `-` for stdin. |

## Options

| Flag | Description | Default |
|------|-------------|---------|
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
| `--list-codes` | List available genetic codes and exit | off |
| `-h, --help` | Show help | |
| `-V, --version` | Show version | |

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, missing file, bad arguments) |
| 2 | Argument parsing error (clap) |

## Environment variables

| Variable | Description |
|----------|-------------|
| `RUST_LOG` | Set log level: `info`, `warn`, `debug` |

## Examples

```bash
# Basic Nei model
eskaks input.fasta

# Li model, 16 threads, CSV output
eskaks input.fasta --model li --workers 16 --format csv -o results

# Vertebrate mitochondrial code, JSON output
eskaks mito.fasta --genetic-code 2 --format json -o mito

# Sliding window with plot
eskaks input.fasta --window-size 100 --window-step 10 --plot -o windows

# From stdin
cat input.fasta | eskaks - --summary -o piped
```
