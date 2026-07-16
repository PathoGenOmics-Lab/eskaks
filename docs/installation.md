# Installation

## Conda (bioconda)

```bash
conda install -c bioconda eskaks
```

Or with mamba:

```bash
mamba install -c bioconda eskaks
```

## From source

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks
make release
```

The binary will be at `target/release/eskaks`. Copy it to your PATH:

```bash
cp target/release/eskaks ~/.local/bin/
```

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70.0

## Performance tip

The `make release` target automatically enables native CPU optimizations (`-C target-cpu=native`). This uses CPU-specific SIMD instructions for maximum performance on your hardware.

## Verify installation

```bash
eskaks --version
eskaks --list-codes
```
