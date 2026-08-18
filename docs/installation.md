---
description: >-
  Install eskaks by building from source with a Rust toolchain, verify the binary,
  and set up shell tab-completion.
---

# Installation

## From a release

Every tagged release attaches a binary for macOS (Apple silicon and Intel), Linux
x86_64 and Windows x86_64. Take the one for your platform from the
[latest release](https://github.com/PathoGenOmics-Lab/eskaks/releases/latest),
unpack it, and put `eskaks` on your `PATH`.

Each release also carries a `SHA256SUMS` file and a signed build provenance
attestation, so you can check that the file you downloaded is the one the release
workflow built, rather than trusting that the download went where you meant it to:

```bash
shasum -a 256 -c SHA256SUMS
gh attestation verify eskaks-<version>-<platform>.tar.gz --owner PathoGenOmics-Lab
```

## From source

eskaks also builds from source with a
[Rust](https://www.rust-lang.org/tools/install) toolchain (Rust ≥ 1.85):

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks
make release          # or: cargo build --release
```

The binary is written to `target/release/eskaks`. Put it on your `PATH`:

```bash
cp target/release/eskaks ~/.local/bin/
```

!!! info "Package managers are not there yet"
    `cargo install eskaks` and a bioconda recipe (`conda install -c bioconda eskaks`)
    are not available. Use a release binary or build from source.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) ≥ 1.85.0, for the source build.
- A C compiler, also for the source build. eskaks itself is Rust, but it reads
  compressed input through `needletail`, which pulls `liblzma-sys` and `zstd-sys`,
  and those build vendored C. A release binary needs neither: it links only the
  platform C library.

## Performance tip

The `make release` target enables native-CPU optimizations (`-C target-cpu=native`):
CPU-specific SIMD instructions for maximum speed on *your* hardware. For a portable
binary you can run on a different CPU, build with `cargo build --release` instead.

## Verify installation

```bash
eskaks --version
eskaks --list-codes
```

## Shell completions

eskaks can print a completion script for your shell, so `<Tab>` completes
subcommands and flags:

```bash
# Bash
eskaks --completions bash | sudo tee /etc/bash_completion.d/eskaks > /dev/null

# Zsh (into a directory on your $fpath, e.g. ~/.zfunc)
eskaks --completions zsh > ~/.zfunc/_eskaks

# Fish
eskaks --completions fish > ~/.config/fish/completions/eskaks.fish
```

`bash`, `zsh`, `fish`, `elvish`, and `powershell` are supported. Restart the
shell (or re-source your profile) afterwards.
