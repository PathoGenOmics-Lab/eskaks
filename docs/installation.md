# Installation

## From source

eskaks currently installs by building from source with a
[Rust](https://www.rust-lang.org/tools/install) toolchain (Rust ≥ 1.70):

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks
make release          # or: cargo build --release
```

The binary is written to `target/release/eskaks`. Put it on your `PATH`:

```bash
cp target/release/eskaks ~/.local/bin/
```

!!! info "Packages are coming"
    `cargo install eskaks`, pre-built binaries, and a bioconda recipe
    (`conda install -c bioconda eskaks`) will ship with the first tagged release.
    Until then, build from source as above.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70.0 — for the source build.

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
