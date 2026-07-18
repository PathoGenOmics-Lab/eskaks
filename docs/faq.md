---
description: >-
  Common eskaks questions and errors — getting started, codon alignment, genetic
  codes, NaN/inf values, VCF/GFF contig mismatches and interpreting the ratios.
---

# FAQ

## I'm new: where do I start?

Run the [getting-started tutorial](tutorial.md). It uses example data that ships in
the `examples/` folder, so you can get a full result in a couple of minutes without
any data of your own. Any unfamiliar term is defined in the [glossary](glossary.md).

## Which mode do I use: `fasta` or `vcf`?

- Comparing **sequences** (one gene from several species/strains, already aligned)?
  → `eskaks fasta` gives you **dN/dS**.
- Comparing **individuals in a population** (you have a reference genome and a VCF of
  variants)? → `eskaks vcf` gives you **pN/pS per gene**.

## What input files do I need?

- **`eskaks fasta`**: one **codon-aligned** FASTA (sequences in frame, all the same
  length). If yours aren't aligned yet, run [MAFFT](https://mafft.cbrc.jp/) +
  [PAL2NAL](http://www.bork.embl.de/pal2nal/) or [MACSE](https://bioweb.supagro.inra.fr/macse/) first.
- **`eskaks vcf`**: three files, a **reference FASTA**, a **GFF3** annotation, and one
  or more **VCF** files. The **contig/chromosome names must match** across all three.

## My `eskaks vcf` output is empty or all `NA`

Almost always a **name mismatch**: the chromosome name in the VCF (e.g. `chr1`) must
be *identical* to the sequence name in the reference FASTA (`>chr1`) and the first
column of the GFF3. eskaks warns when it can't reconcile them, run with `-v` to see
the details. Also check that your GFF3 has `CDS` features and that `--genetic-code`
matches your organism (`11` for bacteria).

## Why is eskaks so much faster?

Three main reasons:

1. **Precomputed lookup tables**: All 4,096 codon pairs have their diffs/sites precomputed at startup. Each pairwise comparison is just table lookups + a final correction formula.

2. **Cache optimization**: The Li model's lookup table (288 KB) fits in L2 cache. Identical codons (~95% in typical alignments) use a separate 1.5 KB table that fits in L1 cache.

3. **Parallel + streaming**: Rayon distributes pairs across CPU cores, and a dedicated writer thread handles I/O without blocking computation.

## Why do I get NaN values?

NaN means the substitution proportion has reached **saturation** (p ≥ 0.749 for Nei-Gojobori, or denominator ≤ 0 for Li). This happens when sequences are too divergent, the correction formula breaks down because there have been so many substitutions that the signal is lost.

**What to do**: This is biologically correct. Very divergent sequences simply can't be reliably compared at the nucleotide level. Consider using protein-level methods instead.

## Which model should I use?

- **Nei-Gojobori**: Simpler, faster, good for quick scans. Use when you want results comparable to most published analyses.
- **Li (1993)**: More accurate because it accounts for the well-known transition/transversion bias. Use when accuracy matters more than speed (though the speed difference is minimal).

## My sequences have internal stop codons

eskaks warns about these automatically. Common causes:
- **Wrong reading frame**: Shift your sequences by 1 or 2 bases
- **Pseudogene**: The gene has been inactivated
- **Frameshift mutation**: An indel has disrupted the reading frame
- **Wrong genetic code**: Try `--genetic-code 2` (mitochondrial) etc.

## Can I use eskaks with non-standard genetic codes?

Yes! Use `--genetic-code <N>` with any of the 20 supported NCBI tables. Run `eskaks --list-codes` to see all options.

## How do I enable tab-completion?

Print a completion script for your shell and install it, e.g. for Bash:

```bash
eskaks --completions bash | sudo tee /etc/bash_completion.d/eskaks > /dev/null
```

`zsh`, `fish`, `elvish`, and `powershell` are also supported — see [Installation](./installation.md#shell-completions).

## How do I control what eskaks prints?

Each run ends with a short "Done" confirmation (counts, model, output files) on stderr. Use `--quiet` to show only errors, `-v`/`-vv` for progress and debug detail, or `--summary` for the full statistics block. `RUST_LOG` overrides the log level entirely.

## How do I cite eskaks?

See the [CITATION.cff](https://github.com/PathoGenOmics-Lab/eskaks/blob/main/CITATION.cff) file, or use GitHub's "Cite this repository" button.

## What is the difference between dN/dS and pN/pS?

**dN/dS** measures fixed substitutions between diverged sequences (from aligned FASTA). It applies corrections for multiple substitutions (Jukes-Cantor or Kimura).

**pN/pS** measures polymorphism within a population (from VCF). It counts raw variant proportions without multiple-hit correction. Use `eskaks fasta` for dN/dS and `eskaks vcf` for pN/pS.

## Can I use eskaks as a library?

Yes. eskaks exposes a `lib.rs` with public modules. Add it as a dependency in your `Cargo.toml`:

```toml
[dependencies]
eskaks = { git = "https://github.com/PathoGenOmics-Lab/eskaks.git" }
```

Then use `eskaks::models::nei::NeiTables` or `eskaks::models::li::LiTables` directly.
