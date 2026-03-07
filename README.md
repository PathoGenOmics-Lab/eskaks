
<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/eskaks">
    <img src="https://github.com/PathoGenOmics-Lab/eskaks/blob/main/.github/logos/esKaKs.png" height="350" alt="eskaks">
  </a>
</p>
</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

**eskaks** is a high-performance command-line tool written in Rust for calculating dN/dS (Ka/Ks) ratios between pairs of aligned coding sequences. It supports both the Nei-Gojobori (1986) and Li (1993) models for estimating synonymous (dS) and nonsynonymous (dN) substitution rates. The tool leverages parallel processing via Rayon, sequence deduplication, and precomputed lookup tables for efficient computation on large datasets.

## Features

- **Models**: Choose between the Nei-Gojobori (default) or Li (1993) model for dN/dS calculation.
- **Parallel Processing**: Configurable multi-threaded computation via `--workers` using Rayon.
- **Sequence Deduplication**: Automatically detects and groups identical sequences to avoid redundant pairwise computations.
- **Precomputed Lookup Tables**: The Li model uses flat-array lookup tables (64x64 codon pairs) for cache-friendly, zero-allocation pair evaluation.
- **Streaming Output**: Pairwise results are written via a dedicated writer thread with buffered I/O to minimize memory usage.
- **Group Summaries**:
  - **Group Average**: Compute mean dN/dS between predefined groups of sequences with 95% confidence intervals.
  - **Lineage Summary**: Compute mean dN/dS for each sequence against sequences grouped by lineage.
  - **Configurable Lineage Grouping**: Group by splitting on `_` (default) or by the first character of each sequence name (`--first_letter_lineage`).

## Models

### Nei-Gojobori (1986)

The Nei-Gojobori method is a straightforward counting approach:

1. Counts synonymous and nonsynonymous sites for each codon using a precomputed site table.
2. Classifies codon differences as synonymous (same amino acid) or nonsynonymous (different amino acid).
3. Applies the **Jukes-Cantor correction** to account for multiple substitutions at the same site.
4. Returns `NaN` when substitution proportions reach saturation (p >= 0.749).

**Best for**: Fast exploratory analyses, large datasets, and when simplicity is preferred.

### Li (1993)

The Li method is more sophisticated:

1. Classifies sites into three categories (0-fold, 2-fold, 4-fold degenerate) based on the genetic code.
2. Considers all possible evolutionary pathways for codons differing at 2 or 3 positions, weighted by amino acid similarity.
3. Separately estimates transition and transversion rates for each site category.
4. Applies the **Kimura two-parameter correction** to each category.
5. Uses precomputed 64x64 flat-array lookup tables for zero-allocation pair computation.

**Best for**: More accurate estimates that account for codon usage bias, transition/transversion differences, and site degeneracy.

## Input

A FASTA file containing multiple aligned coding sequences. Requirements:

- All sequences must be the same length.
- Sequence length must be a multiple of 3 (complete codons only).
- Standard DNA/RNA alphabet (A, C, G, T/U); ambiguous bases (N, etc.) are handled gracefully.

## Output

- **Pairwise results** (`<prefix>_pairwise_results.tsv`): Tab-separated file with columns:
  - `Seq1`, `Seq2`, `dN`, `dS`, `dN/dS` (Nei model)
  - `Seq1`, `Seq2`, `dN(Ka)`, `dS(Ks)`, `dN/dS` (Li model)

- **Group average** (`<prefix>_group_avg_dn_ds.tsv`): Mean dN/dS between groups with columns:
  - `Group1`, `Group2`, `NumSeqs1`, `NumSeqs2`, `NumComparisons`, `Mean_dN/dS`, `StdError`, `95%CI`

- **Lineage summary** (`<prefix>_lineage_summary.tsv`): Mean dN/dS per lineage with columns:
  - `Genome`, `Against_Lineage`, `Mean_dN`, `Mean_dS`, `dN/dS_Ratio`

## Installation

### From source (recommended)

```bash
# Clone the repository
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks

# Build with maximum optimizations
RUSTFLAGS="-C target-cpu=native" cargo build --release

# The binary will be at target/release/eskaks
```

### Requirements

- [Rust](https://www.rust-lang.org/tools/install) >= 1.70.0

### Performance tip

For maximum performance on your specific CPU, compile with native target optimizations:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This enables CPU-specific SIMD instructions and other architecture optimizations.

## Usage

```bash
eskaks <input_file> [options]
```

### Required Arguments

- `<input_file>`: Aligned coding sequences in FASTA format.

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output <PREFIX>` | Base name for output files | `output` |
| `-w, --workers <N>` | Number of parallel threads | `4` |
| `--model <nei\|li>` | Model to use for calculation | `nei` |
| `--lineage` | Compute summary results by lineage | off |
| `--group_average` | Compute average dN/dS between groups | off |
| `--first_letter_lineage` | Group sequences by first letter (requires `--lineage`) | off |

## Examples

1. **Basic Nei model calculation:**
   ```bash
   eskaks input.fasta
   ```
   Produces `output_pairwise_results.tsv` with pairwise dN/dS (Nei model, 4 threads).

2. **Li model with 8 threads:**
   ```bash
   eskaks input.fasta --model li --workers 8 -o results
   ```
   Produces `results_pairwise_results.tsv` using the Li (1993) model with 8 parallel threads.

3. **Lineage summary:**
   ```bash
   eskaks input.fasta --lineage -o analysis
   ```
   Produces `analysis_lineage_summary.tsv` with mean dN/dS for each sequence against each lineage group (lineages determined by splitting sequence IDs on `_`).

4. **Group average with first-letter lineage grouping:**
   ```bash
   eskaks input.fasta --lineage --first_letter_lineage
   ```
   Groups sequences by their first character instead of splitting by `_`.

5. **Group average with confidence intervals:**
   ```bash
   eskaks input.fasta --group_average --workers 16
   ```
   Computes mean dN/dS between all group pairs with standard error and 95% confidence intervals.

## Architecture

```
src/
  main.rs          - CLI, FASTA reading, deduplication, parallel dispatch
  codon.rs         - DNA5 encoding and codon index conversion
  output.rs        - Parallel output writers (pairwise, lineage, group)
  models/
    mod.rs         - Model enum and shared types
    nei.rs         - Nei-Gojobori (1986) with Jukes-Cantor correction
    li.rs          - Li (1993) with precomputed flat-array lookup tables
```

## License

This project is licensed under the [GNU Affero General Public License v3.0](LICENSE).

---
<h2 id="contributors" align="center">

Contributors
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
eskaks is developed by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---
