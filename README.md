
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
- **Memory-Efficient Codon Encoding**: FASTA sequences are converted directly to codon indices (L/3 bytes per sequence instead of L), achieving a 3x memory reduction compared to storing raw nucleotides.
- **Sequence Deduplication**: Automatically detects and groups identical sequences to avoid redundant pairwise computations.
- **Precomputed Lookup Tables**: The Li model uses flat-array lookup tables (64×64 codon pairs, AoS layout) for cache-friendly, zero-allocation pair evaluation, with an L1-cache-resident fast path for identical codons (~95% of comparisons in typical alignments).
- **Streaming Output**: Results are streamed via dedicated writer threads with buffered I/O and lazy per-row caching using generation counters, keeping memory usage at O(U) per thread instead of O(U²) total.
- **Group Summaries**:
  - **Group Average**: Compute mean dN/dS between predefined groups of sequences with 95% confidence intervals.
  - **Lineage Summary**: Compute mean dN/dS for each sequence against sequences grouped by lineage.
  - **Configurable Lineage Grouping**: Group by splitting on `_` (default) or by the first character of each sequence name (`--first-letter-lineage`).

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
5. Uses precomputed 64×64 flat-array lookup tables (~288 KB, AoS layout) for zero-allocation pair computation, with a compact 1.5 KB `same_l` table for identical-codon fast path.

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
| `--lineage` | Compute summary results by lineage (mutually exclusive with `--group-average`) | off |
| `--group-average` | Compute average dN/dS between groups (mutually exclusive with `--lineage`) | off |
| `--first-letter-lineage` | Group sequences by first letter (requires `--lineage`) | off |

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
   eskaks input.fasta --lineage --first-letter-lineage
   ```
   Groups sequences by their first character instead of splitting by `_`.

5. **Group average with confidence intervals:**
   ```bash
   eskaks input.fasta --group-average --workers 16
   ```
   Computes mean dN/dS between all group pairs with standard error and 95% confidence intervals.

## Architecture

```
src/
  main.rs          - CLI, FASTA→codon index reading, deduplication, parallel dispatch
  codon.rs         - DNA5 encoding, codon index conversion, group key extraction
  output.rs        - Streaming output writers with for_each_init + generation counters
  models/
    mod.rs         - Model enum and shared types (DsDn, Z_95_CONFIDENCE)
    nei.rs         - Nei-Gojobori (1986) with Jukes-Cantor correction
    li.rs          - Li (1993) with AoS lookup tables + same_l fast path
tests/
  integration.rs   - 11 integration tests against the compiled binary
  data/
    synthetic.fasta         - 5 hand-crafted sequences for pairwise/lineage tests
    synthetic_grouped.fasta - 5 sequences with A_/B_ prefixes for group-average tests
```

### Key Performance Patterns

- **3x memory reduction**: Sequences stored as codon indices (L/3 bytes) instead of raw nucleotides (L bytes)
- **Sort-based deduplication**: O(n log n) on compact codon indices; only unique pairs are computed
- **Chunk-of-4 processing**: Both models process 4 codons at a time with bitwise validity checks for branch-free fast paths
- **Generation counters**: Thread-local per-row caches are invalidated in O(1) via wrapping counters instead of O(n) clearing
- **L1-cache fast path (Li)**: Identical codons (~95% in typical alignments) use a 1.5 KB table instead of the full 288 KB lookup
- **Streaming I/O**: Crossbeam channels with 64 KB batched buffers feed a dedicated writer thread

## Benchmarks

eskaks was benchmarked against established dN/dS tools on synthetic datasets of varying sizes. All benchmarks were run on a single machine; eskaks timings include both single-threaded (1t) and multi-threaded (4t) runs.

### Accuracy

Numerical accuracy was validated by comparing pairwise dN and dS values against [KaKs_Calculator](https://github.com/kullrich/kakscalculator2) and [BioPython](https://biopython.org/) on a dataset of 20 sequences (300 codons each, 190 pairs).

| Comparison | Metric | n | Mean |diff| | Max |diff| | R² |
|---|---|---|---|---|---|
| eskaks Li vs KaKs_Calc LPB | dN | 154 | 0.000000 | 0.000001 | 1.000000 |
| eskaks Li vs KaKs_Calc LPB | dS | 175 | 0.000000 | 0.000001 | 1.000000 |
| eskaks Nei vs KaKs_Calc NG | dN | 124 | 0.000150 | 0.000416 | 0.999397 |
| eskaks Nei vs KaKs_Calc NG | dS | 184 | 0.001146 | 0.003315 | 0.995155 |
| eskaks Nei vs BioPython NG86 | dN | 190 | 0.000114 | 0.001169 | 0.998169 |
| eskaks Nei vs BioPython NG86 | dS | 190 | 0.000338 | 0.003554 | 0.996981 |

The Li model achieves near-exact agreement (R² = 1.0) with KaKs_Calculator's LPB method. Small differences in the Nei model are due to minor pathway-counting heuristics and are consistent with the inter-tool variation observed between KaKs_Calculator and BioPython themselves (R² = 0.993-0.996).

### Performance

Wall-clock time in milliseconds for pairwise dN/dS computation:

| Dataset | eskaks Nei (1t) | eskaks Nei (4t) | eskaks Li (1t) | eskaks Li (4t) | KaKs_Calc NG | KaKs_Calc LPB | yn00 | BioPython NG |
|---|---|---|---|---|---|---|---|---|
| Small (20 seq, 300 cod) | 3 ms | 2 ms | 6 ms | 6 ms | 34 ms | 48 ms | 8 ms | 610 ms |
| Medium (100 seq, 3000 cod) | 12 ms | 6 ms | 17 ms | 10 ms | 7,703 ms | 10,860 ms | 697 ms | 111,619 ms |
| Large (500 seq, 3000 cod) | 227 ms | 74 ms | 235 ms | 88 ms | 195,456 ms | 271,807 ms | - | - |

On the medium dataset (100 sequences, 3000 codons), eskaks Nei (4t) is **~1,280x faster** than KaKs_Calculator NG and **~18,600x faster** than BioPython. On the large dataset (500 sequences, 3000 codons), eskaks computes 124,750 pairs in under 100 ms.

### Reproducing Benchmarks

The full benchmark suite is in the `benchmark/` directory:

```bash
# Generate synthetic test datasets
python benchmark/generate_seqs.py

# Run cross-tool benchmarks (requires KaKs_Calculator and BioPython installed)
python benchmark/cross_tool_benchmark.py

# Results are saved to benchmark/cross_tool_results.json
# Plots are saved to benchmark/plots/
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
