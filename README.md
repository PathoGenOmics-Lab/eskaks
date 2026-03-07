
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

**eskaks** is a command-line tool designed to calculate dN/dS (also known as Ka/Ks) between pairs of aligned coding sequences. It supports both the Nei and Li models for calculating synonymous (dS) and nonsynonymous (dN) substitution rates. The tool leverages parallel processing for efficient computation on large datasets and provides flexible options for grouping sequences and summarizing results by groups or lineages.

## Features

- **Models**: Choose between Nei (default) or Li models to calculate dN/dS.
- **Parallelization**: Speed up computation using multiple threads (configurable via `--workers`).
- **Filtering Conditions**: Results are only written if certain conditions are met (e.g., ds > 0.0, dn > 0.0, ps > 0.0, pn > 0.0), ensuring that only valid ratios are reported.
- **Group Summaries**:
  - **Group Average**: Compute average dN/dS between predefined groups of sequences.
  - **Lineage Summary**: Compute average dN/dS for each sequence against sequences grouped by lineage.  
  - **Configurable Lineage Grouping**: Either split lineage by `_` (default) or use the first character of each sequence name (`--first_letter_lineage`).

## Differences Between the Li and Nei Models

The **Nei (1986)** and **Li (1993)** models are two approaches for calculating synonymous (dS) and nonsynonymous (dN) substitution rates between aligned coding sequences. Although both methods aim to estimate dN/dS (Ka/Ks) ratios, they differ in complexity, underlying assumptions, and computational requirements.

1. **Methodological Complexity**:  
   - **Nei Model**: The Nei method is relatively straightforward, relying on direct counting of synonymous and nonsynonymous substitutions and then correcting for the total number of synonymous and nonsynonymous sites. Its simplicity makes it quicker and less resource-intensive, serving as a standard or baseline approach.
   - **Li Model**: The Li approach is more sophisticated. It often involves precomputed tables and a probabilistic framework that accounts for variations in codon usage, nucleotide frequencies, and transition-transversion biases. This complexity can provide more nuanced results but requires more computation and careful implementation.

2. **Corrections and Assumptions**:  
   - **Nei Model**: Typically uses simpler assumptions and fewer corrections, focusing on directly observed differences in codons.
   - **Li Model**: Incorporates more corrections and a more refined probabilistic model. By adjusting for factors like codon bias and different classes of substitutions, it aims to provide a more accurate and biologically realistic estimate of dN/dS.

3. **Computational Cost**:  
   - **Nei Model**: Easier and faster to compute, suitable for large datasets and quick analyses.
   - **Li Model**: More computationally demanding due to pre-calculations, matrix-based probability models, and more complex algorithms. This can offer improved accuracy and insight but at greater computational expense.

In summary, the Nei model is often chosen for its simplicity and speed, serving as a commonly used benchmark method. The Li model, although more complex and computationally heavier, may yield more detailed and refined estimates of dN/dS that capture subtle evolutionary patterns.  

## Input

The input should be a FASTA file containing multiple aligned coding sequences of equal length (no partial codons at the end). Each sequence must be fully aligned with the others.

## Output

- A TSV file with pairwise results: `Seq1`, `Seq2`, `Sd`, `Sn`, `S`, `N`, `ps`, `pn`, `ds`, `dn`, and `dn/ds`.
- Optional additional summary files:
  - `*_group_dn_ds.tsv`: Mean dN/dS between specified groups.
  - `*_lineage_summary.tsv`: Mean dN/dS per lineage grouping.

## Usage

```bash
eskaks <input_file> [options]
```
### Required Arguments:
`<input_file>`: Aligned sequences in FASTA format.

### Common Options:
- `-o, --output <PREFIX>`: Base name for output files (default: output).
- `-w, --workers <N>`: Number of parallel threads (default: 4).
- `--model <nei|li>`: Model to use for calculation (default: nei).
- `--lineage`: Compute summary results by lineage.
- `--group_average`: Compute average dN/dS between groups.
- `--first_letter_lineage`: If used with --lineage, group sequences by the first letter of their names instead of splitting by `_`.

## Examples
1. Basic Nei model calculation:
   ```bash
   eskaks input.fasta
   ```
   This produces `output.tsv` with pairwise dN/dS calculations (Nei model).
2. Li model and multiple workers:
   ```bash
   eskaks input.fasta --model li --workers 8
   ```
   Uses the Li model and spawns 8 parallel threads for calculation.
3. Lineage summary:
   ```bash
   eskaks input.fasta --lineage
   ```
   Besides the standard pairwise file, this also produces `output_lineage_summary.tsv` summarizing mean dN/dS per lineage     
   group.
4. Group average by first letter:
   ```bash
   eskaks input.fasta --lineage --first_letter_lineage
   ```
   Groups sequences by their first character rather than splitting by `_` for the lineage summary.

## Requirements
## Installation

---
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
eskaks is developed with ❤️ by:
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
      <a href="" title="Desing">🎨</a>
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
