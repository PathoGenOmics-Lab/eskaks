# Genetic Codes

eskaks supports 20 NCBI translation tables via `--genetic-code <N>`.

## Available tables

```bash
eskaks --list-codes
```

| ID | Name | Common use |
|----|------|-----------|
| 1 | Standard | **Default**. Most nuclear genes |
| 2 | Vertebrate Mitochondrial | Human, mouse, fish mito |
| 3 | Yeast Mitochondrial | *S. cerevisiae* mito |
| 4 | Mold/Protozoan/Coelenterate Mito | Also Mycoplasma/Spiroplasma |
| 5 | Invertebrate Mitochondrial | *Drosophila*, *C. elegans* mito |
| 6 | Ciliate Nuclear | *Tetrahymena*, *Paramecium* |
| 9 | Echinoderm/Flatworm Mito | Sea urchin, planaria mito |
| 10 | Euplotid Nuclear | *Euplotes* |
| 11 | Bacterial/Archaeal/Plant Plastid | Prokaryotes, chloroplasts |
| 12 | Alternative Yeast Nuclear | *Candida* |
| 13 | Ascidian Mitochondrial | Tunicate mito |
| 14 | Alternative Flatworm Mito | Some flatworm mito |
| 16 | Chlorophycean Mito | Green algae mito |
| 21 | Trematode Mitochondrial | *Schistosoma* mito |
| 22 | *Scenedesmus obliquus* Mito | |
| 23 | Thraustochytrium Mito | |
| 24 | Rhabdopleuridae Mito | |
| 25 | Candidate Division SR1/Gracilibacteria | |
| 26 | Pachysolen tannophilus Nuclear | |
| 33 | Cephalodiscidae Mito (UAA=Tyr) | |

## Usage

```bash
# Standard code (default, same as --genetic-code 1)
eskaks genes.fasta

# Vertebrate mitochondrial
eskaks mito_genes.fasta --genetic-code 2

# Bacterial
eskaks prokaryote_genes.fasta --genetic-code 11
```

## How it works

The genetic code affects:
- **Site classification**: Which positions are synonymous vs nonsynonymous
- **Pathway analysis**: Which intermediate codons are stop codons (excluded from pathways)
- **Degeneracy classification** (Li model): 0-fold, 2-fold, 4-fold at each position

Example: `AGA` encodes Arginine in the standard code but is a **stop codon** in vertebrate mitochondrial code. This changes both the site counts and the pathway analysis for any codon pair involving `AGA`.
