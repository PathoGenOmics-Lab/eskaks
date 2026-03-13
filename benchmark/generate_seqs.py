#!/usr/bin/env python3
"""Generate synthetic aligned coding sequences for benchmarking.

Creates FASTA files with N sequences of length L (codons), introducing
random synonymous and nonsynonymous mutations at a specified rate.
"""
import random
import sys
import os

CODONS_BY_AA = {
    'F': ['TTT','TTC'], 'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'I': ['ATT','ATC','ATA'], 'M': ['ATG'], 'V': ['GTT','GTC','GTA','GTG'],
    'S': ['TCT','TCC','TCA','TCG','AGT','AGC'], 'P': ['CCT','CCC','CCA','CCG'],
    'T': ['ACT','ACC','ACA','ACG'], 'A': ['GCT','GCC','GCA','GCG'],
    'Y': ['TAT','TAC'], 'H': ['CAT','CAC'], 'Q': ['CAA','CAG'],
    'N': ['AAT','AAC'], 'K': ['AAA','AAG'], 'D': ['GAT','GAC'],
    'E': ['GAA','GAG'], 'C': ['TGT','TGC'], 'W': ['TGG'],
    'R': ['CGT','CGC','CGA','CGG','AGA','AGG'], 'G': ['GGT','GGC','GGA','GGG'],
}

CODON_TO_AA = {}
for aa, codons in CODONS_BY_AA.items():
    for c in codons:
        CODON_TO_AA[c] = aa

ALL_CODONS = list(CODON_TO_AA.keys())
BASES = 'ACGT'


def random_coding_sequence(n_codons):
    """Generate a random coding sequence of n_codons codons (no stop codons)."""
    seq = []
    for _ in range(n_codons):
        codon = random.choice(ALL_CODONS)
        seq.append(codon)
    return seq


def mutate_sequence(codons, syn_rate=0.02, nonsyn_rate=0.01):
    """Mutate a codon sequence with given syn and nonsyn rates per codon."""
    result = list(codons)
    for i in range(len(result)):
        r = random.random()
        if r < syn_rate:
            # Synonymous mutation: pick a different codon for the same AA
            aa = CODON_TO_AA[result[i]]
            synonyms = [c for c in CODONS_BY_AA[aa] if c != result[i]]
            if synonyms:
                result[i] = random.choice(synonyms)
        elif r < syn_rate + nonsyn_rate:
            # Nonsynonymous mutation: pick a codon for a different AA
            aa = CODON_TO_AA[result[i]]
            other_codons = [c for c in ALL_CODONS if CODON_TO_AA[c] != aa]
            result[i] = random.choice(other_codons)
    return result


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <n_seqs> <n_codons> <output_dir> [seed]")
        sys.exit(1)

    n_seqs = int(sys.argv[1])
    n_codons = int(sys.argv[2])
    output_dir = sys.argv[3]
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42

    random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    # Generate reference sequence
    ref = random_coding_sequence(n_codons)

    # Generate mutated sequences
    sequences = []
    for i in range(n_seqs):
        mutated = mutate_sequence(ref, syn_rate=0.02, nonsyn_rate=0.01)
        sequences.append((f"seq{i+1}", ''.join(mutated)))

    # Write FASTA
    fasta_path = os.path.join(output_dir, f"bench_{n_seqs}x{n_codons*3}.fasta")
    with open(fasta_path, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")

    # Write AXT format for KaKs_Calculator (all pairs)
    axt_path = os.path.join(output_dir, f"bench_{n_seqs}x{n_codons*3}.axt")
    with open(axt_path, 'w') as f:
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                f.write(f"{sequences[i][0]}-{sequences[j][0]}\n")
                f.write(f"{sequences[i][1]}\n")
                f.write(f"{sequences[j][1]}\n\n")

    # Write PHYLIP format for yn00 (sequential)
    phy_path = os.path.join(output_dir, f"bench_{n_seqs}x{n_codons*3}.phy")
    seq_len = n_codons * 3
    with open(phy_path, 'w') as f:
        f.write(f"  {n_seqs}  {seq_len}\n")
        for name, seq in sequences:
            # yn00 needs name padded to at least 10 chars
            f.write(f"{name:<10}{seq}\n")

    print(f"Generated {n_seqs} sequences of {n_codons*3} bp:")
    print(f"  FASTA: {fasta_path}")
    print(f"  AXT:   {axt_path}")
    print(f"  PHYLIP: {phy_path}")


if __name__ == '__main__':
    main()
