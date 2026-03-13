#!/usr/bin/env python3
"""Benchmark BioPython's dN/dS calculation (Nei-Gojobori model).
Usage: biopython_bench.py <fasta_file> <output_file>
"""
import sys
import time
import warnings
warnings.filterwarnings("ignore")

from Bio import SeqIO
from Bio.codonalign import CodonSeq
from Bio.codonalign.codonseq import cal_dn_ds


def nei_gojobori_dn_ds(seq1_str, seq2_str):
    """Calculate dN/dS using BioPython's Nei-Gojobori implementation."""
    try:
        codon_seq1 = CodonSeq(seq1_str)
        codon_seq2 = CodonSeq(seq2_str)
        dn, ds = cal_dn_ds(codon_seq1, codon_seq2, method='NG86')
        return dn, ds
    except Exception:
        return float('nan'), float('nan')


def main():
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]

    records = list(SeqIO.parse(fasta_file, "fasta"))
    seqs = [(r.id, str(r.seq)) for r in records]
    n = len(seqs)

    start = time.time()
    results = []
    for i in range(n):
        for j in range(i+1, n):
            dn, ds = nei_gojobori_dn_ds(seqs[i][1], seqs[j][1])
            ratio = dn / ds if ds and ds > 0 else float('inf') if dn and dn > 0 else 0.0
            results.append((seqs[i][0], seqs[j][0], dn, ds, ratio))
    elapsed = time.time() - start

    with open(output_file, 'w') as f:
        f.write("Seq1\tSeq2\tdN\tdS\tdN/dS\n")
        for r in results:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]:.6f}\t{r[3]:.6f}\t{r[4]:.6f}\n")

    print(f"BioPython NG86: {n} seqs, {len(results)} pairs in {elapsed:.3f}s")


if __name__ == '__main__':
    main()
