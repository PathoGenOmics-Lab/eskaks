//! CDS reconstruction, codon translation, and N/S site counting.

use super::*;

/// Extract the CDS sequence from the reference, handling strand and multi-exon genes.
pub(crate) fn extract_cds_sequence(gene: &Gene, ref_seq: &[u8]) -> Vec<u8> {
    let mut cds = Vec::with_capacity(gene.length_bp);

    // Assemble the forward-strand bases in ASCENDING genomic order (regardless of
    // strand) so the single reverse-complement below yields the coding sequence in
    // the exact order genomic_to_cds_offset assumes (5'-most coding exon first).
    // gene.exons are stored descending for minus strand, so iterate them reversed;
    // otherwise a multi-exon minus gene ends up with its exons swapped — because
    // revcomp(A ++ B) = revcomp(B) ++ revcomp(A) — and every SNP is mis-codon'd.
    let ordered: Vec<_> = if gene.strand == Strand::Minus {
        gene.exons.iter().rev().collect()
    } else {
        gene.exons.iter().collect()
    };
    for exon in ordered {
        let start_idx = exon.start.saturating_sub(1); // Convert 1-based to 0-based
        let end_idx = exon.end.min(ref_seq.len());
        if start_idx >= ref_seq.len() || start_idx >= end_idx {
            continue;
        }
        cds.extend_from_slice(&ref_seq[start_idx..end_idx]);
    }

    if gene.strand == Strand::Minus {
        // Reverse complement the ascending-genomic buffer → coding sequence.
        cds = reverse_complement(&cds);
    }

    // Apply phase: skip leading bases
    if !gene.exons.is_empty() && gene.exons[0].phase > 0 {
        let skip = gene.exons[0].phase as usize;
        if skip < cds.len() {
            cds = cds[skip..].to_vec();
        }
    }

    // Uppercase
    for b in cds.iter_mut() {
        *b = b.to_ascii_uppercase();
    }

    // A CDS whose length is not a multiple of 3 signals a frame/annotation
    // problem (bad phase, wrong exon bounds, or an origin-spanning gene on a
    // circular genome). The trailing 1-2 bases are dropped by the len/3 loops;
    // flag it rather than silently skewing the site counts.
    if !cds.is_empty() && cds.len() % 3 != 0 {
        warn!(
            "Gene {} CDS length {} bp is not a multiple of 3 (phase/annotation issue?); trailing {} base(s) ignored",
            gene.name, cds.len(), cds.len() % 3
        );
    }

    cds
}

/// Convert a genomic position (1-based) to a CDS offset (0-based).
/// Returns None if the position is not within any exon of the gene.
pub(crate) fn genomic_to_cds_offset(gene: &Gene, pos: usize) -> Option<usize> {
    let mut cds_offset = 0usize;

    // Handle phase offset from first exon
    let phase_offset = if !gene.exons.is_empty() {
        gene.exons[0].phase as usize
    } else {
        0
    };

    match gene.strand {
        Strand::Plus => {
            for exon in &gene.exons {
                if pos >= exon.start && pos <= exon.end {
                    let offset = cds_offset + (pos - exon.start);
                    if offset >= phase_offset {
                        return Some(offset - phase_offset);
                    } else {
                        return None; // Within phase region
                    }
                }
                cds_offset += exon.end - exon.start + 1;
            }
        }
        Strand::Minus => {
            for exon in &gene.exons {
                // For minus strand, exons are sorted descending by start
                if pos >= exon.start && pos <= exon.end {
                    let offset = cds_offset + (exon.end - pos);
                    if offset >= phase_offset {
                        return Some(offset - phase_offset);
                    } else {
                        return None;
                    }
                }
                cds_offset += exon.end - exon.start + 1;
            }
        }
    }

    None
}

/// Count nonsynonymous (N) and synonymous (S) sites for a CDS sequence.
///
/// For each codon, enumerates all 9 possible single-nucleotide changes
/// and classifies each as synonymous or nonsynonymous. Fractional sites
/// are assigned as: S_sites = syn_changes/3, N_sites = nonsyn_changes/3.
pub(crate) fn count_sites(cds: &[u8], gc: &GeneticCode) -> (f64, f64) {
    let mut n_sites = 0.0f64;
    let mut s_sites = 0.0f64;

    let codons = cds.len() / 3;
    for i in 0..codons {
        let codon = [cds[i * 3], cds[i * 3 + 1], cds[i * 3 + 2]];
        let ref_aa = codon_to_aa(&codon, gc);
        if ref_aa.is_none() {
            continue; // Skip ambiguous codons
        }
        let ref_aa = ref_aa.unwrap();

        if ref_aa == b'*' {
            continue; // Skip stop codons
        }

        let mut syn = 0;
        let mut nonsyn = 0;

        for pos in 0..3 {
            for &alt_base in b"ACGT" {
                if alt_base == codon[pos] {
                    continue;
                }
                let mut alt_codon = codon;
                alt_codon[pos] = alt_base;
                if let Some(alt_aa) = codon_to_aa(&alt_codon, gc) {
                    if alt_aa == b'*' {
                        continue; // Exclude changes to stop codons from site counts
                    } else if alt_aa == ref_aa {
                        syn += 1;
                    } else {
                        nonsyn += 1;
                    }
                }
                // Skip ambiguous codons entirely (don't count as either)
            }
        }

        // Each codon contributes exactly 3 sites. With stop-codon changes
        // excluded, redistribute proportionally: S = 3 × syn/(syn+nonsyn)
        let valid_changes = (syn + nonsyn) as f64;
        if valid_changes > 0.0 {
            s_sites += 3.0 * syn as f64 / valid_changes;
            n_sites += 3.0 * nonsyn as f64 / valid_changes;
        }
    }

    (n_sites, s_sites)
}

/// Is the change `from`→`to` a transition (A↔G or C↔T)? Purine↔purine or
/// pyrimidine↔pyrimidine. Strand-symmetric, so it holds on either strand.
#[inline]
pub(crate) fn is_transition(from: u8, to: u8) -> bool {
    matches!(
        (from, to),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// Count N and S sites with mutation-spectrum weighting by a transition/
/// transversion rate ratio `kappa`.
///
/// This is the mutation-spectrum-aware generalisation of [`count_sites`]: each
/// candidate single-nucleotide change is weighted by its relative mutation rate
/// (`kappa` for a transition, `1.0` for a transversion) instead of counting 1.
/// It uses the *same* codon-level normalisation as [`count_sites`] — each codon
/// contributes exactly 3 sites, split between S and N by the rate-weighted
/// synonymous fraction — so at `kappa == 1.0` it reduces to [`count_sites`]
/// bit-for-bit (the weights become all-1 and the sums equal the plain counts).
/// Keeping the normalisation identical means `kappa` is the *only* thing that
/// changes, so the correction is not confounded with a normalisation switch.
///
/// Under a transition-biased spectrum (`kappa > 1`), synonymous changes at
/// 2-fold degenerate sites — reached almost exclusively by transitions — get
/// up-weighted, so a 2-fold synonymous-transition site moves from `1/3` toward
/// `kappa/(kappa+2)`. Fully (4-fold) degenerate positions are synonymous for
/// every change, so they stay `kappa`-invariant. Across the coding genome this
/// generally raises total S and lowers total N (individual codons can move
/// either way depending on whether their transition-reachable changes are
/// mostly synonymous or nonsynonymous).
pub(crate) fn count_sites_weighted(cds: &[u8], gc: &GeneticCode, kappa: f64) -> (f64, f64) {
    let mut n_sites = 0.0f64;
    let mut s_sites = 0.0f64;

    let codons = cds.len() / 3;
    for i in 0..codons {
        let codon = [cds[i * 3], cds[i * 3 + 1], cds[i * 3 + 2]];
        let ref_aa = match codon_to_aa(&codon, gc) {
            Some(aa) if aa != b'*' => aa,
            _ => continue, // Skip ambiguous and stop codons
        };

        // Pool rate-weighted synonymous / total over all three positions of the
        // codon, exactly as count_sites pools raw counts.
        let mut syn_rate = 0.0f64;
        let mut tot_rate = 0.0f64;
        for pos in 0..3 {
            for &alt_base in b"ACGT" {
                if alt_base == codon[pos] {
                    continue;
                }
                let mut alt_codon = codon;
                alt_codon[pos] = alt_base;
                if let Some(alt_aa) = codon_to_aa(&alt_codon, gc) {
                    if alt_aa == b'*' {
                        continue; // Exclude changes to stop codons
                    }
                    let w = if is_transition(codon[pos], alt_base) { kappa } else { 1.0 };
                    tot_rate += w;
                    if alt_aa == ref_aa {
                        syn_rate += w;
                    }
                }
                // Ambiguous alternates are skipped (contribute to neither)
            }
        }

        // Each codon contributes exactly 3 sites, split by the rate-weighted
        // synonymous fraction (reduces to count_sites at kappa == 1.0).
        if tot_rate > 0.0 {
            s_sites += 3.0 * syn_rate / tot_rate;
            n_sites += 3.0 * (tot_rate - syn_rate) / tot_rate;
        }
    }

    (n_sites, s_sites)
}

/// Look up amino acid for a codon using the genetic code table (Li indexing).
/// Returns None for codons with ambiguous bases.
pub(crate) fn codon_to_aa(codon: &[u8; 3], gc: &GeneticCode) -> Option<u8> {
    let b1 = base_to_li(codon[0])?;
    let b2 = base_to_li(codon[1])?;
    let b3 = base_to_li(codon[2])?;
    let idx = 16 * b1 + 4 * b2 + b3;
    Some(gc.aa_table[idx])
}

/// Convert a base to Li index: A=0, C=1, G=2, T=3.
#[inline]
pub(crate) fn base_to_li(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Complement a DNA base.
#[inline]
pub(crate) fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'U' | b'u' => b'A',
        _ => b'N',
    }
}

/// Reverse complement a DNA sequence.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Parse a reference FASTA file into a map of sequence_name -> sequence.
/// Uses simple manual parsing (no needletail dependency for this).
pub fn parse_reference_fasta(path: &std::path::Path) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    use anyhow::Context;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path)
        .with_context(|| format!("Failed to open reference FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                seqs.insert(name, std::mem::take(&mut current_seq));
            }
            // Take the first word as the sequence name
            let name = line.strip_prefix('>').unwrap_or("").split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
            current_seq = Vec::new();
        } else {
            current_seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(name) = current_name {
        seqs.insert(name, current_seq);
    }

    if seqs.is_empty() {
        anyhow::bail!("No sequences found in reference FASTA: {}", path.display());
    }

    info!("Loaded {} reference sequences", seqs.len());
    Ok(seqs)
}

