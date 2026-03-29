//! FASTA reading, validation, deduplication, and filtering.

use anyhow::{bail, Context};
use log::{info, warn};
use needletail::parse_fastx_file;

use crate::codon;
use crate::models::Model;

/// Fully processed input data ready for pairwise computation.
pub struct SequenceData {
    /// Original sequence IDs in input order.
    pub ids: Vec<String>,
    /// Map from sequence index → unique index.
    pub uidx_by_id: Vec<usize>,
    /// Codon indices for unique sequences only.
    pub unique_codon_indices: Vec<Vec<u8>>,
    /// Number of unique sequences.
    pub n_unique: usize,
}

/// Read FASTA, validate, filter, and deduplicate.
pub fn load_sequences(input_file: &str, model: Model, min_codons: usize) -> anyhow::Result<SequenceData> {
    info!("Reading sequences from: {}", input_file);

    let mut all_codon_indices: Vec<Vec<u8>> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    let mut gap_count_total = 0usize;
    let mut seqs_with_gaps = 0usize;

    {
        let mut reader = parse_fastx_file(input_file)
            .with_context(|| format!("Failed to open input file '{}'", input_file))?;
        while let Some(record) = reader.next() {
            let rec = record?;
            let seq = rec.seq();
            let gaps = codon::count_gaps(&seq);
            if gaps > 0 {
                seqs_with_gaps += 1;
                gap_count_total += gaps;
            }
            if seq.len() % 3 != 0 {
                warn!(
                    "Sequence '{}' length {} is not divisible by 3; {} trailing base(s) ignored.",
                    String::from_utf8_lossy(rec.id()),
                    seq.len(),
                    seq.len() % 3
                );
            }
            let codons = codon::fasta_to_codon_indices(&seq, model);
            let id_str = String::from_utf8_lossy(rec.id()).into_owned();
            all_codon_indices.push(codons);
            ids.push(id_str);
        }
    }

    if ids.is_empty() {
        bail!("No sequences found in '{}'", input_file);
    }
    if ids.len() < 2 {
        bail!(
            "Need at least 2 sequences for pairwise comparison, found {}",
            ids.len()
        );
    }
    info!("Found {} total sequences.", ids.len());

    // Warn about sequences with no valid codons
    let empty_seqs: Vec<&str> = ids
        .iter()
        .zip(all_codon_indices.iter())
        .filter(|(_, codons)| codons.iter().all(|&c| c == codon::INVALID_CODON))
        .map(|(id, _)| id.as_str())
        .collect();
    if !empty_seqs.is_empty() {
        warn!(
            "{} sequence(s) have no valid codons: {}",
            empty_seqs.len(),
            if empty_seqs.len() <= 5 {
                empty_seqs.join(", ")
            } else {
                format!(
                    "{}, ... and {} more",
                    empty_seqs[..5].join(", "),
                    empty_seqs.len() - 5
                )
            }
        );
    }

    // Diagnostics
    if seqs_with_gaps > 0 {
        warn!(
            "{} sequence(s) contain {} total gap character(s) ('-' or '.'), treated as ambiguous.",
            seqs_with_gaps, gap_count_total
        );
    }
    if let Some(first_len) = all_codon_indices.first().map(|v| v.len()) {
        let mismatched = all_codon_indices
            .iter()
            .filter(|v| v.len() != first_len)
            .count();
        if mismatched > 0 {
            warn!(
                "{} sequence(s) have different codon lengths than the first ({} codons). Sequences may not be aligned.",
                mismatched, first_len
            );
        }
        info!(
            "Alignment length: {} codons ({} bp).",
            first_len,
            first_len * 3
        );
    }

    // Filter by --min-codons
    if min_codons > 0 {
        let before = ids.len();
        let keep: Vec<usize> = (0..ids.len())
            .filter(|&i| {
                all_codon_indices[i]
                    .iter()
                    .filter(|&&c| c != codon::INVALID_CODON)
                    .count()
                    >= min_codons
            })
            .collect();
        if keep.len() < before {
            warn!(
                "Filtered {} sequence(s) with fewer than {} valid codons.",
                before - keep.len(),
                min_codons
            );
            let new_ids: Vec<String> = keep
                .iter()
                .map(|&i| std::mem::take(&mut ids[i]))
                .collect();
            let new_codons: Vec<Vec<u8>> = keep
                .iter()
                .map(|&i| std::mem::take(&mut all_codon_indices[i]))
                .collect();
            ids = new_ids;
            all_codon_indices = new_codons;
        }
        if ids.is_empty() {
            bail!("All sequences filtered out by --min-codons {}", min_codons);
        }
    }

    // Sort-based deduplication
    let (uidx_by_id, unique_codon_indices, n_unique) = deduplicate(&ids, all_codon_indices);
    info!("Found {} unique sequences.", n_unique);

    Ok(SequenceData {
        ids,
        uidx_by_id,
        unique_codon_indices,
        n_unique,
    })
}

/// Sort-based deduplication on codon indices.
/// Returns (uidx_by_id, unique_codon_indices, n_unique).
fn deduplicate(
    ids: &[String],
    mut all_codon_indices: Vec<Vec<u8>>,
) -> (Vec<usize>, Vec<Vec<u8>>, usize) {
    let mut sorted_indices: Vec<usize> = (0..ids.len()).collect();
    sorted_indices.sort_unstable_by(|&a, &b| all_codon_indices[a].cmp(&all_codon_indices[b]));

    let mut uidx_by_sorted = Vec::with_capacity(ids.len());
    let mut unique_repr_indices: Vec<usize> = Vec::new();
    let mut current_uidx = 0usize;

    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        if pos > 0 && all_codon_indices[orig_idx] != all_codon_indices[sorted_indices[pos - 1]] {
            current_uidx += 1;
        }
        if pos == 0 || all_codon_indices[orig_idx] != all_codon_indices[sorted_indices[pos - 1]] {
            unique_repr_indices.push(orig_idx);
        }
        uidx_by_sorted.push(current_uidx);
    }
    let n_unique = unique_repr_indices.len();

    let mut uidx_by_id: Vec<usize> = vec![0; ids.len()];
    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        uidx_by_id[orig_idx] = uidx_by_sorted[pos];
    }

    let unique_codon_indices: Vec<Vec<u8>> = unique_repr_indices
        .iter()
        .map(|&orig_idx| std::mem::take(&mut all_codon_indices[orig_idx]))
        .collect();

    (uidx_by_id, unique_codon_indices, n_unique)
}
