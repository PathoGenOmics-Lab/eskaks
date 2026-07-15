//! FASTA reading, validation, deduplication, and filtering.

use anyhow::{bail, Context};
use log::{info, warn};
use needletail::{parse_fastx_file, parse_fastx_reader};

use crate::codon;
use crate::models::Model;

/// Fully processed input data ready for pairwise computation.
#[derive(Debug)]
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

/// Returns true if the input should be read from stdin.
fn is_stdin(path: &str) -> bool {
    path == "-" || path == "/dev/stdin"
}

/// Read FASTA, validate, filter, and deduplicate.
/// `stop_codons` is an optional set of codon indices that are stop codons for the active genetic code.
pub fn load_sequences(
    input_file: &str,
    model: Model,
    min_codons: usize,
    stop_codons: Option<&[u8]>,
) -> anyhow::Result<SequenceData> {
    let source = if is_stdin(input_file) { "stdin" } else { input_file };
    info!("Reading sequences from: {}", source);

    let mut all_codon_indices: Vec<Vec<u8>> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    let mut gap_count_total = 0usize;
    let mut seqs_with_gaps = 0usize;

    // Buffer for stdin data (must outlive the reader)
    let stdin_buf: Vec<u8> = if is_stdin(input_file) {
        let mut buf = Vec::new();
        std::io::Read::read_to_end(&mut std::io::stdin().lock(), &mut buf)
            .context("Failed to read from stdin")?;
        buf
    } else {
        Vec::new()
    };

    {
        let mut reader: Box<dyn needletail::FastxReader> = if is_stdin(input_file) {
            parse_fastx_reader(&stdin_buf[..])
                .context("Failed to parse FASTA from stdin")?
        } else {
            parse_fastx_file(input_file)
                .with_context(|| format!("Failed to open input file '{}'", input_file))?
        };
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

            // Check for internal stop codons (last codon excluded — it's expected)
            if let Some(stops) = stop_codons {
                let id = String::from_utf8_lossy(rec.id());
                let last = codons.len().saturating_sub(1);
                for (pos, &c) in codons.iter().enumerate() {
                    if pos < last && c != codon::INVALID_CODON && stops.contains(&c) {
                        warn!(
                            "Sequence '{}' has internal stop codon at codon position {} (0-based). \
                             This may indicate a frameshift, pseudogene, or wrong reading frame.",
                            id, pos
                        );
                        break; // one warning per sequence is enough
                    }
                }
            }

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
        // Pairwise comparison needs at least two sequences; --min-codons can drop
        // below that without emptying the set, so check the count, not just emptiness.
        if ids.len() < 2 {
            bail!(
                "Only {} sequence(s) remain after --min-codons {} (need at least 2 for pairwise comparison)",
                ids.len(), min_codons
            );
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn write_temp_fasta(name: &str, content: &str) -> String {
        let path = format!("/tmp/eskaks_input_test_{}.fasta", name);
        fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn load_two_identical_sequences() {
        let path = write_temp_fasta("ident", ">s1\nATGGCTGCT\n>s2\nATGGCTGCT\n");
        let data = load_sequences(&path, Model::Nei, 0, None).unwrap();
        assert_eq!(data.ids.len(), 2);
        assert_eq!(data.n_unique, 1, "identical seqs should deduplicate to 1");
        assert_eq!(data.uidx_by_id[0], data.uidx_by_id[1]);
        fs::remove_file(&path).ok();
    }

    #[test]
    fn load_two_different_sequences() {
        let path = write_temp_fasta("diff", ">s1\nATGGCTGCT\n>s2\nATGATTGCT\n");
        let data = load_sequences(&path, Model::Nei, 0, None).unwrap();
        assert_eq!(data.ids.len(), 2);
        assert_eq!(data.n_unique, 2);
        assert_ne!(data.uidx_by_id[0], data.uidx_by_id[1]);
        fs::remove_file(&path).ok();
    }

    #[test]
    fn load_filters_by_min_codons() {
        // seq1: 3 valid codons, seq2: 1 valid + 2 N-codons, seq3: 3 valid
        let path = write_temp_fasta(
            "mincodon",
            ">s1\nATGGCTGCT\n>s2\nATGNNNNNN\n>s3\nATGATTGCT\n",
        );
        let data = load_sequences(&path, Model::Nei, 2, None).unwrap();
        // s2 has only 1 valid codon → filtered out
        assert_eq!(data.ids.len(), 2);
        assert!(data.ids.contains(&"s1".to_string()));
        assert!(data.ids.contains(&"s3".to_string()));
        fs::remove_file(&path).ok();
    }

    #[test]
    fn load_empty_file_errors() {
        let path = write_temp_fasta("empty", "");
        let result = load_sequences(&path, Model::Nei, 0, None);
        assert!(result.is_err());
        fs::remove_file(&path).ok();
    }

    #[test]
    fn load_single_sequence_errors() {
        let path = write_temp_fasta("single", ">s1\nATGGCTGCT\n");
        let result = load_sequences(&path, Model::Nei, 0, None);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("at least 2"), "error: {}", err);
        fs::remove_file(&path).ok();
    }

    #[test]
    fn load_nonexistent_file_errors() {
        let result = load_sequences("/tmp/eskaks_nonexistent_xyz.fasta", Model::Nei, 0, None);
        assert!(result.is_err());
    }

    #[test]
    fn dedup_preserves_order() {
        // 4 seqs: A, B, A, C → uidx: [0, 1, 0, 2], n_unique=3
        let ids = vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()];
        let codons = vec![
            vec![1u8, 2, 3],
            vec![4, 5, 6],
            vec![1, 2, 3], // same as s1
            vec![7, 8, 9],
        ];
        let (uidx, unique, n) = deduplicate(&ids, codons);
        assert_eq!(n, 3);
        assert_eq!(uidx[0], uidx[2], "s1 and s3 should share unique index");
        assert_ne!(uidx[0], uidx[1]);
        assert_ne!(uidx[1], uidx[3]);
        assert_eq!(unique.len(), 3);
    }
}
