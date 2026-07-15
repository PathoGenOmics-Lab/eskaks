//! GFF3 file parser for CDS feature extraction.
//!
//! Parses GFF3 annotation files to extract CDS (coding sequence) features,
//! groups multi-exon genes by their Parent or gene_id attribute, and returns
//! sorted CDS regions per gene.

use anyhow::Context;
use log::warn;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Strand orientation of a genomic feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

/// A single CDS exon from the GFF3 file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CdsExon {
    /// Chromosome/contig name
    pub seqid: String,
    /// 1-based start position (inclusive)
    pub start: usize,
    /// 1-based end position (inclusive)
    pub end: usize,
    /// Strand
    pub strand: Strand,
    /// Phase (0, 1, or 2) — number of bases to skip at the start for reading frame
    pub phase: u8,
}

/// A gene composed of one or more CDS exons.
#[derive(Debug, Clone)]
pub struct Gene {
    /// Gene name (from Name=, gene=, or ID= attribute)
    pub name: String,
    /// Chromosome/contig
    pub seqid: String,
    /// Strand
    pub strand: Strand,
    /// CDS exons sorted by position (ascending for +, descending for -)
    pub exons: Vec<CdsExon>,
    /// Total CDS length in bp
    pub length_bp: usize,
    /// Genomic start of the first exon (for plotting/ordering)
    pub start: usize,
}

/// Parse a GFF3 file and return a list of genes with their CDS regions.
///
/// Groups CDS entries by their Parent or gene_id attribute.
/// Returns genes sorted by chromosome and start position.
pub fn parse_gff3(path: &Path) -> anyhow::Result<Vec<Gene>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open GFF3 file: {}", path.display()))?;
    let reader = BufReader::new(file);

    // Map: gene_id -> (name, seqid, strand, Vec<CdsExon>)
    let mut gene_map: BTreeMap<String, (String, String, Strand, Vec<CdsExon>)> = BTreeMap::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {} of GFF3", line_no + 1))?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Only process CDS features
        if fields[2] != "CDS" {
            continue;
        }

        let seqid = fields[0].to_string();
        let start: usize = fields[3]
            .parse()
            .with_context(|| format!("Invalid start at GFF3 line {}", line_no + 1))?;
        let end: usize = fields[4]
            .parse()
            .with_context(|| format!("Invalid end at GFF3 line {}", line_no + 1))?;
        // GFF3 requires start <= end (strand is given separately). A malformed line
        // with end < start would underflow `end - start` (panic in debug, a huge
        // wrapped length → capacity-overflow panic in release), so skip it.
        if end < start {
            warn!("CDS end {} < start {} at GFF3 line {}, skipping", end, start, line_no + 1);
            continue;
        }
        let strand = match fields[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => {
                warn!("Unknown strand '{}' at GFF3 line {}, skipping", fields[6], line_no + 1);
                continue;
            }
        };
        let phase: u8 = match fields[7] {
            "0" => 0,
            "1" => 1,
            "2" => 2,
            "." => 0, // Default phase
            _ => {
                warn!("Invalid phase '{}' at GFF3 line {}, using 0", fields[7], line_no + 1);
                0
            }
        };

        let attributes = fields[8];
        let gene_id = extract_gene_id(attributes);
        let gene_name = extract_gene_name(attributes).unwrap_or_else(|| gene_id.clone());

        let entry = gene_map
            .entry(gene_id)
            .or_insert_with(|| (gene_name.clone(), seqid.clone(), strand, Vec::new()));

        entry.3.push(CdsExon {
            seqid,
            start,
            end,
            strand,
            phase,
        });
    }

    // Convert map to sorted gene list
    let mut genes: Vec<Gene> = gene_map
        .into_iter()
        .filter_map(|(_id, (name, seqid, strand, mut exons))| {
            if exons.is_empty() {
                return None;
            }

            // Sort exons by position
            match strand {
                Strand::Plus => exons.sort_by_key(|e| e.start),
                Strand::Minus => exons.sort_by(|a, b| b.start.cmp(&a.start)),
            }

            let length_bp: usize = exons.iter().map(|e| e.end - e.start + 1).sum();
            let start = exons.iter().map(|e| e.start).min().unwrap();

            let display_name = name;
            Some(Gene {
                name: display_name,
                seqid,
                strand,
                exons,
                length_bp,
                start,
            })
        })
        .collect();

    // Sort by chromosome and start position
    genes.sort_by(|a, b| a.seqid.cmp(&b.seqid).then(a.start.cmp(&b.start)));

    if genes.is_empty() {
        anyhow::bail!("No CDS features found in GFF3 file: {}", path.display());
    }

    Ok(genes)
}

/// Extract the gene identifier from GFF3 attributes.
/// Tries Parent first, then gene_id, then ID.
fn extract_gene_id(attributes: &str) -> String {
    // Try Parent= (most common for CDS features in GFF3)
    if let Some(val) = extract_attribute(attributes, "Parent") {
        return val;
    }
    // Try gene_id= (GTF-style)
    if let Some(val) = extract_attribute(attributes, "gene_id") {
        return val;
    }
    // Fall back to ID=
    if let Some(val) = extract_attribute(attributes, "ID") {
        return val;
    }
    // Last resort: use the whole attributes string
    attributes.to_string()
}

/// Extract a human-readable gene name from GFF3 attributes.
fn extract_gene_name(attributes: &str) -> Option<String> {
    // Try gene= first (common in prokaryotic GFF3)
    if let Some(val) = extract_attribute(attributes, "gene") {
        return Some(val);
    }
    // Try Name=
    if let Some(val) = extract_attribute(attributes, "Name") {
        return Some(val);
    }
    // Try locus_tag=
    if let Some(val) = extract_attribute(attributes, "locus_tag") {
        return Some(val);
    }
    None
}

/// Extract a single attribute value from a GFF3 attributes string.
fn extract_attribute(attributes: &str, key: &str) -> Option<String> {
    for entry in attributes.split(';') {
        let entry = entry.trim();
        if let Some(rest) = entry.strip_prefix(key) {
            if let Some(val) = rest.strip_prefix('=') {
                // URL-decode common escapes
                return Some(url_decode(val));
            }
        }
    }
    None
}

/// Simple URL decoding for GFF3 attribute values.
fn url_decode(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                result.push(byte as char);
            } else {
                result.push('%');
                result.push_str(&hex);
            }
        } else {
            result.push(c);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_gff(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn parse_simple_gff3() {
        let gff = "\
##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=geneA
chr1\t.\tCDS\t100\t200\t.\t+\t0\tParent=gene1;gene=geneA\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "geneA");
        assert_eq!(genes[0].exons.len(), 1);
        assert_eq!(genes[0].exons[0].start, 100);
        assert_eq!(genes[0].exons[0].end, 200);
        assert_eq!(genes[0].strand, Strand::Plus);
    }

    #[test]
    fn cds_with_end_before_start_is_skipped_not_panicking() {
        // Regression: a malformed CDS (end < start) used to underflow `end - start`
        // (panic in debug, capacity-overflow panic in release). It must be skipped.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t200\t100\t.\t+\t0\tParent=g1;gene=bad
chr1\t.\tCDS\t300\t360\t.\t+\t0\tParent=g2;gene=ok\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        // The bad line is dropped; the valid gene survives.
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "ok");
    }

    #[test]
    fn parse_multi_exon_gene() {
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tParent=gene1;gene=geneA
chr1\t.\tCDS\t200\t250\t.\t+\t0\tParent=gene1;gene=geneA\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].exons.len(), 2);
        assert_eq!(genes[0].length_bp, 102); // (150-100+1) + (250-200+1)
    }

    #[test]
    fn parse_minus_strand() {
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t200\t.\t-\t0\tParent=gene1;gene=geneA\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes[0].strand, Strand::Minus);
    }
}
