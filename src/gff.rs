//! GFF3 file parser for CDS feature extraction.
//!
//! Parses GFF3 annotation files to extract CDS (coding sequence) features,
//! groups multi-exon genes by their Parent or gene_id attribute, and returns
//! sorted CDS regions per gene.

use anyhow::Context;
use log::warn;
use std::collections::BTreeMap;
use std::io::BufRead;
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
    let reader = crate::input::open_text(path, "GFF3 file")?;

    // Map: (seqid, gene_id) -> (name, seqid, strand, Vec<CdsExon>). Keying on the pair,
    // not the id alone, keeps two genes that (against the GFF3 spec) reuse an id on
    // different contigs separate — otherwise their exons pool into one gene whose CDS is
    // assembled from a single contig and the other gene silently vanishes.
    type GeneGroup = (String, String, Strand, Vec<CdsExon>);
    let mut gene_map: BTreeMap<(String, String), GeneGroup> = BTreeMap::new();
    // Format diagnostics: did we ever see a proper 9-column GFF3 record?
    let mut any_structured = false;
    let mut first_line: Option<String> = None;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {} of GFF3", line_no + 1))?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if first_line.is_none() {
            first_line = Some(line.chars().take(60).collect());
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        any_structured = true;

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
        // A single exon wider than any real chromosome is almost certainly a malformed
        // coordinate; keeping it would blow up `length_bp` and the CDS `Vec::with_capacity`
        // (out-of-memory / capacity-overflow panic). Skip it instead of crashing.
        const MAX_CDS_SPAN: usize = 1 << 30; // ~1.07 Gb
        if end - start >= MAX_CDS_SPAN {
            warn!("CDS span {}..{} at GFF3 line {} is implausibly large, skipping", start, end, line_no + 1);
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
            .entry((seqid.clone(), gene_id))
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
                Strand::Minus => exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
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
        // Distinguish "valid GFF3 with no CDS rows" from "this isn't a GFF3 at all".
        if !any_structured && first_line.is_some() {
            let hint = match first_line.as_deref() {
                Some(l) if l.starts_with('>') => " (the first line starts with '>': is this a FASTA?)",
                _ => "",
            };
            anyhow::bail!(
                "{}: no tab-separated 9-column GFF3 records found — is this really a GFF3?{}",
                path.display(), hint
            );
        }
        anyhow::bail!(
            "No CDS features found in GFF3 file: {} (the file parsed as GFF3 but has no 'CDS' rows — \
             check the feature type column).",
            path.display()
        );
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

/// URL-decode a GFF3 attribute value. Decoded bytes are accumulated and converted
/// as UTF-8 at the end, so a multi-byte escape such as `%C3%A9` reassembles to `é`
/// instead of the two Latin-1 characters `Ã©` a per-byte `byte as char` would give.
fn url_decode(s: &str) -> String {
    let mut bytes: Vec<u8> = Vec::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                bytes.push(byte);
            } else {
                bytes.push(b'%');
                bytes.extend_from_slice(hex.as_bytes());
            }
        } else {
            let mut buf = [0u8; 4];
            bytes.extend_from_slice(c.encode_utf8(&mut buf).as_bytes());
        }
    }
    String::from_utf8_lossy(&bytes).into_owned()
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
    fn colliding_gene_id_on_different_contigs_stays_separate() {
        // Regression: a gene id reused on two contigs (against the GFF3 spec) must not
        // merge into one gene — which would assemble one gene's CDS from a single contig
        // and silently drop the other. Keying on (seqid, id) keeps them separate.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t1\t15\t.\t+\t0\tParent=g1;gene=alpha
chr2\t.\tCDS\t1\t15\t.\t+\t0\tParent=g1;gene=beta\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 2, "colliding ids on different contigs must stay separate");
        let seqids: Vec<&str> = genes.iter().map(|g| g.seqid.as_str()).collect();
        assert!(seqids.contains(&"chr1") && seqids.contains(&"chr2"), "one gene per contig");
        for g in &genes {
            assert_eq!(g.length_bp, 15, "each gene keeps its own single 15 bp CDS, not a doubled one");
        }
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

    #[test]
    fn url_decode_reassembles_utf8() {
        // Regression: %C3%A9 is the two-byte UTF-8 for 'é' and must decode to 'é',
        // not the two Latin-1 chars 'Ã©' a per-byte `byte as char` would produce.
        assert_eq!(url_decode("caf%C3%A9"), "café");
        assert_eq!(url_decode("a%20b"), "a b");
        assert_eq!(url_decode("plain"), "plain");
        assert_eq!(url_decode("bad%ZZ"), "bad%ZZ"); // invalid escape kept literally
    }

    // ---- gene-id / gene-name attribute fallbacks -------------------------

    #[test]
    fn gene_id_falls_back_to_gene_id_then_id() {
        // No Parent=: extract_gene_id tries gene_id=, then ID=. gene_id wins here.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tgene_id=myGene;Name=NN
chr1\t.\tCDS\t400\t450\t.\t+\t0\tID=cdsOnly;locus_tag=LT1\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 2);
        // Sorted by start: myGene (100) then the ID-keyed gene (400).
        assert_eq!(genes[0].name, "NN", "Name= should supply the display name");
        // Second gene has no Parent/gene_id: ID= supplies the group id, locus_tag= the name.
        assert_eq!(genes[1].name, "LT1", "locus_tag= should be used when gene=/Name= absent");
    }

    #[test]
    fn gene_name_falls_back_to_id_when_no_name_attrs() {
        // Only Parent= present: no gene=/Name=/locus_tag= => name defaults to the id.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tParent=g9\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes[0].name, "g9");
    }

    // ---- malformed / non-CDS line handling -------------------------------

    #[test]
    fn unknown_strand_line_is_skipped() {
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t?\t0\tParent=bad
chr1\t.\tCDS\t400\t450\t.\t-\t0\tParent=ok\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "ok");
        assert_eq!(genes[0].strand, Strand::Minus);
    }

    #[test]
    fn invalid_or_dot_phase_defaults_to_zero() {
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t9\tParent=g1;gene=badPhase
chr1\t.\tCDS\t400\t450\t.\t+\t.\tParent=g2;gene=dotPhase
chr1\t.\tCDS\t700\t750\t.\t+\t2\tParent=g3;gene=twoPhase\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 3);
        assert_eq!(genes[0].exons[0].phase, 0, "invalid phase '9' => 0");
        assert_eq!(genes[1].exons[0].phase, 0, "phase '.' => 0");
        assert_eq!(genes[2].exons[0].phase, 2, "phase '2' preserved");
    }

    #[test]
    fn line_with_fewer_than_nine_fields_is_skipped() {
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+
chr1\t.\tCDS\t400\t450\t.\t+\t0\tParent=ok\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "ok");
    }

    #[test]
    fn no_cds_features_is_error() {
        // Only gene/mRNA rows, no CDS => nothing to analyse => error.
        let gff = "\
##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1
chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=mrna1;Parent=gene1\n";
        let f = write_temp_gff(gff);
        assert!(parse_gff3(f.path()).is_err(), "a GFF3 with no CDS must error");
    }

    #[test]
    fn minus_strand_orders_exons_descending() {
        // For a minus-strand gene exons are stored 3'->5' (descending genomic start),
        // but Gene.start (for plotting) is the minimum genomic coordinate.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t-\t0\tParent=g1;gene=m
chr1\t.\tCDS\t200\t250\t.\t-\t0\tParent=g1;gene=m\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].exons.len(), 2);
        assert_eq!(genes[0].exons[0].start, 200, "first stored exon is the highest-coordinate one");
        assert_eq!(genes[0].exons[1].start, 100);
        assert_eq!(genes[0].start, 100, "Gene.start is the minimum genomic coordinate");
        assert_eq!(genes[0].length_bp, 102);
    }

    #[test]
    fn implausibly_large_cds_span_is_skipped() {
        // A ~10 Gb CDS would blow up length_bp / the CDS Vec::with_capacity; it must be
        // skipped with a warning, not panic/OOM. The valid gene survives.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t1\t9999999999\t.\t+\t0\tParent=g1;gene=huge
chr1\t.\tCDS\t100\t150\t.\t+\t0\tParent=g2;gene=ok\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "ok");
    }

    #[test]
    fn fasta_content_as_gff3_reports_wrong_format() {
        let gff = ">chr1\nATGGCTGCTAAA\n";
        let f = write_temp_gff(gff);
        let err = parse_gff3(f.path()).unwrap_err().to_string();
        assert!(err.contains("GFF3") && err.contains("FASTA"), "err: {}", err);
    }
}
