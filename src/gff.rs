//! GFF3 file parser for CDS feature extraction.
//!
//! Parses GFF3 annotation files to extract CDS (coding sequence) features,
//! groups multi-exon genes by their Parent or gene_id attribute, and returns
//! sorted CDS regions per gene.

use anyhow::Context;
use log::warn;
use std::collections::{BTreeMap, BTreeSet, HashMap};
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

/// A CDS row buffered until the whole file has been read.
struct RawCds {
    exon: CdsExon,
    key: CdsKey,
    /// Display name from gene=/Name=/locus_tag=, if the row carried one.
    name: Option<String>,
}

/// How a CDS row identifies the feature it belongs to.
enum CdsKey {
    /// One or more `Parent=` ids. GFF3 allows a comma-separated list when one exon is
    /// shared by several transcripts. A Parent is either a transcript row (then the gene
    /// is that row's own Parent) or the gene itself (prokaryotic GFF3 has no mRNA level).
    Parents(Vec<String>),
    /// A `gene_id=` or `ID=` value: an id of the row itself, never a pointer to a
    /// transcript row, so it is used as the gene id directly.
    Direct(String),
}

/// One gene under assembly: its CDS exons kept separately per transcript, so several
/// isoforms of the same gene can be collapsed to one representative below.
struct GeneGroup {
    name: Option<String>,
    strand: Strand,
    transcripts: BTreeMap<String, Vec<CdsExon>>,
}

/// Parse a GFF3 file and return a list of genes with their CDS regions.
///
/// Groups CDS entries by the gene they belong to: the `Parent=` attribute, resolved
/// through any `mRNA`/`transcript` row so that several isoforms of one gene collapse
/// into a single gene (the longest CDS wins), else `gene_id=`, else `ID=`.
/// Returns genes sorted by chromosome and start position.
pub fn parse_gff3(path: &Path) -> anyhow::Result<Vec<Gene>> {
    let reader = crate::input::open_text(path, "GFF3 file")?;

    // CDS rows are buffered and only resolved into genes at EOF: a CDS may reference an
    // mRNA row that appears further down the file, and that mRNA row is what says which
    // gene the transcript belongs to.
    let mut cds_rows: Vec<RawCds> = Vec::new();
    // (seqid, transcript ID) -> gene id, harvested from mRNA/transcript rows.
    let mut transcript_gene: HashMap<(String, String), String> = HashMap::new();
    // Format diagnostics: did we ever see a proper 9-column GFF3 record?
    let mut any_structured = false;
    let mut first_line: Option<String> = None;
    // CDS rows carrying no Parent=/gene_id=/ID= at all: they name no feature, so they
    // are skipped rather than keyed on their own attribute text.
    let mut unidentified_cds: usize = 0;

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

        let attributes = fields[8];
        // GTF writes `key "value"; ...` where GFF3 writes `key=value;...`. Read as GFF3,
        // a GTF row matches no Parent=/gene_id=/ID= at all, and the parser used to fall
        // back to the raw attribute text as the gene id: every exon became its own "gene"
        // named after its own attribute string, exit status 0, no warning. Refuse loudly
        // instead of inventing genes.
        if looks_like_gtf(attributes) {
            anyhow::bail!(gtf_not_gff3_error(path, line_no + 1, attributes));
        }

        // Remember transcript -> gene links, so the CDS rows of several isoforms of one
        // gene resolve to that single gene instead of to one "gene" per transcript.
        if is_transcript_type(fields[2]) {
            if let (Some(id), Some(parent)) = (
                extract_attribute(attributes, "ID"),
                extract_attribute(attributes, "Parent"),
            ) {
                if !id.is_empty() && !parent.is_empty() {
                    transcript_gene.insert((fields[0].to_string(), id), parent);
                }
            }
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

        let Some(key) = extract_cds_key(attributes) else {
            // No Parent=, no gene_id=, no ID=: this row names no feature. Keying it on
            // its raw attribute text (the old behaviour) either merges unrelated rows
            // that happen to share a comment or splits one gene into many; both invent
            // a gene table out of nothing, so skip the row and say so.
            unidentified_cds += 1;
            if unidentified_cds <= 3 {
                warn!(
                    "CDS at GFF3 line {} has no Parent=, gene_id= or ID= attribute, skipping (nothing says which gene it belongs to)",
                    line_no + 1
                );
            }
            continue;
        };

        cds_rows.push(RawCds {
            exon: CdsExon {
                seqid,
                start,
                end,
                strand,
                phase,
            },
            key,
            name: extract_gene_name(attributes),
        });
    }

    if unidentified_cds > 3 {
        warn!(
            "{} CDS rows in total had no Parent=, gene_id= or ID= attribute and were skipped",
            unidentified_cds
        );
    }

    // Map: (seqid, gene id) -> gene under assembly. Keying on the pair, not the id alone,
    // keeps two genes that (against the GFF3 spec) reuse an id on different contigs
    // separate: otherwise their exons pool into one gene whose CDS is assembled from a
    // single contig and the other gene silently vanishes.
    let mut gene_map: BTreeMap<(String, String), GeneGroup> = BTreeMap::new();
    for raw in cds_rows {
        let seqid = raw.exon.seqid.clone();
        let (gene_id, tx_ids) = resolve_gene_and_transcripts(&seqid, &raw.key, &transcript_gene);
        let entry = gene_map
            .entry((seqid, gene_id))
            .or_insert_with(|| GeneGroup {
                name: raw.name.clone(),
                strand: raw.exon.strand,
                transcripts: BTreeMap::new(),
            });
        if entry.name.is_none() {
            entry.name = raw.name.clone();
        }
        for tx in tx_ids {
            entry.transcripts.entry(tx).or_default().push(raw.exon.clone());
        }
    }

    // Multi-isoform bookkeeping, reported once at the end instead of per gene.
    let mut collapsed_genes = 0usize;
    let mut dropped_transcripts = 0usize;
    let mut collapse_examples: Vec<String> = Vec::new();

    // Convert map to sorted gene list
    let mut genes: Vec<Gene> = gene_map
        .into_iter()
        .filter_map(|((seqid, gene_id), group)| {
            let n_transcripts = group.transcripts.len();
            // One gene, several transcripts: keep one representative, the longest CDS,
            // which is the usual convention. Emitting one row per isoform would not just
            // duplicate the gene in the table, it would enlarge the multiple-testing
            // family and shift every q-value in the genome. BTreeMap iterates in id
            // order and the comparison is strict, so ties go to the first id: the choice
            // is deterministic for a given file.
            let mut chosen: Option<(String, Vec<CdsExon>)> = None;
            let mut best_len = 0usize;
            for (tx_id, exons) in group.transcripts {
                let len: usize = exons.iter().map(|e| e.end - e.start + 1).sum();
                if chosen.is_none() || len > best_len {
                    best_len = len;
                    chosen = Some((tx_id, exons));
                }
            }
            let (tx_id, mut exons) = chosen?;
            if exons.is_empty() {
                return None;
            }
            if n_transcripts > 1 {
                collapsed_genes += 1;
                dropped_transcripts += n_transcripts - 1;
                if collapse_examples.len() < 3 {
                    collapse_examples.push(format!(
                        "{} ({} transcripts, kept {} at {} bp)",
                        gene_id, n_transcripts, tx_id, best_len
                    ));
                }
            }

            // Sort exons by position
            match group.strand {
                Strand::Plus => exons.sort_by_key(|e| e.start),
                Strand::Minus => exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
            }

            let length_bp: usize = exons.iter().map(|e| e.end - e.start + 1).sum();
            let start = exons.iter().map(|e| e.start).min().unwrap();

            let display_name = group.name.unwrap_or(gene_id);
            Some(Gene {
                name: display_name,
                seqid,
                strand: group.strand,
                exons,
                length_bp,
                start,
            })
        })
        .collect();

    if collapsed_genes > 0 {
        warn!(
            "{} gene(s) in {} have several transcripts; kept the longest CDS of each and dropped {} alternative transcript(s), so every gene is counted and tested exactly once (e.g. {}). Convert to a single-isoform GFF3 first if you need a specific transcript.",
            collapsed_genes,
            path.display(),
            dropped_transcripts,
            collapse_examples.join("; ")
        );
    }

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
        if unidentified_cds > 0 {
            anyhow::bail!(
                "{}: all {} CDS row(s) lack a Parent=, gene_id= or ID= attribute, so no gene \
                 could be identified. eskaks groups CDS exons into genes by Parent= (resolved \
                 through the mRNA row when there is one), else gene_id=, else ID=; add one of \
                 them to the attributes column.",
                path.display(),
                unidentified_cds
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

/// Feature types that sit between a gene and its CDS rows, i.e. whose `ID` a CDS names
/// as its `Parent` and whose own `Parent` is the gene.
fn is_transcript_type(feature_type: &str) -> bool {
    matches!(
        feature_type.to_ascii_lowercase().as_str(),
        "mrna" | "transcript" | "primary_transcript" | "pseudogenic_transcript"
    )
}

/// Extract the identifier a CDS row uses to say what it belongs to.
/// Tries Parent first, then gene_id, then ID. Returns None when the row carries none of
/// them: such a row names no feature, and inventing one from its attribute text is how a
/// GTF used to come out as one "gene" per exon.
fn extract_cds_key(attributes: &str) -> Option<CdsKey> {
    // Parent= is the usual one for CDS features in GFF3. GFF3 allows a comma-separated
    // list of parents (one exon shared by two isoforms); a comma inside a single value
    // must be written %2C, so splitting before URL-decoding is safe.
    if let Some(raw) = extract_attribute_raw(attributes, "Parent") {
        let ids: Vec<String> = raw
            .split(',')
            .map(|p| url_decode(p.trim()))
            .filter(|p| !p.is_empty())
            .collect();
        if !ids.is_empty() {
            return Some(CdsKey::Parents(ids));
        }
    }
    // gene_id= (a GFF3 file that spells the gene id this way, not a GTF: GTF is rejected
    // before we get here because it writes `gene_id "x"` with no '=').
    for key in ["gene_id", "ID"] {
        if let Some(val) = extract_attribute(attributes, key) {
            if !val.is_empty() {
                return Some(CdsKey::Direct(val));
            }
        }
    }
    None
}

/// Resolve a CDS row's key to the gene it belongs to and the transcript(s) it is part of.
///
/// A `Parent` is either a transcript row, in which case the gene is that row's own
/// `Parent`, or the gene itself (prokaryotic GFF3 has no mRNA level, and then transcript
/// and gene id coincide).
fn resolve_gene_and_transcripts(
    seqid: &str,
    key: &CdsKey,
    transcript_gene: &HashMap<(String, String), String>,
) -> (String, Vec<String>) {
    match key {
        CdsKey::Direct(id) => (id.clone(), vec![id.clone()]),
        CdsKey::Parents(parents) => {
            let genes: BTreeSet<String> = parents
                .iter()
                .map(|p| {
                    transcript_gene
                        .get(&(seqid.to_string(), p.clone()))
                        .cloned()
                        .unwrap_or_else(|| p.clone())
                })
                .collect();
            match genes.len() {
                1 => (
                    genes.into_iter().next().expect("one gene"),
                    parents.clone(),
                ),
                _ => {
                    // An exon claimed by transcripts of different genes: there is no single
                    // owner, so keep it as its own entry rather than silently giving it to one.
                    let joined = parents.join(",");
                    warn!(
                        "CDS with Parent={} is shared by {} different genes; keeping it as a separate entry",
                        joined,
                        genes.len()
                    );
                    (joined.clone(), vec![joined])
                }
            }
        }
    }
}

/// Does this attributes column use the GTF spelling `key "value"` instead of the GFF3
/// spelling `key=value`?
fn looks_like_gtf(attributes: &str) -> bool {
    // Cheap reject first: GFF3 values very rarely contain a double quote, GTF always does.
    if !attributes.contains('"') {
        return false;
    }
    attributes.split(';').any(|entry| {
        let entry = entry.trim();
        let Some((key, value)) = entry.split_once(|c: char| c.is_whitespace()) else {
            return false;
        };
        let value = value.trim();
        // The key must be a bare identifier (a GFF3 `Note=text with "quotes"` splits into a
        // key containing '=' and is rejected here), and the value fully quoted.
        !key.is_empty()
            && key.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
            && value.len() >= 2
            && value.starts_with('"')
            && value.ends_with('"')
    })
}

/// The error raised when the annotation turns out to be a GTF.
fn gtf_not_gff3_error(path: &Path, line_no: usize, attributes: &str) -> String {
    let sample: String = attributes.chars().take(70).collect();
    let p = path.display();
    format!(
        "{p}: this is a GTF, not a GFF3. Line {line_no} writes its attributes the GTF way, \
         `key \"value\"` instead of `key=value`:\n  {sample}\n\
         eskaks reads GFF3 only, so convert the annotation first, for example:\n  \
         gffread {p} -o annotation.gff3\n  \
         (or: agat_convert_sp_gxf2gxf.pl -g {p} -o annotation.gff3)\n\
         Parsing the GTF as it stands would find no Parent=, gene_id= or ID= on any row and \
         would break every gene into one entry per exon. Converting is not cosmetic either: a \
         GTF CDS excludes the stop codon, which is a separate stop_codon row, and a converter \
         puts the coding sequence back together."
    )
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

/// Extract a single attribute value from a GFF3 attributes string, URL-decoded.
fn extract_attribute(attributes: &str, key: &str) -> Option<String> {
    extract_attribute_raw(attributes, key).map(url_decode)
}

/// Extract a single attribute value without decoding it, for callers that must split the
/// raw text first (a `Parent` list, where an encoded `%2C` is not a separator).
fn extract_attribute_raw<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    for entry in attributes.split(';') {
        let entry = entry.trim();
        if let Some(rest) = entry.strip_prefix(key) {
            if let Some(val) = rest.strip_prefix('=') {
                return Some(val);
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
            // Decode only a COMPLETE two-hex-digit escape; a truncated tail (`%4`, a bare
            // `%`) or a non-hex pair (`%ZZ`) is kept literal instead of decoding to a
            // stray control byte. `from_str_radix` alone would accept a single digit.
            match u8::from_str_radix(&hex, 16) {
                Ok(byte) if hex.len() == 2 => bytes.push(byte),
                _ => {
                    bytes.push(b'%');
                    bytes.extend_from_slice(hex.as_bytes());
                }
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
        // A truncated single-hex-digit tail must stay literal too, not decode to a stray
        // control byte (from_str_radix accepts one digit, so the length must be checked).
        assert_eq!(url_decode("z%4"), "z%4");
        assert_eq!(url_decode("end%A"), "end%A");
        assert_eq!(url_decode("tail%"), "tail%");
        assert_eq!(url_decode("mid%4z"), "mid%4z"); // '%' + '4' + 'z': non-hex pair, literal
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

    // ---- GTF is refused, never silently shattered ------------------------

    /// A real two-exon GTF gene as Ensembl/GENCODE/SnpEff hand it out.
    const TWO_EXON_GTF: &str = "\
#!genome-build toy
chr1\teskaks\tgene\t1\t399\t.\t+\t.\tgene_id \"gA\";
chr1\teskaks\ttranscript\t1\t399\t.\t+\t.\tgene_id \"gA\"; transcript_id \"tA\";
chr1\teskaks\tCDS\t1\t198\t.\t+\t0\tgene_id \"gA\"; transcript_id \"tA\"; exon_number \"1\";
chr1\teskaks\tCDS\t202\t399\t.\t+\t0\tgene_id \"gA\"; transcript_id \"tA\"; exon_number \"2\";
chr1\teskaks\tstop_codon\t397\t399\t.\t+\t0\tgene_id \"gA\"; transcript_id \"tA\";
";

    #[test]
    fn gtf_is_rejected_with_a_conversion_hint() {
        let f = write_temp_gff(TWO_EXON_GTF);
        let err = parse_gff3(f.path()).unwrap_err().to_string();
        assert!(err.contains("GTF"), "the error must name the format: {}", err);
        assert!(
            err.contains("gffread") || err.contains("agat"),
            "the error must tell the user how to convert: {}",
            err
        );
    }

    #[test]
    fn gtf_never_produces_genes_named_after_their_own_attributes() {
        // Regression: extract_gene_id used to fall back to the whole attribute string, so
        // this exact file parsed with exit status 0 into two "genes" called
        // `gene_id "gA"; transcript_id "tA"; exon_number "1";` and the same with
        // exon_number "2". Any Ok(..) here means that silent path is back.
        let f = write_temp_gff(TWO_EXON_GTF);
        match parse_gff3(f.path()) {
            Ok(genes) => panic!(
                "a GTF must never parse into a gene table, got {:?}",
                genes.iter().map(|g| g.name.as_str()).collect::<Vec<_>>()
            ),
            Err(e) => assert!(!e.to_string().contains("no 'CDS' rows"), "wrong diagnosis: {}", e),
        }
    }

    #[test]
    fn gff3_attributes_with_quotes_or_spaces_are_not_mistaken_for_gtf() {
        // The GTF detector must not fire on legitimate GFF3 free-text values, which may
        // contain spaces and even quotes.
        assert!(!looks_like_gtf("Parent=g1;product=hypothetical protein"));
        assert!(!looks_like_gtf("Parent=g1;Note=the so called \"core\" region"));
        assert!(!looks_like_gtf("ID=cds1;Dbxref=GO:0003677"));
        assert!(!looks_like_gtf("."));
        assert!(!looks_like_gtf(""));
        // And it must fire on the GTF spellings, with or without the trailing semicolon.
        assert!(looks_like_gtf("gene_id \"gA\"; transcript_id \"tA\";"));
        assert!(looks_like_gtf("gene_id \"gA\""));
    }

    #[test]
    fn a_gff3_carrying_quoted_free_text_still_parses() {
        // End-to-end guard for the detector: this is a GFF3, quotes and all.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tParent=g1;gene=geneA;Note=the so called \"core\" region\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "geneA");
    }

    // ---- multi-isoform GFF3 collapses to one gene ------------------------

    #[test]
    fn multi_isoform_gff3_collapses_to_one_row_per_gene() {
        // NCBI-style three-level GFF3: gene -> mRNA -> CDS, two isoforms of gene01.
        // Each CDS names a different Parent, so the gene used to appear twice: a duplicated
        // row in the table AND an extra member of the multiple-testing family, which moves
        // every q-value in the genome. It must collapse to one gene, the longest CDS.
        let gff = "\
##gff-version 3
chr1\t.\tgene\t1\t399\t.\t+\t.\tID=g1;Name=gene01
chr1\t.\tmRNA\t1\t399\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\t.\tmRNA\t1\t198\t.\t+\t.\tID=g1.t2;Parent=g1
chr1\t.\tCDS\t1\t198\t.\t+\t0\tID=c1;Parent=g1.t1;gene=gene01
chr1\t.\tCDS\t202\t399\t.\t+\t0\tID=c2;Parent=g1.t1;gene=gene01
chr1\t.\tCDS\t1\t198\t.\t+\t0\tID=c3;Parent=g1.t2;gene=gene01
chr1\t.\tgene\t441\t695\t.\t+\t.\tID=g2;Name=gene02
chr1\t.\tCDS\t441\t695\t.\t+\t0\tID=c4;Parent=g2;gene=gene02\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 2, "one row per gene, not one per transcript: {:?}", genes.iter().map(|g| (g.name.as_str(), g.length_bp)).collect::<Vec<_>>());
        assert_eq!(genes[0].name, "gene01");
        assert_eq!(genes[0].exons.len(), 2, "the two-exon isoform is the one kept");
        assert_eq!(genes[0].length_bp, 396, "longest CDS wins (198 + 198), not the 198 bp isoform");
        assert_eq!(genes[1].name, "gene02");
        assert_eq!(genes[1].length_bp, 255);
    }

    #[test]
    fn longest_transcript_wins_regardless_of_file_order() {
        // The longer isoform is listed second here, so "keep the longest" cannot be
        // satisfied by accident by keeping the first one seen.
        let gff = "\
##gff-version 3
chr1\t.\tmRNA\t1\t99\t.\t+\t.\tID=t1;Parent=gX
chr1\t.\tmRNA\t1\t300\t.\t+\t.\tID=t2;Parent=gX
chr1\t.\tCDS\t1\t99\t.\t+\t0\tParent=t1
chr1\t.\tCDS\t1\t300\t.\t+\t0\tParent=t2\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].length_bp, 300);
        assert_eq!(genes[0].name, "gX", "with no gene=/Name=, the gene id names the gene, not the transcript id");
    }

    #[test]
    fn equal_length_isoforms_collapse_deterministically() {
        // A tie must still yield exactly one gene, and the same one on every run.
        let gff = "\
##gff-version 3
chr1\t.\tmRNA\t1\t99\t.\t+\t.\tID=t2;Parent=gX
chr1\t.\tmRNA\t201\t299\t.\t+\t.\tID=t1;Parent=gX
chr1\t.\tCDS\t1\t99\t.\t+\t0\tParent=t2
chr1\t.\tCDS\t201\t299\t.\t+\t0\tParent=t1\n";
        let f = write_temp_gff(gff);
        let first = parse_gff3(f.path()).unwrap();
        let second = parse_gff3(f.path()).unwrap();
        assert_eq!(first.len(), 1);
        assert_eq!(first[0].length_bp, 99);
        assert_eq!(first[0].start, second[0].start, "the tie-break must be stable");
        assert_eq!(first[0].start, 201, "ties go to the first transcript id in order (t1)");
    }

    #[test]
    fn single_isoform_through_an_mrna_row_keys_on_the_gene() {
        // The common eukaryotic layout with one transcript: the gene id, not the mRNA id,
        // identifies the gene, so a later isoform would join it instead of forming a new one.
        let gff = "\
##gff-version 3
chr1\t.\tgene\t1\t300\t.\t+\t.\tID=gene-LOC1
chr1\t.\tmRNA\t1\t300\t.\t+\t.\tID=rna-XM_1;Parent=gene-LOC1
chr1\t.\tCDS\t1\t300\t.\t+\t0\tID=cds-XP_1;Parent=rna-XM_1\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "gene-LOC1");
    }

    #[test]
    fn exon_shared_by_two_isoforms_counts_for_the_transcript_that_is_kept() {
        // GFF3 allows a comma-separated Parent list for an exon shared by two isoforms.
        // The shared exon belongs to both, so the longest-CDS choice must see it in both.
        let gff = "\
##gff-version 3
chr1\t.\tmRNA\t1\t400\t.\t+\t.\tID=t1;Parent=gX
chr1\t.\tmRNA\t1\t400\t.\t+\t.\tID=t2;Parent=gX
chr1\t.\tCDS\t1\t100\t.\t+\t0\tParent=t1,t2
chr1\t.\tCDS\t301\t400\t.\t+\t0\tParent=t2\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1, "one gene, not one per parent list");
        assert_eq!(genes[0].length_bp, 200, "t2 = shared exon + its own exon");
        assert_eq!(genes[0].exons.len(), 2);
    }

    // ---- a CDS with no identity at all -----------------------------------

    #[test]
    fn cds_without_any_id_attribute_is_skipped_not_keyed_on_its_raw_text() {
        // Regression: the old fallback turned the attribute column itself into the gene
        // id, so this row became a gene called `Note=orphan;product=hypothetical protein`.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tNote=orphan;product=hypothetical protein
chr1\t.\tCDS\t400\t450\t.\t+\t0\tParent=g2;gene=ok\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1, "the identity-less row must not become a gene");
        assert_eq!(genes[0].name, "ok");
        assert!(
            !genes.iter().any(|g| g.name.contains("Note=") || g.name.contains("product=")),
            "no gene may be named after raw attribute text"
        );
    }

    #[test]
    fn a_file_of_identity_less_cds_rows_is_an_error_not_one_fake_gene() {
        // All rows share the same attribute text, so the old fallback pooled every exon
        // into a single "gene" named after it. There is no gene here: say so.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tsource=manual
chr1\t.\tCDS\t400\t450\t.\t+\t0\tsource=manual\n";
        let f = write_temp_gff(gff);
        let err = parse_gff3(f.path()).unwrap_err().to_string();
        assert!(err.contains("Parent="), "the error must say what is missing: {}", err);
    }

    #[test]
    fn empty_parent_value_falls_through_to_id() {
        // `Parent=` with nothing after it names no feature; ID= still does.
        let gff = "\
##gff-version 3
chr1\t.\tCDS\t100\t150\t.\t+\t0\tParent=;ID=cds1\n";
        let f = write_temp_gff(gff);
        let genes = parse_gff3(f.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].name, "cds1");
    }
}
