//! Core pN/pS analysis from VCF data.
//!
//! For each gene, reconstructs reference and alternate codons from SNPs,
//! classifies mutations as synonymous or nonsynonymous using the genetic
//! code table, and computes pN/pS ratios.

use crate::gff::{Gene, Strand};
use crate::genetic_code::GeneticCode;
use crate::vcf::VcfSnp;
use log::{info, warn};
use std::collections::HashMap;

/// Result of pN/pS analysis for a single gene.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GenePnPs {
    /// Gene name
    pub name: String,
    /// Total CDS length in bp
    pub length_bp: usize,
    /// Number of nonsynonymous sites (fractional)
    pub n_sites: f64,
    /// Number of synonymous sites (fractional)
    pub s_sites: f64,
    /// Proportion of nonsynonymous substitutions (nonsyn_snps / N_sites)
    pub pn: f64,
    /// Proportion of synonymous substitutions (syn_snps / S_sites)
    pub ps: f64,
    /// pN/pS ratio
    pub pn_ps: f64,
    /// Count of nonsynonymous SNPs (AF-weighted if --af-weighted)
    pub nonsyn_snps: f64,
    /// Count of synonymous SNPs (AF-weighted if --af-weighted)
    pub syn_snps: f64,
    /// Total SNPs in this gene (AF-weighted if --af-weighted)
    pub total_snps: f64,
    /// Genomic start position (for plotting)
    pub genome_start: usize,
    /// Chromosome
    pub chrom: String,
}

/// Compute pN/pS for all genes given a reference sequence, gene annotations, and SNPs.
///
/// If `af_weighted` is true, each SNP contributes its allele frequency to the
/// syn/nonsyn count instead of 1.0 (πN/πS instead of pN/pS).
pub fn compute_pn_ps(
    reference: &HashMap<String, Vec<u8>>,
    genes: &[Gene],
    snps: &[VcfSnp],
    gc: &GeneticCode,
    af_weighted: bool,
) -> Vec<GenePnPs> {
    // Index SNPs by chromosome and position for fast lookup
    let mut snp_map: HashMap<&str, Vec<&VcfSnp>> = HashMap::new();
    for snp in snps {
        snp_map.entry(snp.chrom.as_str()).or_default().push(snp);
    }
    // Sort SNPs by position within each chromosome
    for snps_list in snp_map.values_mut() {
        snps_list.sort_by_key(|s| s.pos);
    }

    let mut results = Vec::with_capacity(genes.len());

    for gene in genes {
        let ref_seq = match reference.get(&gene.seqid) {
            Some(seq) => seq,
            None => {
                warn!("Reference sequence not found for {}, skipping gene {}", gene.seqid, gene.name);
                continue;
            }
        };

        let chrom_snps = snp_map.get(gene.seqid.as_str());

        // Extract the full CDS sequence from the reference
        let cds_seq = extract_cds_sequence(gene, ref_seq);
        if cds_seq.len() < 3 {
            warn!("Gene {} CDS too short ({} bp), skipping", gene.name, cds_seq.len());
            continue;
        }

        // Count S and N sites from reference codons
        let (n_sites, s_sites) = count_sites(&cds_seq, gc);

        // Find SNPs that fall within this gene's CDS regions
        // Counts are f64 to support AF-weighted mode (πN/πS)
        let mut nonsyn_count = 0.0f64;
        let mut syn_count = 0.0f64;

        if let Some(snps_list) = chrom_snps {
            for snp in snps_list.iter() {
                // Check if this SNP falls within any exon of this gene
                if let Some(cds_offset) = genomic_to_cds_offset(gene, snp.pos) {
                    let codon_idx = cds_offset / 3;
                    let pos_in_codon = cds_offset % 3;
                    let codon_start = codon_idx * 3;

                    if codon_start + 3 > cds_seq.len() {
                        continue;
                    }

                    // Get reference codon from the extracted CDS
                    let ref_codon = [
                        cds_seq[codon_start],
                        cds_seq[codon_start + 1],
                        cds_seq[codon_start + 2],
                    ];

                    // Verify VCF REF allele matches the reference sequence
                    let expected_ref = if gene.strand == Strand::Minus {
                        complement(cds_seq[codon_start + pos_in_codon])
                    } else {
                        cds_seq[codon_start + pos_in_codon]
                    };
                    if snp.ref_allele != expected_ref {
                        warn!(
                            "VCF REF mismatch at {}:{} — VCF says {}, reference has {}. Skipping.",
                            snp.chrom, snp.pos,
                            snp.ref_allele as char, expected_ref as char
                        );
                        continue;
                    }

                    // For each ALT allele at this position
                    for (alt_idx, alt_base) in snp.alt_alleles.iter().enumerate() {
                        // Weight: AF if weighted mode, 1.0 otherwise
                        let weight = if af_weighted {
                            snp.alt_freqs.get(alt_idx).copied().unwrap_or(1.0)
                        } else {
                            1.0
                        };

                        // Build alternate codon
                        let mut alt_codon = ref_codon;
                        let alt_in_cds = if gene.strand == Strand::Minus {
                            complement(*alt_base)
                        } else {
                            *alt_base
                        };
                        alt_codon[pos_in_codon] = alt_in_cds;

                        // Look up amino acids
                        let ref_aa = codon_to_aa(&ref_codon, gc);
                        let alt_aa = codon_to_aa(&alt_codon, gc);

                        match (ref_aa, alt_aa) {
                            (Some(r), Some(a)) if r != b'*' && a != b'*' => {
                                if r == a {
                                    syn_count += weight;
                                } else {
                                    nonsyn_count += weight;
                                }
                            }
                            // Skip: ambiguous codons or mutations to/from stop
                            _ => continue,
                        }
                    }
                }
            }
        }

        let total_snps = nonsyn_count + syn_count;
        let pn = if n_sites > 0.0 {
            nonsyn_count / n_sites
        } else {
            0.0
        };
        let ps = if s_sites > 0.0 {
            syn_count / s_sites
        } else {
            0.0
        };
        let pn_ps = if ps > 0.0 {
            pn / ps
        } else if pn > 0.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };

        results.push(GenePnPs {
            name: gene.name.clone(),
            length_bp: gene.length_bp,
            n_sites,
            s_sites,
            pn,
            ps,
            pn_ps,
            nonsyn_snps: nonsyn_count,
            syn_snps: syn_count,
            total_snps,
            genome_start: gene.start,
            chrom: gene.seqid.clone(),
        });
    }

    info!("Computed pN/pS for {} genes", results.len());
    results
}

/// Extract the CDS sequence from the reference, handling strand and multi-exon genes.
fn extract_cds_sequence(gene: &Gene, ref_seq: &[u8]) -> Vec<u8> {
    let mut cds = Vec::with_capacity(gene.length_bp);

    // Exons are already sorted in coding order (ascending for +, descending for -)
    for exon in &gene.exons {
        let start_idx = exon.start.saturating_sub(1); // Convert 1-based to 0-based
        let end_idx = exon.end.min(ref_seq.len());
        if start_idx >= ref_seq.len() || start_idx >= end_idx {
            continue;
        }
        cds.extend_from_slice(&ref_seq[start_idx..end_idx]);
    }

    if gene.strand == Strand::Minus {
        // Reverse complement
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

    cds
}

/// Convert a genomic position (1-based) to a CDS offset (0-based).
/// Returns None if the position is not within any exon of the gene.
fn genomic_to_cds_offset(gene: &Gene, pos: usize) -> Option<usize> {
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
fn count_sites(cds: &[u8], gc: &GeneticCode) -> (f64, f64) {
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

/// Look up amino acid for a codon using the genetic code table (Li indexing).
/// Returns None for codons with ambiguous bases.
fn codon_to_aa(codon: &[u8; 3], gc: &GeneticCode) -> Option<u8> {
    let b1 = base_to_li(codon[0])?;
    let b2 = base_to_li(codon[1])?;
    let b3 = base_to_li(codon[2])?;
    let idx = 16 * b1 + 4 * b2 + b3;
    Some(gc.aa_table[idx])
}

/// Convert a base to Li index: A=0, C=1, G=2, T=3.
#[inline]
fn base_to_li(b: u8) -> Option<usize> {
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
fn complement(b: u8) -> u8 {
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
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

/// Write pN/pS results to a file.
pub fn write_results(
    results: &[GenePnPs],
    prefix: &str,
    format: &crate::models::OutputFormat,
) -> anyhow::Result<String> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let ext = format.extension();
    let output_path = format!("{}_pnps.{}", prefix, ext);

    match format {
        crate::models::OutputFormat::Json => {
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(file, "[")?;
            for (i, r) in results.iter().enumerate() {
                let comma = if i + 1 < results.len() { "," } else { "" };
                writeln!(
                    file,
                    "  {{\"gene\":\"{}\",\"length_bp\":{},\"N_sites\":{:.4},\"S_sites\":{:.4},\"pN\":{},\"pS\":{},\"pN_pS\":{},\"nonsyn_snps\":{:.4},\"syn_snps\":{:.4},\"total_snps\":{:.4}}}{}",
                    r.name, r.length_bp, r.n_sites, r.s_sites,
                    format_json_f64(r.pn), format_json_f64(r.ps), format_json_f64(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps, comma
                )?;
            }
            writeln!(file, "]")?;
        }
        _ => {
            let sep = format.separator();
            let mut file = BufWriter::new(File::create(&output_path)?);
            writeln!(
                file,
                "Gene{s}Length_bp{s}N_sites{s}S_sites{s}pN{s}pS{s}pN/pS{s}Nonsyn_SNPs{s}Syn_SNPs{s}Total_SNPs",
                s = sep
            )?;
            for r in results {
                writeln!(
                    file,
                    "{}{s}{}{s}{:.4}{s}{:.4}{s}{:.6}{s}{:.6}{s}{}{s}{:.4}{s}{:.4}{s}{:.4}",
                    r.name, r.length_bp, r.n_sites, r.s_sites,
                    r.pn, r.ps, format_pnps(r.pn_ps),
                    r.nonsyn_snps, r.syn_snps, r.total_snps,
                    s = sep
                )?;
            }
        }
    }

    Ok(output_path)
}

fn format_pnps(v: f64) -> String {
    if v.is_nan() {
        "NaN".to_string()
    } else if v.is_infinite() {
        "inf".to_string()
    } else {
        format!("{:.6}", v)
    }
}

fn format_json_f64(v: f64) -> String {
    if v.is_nan() || v.is_infinite() {
        "null".to_string()
    } else if v == 0.0 {
        "0.000000".to_string()
    } else {
        format!("{:.6}", v)
    }
}

/// Generate an SVG Manhattan-style plot of pN/pS per gene along the genome.
pub fn write_pnps_plot(results: &[GenePnPs], prefix: &str) -> anyhow::Result<String> {
    use std::fmt::Write as FmtWrite;
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Color constants — passed as format arguments so they never appear
    // literally inside r#"..."# delimiters (which `"#` would terminate).
    const C_GRID: &str = "#e0e0e0";
    const C_POSITIVE: &str = "#d94a4a";
    const C_PURIFYING: &str = "#4a90d9";
    const C_AXIS: &str = "#333333";

    let plot_path = format!("{}_pnps_manhattan.svg", prefix);

    // Filter out genes with no SNPs or NaN/infinite pN/pS
    let plot_data: Vec<&GenePnPs> = results
        .iter()
        .filter(|r| r.total_snps > 0.0 && r.pn_ps.is_finite())
        .collect();

    if plot_data.is_empty() {
        info!("No valid data points for pN/pS plot");
        return Ok(plot_path);
    }

    let width = 900.0f64;
    let height = 500.0f64;
    let margin_top = 50.0f64;
    let margin_right = 40.0f64;
    let margin_bottom = 80.0f64;
    let margin_left = 80.0f64;
    let plot_w = width - margin_left - margin_right;
    let plot_h = height - margin_top - margin_bottom;

    let max_pos = plot_data.iter().map(|r| r.genome_start).max().unwrap_or(1) as f64;
    let min_pos = plot_data.iter().map(|r| r.genome_start).min().unwrap_or(0) as f64;
    let pos_range = (max_pos - min_pos).max(1.0);

    let max_y = plot_data
        .iter()
        .map(|r| r.pn_ps)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.5)
        * 1.1;

    let to_x = |pos: usize| margin_left + ((pos as f64 - min_pos) / pos_range) * plot_w;
    let to_y = |v: f64| margin_top + plot_h * (1.0 - v / max_y);

    let mut svg = String::with_capacity(4096);

    // Header
    svg.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" width=\"{w}\" height=\"{h}\">",
        w = width, h = height
    );
    let _ = writeln!(svg, "<style>");
    let _ = writeln!(svg, "  text {{ font-family: sans-serif; fill: {c}; }}", c = C_AXIS);
    let _ = writeln!(svg, "  .title {{ font-size: 16px; font-weight: bold; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .axis-label {{ font-size: 12px; text-anchor: middle; }}");
    let _ = writeln!(svg, "  .tick-label {{ font-size: 10px; }}");
    let _ = writeln!(svg, "</style>");
    let _ = writeln!(svg, "<rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>", w = width, h = height);
    let _ = writeln!(
        svg,
        "<text x=\"{cx}\" y=\"30\" class=\"title\">pN/pS per Gene (Manhattan Plot)</text>",
        cx = width / 2.0
    );

    // Y-axis grid lines and labels
    let num_y_ticks = 5;
    for i in 0..=num_y_ticks {
        let frac = i as f64 / num_y_ticks as f64;
        let val = max_y * frac;
        let y = to_y(val);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"0.5\"/>",
            margin_left, y, margin_left + plot_w, y, C_GRID
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x}\" y=\"{y:.1}\" class=\"tick-label\" text-anchor=\"end\" dominant-baseline=\"middle\">{val:.2}</text>",
            x = margin_left - 8.0, y = y, val = val
        );
    }

    // Neutral line at pN/pS = 1.0
    if 1.0 <= max_y {
        let y1 = to_y(1.0);
        let _ = writeln!(
            svg,
            "<line x1=\"{}\" y1=\"{:.1}\" x2=\"{}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1\" stroke-dasharray=\"6,3\"/>",
            margin_left, y1, margin_left + plot_w, y1, C_POSITIVE
        );
        let _ = writeln!(
            svg,
            "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"tick-label\" fill=\"{c}\">pN/pS = 1</text>",
            x = margin_left + plot_w + 3.0, y = y1 + 3.0, c = C_POSITIVE
        );
    }

    // Data points
    for r in &plot_data {
        let x = to_x(r.genome_start);
        let y = to_y(r.pn_ps);
        let color = if r.pn_ps < 1.0 { C_PURIFYING } else { C_POSITIVE };
        let radius = r.total_snps.sqrt().clamp(2.0, 8.0);
        let _ = writeln!(
            svg,
            "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{r:.1}\" fill=\"{c}\" opacity=\"0.7\">",
            x = x, y = y, r = radius, c = color
        );
        let _ = writeln!(
            svg,
            "  <title>{name}: pN/pS={ratio:.4} ({s}S/{n}N SNPs)</title>",
            name = r.name, ratio = r.pn_ps, s = r.syn_snps, n = r.nonsyn_snps
        );
        svg.push_str("</circle>\n");
    }

    // Axes
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top, margin_left, margin_top + plot_h, C_AXIS
    );
    let _ = writeln!(
        svg,
        "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        margin_left, margin_top + plot_h, margin_left + plot_w, margin_top + plot_h, C_AXIS
    );

    // Axis labels
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"{y:.1}\" class=\"axis-label\">Genome Position</text>",
        x = margin_left + plot_w / 2.0, y = height - 10.0
    );
    let _ = writeln!(
        svg,
        "<text x=\"15\" y=\"{y:.1}\" class=\"axis-label\" transform=\"rotate(-90,15,{y:.1})\">pN/pS</text>",
        y = margin_top + plot_h / 2.0
    );

    // Legend
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"40\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_PURIFYING
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"50\" class=\"tick-label\">Purifying (pN/pS &lt; 1)</text>",
        x = width - 184.0
    );
    let _ = writeln!(
        svg,
        "<rect x=\"{x:.1}\" y=\"56\" width=\"12\" height=\"12\" fill=\"{c}\"/>",
        x = width - 200.0, c = C_POSITIVE
    );
    let _ = writeln!(
        svg,
        "<text x=\"{x:.1}\" y=\"66\" class=\"tick-label\">Positive (pN/pS &ge; 1)</text>",
        x = width - 184.0
    );

    svg.push_str("</svg>\n");

    let mut file = BufWriter::new(File::create(&plot_path)?);
    file.write_all(svg.as_bytes())?;

    Ok(plot_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_gc() -> &'static GeneticCode {
        crate::genetic_code::get_table(1).unwrap()
    }

    #[test]
    fn test_codon_to_aa() {
        let gc = make_gc();
        // ATG = Met
        assert_eq!(codon_to_aa(&[b'A', b'T', b'G'], gc), Some(b'M'));
        // TAA = Stop
        assert_eq!(codon_to_aa(&[b'T', b'A', b'A'], gc), Some(b'*'));
        // GCT = Ala
        assert_eq!(codon_to_aa(&[b'G', b'C', b'T'], gc), Some(b'A'));
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"AATTCC"), b"GGAATT");
    }

    #[test]
    fn test_count_sites() {
        let gc = make_gc();
        // ATG GCT = Met Ala (2 codons)
        let cds = b"ATGGCT";
        let (n, s) = count_sites(cds, gc);
        // Each codon contributes exactly 3 sites total (S + N = 3 per codon, 6 total)
        let total = n + s;
        assert!((total - 6.0).abs() < 0.01, "Total sites should be 6.0, got {}", total);
        // ATG (Met): 0 syn, 8 nonsyn (1 change → stop TAA excluded) → N=3*8/8=3, S=0
        // GCT (Ala): 3 syn (GCA,GCC,GCG), 6 nonsyn, 0 stops → N=3*6/9=2, S=3*3/9=1
        assert!(n > 4.5 && n < 5.5, "N_sites: expected ~5.0, got {}", n);
        assert!(s > 0.5 && s < 1.5, "S_sites: expected ~1.0, got {}", s);
    }

    #[test]
    fn test_genomic_to_cds_offset_plus_strand() {
        let gene = Gene {
            name: "test".to_string(),
            seqid: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![crate::gff::CdsExon {
                seqid: "chr1".to_string(),
                start: 100,
                end: 109,
                strand: Strand::Plus,
                phase: 0,
            }],
            length_bp: 10,
            start: 100,
        };

        assert_eq!(genomic_to_cds_offset(&gene, 100), Some(0));
        assert_eq!(genomic_to_cds_offset(&gene, 105), Some(5));
        assert_eq!(genomic_to_cds_offset(&gene, 109), Some(9));
        assert_eq!(genomic_to_cds_offset(&gene, 110), None);
        assert_eq!(genomic_to_cds_offset(&gene, 99), None);
    }

    #[test]
    fn test_genomic_to_cds_offset_minus_strand() {
        let gene = Gene {
            name: "test".to_string(),
            seqid: "chr1".to_string(),
            strand: Strand::Minus,
            exons: vec![crate::gff::CdsExon {
                seqid: "chr1".to_string(),
                start: 100,
                end: 109,
                strand: Strand::Minus,
                phase: 0,
            }],
            length_bp: 10,
            start: 100,
        };

        // Minus strand: position 109 maps to CDS offset 0
        assert_eq!(genomic_to_cds_offset(&gene, 109), Some(0));
        assert_eq!(genomic_to_cds_offset(&gene, 100), Some(9));
    }
}
