//! Command-line argument parsing and validation.

use clap::{Parser, Subcommand};

use crate::models::{Model, OutputFormat};

/// Calculates dN/dS for sequences using Nei-Gojobori or Li (1993) models.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// List all available genetic code tables and exit
    #[arg(long)]
    pub list_codes: bool,

    #[command(subcommand)]
    pub command: Option<SubCmd>,
}

#[derive(Subcommand, Debug)]
pub enum SubCmd {
    /// Compute pairwise dN/dS from codon-aligned FASTA sequences
    Fasta(FastaArgs),
    /// Compute pN/pS per gene from a VCF file, reference FASTA, and GFF3 annotation
    Vcf(VcfArgs),
}

/// Arguments for the FASTA subcommand (original behavior).
#[derive(Parser, Debug)]
pub struct FastaArgs {
    /// Input file with aligned sequences in FASTA format
    pub input_file: String,

    /// Base name for output files
    #[arg(short, long, default_value = "output")]
    pub output: String,

    /// Number of parallel threads
    #[arg(short, long, default_value_t = 4)]
    pub workers: usize,

    /// Compute mean dN and dS grouped by lineage against all others
    #[arg(long, group = "output_mode")]
    pub lineage: bool,

    /// Compute mean dN/dS between predefined groups
    #[arg(long, group = "output_mode")]
    pub group_average: bool,

    /// Group by the first letter of the sequence ID instead of splitting on '_'
    /// (requires --lineage or --group-average)
    #[arg(long, requires = "output_mode")]
    pub first_letter_lineage: bool,

    /// Model to use for calculation
    #[arg(long, value_enum, default_value_t = Model::Nei)]
    pub model: Model,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    pub format: OutputFormat,

    /// Minimum valid codons per sequence (sequences below this are filtered out)
    #[arg(long, default_value_t = 0)]
    pub min_codons: usize,

    /// Sliding window size in codons (pairwise mode only)
    #[arg(long)]
    pub window_size: Option<usize>,

    /// Sliding window step in codons (default: 1)
    #[arg(long, default_value_t = 1)]
    pub window_step: usize,

    /// Print statistical summary to stderr after computation
    #[arg(long)]
    pub summary: bool,

    /// Generate SVG plot file
    #[arg(long)]
    pub plot: bool,

    /// NCBI genetic code table number (default: 1 = Standard).
    /// Use --list-codes to see all available tables.
    /// Common alternatives: 2 (Vertebrate Mito), 4 (Mycoplasma), 11 (Bacterial)
    #[arg(long, default_value_t = 1)]
    pub genetic_code: u8,
}

/// Arguments for the VCF subcommand.
#[derive(Parser, Debug)]
pub struct VcfArgs {
    /// Reference FASTA file
    #[arg(long = "ref")]
    pub reference: String,

    /// GFF3 annotation file
    #[arg(long)]
    pub gff: String,

    /// VCF file(s) with variants. One per sample.
    /// Use multiple times: --vcf s1.vcf --vcf s2.vcf
    /// Or provide a file with one VCF path per line: --vcf-list samples.txt
    #[arg(long)]
    pub vcf: Vec<String>,

    /// File containing one VCF path per line (alternative to multiple --vcf flags)
    #[arg(long)]
    pub vcf_list: Option<String>,

    /// Base name for output files
    #[arg(short, long, default_value = "output")]
    pub output: String,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    pub format: OutputFormat,

    /// NCBI genetic code table number (default: 1 = Standard).
    /// Use --list-codes to see all available tables.
    /// Common alternatives: 2 (Vertebrate Mito), 4 (Mycoplasma), 11 (Bacterial)
    #[arg(long, default_value_t = 1)]
    pub genetic_code: u8,

    /// Only include variants with FILTER=PASS (or '.')
    #[arg(long)]
    pub pass_only: bool,

    /// Minimum allele frequency threshold (0.0-1.0)
    #[arg(long)]
    pub min_af: Option<f64>,

    /// Maximum allele frequency threshold (0.0-1.0).
    /// Use --max-af 0.99 to exclude fixed variants (AF=1.0) and keep only
    /// segregating polymorphisms for pN/pS analysis.
    #[arg(long)]
    pub max_af: Option<f64>,

    /// Minimum read depth (from INFO/DP)
    #[arg(long)]
    pub min_depth: Option<u32>,

    /// Weight SNP counts by allele frequency instead of counting each SNP as 1.
    /// With this flag, a SNP at AF=0.3 contributes 0.3 to the syn/nonsyn count
    /// instead of 1.0. This computes πN/πS (nucleotide diversity ratio) rather
    /// than simple pN/pS.
    #[arg(long)]
    pub af_weighted: bool,

    /// Generate SVG Manhattan-style plot of pN/pS per gene
    #[arg(long)]
    pub plot: bool,
}
