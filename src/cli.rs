//! Command-line argument parsing and validation.

use clap::Parser;

use crate::models::{Model, OutputFormat};

/// Calculates dN/dS for sequences using Nei-Gojobori or Li (1993) models.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Input file with aligned sequences in FASTA format
    #[arg(required_unless_present = "list_codes")]
    pub input_file: Option<String>,

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

    /// List all available genetic code tables and exit
    #[arg(long)]
    pub list_codes: bool,
}
