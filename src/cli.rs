//! Command-line argument parsing and validation.

use clap::{Parser, Subcommand};

use crate::models::{Model, OutputFormat};

const FASTA_AFTER_HELP: &str = "\
EXAMPLES:
  # Pairwise dN/dS from a codon alignment  (writes results_pairwise_results.tsv)
  eskaks fasta alignment.fasta -o results

  # Per-lineage means plus an interactive HTML report
  eskaks fasta alignment.fasta --lineage --report -o results

  # Li/LPB93 model, a sliding window, and a printed summary
  eskaks fasta alignment.fasta --model li --window-size 30 --summary

  # Bacterial genetic code, bootstrap CIs, JSON output
  eskaks fasta alignment.fasta --genetic-code 11 --bootstrap 1000 --format json

  IDs are split on '_' for --lineage / --group-average (or use --first-letter-lineage).";

const VCF_AFTER_HELP: &str = "\
EXAMPLES:
  # Per-gene pN/pS from a VCF, reference, and annotation
  eskaks vcf --ref genome.fa --gff genes.gff3 --vcf sample.vcf -o results

  # M. tuberculosis-style: bacterial code, ts/tv weighting, clonal correction, report
  eskaks vcf --ref ref.fa --gff genes.gff3 --vcf calls.vcf \\
      --genetic-code 11 --kappa 2 --genomic-control --exclude-repetitive --report

  # Many samples (one VCF path per line) plus a McDonald-Kreitman test
  eskaks vcf --ref ref.fa --gff genes.gff3 --vcf-list samples.txt --mk -o pop";

const TOP_AFTER_HELP: &str = "\
EXAMPLES:
  eskaks fasta alignment.fasta -o results          # pairwise dN/dS from a codon alignment
  eskaks vcf --ref g.fa --gff g.gff3 --vcf v.vcf   # per-gene pN/pS from variants

Run `eskaks help fasta` or `eskaks help vcf` for all options and more examples.
Docs & guide: https://pathogenomics-lab.github.io/eskaks/";

/// Calculates dN/dS for sequences using Nei-Gojobori or Li (1993) models.
#[derive(Parser, Debug)]
#[command(version = env!("ESKAKS_VERSION"), about, long_about = None, after_help = TOP_AFTER_HELP)]
pub struct Args {
    /// List all available genetic code tables and exit
    #[arg(long)]
    pub list_codes: bool,

    /// Run a quick demo on bundled example data (no input files needed) and exit
    #[arg(long)]
    pub demo: bool,

    /// Generate a shell completion script for the given shell and exit
    /// (e.g. `eskaks --completions bash > /etc/bash_completion.d/eskaks`).
    #[arg(long, value_enum, value_name = "SHELL")]
    pub completions: Option<clap_complete::Shell>,

    /// Increase log verbosity: -v shows info, -vv shows debug.
    /// Data-quality warnings are shown by default; RUST_LOG overrides this.
    #[arg(short = 'v', long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Silence all logs except errors.
    #[arg(short = 'q', long, global = true)]
    pub quiet: bool,

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
#[command(after_help = FASTA_AFTER_HELP)]
pub struct FastaArgs {
    /// Input file with aligned sequences in FASTA format
    pub input_file: String,

    /// Number of parallel threads (0 = use all available cores)
    #[arg(short, long, default_value_t = 4)]
    pub workers: usize,

    /// Base name for output files
    #[arg(short, long, default_value = "output", help_heading = "Output options")]
    pub output: String,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv, help_heading = "Output options")]
    pub format: OutputFormat,

    /// Print statistical summary to stderr after computation
    #[arg(long, help_heading = "Output options")]
    pub summary: bool,

    /// Generate SVG plot file
    #[arg(long, help_heading = "Output options")]
    pub plot: bool,

    /// Write a self-contained interactive HTML report (<output>_report.html)
    /// with the dynamic visualizations for whichever mode was run (pairwise
    /// distribution, lineage/group scatter with per-group means). No CDN needed.
    #[arg(long, help_heading = "Output options")]
    pub report: bool,

    /// Model to use for calculation
    #[arg(long, value_enum, default_value_t = Model::Nei, help_heading = "Analysis options")]
    pub model: Model,

    /// NCBI genetic code table number (default: 1 = Standard).
    /// Use --list-codes to see all available tables.
    /// Common alternatives: 2 (Vertebrate Mito), 4 (Mycoplasma), 11 (Bacterial)
    #[arg(long, default_value_t = 1, help_heading = "Analysis options")]
    pub genetic_code: u8,

    /// Minimum valid codons per sequence (sequences below this are filtered out)
    #[arg(long, default_value_t = 0, help_heading = "Analysis options")]
    pub min_codons: usize,

    /// Compute mean dN and dS grouped by lineage against all others
    #[arg(long, group = "output_mode", help_heading = "Analysis options")]
    pub lineage: bool,

    /// Compute mean dN/dS between predefined groups
    #[arg(long, group = "output_mode", help_heading = "Analysis options")]
    pub group_average: bool,

    /// Group by the first letter of the sequence ID instead of splitting on '_'
    /// (requires --lineage or --group-average)
    #[arg(long, requires = "output_mode", help_heading = "Analysis options")]
    pub first_letter_lineage: bool,

    /// Sliding window size in codons (pairwise mode only)
    #[arg(long, help_heading = "Analysis options")]
    pub window_size: Option<usize>,

    /// Sliding window step in codons (default: 1)
    #[arg(long, default_value_t = 1, help_heading = "Analysis options")]
    pub window_step: usize,

    /// Also write a per-pair Nei-Gojobori neutrality test
    /// (<output>_pairwise_tests): dN, dS, standard errors, Z, and a two-sided
    /// p-value. Analytic variances are Nei-only (NaN for --model li).
    #[arg(long, help_heading = "Statistics options")]
    pub neutrality: bool,

    /// Bootstrap replicates for per-pair 95% CIs on dN, dS, and dN/dS
    /// (resamples codon columns; works for both models). Writes
    /// <output>_pairwise_bootstrap. 0 disables it.
    #[arg(long, default_value_t = 0, help_heading = "Statistics options")]
    pub bootstrap: usize,

    /// Seed for reproducible bootstrap resampling (only used with --bootstrap).
    #[arg(long, default_value_t = 42, help_heading = "Statistics options")]
    pub seed: u64,
}

/// Arguments for the VCF subcommand.
#[derive(Parser, Debug)]
#[command(after_help = VCF_AFTER_HELP, allow_negative_numbers = true)]
pub struct VcfArgs {
    /// Reference FASTA file
    #[arg(long = "ref", help_heading = "Input options")]
    pub reference: String,

    /// GFF3 annotation file
    #[arg(long, help_heading = "Input options")]
    pub gff: String,

    /// VCF file(s) with variants. One per sample.
    /// Use multiple times: --vcf s1.vcf --vcf s2.vcf
    /// Or provide a file with one VCF path per line: --vcf-list samples.txt
    #[arg(long, help_heading = "Input options")]
    pub vcf: Vec<String>,

    /// File containing one VCF path per line (alternative to multiple --vcf flags)
    #[arg(long, help_heading = "Input options")]
    pub vcf_list: Option<String>,

    /// Number of parallel threads for the per-gene pN/pS computation (0 = all cores)
    #[arg(short, long, default_value_t = 4)]
    pub workers: usize,

    /// Base name for output files
    #[arg(short, long, default_value = "output", help_heading = "Output options")]
    pub output: String,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv, help_heading = "Output options")]
    pub format: OutputFormat,

    /// Print the pN/pS summary even under --quiet. The summary is shown by default;
    /// this flag is accepted for symmetry with `eskaks fasta` and forces it when quiet.
    #[arg(long, help_heading = "Output options")]
    pub summary: bool,

    /// NCBI genetic code table number (default: 1 = Standard).
    /// Use --list-codes to see all available tables.
    /// Common alternatives: 2 (Vertebrate Mito), 4 (Mycoplasma), 11 (Bacterial)
    #[arg(long, default_value_t = 1, help_heading = "Analysis options")]
    pub genetic_code: u8,

    /// Only include variants with FILTER=PASS (or '.')
    #[arg(long, help_heading = "Filtering options")]
    pub pass_only: bool,

    /// Minimum allele frequency threshold (0.0-1.0)
    #[arg(long, help_heading = "Filtering options")]
    pub min_af: Option<f64>,

    /// Maximum allele frequency threshold (0.0-1.0).
    /// Use --max-af 0.99 to exclude fixed variants (AF=1.0) and keep only
    /// segregating polymorphisms for pN/pS analysis.
    #[arg(long, help_heading = "Filtering options")]
    pub max_af: Option<f64>,

    /// Minimum read depth (from INFO/DP)
    #[arg(long, help_heading = "Filtering options")]
    pub min_depth: Option<u32>,

    /// Drop genes with fewer than this many SNPs from the per-gene table,
    /// plot, and neutrality test (their pN/pS is unreliable). The genome-wide
    /// pooled estimate still uses all genes.
    #[arg(long, default_value_t = 0, help_heading = "Filtering options")]
    pub min_snps: usize,

    /// Exclude repetitive / hard-to-map genes (PE/PPE/PGRS, IS elements,
    /// maturases) from the genome-wide pooled estimate and the neutrality-test
    /// family, since their SNP calls are frequently mapping artefacts. They still
    /// appear in the per-gene table, flagged.
    #[arg(long, help_heading = "Filtering options")]
    pub exclude_repetitive: bool,

    /// Weight SNP counts by allele frequency instead of counting each SNP as 1.
    /// With this flag, a SNP at AF=0.3 contributes 0.3 to the syn/nonsyn count
    /// instead of 1.0. This computes πN/πS (nucleotide diversity ratio) rather
    /// than simple pN/pS.
    #[arg(long, help_heading = "Analysis options")]
    pub af_weighted: bool,

    /// Transition/transversion rate ratio (kappa) for mutation-spectrum-aware
    /// counting of synonymous vs nonsynonymous SITES. kappa > 1 up-weights
    /// transitions (A<->G, C<->T), correcting the equal-rates (Nei-Gojobori)
    /// bias that understates pN/pS under a transition-skewed spectrum such as
    /// M. tuberculosis. The default 1.0 is the classic equal-rates counting.
    #[arg(long, default_value_t = 1.0, help_heading = "Analysis options")]
    pub kappa: f64,

    /// Apply a genomic-control correction to the neutrality test: divide each
    /// gene's chi-square by the median-based inflation factor lambda (never below
    /// one) and report corrected p/q columns. Recommended for clonal organisms
    /// such as M. tuberculosis, where genome-wide linkage inflates significance.
    #[arg(long, help_heading = "Analysis options")]
    pub genomic_control: bool,

    /// False-discovery-rate threshold for calling genes significant in the
    /// per-gene neutrality test (used for the summary count and plot highlight).
    #[arg(long, default_value_t = 0.05, help_heading = "Statistics options")]
    pub fdr: f64,

    /// Bootstrap replicates for a 95% CI on the genome-wide pooled pN/pS
    /// (resamples genes with replacement). 0 disables it.
    #[arg(long, default_value_t = 0, help_heading = "Statistics options")]
    pub bootstrap: usize,

    /// Seed for the bootstrap resampling, for reproducible CIs (only used with --bootstrap).
    #[arg(long, default_value_t = 42, help_heading = "Statistics options")]
    pub seed: u64,

    /// Also run a per-gene McDonald-Kreitman test (writes <output>_mk.<ext>):
    /// 2x2 fixed/polymorphic table, Neutrality Index, alpha, and Fisher exact p.
    #[arg(long, help_heading = "Statistics options")]
    pub mk: bool,

    /// Write the per-variant table (<output>_variants.<ext>): one row per coding
    /// SNP with genomic position, base and amino-acid change (e.g. S315T), allele
    /// frequency, and effect (synonymous/missense/nonsense/stop_loss). The
    /// mutation-level key for joining to a resistance catalogue.
    #[arg(long, help_heading = "Output options")]
    pub variants: bool,

    /// Add Codon_Shared and Codon_Change columns to the --variants table. Codon_Shared is
    /// true when a sample carrying this ALT also carries another SNP in the same codon,
    /// so the row's amino-acid change is the JOINT (multi-nucleotide) one; false when no
    /// sample does; NA when the input has no per-sample genotypes to check, in which case
    /// the row was scored against the reference codon alone and may be wrong.
    /// Codon_Change is the codon change itself (CTT>TTA, coding strand). Opt-in so the
    /// default table keeps its column layout.
    #[arg(long, requires = "variants", help_heading = "Output options")]
    pub shared_codons: bool,

    /// Write per-gene population-diversity statistics (<output>_diversity.<ext>):
    /// per-site πN and πS (nucleotide diversity) and their ratio, Watterson's θ,
    /// and Tajima's D (the SFS neutrality test). Needs the sample size, so it
    /// requires a multi-sample VCF or several single-sample VCFs (--vcf-list).
    /// Assumes haploid genotypes (M. tuberculosis); only segregating sites count
    /// (variants fixed within the sample are excluded from π/θ/D).
    #[arg(long, help_heading = "Statistics options")]
    pub diversity: bool,

    /// Rank codons by how many DISTINCT nonsynonymous alleles they carry
    /// (writes <output>_codons.<ext>). Tests allelic multiplicity, not carrier
    /// counts and not a per-codon pN/pS: a codon has at most 9 possible alleles,
    /// so no per-codon pN/pS can reach even p = 0.05. Recurrence at one residue is
    /// a hallmark of drug resistance, but a mutational or mapping hotspot looks
    /// identical, so pair it with --exclude-repetitive.
    #[arg(long, help_heading = "Statistics options")]
    pub codon_scan: bool,

    /// Newick tree over the cohort's samples. Adds independent-origin counting to
    /// --codon-scan (columns Nonsyn_Origins, Max_Allele_Origins, Exp_Nonsyn_Origins,
    /// P_Origins, Q_Origins) and to --variants (column Origins), which is what lets a
    /// SINGLE allele that arose many times independently be seen as convergent. Tip
    /// names must match the VCF's sample names exactly (for one VCF per sample: the
    /// file's own sample column, else its basename); any mismatch is a hard error.
    ///
    /// Read the caveats. Origins are counted by Fitch parsimony, so: a genotyping
    /// error inflates them (mitigated, not removed, by --min-origin-support); missing
    /// data is invisible, since a carrier set cannot tell REF from a no-call, so heavy
    /// missingness inflates them too; a mutational or mapping hotspot is
    /// indistinguishable from selection, so pair this with --exclude-repetitive;
    /// RECOMBINATION reads as convergence, which is near-zero in M. tuberculosis but
    /// not in every organism; an unresolved tree inflates origins; and the counts are a
    /// property of the sampling frame, not of nature, so they are not portable between
    /// collections. The ranking is what is trustworthy, not the literal p-values.
    #[arg(long, value_name = "FILE", help_heading = "Statistics options")]
    pub tree: Option<String>,

    /// Minimum carriers in the clade an origin subtends for that origin to count
    /// (requires --tree). The default of 2 is deliberate and load-bearing: raw
    /// parsimony counts every isolated false call as its own origin, and three of them
    /// on top of one real clade allele are enough to fake genome-wide significance.
    /// Requiring two carriers per origin removes that failure mode and demotes calling
    /// artefacts below real signals. It biases conservative: a mutation that genuinely
    /// arose once in one sampled genome goes uncounted. Use 1 only for a cohort you
    /// trust completely.
    #[arg(long, default_value_t = 2, requires = "tree", help_heading = "Statistics options")]
    pub min_origin_support: u32,

    /// Allele frequency at/above which a variant is treated as "fixed"
    /// (divergence) rather than polymorphic in the McDonald-Kreitman test.
    #[arg(long, default_value_t = 0.99, help_heading = "Statistics options")]
    pub mk_fixed_af: f64,

    /// Generate SVG Manhattan-style plot of pN/pS per gene
    #[arg(long, help_heading = "Output options")]
    pub plot: bool,

    /// Write a self-contained interactive HTML report (<output>_report.html)
    /// with summary cards, an interactive Manhattan plot, and a sortable,
    /// filterable per-gene table. No internet/CDN needed.
    #[arg(long, help_heading = "Output options")]
    pub report: bool,

    /// Per-gene divergence dN/dS table for the report's polymorphism-vs-divergence
    /// panel (a 2-column TSV/CSV: gene<TAB>dN/dS; header optional). Genes are
    /// matched to the pN/pS results by name.
    #[arg(long, help_heading = "Output options")]
    pub divergence: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        // clap's own invariant checker: catches duplicate/conflicting args, bad
        // help_heading / group usage, and malformed defaults across the whole CLI
        // (top-level flags plus both subcommands).
        Args::command().debug_assert();
    }
}
