mod codon;
mod genetic_code;
mod models;
mod output;
mod plot;
mod stats;

use clap::Parser;
use log::{info, warn};
use needletail::parse_fastx_file;

use models::{DsDn, Model, OutputFormat};
use stats::SummaryStats;

/// Calculates dN/dS for sequences using Nei-Gojobori or Li (1993) models.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input file with aligned sequences in FASTA format
    #[arg(required_unless_present = "list_codes")]
    input_file: Option<String>,

    /// Base name for output files
    #[arg(short, long, default_value = "output")]
    output: String,

    /// Number of parallel threads
    #[arg(short, long, default_value_t = 4)]
    workers: usize,

    /// Compute mean dN and dS grouped by lineage against all others
    #[arg(long, group = "output_mode")]
    lineage: bool,

    /// Compute mean dN/dS between predefined groups
    #[arg(long, group = "output_mode")]
    group_average: bool,

    /// Group by the first letter of the sequence ID instead of splitting on '_'
    /// (requires --lineage or --group-average)
    #[arg(long, requires = "output_mode")]
    first_letter_lineage: bool,

    /// Model to use for calculation
    #[arg(long, value_enum, default_value_t = Model::Nei)]
    model: Model,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    format: OutputFormat,

    /// Minimum valid codons per sequence (sequences below this are filtered out)
    #[arg(long, default_value_t = 0)]
    min_codons: usize,

    /// Sliding window size in codons (pairwise mode only)
    #[arg(long)]
    window_size: Option<usize>,

    /// Sliding window step in codons (default: 1)
    #[arg(long, default_value_t = 1)]
    window_step: usize,

    /// Print statistical summary to stderr after computation
    #[arg(long)]
    summary: bool,

    /// Generate SVG plot file
    #[arg(long)]
    plot: bool,

    /// NCBI genetic code table number (default: 1 = Standard).
    /// Use --list-codes to see all available tables.
    /// Common alternatives: 2 (Vertebrate Mito), 4 (Mycoplasma), 11 (Bacterial)
    #[arg(long, default_value_t = 1)]
    genetic_code: u8,

    /// List all available genetic code tables and exit
    #[arg(long)]
    list_codes: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();

    // Handle --list-codes
    if args.list_codes {
        eprintln!("Available NCBI genetic code tables:");
        for (id, name) in genetic_code::list_tables() {
            eprintln!("  {:>2}  {}", id, name);
        }
        return Ok(());
    }

    // Validate genetic code
    let gc = genetic_code::get_table(args.genetic_code)
        .unwrap_or_else(|| {
            eprintln!("Error: unknown genetic code table {}. Use --list-codes to see available tables.", args.genetic_code);
            std::process::exit(1);
        });
    if args.genetic_code != 1 {
        info!("Using genetic code table {}: {}", gc.id, gc.name);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    // --- Read sequences: convert FASTA bytes directly to codon indices ---
    let input_file = args.input_file.as_ref().unwrap_or_else(|| {
        eprintln!("Error: no input file specified.");
        std::process::exit(1);
    });

    info!("Reading sequences from: {}", input_file);
    let mut all_codon_indices: Vec<Vec<u8>> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    let mut gap_count_total = 0usize;
    let mut seqs_with_gaps = 0usize;
    {
        let mut reader = parse_fastx_file(input_file)?;
        while let Some(record) = reader.next() {
            let rec = record?;
            let seq = rec.seq();
            let gaps = codon::count_gaps(&seq);
            if gaps > 0 {
                seqs_with_gaps += 1;
                gap_count_total += gaps;
            }
            if seq.len() % 3 != 0 {
                warn!("Sequence '{}' length {} is not divisible by 3; {} trailing base(s) ignored.",
                    String::from_utf8_lossy(rec.id()), seq.len(), seq.len() % 3);
            }
            let codons = codon::fasta_to_codon_indices(&seq, args.model);
            let id_str = String::from_utf8_lossy(rec.id()).into_owned();
            all_codon_indices.push(codons);
            ids.push(id_str);
        }
    }
    if ids.is_empty() {
        eprintln!("No sequences found in the input file.");
        std::process::exit(1);
    }
    info!("Found {} total sequences.", ids.len());

    // --- Diagnostics: gaps and alignment ---
    if seqs_with_gaps > 0 {
        warn!("{} sequence(s) contain {} total gap character(s) ('-' or '.'), treated as ambiguous.",
            seqs_with_gaps, gap_count_total);
    }
    if let Some(first_len) = all_codon_indices.first().map(|v| v.len()) {
        let mismatched = all_codon_indices.iter().filter(|v| v.len() != first_len).count();
        if mismatched > 0 {
            warn!("{} sequence(s) have different codon lengths than the first ({} codons). Sequences may not be aligned.",
                mismatched, first_len);
        }
        info!("Alignment length: {} codons ({} bp).", first_len, first_len * 3);
    }

    // --- Filter by --min-codons ---
    if args.min_codons > 0 {
        let before = ids.len();
        let mut keep = Vec::with_capacity(ids.len());
        for i in 0..ids.len() {
            let valid = all_codon_indices[i].iter().filter(|&&c| c != codon::INVALID_CODON).count();
            if valid >= args.min_codons {
                keep.push(i);
            }
        }
        if keep.len() < before {
            let filtered = before - keep.len();
            warn!("Filtered {} sequence(s) with fewer than {} valid codons.", filtered, args.min_codons);
            let new_ids: Vec<String> = keep.iter().map(|&i| std::mem::take(&mut ids[i])).collect();
            let new_codons: Vec<Vec<u8>> = keep.iter().map(|&i| std::mem::take(&mut all_codon_indices[i])).collect();
            ids = new_ids;
            all_codon_indices = new_codons;
        }
        if ids.is_empty() {
            eprintln!("All sequences filtered out by --min-codons {}.", args.min_codons);
            std::process::exit(1);
        }
    }

    // --- Sort-based deduplication on codon indices ---
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
    let n_u = unique_repr_indices.len();

    let mut uidx_by_id: Vec<usize> = vec![0; ids.len()];
    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        uidx_by_id[orig_idx] = uidx_by_sorted[pos];
    }
    drop(sorted_indices);
    drop(uidx_by_sorted);
    info!("Found {} unique sequences.", n_u);

    // --- Extract only unique codon indices, drop the rest ---
    let unique_codon_indices: Vec<Vec<u8>> = {
        let mut unique_vecs: Vec<Vec<u8>> = Vec::with_capacity(n_u);
        for &orig_idx in &unique_repr_indices {
            unique_vecs.push(std::mem::take(&mut all_codon_indices[orig_idx]));
        }
        drop(all_codon_indices);
        unique_vecs
    };

    // --- Precomputation for model lookup tables ---
    let li_tables = if args.model == Model::Li {
        info!("Precomputing lookup tables for Li (1993) model...");
        let tables = models::li::LiTables::with_genetic_code(&gc.aa_table);
        info!("Li precomputation finished.");
        Some(tables)
    } else {
        None
    };

    let nei_tables = if args.model == Model::Nei {
        info!("Precomputing lookup tables for Nei-Gojobori (1986) model...");
        let tables = models::nei::NeiTables::with_genetic_code(gc);
        info!("Nei precomputation finished.");
        Some(tables)
    } else {
        None
    };

    // --- Summary stats (only allocated when --summary or --plot is active) ---
    let summary_stats = if args.summary || args.plot {
        Some(SummaryStats::new())
    } else {
        None
    };

    let compute_pair = |u_i: usize, u_j: usize| -> DsDn {
        if u_i == u_j {
            return DsDn { dn: 0.0, ds: 0.0 };
        }
        let (dn, ds) = match (&li_tables, &nei_tables) {
            (Some(tables), _) => tables.compute_pair(&unique_codon_indices[u_i], &unique_codon_indices[u_j]),
            (_, Some(tables)) => tables.compute_pair(&unique_codon_indices[u_i], &unique_codon_indices[u_j]),
            _ => unreachable!("one of li_tables or nei_tables must be Some"),
        };
        DsDn { dn, ds }
    };

    let compute_pair_slices = |s1: &[u8], s2: &[u8]| -> DsDn {
        let (dn, ds) = match (&li_tables, &nei_tables) {
            (Some(tables), _) => tables.compute_pair(s1, s2),
            (_, Some(tables)) => tables.compute_pair(s1, s2),
            _ => unreachable!("one of li_tables or nei_tables must be Some"),
        };
        DsDn { dn, ds }
    };

    let ext = args.format.extension();
    let sep = args.format.separator();
    let stats_ref = summary_stats.as_ref();

    if args.group_average {
        if args.window_size.is_some() {
            eprintln!("--window-size cannot be used with --group-average.");
            std::process::exit(1);
        }
        info!("Computing group average dN/dS...");
        let plot_data = output::write_group_average(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, args.first_letter_lineage, sep, ext, stats_ref)?;
        info!("Results saved to {}_group_avg_dn_ds.{}", args.output, ext);

        if args.plot && !plot_data.is_empty() {
            let plot_path = format!("{}_group_dnds.svg", args.output);
            plot::group_bar_svg(&plot_data, &plot_path)?;
            info!("Plot saved to {}", plot_path);
        }
    } else if args.lineage {
        if args.window_size.is_some() {
            eprintln!("--window-size cannot be used with --lineage.");
            std::process::exit(1);
        }
        info!("Computing dN/dS lineage summary...");
        let mut lineage_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
        let mut lineage_names: Vec<String> = Vec::new();
        let lineage_indices: Vec<usize> = ids.iter().map(|id| {
            let key = codon::extract_group_key(id, args.first_letter_lineage);
            let next_idx = lineage_names.len();
            *lineage_map.entry(key).or_insert_with(|| {
                lineage_names.push(key.to_string());
                next_idx
            })
        }).collect();
        let plot_data = output::write_lineage(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, &lineage_indices, &lineage_names, sep, ext, stats_ref)?;
        info!("Lineage summary saved to {}_lineage_summary.{}", args.output, ext);

        if args.plot && !plot_data.is_empty() {
            let lineage_plot_data: Vec<plot::LineagePlotData> = plot_data.into_iter()
                .map(|(genome, lineage, ratio)| plot::LineagePlotData { genome, lineage, ratio })
                .collect();
            let plot_path = format!("{}_lineage_dnds.svg", args.output);
            plot::lineage_bar_svg(&lineage_plot_data, &plot_path)?;
            info!("Plot saved to {}", plot_path);
        }
    } else if let Some(win_size) = args.window_size {
        let seq_len = unique_codon_indices.first().map(|v| v.len()).unwrap_or(0);
        if seq_len == 0 {
            eprintln!("Cannot use --window-size with empty sequences.");
            std::process::exit(1);
        }
        // Window mode requires uniform sequence lengths; misaligned sequences
        // would cause out-of-bounds panics when slicing window ranges.
        let misaligned = unique_codon_indices.iter().any(|v| v.len() != seq_len);
        if misaligned {
            eprintln!("--window-size requires all sequences to have equal length. Sequences are not aligned.");
            std::process::exit(1);
        }
        if win_size == 0 || win_size > seq_len {
            eprintln!("--window-size must be between 1 and {} (sequence length in codons).", seq_len);
            std::process::exit(1);
        }
        if args.window_step == 0 {
            eprintln!("--window-step must be at least 1.");
            std::process::exit(1);
        }
        let num_windows = (seq_len - win_size) / args.window_step + 1;
        let window_stats = if args.plot {
            Some(stats::WindowStats::new(num_windows))
        } else {
            None
        };
        info!("Generating sliding window pairwise results (window={}, step={})...", win_size, args.window_step);
        output::write_pairwise_windows(&ids, &uidx_by_id, &unique_codon_indices, &compute_pair_slices, &args.output, args.model, sep, ext, win_size, args.window_step, stats_ref, window_stats.as_ref())?;
        info!("Results saved to {}_pairwise_windows.{}", args.output, ext);

        if let Some(ws) = &window_stats {
            let plot_path = format!("{}_window_plot.svg", args.output);
            plot::window_plot_svg(ws, &plot_path, win_size, args.window_step)?;
            info!("Plot saved to {}", plot_path);
        }
    } else {
        info!("Generating pairwise results...");
        output::write_pairwise(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, args.model, sep, ext, stats_ref)?;
        info!("Results saved to {}_pairwise_results.{}", args.output, ext);

        if args.plot {
            if let Some(ref stats) = summary_stats {
                let plot_path = format!("{}_dnds_histogram.svg", args.output);
                plot::histogram_svg(stats, &plot_path)?;
                info!("Plot saved to {}", plot_path);
            }
        }
    }

    // --- Print summary if requested ---
    if let Some(ref stats) = summary_stats {
        if args.summary {
            stats.print_summary();
        }
    }

    info!("All processes completed successfully.");
    Ok(())
}
