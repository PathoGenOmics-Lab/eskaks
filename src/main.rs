mod codon;
mod models;
mod output;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{error, info, warn};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::sync::Arc;

use codon::{dna5_to_codon_indices, extract_group_key, seq_to_dna5};
use models::{DsDn, Model};

/// Calculates dN/dS for sequences using Nei-Gojobori or Li (1993) models.
#[derive(Parser, Debug)]
#[command(version = "1.0.0", about, long_about = None)]
struct Args {
    /// Input file with aligned sequences in FASTA format
    #[arg(required = true)]
    input_file: String,

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
    #[arg(long)]
    first_letter_lineage: bool,

    /// Model to use for calculation
    #[arg(long, value_enum, default_value_t = Model::Nei)]
    model: Model,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();

    if args.first_letter_lineage && !args.lineage && !args.group_average {
        error!("--first_letter_lineage requires --lineage or --group_average.");
        std::process::exit(1);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    // --- Read sequences with sort-based deduplication ---
    info!("Reading sequences from: {}", args.input_file);
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} {msg}")
        .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "));
    spinner.set_message(format!("Reading sequences from {}...", args.input_file));
    spinner.enable_steady_tick(80);

    let mut entries: Vec<(Vec<u8>, String)> = Vec::new();
    let mut reader = parse_fastx_file(&args.input_file)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq_dna5 = seq_to_dna5(&rec.seq());
        let id_str = String::from_utf8_lossy(rec.id()).into_owned();
        entries.push((seq_dna5, id_str));
    }
    spinner.finish_and_clear();

    if entries.is_empty() {
        error!("No sequences found in the input file.");
        std::process::exit(1);
    }
    info!("Found {} total sequences.", entries.len());

    // Validate sequence lengths
    let first_len = entries[0].0.len();
    let all_same_length = entries.iter().all(|(seq, _)| seq.len() == first_len);
    if !all_same_length {
        error!("Not all sequences have the same length. Aligned sequences are required.");
        std::process::exit(1);
    }
    if !first_len.is_multiple_of(3) {
        warn!("Sequence length ({}) is not a multiple of 3. Trailing bases will be ignored.", first_len);
    }

    // Sort by sequence for deduplication without HashMap
    let mut sorted_indices: Vec<usize> = (0..entries.len()).collect();
    sorted_indices.sort_unstable_by(|&a, &b| entries[a].0.cmp(&entries[b].0));

    // Assign unique indices by grouping identical sequences
    let mut uidx_by_sorted = Vec::with_capacity(entries.len());
    let mut unique_repr_indices: Vec<usize> = Vec::new();
    let mut current_uidx = 0usize;
    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        if pos > 0 && entries[orig_idx].0 != entries[sorted_indices[pos - 1]].0 {
            current_uidx += 1;
        }
        if pos == 0 || entries[orig_idx].0 != entries[sorted_indices[pos - 1]].0 {
            unique_repr_indices.push(orig_idx);
        }
        uidx_by_sorted.push(current_uidx);
    }
    let n_u = unique_repr_indices.len();

    // Map original index -> unique index
    let mut uidx_by_id: Vec<usize> = vec![0; entries.len()];
    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        uidx_by_id[orig_idx] = uidx_by_sorted[pos];
    }

    let ids: Vec<String> = entries.iter().map(|(_, id)| id.clone()).collect();
    info!("Found {} unique sequences.", n_u);

    // --- Encode unique sequences to codon indices ---
    info!("Encoding unique sequences to codon indices...");
    let unique_codon_indices: Arc<Vec<Vec<u8>>> = Arc::new(
        unique_repr_indices.par_iter()
            .map(|&orig_idx| dna5_to_codon_indices(&entries[orig_idx].0, args.model))
            .collect()
    );
    drop(entries);

    // --- Precomputation for Li model ---
    let li_tables = if args.model == Model::Li {
        info!("Precomputing lookup tables for Li (1993) model...");
        let tables = models::li::LiTables::new();
        info!("Li precomputation finished.");
        Some(tables)
    } else {
        None
    };

    // --- Parallel computation of unique pairs ---
    let total_unique_pairs = if n_u > 0 { n_u * (n_u - 1) / 2 } else { 0 };
    info!("Computing dN/dS for {} unique sequence pairs...", total_unique_pairs);
    let mut pair_results: Vec<DsDn> = vec![DsDn { dn: f64::NAN, ds: f64::NAN }; total_unique_pairs];

    let pb_calc = ProgressBar::new(n_u as u64);
    pb_calc.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Computing unique pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let row_slices: Vec<&mut [DsDn]> = {
        let mut remaining = pair_results.as_mut_slice();
        let mut slices = Vec::with_capacity(n_u);
        for i in 0..n_u {
            let row_len = n_u - i - 1;
            let (row, rest) = remaining.split_at_mut(row_len);
            slices.push(row);
            remaining = rest;
        }
        slices
    };

    row_slices.into_par_iter().enumerate().for_each(|(i, row_slice)| {
        let cod_i = &unique_codon_indices[i];
        for (j_offset, j) in ((i + 1)..n_u).enumerate() {
            let cod_j = &unique_codon_indices[j];
            let (dn, ds) = match &li_tables {
                Some(tables) => tables.compute_pair(cod_i, cod_j),
                None => models::nei::calculate_syn_nonsyn_from_indices(cod_i, cod_j),
            };
            row_slice[j_offset] = DsDn { dn, ds };
        }
        pb_calc.inc(1);
    });
    pb_calc.finish_with_message("Unique pair computation completed.");

    // --- Output generation ---
    let get_result = |u_i: usize, u_j: usize| -> DsDn {
        if u_i == u_j {
            DsDn { dn: 0.0, ds: 0.0 }
        } else {
            let (i, j) = if u_i < u_j { (u_i, u_j) } else { (u_j, u_i) };
            let flat_idx = (i * n_u + j) - ((i + 1) * (i + 2)) / 2;
            pair_results[flat_idx]
        }
    };

    if args.group_average {
        info!("Computing group average dN/dS...");
        output::write_group_average(&ids, &uidx_by_id, get_result, &args.output, args.first_letter_lineage)?;
        info!("Results saved to {}_group_avg_dn_ds.tsv", args.output);
    } else if args.lineage {
        info!("Computing dN/dS lineage summary...");
        // Pre-compute lineage indices to avoid allocations in the hot loop
        let mut lineage_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
        let mut lineage_names: Vec<String> = Vec::new();
        let lineage_indices: Vec<usize> = ids.iter().map(|id| {
            let key = extract_group_key(id, args.first_letter_lineage);
            let next_idx = lineage_names.len();
            *lineage_map.entry(key).or_insert_with(|| {
                lineage_names.push(key.to_string());
                next_idx
            })
        }).collect();
        output::write_lineage(&ids, &uidx_by_id, get_result, &args.output, &lineage_indices, &lineage_names)?;
        info!("Lineage summary saved to {}_lineage_summary.tsv", args.output);
    } else {
        info!("Generating pairwise results...");
        output::write_pairwise(&ids, &uidx_by_id, get_result, &args.output, args.model)?;
        info!("Results saved to {}_pairwise_results.tsv", args.output);
    }

    info!("All processes completed successfully.");
    Ok(())
}
