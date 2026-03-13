mod codon;
mod models;
mod output;

use clap::Parser;
use log::info;
use needletail::parse_fastx_file;

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

    /// In lineage mode, group by the first letter of the sequence ID
    #[arg(long, requires = "lineage")]
    first_letter_lineage: bool,

    /// Model to use for calculation
    #[arg(long, value_enum, default_value_t = Model::Nei)]
    model: Model,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    // --- Read sequences: convert FASTA bytes directly to codon indices ---
    // We store codon indices (L/3 bytes each) instead of raw DNA5 (L bytes each),
    // reducing memory by 3x and avoiding the intermediate DNA5 allocation.
    info!("Reading sequences from: {}", args.input_file);
    let mut all_codon_indices: Vec<Vec<u8>> = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    {
        let mut reader = parse_fastx_file(&args.input_file)?;
        while let Some(record) = reader.next() {
            let rec = record?;
            let codons = codon::fasta_to_codon_indices(&rec.seq(), args.model);
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

    // --- Sort-based deduplication on codon indices (L/3 bytes, not L bytes) ---
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
        let tables = models::li::LiTables::new();
        info!("Li precomputation finished.");
        Some(tables)
    } else {
        None
    };

    let nei_tables = if args.model == Model::Nei {
        info!("Precomputing lookup tables for Nei-Gojobori (1986) model...");
        let tables = models::nei::NeiTables::new();
        info!("Nei precomputation finished.");
        Some(tables)
    } else {
        None
    };

    // --- Output generation (streaming: compute pairs on-the-fly with per-row caching) ---
    // Instead of precomputing and storing all U*(U-1)/2 pair results (O(U²) memory),
    // each output function computes pairs on demand using lazy per-row caching.
    // Memory: O(U) per thread instead of O(U²) total.
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

    if args.group_average {
        info!("Computing group average dN/dS...");
        output::write_group_average(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, args.first_letter_lineage)?;
        info!("Results saved to {}_group_avg_dn_ds.tsv", args.output);
    } else if args.lineage {
        info!("Computing dN/dS lineage summary...");
        let mut lineage_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
        let mut lineage_names: Vec<String> = Vec::new();
        let lineage_indices: Vec<usize> = ids.iter().map(|id| {
            let key = if args.first_letter_lineage {
                &id[..id.chars().next().map(|c| c.len_utf8()).unwrap_or(0)]
            } else {
                id.split('_').next().unwrap_or(id)
            };
            let next_idx = lineage_names.len();
            *lineage_map.entry(key).or_insert_with(|| {
                lineage_names.push(key.to_string());
                next_idx
            })
        }).collect();
        output::write_lineage(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, &lineage_indices, &lineage_names)?;
        info!("Lineage summary saved to {}_lineage_summary.tsv", args.output);
    } else {
        info!("Generating pairwise results...");
        output::write_pairwise(&ids, &uidx_by_id, n_u, &compute_pair, &args.output, args.model)?;
        info!("Results saved to {}_pairwise_results.tsv", args.output);
    }

    info!("All processes completed successfully.");
    Ok(())
}