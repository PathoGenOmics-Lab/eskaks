mod codon;
mod models;
mod output;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::Arc;

use codon::{dna5_to_codon_indices, seq_to_dna5};
use models::{DsDn, Model};

/// Calcula dN/dS para secuencias usando los modelos Nei-Gojobori o Li (1993).
#[derive(Parser, Debug)]
#[command(version = "1.0.0", about, long_about = None)]
struct Args {
    /// Archivo de entrada con secuencias en formato FASTA
    #[arg(required = true)]
    input_file: String,

    /// Nombre base para los archivos de salida
    #[arg(short, long, default_value = "output")]
    output: String,

    /// Numero de procesos paralelos
    #[arg(short, long, default_value_t = 4)]
    workers: usize,

    /// Calcula la media de dN y dS agrupada por linaje contra todos los demas
    #[arg(long, group = "output_mode")]
    lineage: bool,

    /// Calcula la media de dN/dS entre grupos predefinidos
    #[arg(long, group = "output_mode")]
    group_average: bool,

    /// En el modo linaje, agrupa por la primera letra del ID
    #[arg(long, requires = "lineage")]
    first_letter_lineage: bool,

    /// Modelo a utilizar
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

    // --- Lectura y deduplicacion ---
    info!("Leyendo y deduplicando secuencias desde: {}", args.input_file);
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut seq_to_ids: FxHashMap<Vec<u8>, Vec<Vec<u8>>> = FxHashMap::default();
    let mut reader = parse_fastx_file(&args.input_file)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq_dna5 = seq_to_dna5(&rec.seq());
        let id = rec.id().to_vec();
        ids.push(id.clone());
        seq_to_ids.entry(seq_dna5).or_default().push(id);
    }
    info!("Se encontraron {} secuencias en total.", ids.len());
    if ids.is_empty() {
        eprintln!("No se encontraron secuencias en el archivo de entrada.");
        std::process::exit(1);
    }

    let mut unique_seqs_dna5: Vec<Vec<u8>> = seq_to_ids.keys().cloned().collect();
    unique_seqs_dna5.sort_unstable();
    let n_u = unique_seqs_dna5.len();
    info!("Se encontraron {} secuencias unicas.", n_u);

    // --- Mapeo de IDs a indices unicos ---
    let id_to_uidx: FxHashMap<Vec<u8>, usize> = unique_seqs_dna5.iter().enumerate()
        .flat_map(|(u_idx, seq_dna5)| {
            seq_to_ids.get(seq_dna5)
                .expect("BUG: unique sequence not found in seq_to_ids map")
                .iter()
                .map(move |id| (id.clone(), u_idx))
        })
        .collect();

    // --- Codificacion de secuencias unicas ---
    info!("Codificando secuencias unicas a indices de codones...");
    let unique_codon_indices: Arc<Vec<Vec<u8>>> = Arc::new(
        unique_seqs_dna5.par_iter()
            .map(|s_dna5| dna5_to_codon_indices(s_dna5, args.model))
            .collect()
    );
    drop(unique_seqs_dna5);
    drop(seq_to_ids);

    // --- Precomputacion para modelo Li ---
    let li_tables = if args.model == Model::Li {
        info!("Precalculando tablas para el modelo Li (1993)...");
        let tables = models::li::LiTables::new();
        info!("Precalculacion para Li finalizada.");
        Some(tables)
    } else {
        None
    };

    // --- Calculo en paralelo de pares unicos ---
    let total_unique_pairs = if n_u > 0 { n_u * (n_u - 1) / 2 } else { 0 };
    info!("Calculando dN/dS para {} pares de secuencias unicas...", total_unique_pairs);
    let mut pair_results: Vec<DsDn> = vec![DsDn { dn: f64::NAN, ds: f64::NAN }; total_unique_pairs];

    let pb_calc = ProgressBar::new(n_u as u64);
    pb_calc.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Calculando pares unicos: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    // Particionar en slices mutables no superpuestos (uno por fila)
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
    pb_calc.finish_with_message("Calculo de pares unicos completado.");

    // --- Generacion de salida ---
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
        info!("Calculando promedios de grupo dN/dS...");
        output::write_group_average(&ids, &id_to_uidx, &get_result, &args.output)?;
        info!("Resultados guardados en {}_group_avg_dn_ds.tsv", args.output);
    } else if args.lineage {
        info!("Calculando resumen de dN/dS por linaje...");
        output::write_lineage(&ids, &id_to_uidx, &get_result, &args.output, args.first_letter_lineage)?;
        info!("Resumen de linaje guardado en {}_lineage_summary.tsv", args.output);
    } else {
        info!("Generando informe de resultados por pares...");
        output::write_pairwise(&ids, &id_to_uidx, &get_result, &args.output, args.model)?;
        info!("Resultados guardados en {}_pairwise_results.tsv", args.output);
    }

    info!("Todos los procesos han finalizado con exito.");
    Ok(())
}
