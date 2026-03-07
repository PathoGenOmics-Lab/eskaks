mod codon;
mod models;
mod output;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use needletail::parse_fastx_file;
use rayon::prelude::*;
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

    // --- Lectura con deduplicacion basada en sort ---
    info!("Leyendo secuencias desde: {}", args.input_file);
    let mut entries: Vec<(Vec<u8>, String)> = Vec::new(); // (dna5_seq, id_string)
    let mut reader = parse_fastx_file(&args.input_file)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq_dna5 = seq_to_dna5(&rec.seq());
        let id_str = String::from_utf8_lossy(rec.id()).into_owned();
        entries.push((seq_dna5, id_str));
    }
    if entries.is_empty() {
        eprintln!("No se encontraron secuencias en el archivo de entrada.");
        std::process::exit(1);
    }
    info!("Se encontraron {} secuencias en total.", entries.len());

    // Ordenar por secuencia para deduplicar sin HashMap
    let original_order: Vec<usize> = (0..entries.len()).collect();
    let mut sorted_indices: Vec<usize> = original_order.clone();
    sorted_indices.sort_unstable_by(|&a, &b| entries[a].0.cmp(&entries[b].0));

    // Asignar indices unicos agrupando secuencias iguales
    let mut uidx_by_sorted = Vec::with_capacity(entries.len());
    let mut unique_repr_indices: Vec<usize> = Vec::new(); // indice del primer representante de cada grupo
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

    // Mapear indice original -> uidx
    let mut uidx_by_id: Vec<usize> = vec![0; entries.len()];
    for (pos, &orig_idx) in sorted_indices.iter().enumerate() {
        uidx_by_id[orig_idx] = uidx_by_sorted[pos];
    }

    // Extraer IDs como Vec<String> (en orden original)
    let ids: Vec<String> = entries.iter().map(|(_, id)| id.clone()).collect();
    info!("Se encontraron {} secuencias unicas.", n_u);

    // --- Codificacion de secuencias unicas ---
    info!("Codificando secuencias unicas a indices de codones...");
    let unique_codon_indices: Arc<Vec<Vec<u8>>> = Arc::new(
        unique_repr_indices.par_iter()
            .map(|&orig_idx| dna5_to_codon_indices(&entries[orig_idx].0, args.model))
            .collect()
    );
    drop(entries);

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
        output::write_group_average(&ids, &uidx_by_id, &get_result, &args.output)?;
        info!("Resultados guardados en {}_group_avg_dn_ds.tsv", args.output);
    } else if args.lineage {
        info!("Calculando resumen de dN/dS por linaje...");
        // Pre-computar indices de lineage para evitar allocaciones en el hot loop
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
        output::write_lineage(&ids, &uidx_by_id, &get_result, &args.output, &lineage_indices, &lineage_names)?;
        info!("Resumen de linaje guardado en {}_lineage_summary.tsv", args.output);
    } else {
        info!("Generando informe de resultados por pares...");
        output::write_pairwise(&ids, &uidx_by_id, &get_result, &args.output, args.model)?;
        info!("Resultados guardados en {}_pairwise_results.tsv", args.output);
    }

    info!("Todos los procesos han finalizado con exito.");
    Ok(())
}
