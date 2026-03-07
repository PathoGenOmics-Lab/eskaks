use crate::models::{DsDn, Model, Z_95_CONFIDENCE};
use crossbeam::channel::unbounded;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::thread;

/// Escribe resultados pairwise completos usando un hilo escritor dedicado.
pub fn write_pairwise(
    ids: &[Vec<u8>],
    id_to_uidx: &FxHashMap<Vec<u8>, usize>,
    get_result: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
    model: Model,
) -> Result<(), Box<dyn std::error::Error>> {
    let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
    let pb_write = ProgressBar::new(total_pairs_to_write as u64);
    pb_write.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Escribiendo pares: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = format!("{}_pairwise_results.tsv", output_prefix);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(File::create(output_path)?);
            let header = match model {
                Model::Li => "Seq1\tSeq2\tdN(Ka)\tdS(Ks)\tdN/dS\n",
                Model::Nei => "Seq1\tSeq2\tdN\tdS\tdN/dS\n",
            };
            out_file.write_all(header.as_bytes())?;
            for line_batch in rx {
                out_file.write_all(line_batch.as_bytes())?;
            }
            Ok(())
        }
    });

    (0..ids.len()).into_par_iter().for_each_with(tx, |s, i| {
        let mut local_buffer = String::with_capacity(1024 * 8);
        let id_i_bytes = &ids[i];
        let u_i = id_to_uidx[id_i_bytes];

        for j in (i + 1)..ids.len() {
            let id_j_bytes = &ids[j];
            let u_j = id_to_uidx[id_j_bytes];

            let result = get_result(u_i, u_j);
            let ratio = if result.ds == 0.0 {
                if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                result.dn / result.ds
            };

            let _ = write!(local_buffer, "{}\t{}\t{:.6}\t{:.6}\t{:.6}\n",
                String::from_utf8_lossy(id_i_bytes),
                String::from_utf8_lossy(id_j_bytes),
                result.dn, result.ds, ratio);

            if local_buffer.len() > 1024 * 4 {
                s.send(std::mem::take(&mut local_buffer))
                    .expect("Writer thread channel closed unexpectedly");
            }
        }
        if !local_buffer.is_empty() {
            s.send(local_buffer).expect("Writer thread channel closed unexpectedly");
        }
        pb_write.inc((ids.len() - 1 - i) as u64);
    });

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb_write.finish_with_message("Escritura de pares completada.");
    Ok(())
}

/// Escribe resumen de dN/dS por linaje.
pub fn write_lineage(
    ids: &[Vec<u8>],
    id_to_uidx: &FxHashMap<Vec<u8>, usize>,
    get_result: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
    first_letter_lineage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let lineage_results: Vec<String> = (0..ids.len())
        .into_par_iter()
        .map(|i| {
            let id_i_bytes = &ids[i];
            let u_i = id_to_uidx[id_i_bytes];
            let mut local_aggr: FxHashMap<Vec<u8>, (f64, f64, usize)> = FxHashMap::default();

            for j in 0..ids.len() {
                if i == j { continue; }
                let id_j_bytes = &ids[j];
                let u_j = id_to_uidx[id_j_bytes];
                let result = get_result(u_i, u_j);

                if result.dn.is_finite() && result.ds.is_finite() {
                    let lineage_key = if first_letter_lineage {
                        id_j_bytes.iter().next().map(|&b| vec![b]).unwrap_or_default()
                    } else {
                        id_j_bytes.split(|&b| b == b'_').next().unwrap_or(id_j_bytes).to_vec()
                    };
                    let entry = local_aggr.entry(lineage_key).or_default();
                    entry.0 += result.dn;
                    entry.1 += result.ds;
                    entry.2 += 1;
                }
            }

            let mut output_lines = String::new();
            for (lineage_key, (sum_dn, sum_ds, count)) in local_aggr {
                let mean_dn = sum_dn / count as f64;
                let mean_ds = sum_ds / count as f64;
                let ratio = if mean_ds == 0.0 {
                    if mean_dn == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    mean_dn / mean_ds
                };
                let _ = write!(output_lines, "{}\t{}\t{:.6}\t{:.6}\t{:.6}\n",
                    String::from_utf8_lossy(id_i_bytes),
                    String::from_utf8_lossy(&lineage_key),
                    mean_dn, mean_ds, ratio);
            }
            output_lines
        })
        .collect();

    let output_path = format!("{}_lineage_summary.tsv", output_prefix);
    let mut out_file = BufWriter::new(File::create(&output_path)?);
    writeln!(out_file, "Genome\tAgainst_Lineage\tMean_dN\tMean_dS\tdN/dS_Ratio")?;
    for block in lineage_results {
        out_file.write_all(block.as_bytes())?;
    }
    Ok(())
}

/// Escribe promedios de dN/dS agrupados.
pub fn write_group_average(
    ids: &[Vec<u8>],
    id_to_uidx: &FxHashMap<Vec<u8>, usize>,
    get_result: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut group_map: FxHashMap<Vec<u8>, Vec<Vec<u8>>> = FxHashMap::default();
    for id_bytes in ids {
        let group_key = id_bytes.split(|&b| b == b'_').next().unwrap_or(id_bytes).to_vec();
        group_map.entry(group_key).or_default().push(id_bytes.clone());
    }

    let group_keys: Vec<Vec<u8>> = group_map.keys().cloned().collect();
    let mut group_pairs_to_process = Vec::new();
    for (i, g1_key) in group_keys.iter().enumerate() {
        for g2_key in group_keys.iter().skip(i) {
            group_pairs_to_process.push((g1_key.clone(), g2_key.clone()));
        }
    }

    let pb_group = ProgressBar::new(group_pairs_to_process.len() as u64);
    pb_group.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Grupos: {pos}/{len} ({eta})")
            .progress_chars("#>-"),
    );

    let group_results: Vec<String> = group_pairs_to_process
        .par_iter()
        .map(|(g1_key, g2_key)| {
            let ids1 = &group_map[g1_key];
            let ids2 = &group_map[g2_key];
            let mut pair_dn_ds_ratios = Vec::new();

            for (i, id1) in ids1.iter().enumerate() {
                let u_idx1 = id_to_uidx[id1];
                let iter_ids2: Box<dyn Iterator<Item = &Vec<u8>>> = if g1_key == g2_key {
                    Box::new(ids2.iter().skip(i + 1))
                } else {
                    Box::new(ids2.iter())
                };
                for id2 in iter_ids2 {
                    let u_idx2 = id_to_uidx[id2];
                    let result = get_result(u_idx1, u_idx2);
                    if result.dn.is_finite() && result.ds.is_finite() && result.ds > 0.0 {
                        pair_dn_ds_ratios.push(result.dn / result.ds);
                    }
                }
            }

            pb_group.inc(1);
            if pair_dn_ds_ratios.is_empty() {
                format!("{}\t{}\t{}\t{}\t0\tNaN\tNaN\t[NaN, NaN]\n",
                    String::from_utf8_lossy(g1_key), String::from_utf8_lossy(g2_key),
                    ids1.len(), ids2.len())
            } else {
                let num_comparisons = pair_dn_ds_ratios.len();
                let mean_dn_ds: f64 = pair_dn_ds_ratios.iter().sum::<f64>() / num_comparisons as f64;
                let variance: f64 = if num_comparisons > 1 {
                    pair_dn_ds_ratios.iter().map(|&val| (val - mean_dn_ds).powi(2)).sum::<f64>() / (num_comparisons - 1) as f64
                } else { 0.0 };
                let standard_error = if num_comparisons > 0 {
                    (variance / num_comparisons as f64).sqrt()
                } else { 0.0 };
                let ci_half_width = Z_95_CONFIDENCE * standard_error;
                let ci_lower = mean_dn_ds - ci_half_width;
                let ci_upper = mean_dn_ds + ci_half_width;
                format!("{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t[{:.6}, {:.6}]\n",
                    String::from_utf8_lossy(g1_key), String::from_utf8_lossy(g2_key),
                    ids1.len(), ids2.len(), num_comparisons,
                    mean_dn_ds, standard_error, ci_lower, ci_upper)
            }
        }).collect();

    let output_path = format!("{}_group_avg_dn_ds.tsv", output_prefix);
    let mut out_file = BufWriter::new(File::create(&output_path)?);
    writeln!(out_file, "Group1\tGroup2\tNumSeqs1\tNumSeqs2\tNumComparisons\tMean_dN/dS\tStdError\t95%CI")?;
    for line in group_results {
        if !line.is_empty() { out_file.write_all(line.as_bytes())?; }
    }
    pb_group.finish_with_message("Calculo de promedios de grupo completado.");
    Ok(())
}
