use crate::codon::extract_group_key;
use crate::models::{DsDn, Model, Z_95_CONFIDENCE};
use crossbeam::channel::unbounded;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

/// Writes pairwise results using a dedicated writer thread.
/// Computes pairs on-the-fly with lazy per-row caching (O(U) memory per thread).
pub fn write_pairwise(
    ids: &[String],
    uidx_by_id: &[usize],
    n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
    model: Model,
    sep: char,
    ext: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
    let pb_write = ProgressBar::new(total_pairs_to_write as u64);
    pb_write.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Computing & writing pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let nan_count = AtomicUsize::new(0);
    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = format!("{}_pairwise_results.{}", output_prefix, ext);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            let header = match model {
                Model::Li => format!("Seq1{s}Seq2{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
                Model::Nei => format!("Seq1{s}Seq2{s}dN{s}dS{s}dN/dS\n", s = sep),
            };
            out_file.write_all(header.as_bytes())?;
            for line_batch in rx {
                out_file.write_all(line_batch.as_bytes())?;
            }
            Ok(())
        }
    });

    (0..ids.len()).into_par_iter().for_each_init(
        || {
            (
                tx.clone(),
                vec![DsDn { dn: 0.0, ds: 0.0 }; n_u],
                vec![0u32; n_u],
                0u32,
                String::with_capacity(1024 * 64),
            )
        },
        |(sender, row_cache, gen_map, cur_gen, local_buffer), i| {
            *cur_gen = cur_gen.wrapping_add(1);
            let gen = *cur_gen;
            let u_i = uidx_by_id[i];
            row_cache[u_i] = DsDn { dn: 0.0, ds: 0.0 };
            gen_map[u_i] = gen;

            for j in (i + 1)..ids.len() {
                let u_j = uidx_by_id[j];
                if gen_map[u_j] != gen {
                    row_cache[u_j] = compute_pair(u_i, u_j);
                    gen_map[u_j] = gen;
                }
                let result = row_cache[u_j];
                if !result.dn.is_finite() || !result.ds.is_finite() {
                    nan_count.fetch_add(1, Ordering::Relaxed);
                }
                let ratio = if result.ds == 0.0 {
                    if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    result.dn / result.ds
                };

                let _ = writeln!(local_buffer, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                    &ids[i], &ids[j], result.dn, result.ds, ratio, s = sep);

                if local_buffer.len() > 1024 * 32 {
                    sender.send(std::mem::take(local_buffer))
                        .expect("Writer thread channel closed unexpectedly");
                }
            }
            if !local_buffer.is_empty() {
                sender.send(std::mem::take(local_buffer))
                    .expect("Writer thread channel closed unexpectedly");
            }
            pb_write.inc((ids.len() - 1 - i) as u64);
        },
    );
    drop(tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb_write.finish_with_message("Pairwise computation & writing completed.");

    let nan_total = nan_count.load(Ordering::Relaxed);
    if nan_total > 0 {
        info!("{} of {} pairs ({:.1}%) returned NaN due to saturation.",
            nan_total, total_pairs_to_write,
            100.0 * nan_total as f64 / total_pairs_to_write as f64);
    }
    Ok(())
}

/// Writes dN/dS summary by lineage using a dedicated writer thread.
/// Computes pairs on-the-fly with lazy per-row caching.
pub fn write_lineage(
    ids: &[String],
    uidx_by_id: &[usize],
    n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
    lineage_indices: &[usize],
    lineage_names: &[String],
    sep: char,
    ext: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let num_lineages = lineage_names.len();
    let output_path = format!("{}_lineage_summary.{}", output_prefix, ext);

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = output_path.clone();
        let header = format!("Genome{s}Against_Lineage{s}Mean_dN{s}Mean_dS{s}dN/dS_Ratio\n", s = sep);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            out_file.write_all(header.as_bytes())?;
            for block in rx {
                out_file.write_all(block.as_bytes())?;
            }
            Ok(())
        }
    });

    (0..ids.len())
        .into_par_iter()
        .for_each_init(
            || {
                (
                    tx.clone(),
                    vec![DsDn { dn: 0.0, ds: 0.0 }; n_u],
                    vec![0u32; n_u],
                    0u32,
                    vec![(0.0f64, 0.0f64, 0usize); num_lineages],
                )
            },
            |(sender, row_cache, gen_map, cur_gen, local_aggr), i| {
                *cur_gen = cur_gen.wrapping_add(1);
                let gen = *cur_gen;
                let u_i = uidx_by_id[i];
                row_cache[u_i] = DsDn { dn: 0.0, ds: 0.0 };
                gen_map[u_i] = gen;

                for a in local_aggr.iter_mut() { *a = (0.0, 0.0, 0); }

                for j in 0..ids.len() {
                    if i == j { continue; }
                    let u_j = uidx_by_id[j];
                    if gen_map[u_j] != gen {
                        row_cache[u_j] = compute_pair(u_i, u_j);
                        gen_map[u_j] = gen;
                    }
                    let result = row_cache[u_j];

                    if result.dn.is_finite() && result.ds.is_finite() {
                        let lin_idx = lineage_indices[j];
                        local_aggr[lin_idx].0 += result.dn;
                        local_aggr[lin_idx].1 += result.ds;
                        local_aggr[lin_idx].2 += 1;
                    }
                }

                let mut block = String::new();
                for (lin_idx, &(sum_dn, sum_ds, count)) in local_aggr.iter().enumerate() {
                    if count == 0 { continue; }
                    let mean_dn = sum_dn / count as f64;
                    let mean_ds = sum_ds / count as f64;
                    let ratio = if mean_ds == 0.0 {
                        if mean_dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        mean_dn / mean_ds
                    };
                    let _ = writeln!(block, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        &ids[i], &lineage_names[lin_idx], mean_dn, mean_ds, ratio, s = sep);
                }
                if !block.is_empty() {
                    sender.send(block).expect("Writer thread channel closed unexpectedly");
                }
            },
        );
    drop(tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    Ok(())
}

/// Writes grouped dN/dS averages using a dedicated writer thread.
/// Computes pairs on-the-fly (no pre-stored results needed).
pub fn write_group_average(
    ids: &[String],
    uidx_by_id: &[usize],
    _n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    output_prefix: &str,
    first_letter_lineage: bool,
    sep: char,
    ext: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut group_map: rustc_hash::FxHashMap<&str, usize> = rustc_hash::FxHashMap::default();
    let mut group_names: Vec<String> = Vec::new();
    let group_by_id: Vec<usize> = ids.iter().map(|id| {
        let key = extract_group_key(id, first_letter_lineage);
        let next_idx = group_names.len();
        *group_map.entry(key).or_insert_with(|| {
            group_names.push(key.to_string());
            next_idx
        })
    }).collect();

    let num_groups = group_names.len();
    let mut group_members: Vec<Vec<usize>> = vec![Vec::new(); num_groups];
    for (id_idx, &grp) in group_by_id.iter().enumerate() {
        group_members[grp].push(id_idx);
    }

    let mut group_pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..num_groups {
        for j in i..num_groups {
            group_pairs.push((i, j));
        }
    }

    let pb_group = ProgressBar::new(group_pairs.len() as u64);
    pb_group.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Groups: {pos}/{len} ({eta})")
            .progress_chars("#>-"),
    );

    let output_path = format!("{}_group_avg_dn_ds.{}", output_prefix, ext);

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = output_path.clone();
        let header = format!("Group1{s}Group2{s}NumSeqs1{s}NumSeqs2{s}NumComparisons{s}Mean_dN/dS{s}StdError{s}95%CI\n", s = sep);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            out_file.write_all(header.as_bytes())?;
            for line in rx {
                out_file.write_all(line.as_bytes())?;
            }
            Ok(())
        }
    });

    group_pairs.into_par_iter().for_each_with(tx, |s, (g1, g2)| {
        let members1 = &group_members[g1];
        let members2 = &group_members[g2];
        let mut pair_dn_ds_ratios = Vec::new();

        for (pos, &id1_idx) in members1.iter().enumerate() {
            let u1 = uidx_by_id[id1_idx];
            let start = if g1 == g2 { pos + 1 } else { 0 };
            for &id2_idx in &members2[start..] {
                let u2 = uidx_by_id[id2_idx];
                let result = compute_pair(u1, u2);
                if result.dn.is_finite() && result.ds.is_finite() {
                    let ratio = if result.ds == 0.0 {
                        if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        result.dn / result.ds
                    };
                    pair_dn_ds_ratios.push(ratio);
                }
            }
        }

        let line = if pair_dn_ds_ratios.is_empty() {
            format!("{}{s}{}{s}{}{s}{}{s}0{s}NaN{s}NaN{s}[NaN, NaN]\n",
                &group_names[g1], &group_names[g2],
                members1.len(), members2.len(), s = sep)
        } else {
            let n = pair_dn_ds_ratios.len();
            let mean: f64 = pair_dn_ds_ratios.iter().sum::<f64>() / n as f64;
            if n == 1 {
                format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}N/A{s}N/A\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n, mean, s = sep)
            } else {
                let variance: f64 =
                    pair_dn_ds_ratios.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
                let se = (variance / n as f64).sqrt();
                let ci_hw = Z_95_CONFIDENCE * se;
                format!("{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}[{:.6}, {:.6}]\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n,
                    mean, se, mean - ci_hw, mean + ci_hw, s = sep)
            }
        };

        s.send(line).expect("Writer thread channel closed unexpectedly");
        pb_group.inc(1);
    });

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb_group.finish_with_message("Group average computation completed.");
    Ok(())
}

/// Writes pairwise sliding window results using a dedicated writer thread.
pub fn write_pairwise_windows(
    ids: &[String],
    uidx_by_id: &[usize],
    unique_codon_indices: &[Vec<u8>],
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
    output_prefix: &str,
    model: Model,
    sep: char,
    ext: &str,
    window_size: usize,
    window_step: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let seq_len = unique_codon_indices.first().map(|v| v.len()).unwrap_or(0);
    let num_windows = if seq_len >= window_size {
        (seq_len - window_size) / window_step + 1
    } else {
        0
    };
    let total_pairs = ids.len() * (ids.len() - 1) / 2;
    let pb = ProgressBar::new((total_pairs * num_windows) as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Window pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let nan_count = AtomicUsize::new(0);
    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = format!("{}_pairwise_windows.{}", output_prefix, ext);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            let header = match model {
                Model::Li => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
                Model::Nei => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN{s}dS{s}dN/dS\n", s = sep),
            };
            out_file.write_all(header.as_bytes())?;
            for batch in rx {
                out_file.write_all(batch.as_bytes())?;
            }
            Ok(())
        }
    });

    (0..ids.len()).into_par_iter().for_each_init(
        || {
            (
                tx.clone(),
                String::with_capacity(1024 * 64),
            )
        },
        |(sender, local_buffer), i| {
            let u_i = uidx_by_id[i];
            for j in (i + 1)..ids.len() {
                let u_j = uidx_by_id[j];
                let s1 = &unique_codon_indices[u_i];
                let s2 = &unique_codon_indices[u_j];
                for w in 0..num_windows {
                    let start = w * window_step;
                    let end = start + window_size;
                    let result = if u_i == u_j {
                        DsDn { dn: 0.0, ds: 0.0 }
                    } else {
                        compute_pair_slices(&s1[start..end], &s2[start..end])
                    };
                    if !result.dn.is_finite() || !result.ds.is_finite() {
                        nan_count.fetch_add(1, Ordering::Relaxed);
                    }
                    let ratio = if result.ds == 0.0 {
                        if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        result.dn / result.ds
                    };
                    // Window positions in 1-based codon coordinates
                    let _ = writeln!(local_buffer, "{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        &ids[i], &ids[j], start + 1, end, result.dn, result.ds, ratio, s = sep);

                    if local_buffer.len() > 1024 * 32 {
                        sender.send(std::mem::take(local_buffer))
                            .expect("Writer thread channel closed unexpectedly");
                    }
                }
                pb.inc(num_windows as u64);
            }
            if !local_buffer.is_empty() {
                sender.send(std::mem::take(local_buffer))
                    .expect("Writer thread channel closed unexpectedly");
            }
        },
    );
    drop(tx);

    writer_thread.join()
        .expect("Writer thread panicked")
        .expect("Writer thread encountered an I/O error");
    pb.finish_with_message("Sliding window computation completed.");

    let nan_total = nan_count.load(Ordering::Relaxed);
    let total_computations = total_pairs * num_windows;
    if nan_total > 0 {
        info!("{} of {} window pairs ({:.1}%) returned NaN due to saturation.",
            nan_total, total_computations,
            100.0 * nan_total as f64 / total_computations as f64);
    }
    Ok(())
}
