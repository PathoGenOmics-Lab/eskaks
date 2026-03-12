use crate::codon::extract_group_key;
use crate::models::{DsDn, Model, Z_95_CONFIDENCE};
use crossbeam::channel::unbounded;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
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
) -> Result<(), Box<dyn std::error::Error>> {
    let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
    let pb_write = ProgressBar::new(total_pairs_to_write as u64);
    pb_write.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Computing & writing pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = format!("{}_pairwise_results.tsv", output_prefix);
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
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
        let u_i = uidx_by_id[i];
        // Lazy per-row cache: only compute unique pairs as needed
        let mut row_cache: Vec<DsDn> = vec![DsDn { dn: 0.0, ds: 0.0 }; n_u];
        let mut computed: Vec<bool> = vec![false; n_u];
        computed[u_i] = true; // self-pair is (0.0, 0.0) by default

        let mut local_buffer = String::with_capacity(1024 * 8);

        for j in (i + 1)..ids.len() {
            let u_j = uidx_by_id[j];
            if !computed[u_j] {
                row_cache[u_j] = compute_pair(u_i, u_j);
                computed[u_j] = true;
            }
            let result = row_cache[u_j];
            let ratio = if result.ds == 0.0 {
                if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                result.dn / result.ds
            };

            let _ = writeln!(local_buffer, "{}\t{}\t{:.6}\t{:.6}\t{:.6}",
                &ids[i], &ids[j], result.dn, result.ds, ratio);

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
    pb_write.finish_with_message("Pairwise computation & writing completed.");
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
) -> Result<(), Box<dyn std::error::Error>> {
    let num_lineages = lineage_names.len();
    let output_path = format!("{}_lineage_summary.tsv", output_prefix);

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = output_path.clone();
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            out_file.write_all(b"Genome\tAgainst_Lineage\tMean_dN\tMean_dS\tdN/dS_Ratio\n")?;
            for block in rx {
                out_file.write_all(block.as_bytes())?;
            }
            Ok(())
        }
    });

    (0..ids.len())
        .into_par_iter()
        .for_each_with(tx, |s, i| {
            let u_i = uidx_by_id[i];
            // Lazy per-row cache
            let mut row_cache: Vec<DsDn> = vec![DsDn { dn: 0.0, ds: 0.0 }; n_u];
            let mut computed: Vec<bool> = vec![false; n_u];
            computed[u_i] = true;

            let mut local_aggr: Vec<(f64, f64, usize)> = vec![(0.0, 0.0, 0); num_lineages];

            for j in 0..ids.len() {
                if i == j { continue; }
                let u_j = uidx_by_id[j];
                if !computed[u_j] {
                    row_cache[u_j] = compute_pair(u_i, u_j);
                    computed[u_j] = true;
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
                let _ = writeln!(block, "{}\t{}\t{:.6}\t{:.6}\t{:.6}",
                    &ids[i], &lineage_names[lin_idx], mean_dn, mean_ds, ratio);
            }
            if !block.is_empty() {
                s.send(block).expect("Writer thread channel closed unexpectedly");
            }
        });

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

    let output_path = format!("{}_group_avg_dn_ds.tsv", output_prefix);

    let (tx, rx) = unbounded::<String>();

    let writer_thread = thread::spawn({
        let output_path = output_path.clone();
        move || -> Result<(), std::io::Error> {
            let mut out_file = BufWriter::new(
                File::create(&output_path)
                    .map_err(|e| std::io::Error::new(e.kind(), format!("Cannot create '{}': {}", output_path, e)))?
            );
            writeln!(out_file, "Group1\tGroup2\tNumSeqs1\tNumSeqs2\tNumComparisons\tMean_dN/dS\tStdError\t95%CI")?;
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
            format!("{}\t{}\t{}\t{}\t0\tNaN\tNaN\t[NaN, NaN]\n",
                &group_names[g1], &group_names[g2],
                members1.len(), members2.len())
        } else {
            let n = pair_dn_ds_ratios.len();
            let mean: f64 = pair_dn_ds_ratios.iter().sum::<f64>() / n as f64;
            if n == 1 {
                format!("{}\t{}\t{}\t{}\t{}\t{:.6}\tN/A\tN/A\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n, mean)
            } else {
                let variance: f64 =
                    pair_dn_ds_ratios.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
                let se = (variance / n as f64).sqrt();
                let ci_hw = Z_95_CONFIDENCE * se;
                format!("{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t[{:.6}, {:.6}]\n",
                    &group_names[g1], &group_names[g2],
                    members1.len(), members2.len(), n,
                    mean, se, mean - ci_hw, mean + ci_hw)
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
