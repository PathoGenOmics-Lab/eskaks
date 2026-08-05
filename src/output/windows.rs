//! Sliding-window pairwise table writer.

use super::*;

/// Writes pairwise sliding window results using a dedicated writer thread.
#[allow(clippy::too_many_arguments)]
pub fn write_pairwise_windows(
    ids: &[String],
    uidx_by_id: &[usize],
    unique_codon_indices: &[Vec<u8>],
    compute_pair_slices: impl Fn(&[u8], &[u8]) -> DsDn + Sync,
    window_size: usize,
    window_step: usize,
    window_stats: Option<&WindowStats>,
    cfg: &OutputConfig,
) -> anyhow::Result<()> {
    let output_prefix = cfg.prefix;
    let model = cfg.model;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
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
    let is_json = ext == "json";
    let header = match model {
        Model::Li => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
        Model::Nei => format!("Seq1{s}Seq2{s}Window_Start{s}Window_End{s}dN{s}dS{s}dN/dS\n", s = sep),
    };
    let path = format!("{}_pairwise_windows.{}", output_prefix, ext);
    let (tx, writer_thread) = if is_json {
        spawn_ordered_json_writer(path)?
    } else {
        spawn_ordered_writer(path, header)?
    };

    (0..ids.len()).into_par_iter().for_each_init(
        || {
            (
                tx.clone(),
                String::with_capacity(1024 * 64),
                FloatAccum::new(),
            )
        },
        |(sender, local_buffer, local_stats), i| {
            let u_i = uidx_by_id[i];
            for j in (i + 1)..ids.len() {
                let u_j = uidx_by_id[j];
                let s1 = &unique_codon_indices[u_i];
                let s2 = &unique_codon_indices[u_j];
                for w in 0..num_windows {
                    let start = w * window_step;
                    let end = start + window_size;
                    let result = if u_i == u_j {
                        // Identical sequences deduplicate to the same unique index. A
                        // window that has valid codons has zero divergence (0.0), but an
                        // all-N / all-gap window has no comparable codons and is undefined
                        // (NaN), not a spurious 0.0 that reads as strong purifying selection.
                        // This must be checked per window: even for identical sequences one
                        // window can be all-invalid while others are valid.
                        if s1[start..end].iter().all(|&c| c == crate::codon::INVALID_CODON) {
                            DsDn { dn: f64::NAN, ds: f64::NAN }
                        } else {
                            DsDn { dn: 0.0, ds: 0.0 }
                        }
                    } else {
                        compute_pair_slices(&s1[start..end], &s2[start..end])
                    };
                    if !result.dn.is_finite() || !result.ds.is_finite() {
                        nan_count.fetch_add(1, Ordering::Relaxed);
                    }
                    // An undefined dN or dS makes the ratio undefined too: a NaN numerator
                    // must not slip through the ds==0 branch and print as +inf.
                    let ratio = if result.dn.is_nan() || result.ds.is_nan() {
                        f64::NAN
                    } else if result.ds == 0.0 {
                        if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        result.dn / result.ds
                    };

                    if let Some(stats) = summary {
                        stats.record_pair_atomic(result.dn, result.ds, ratio);
                        local_stats.record(result.dn, result.ds, ratio);
                    }
                    if let Some(ws) = window_stats {
                        ws.record(w, ratio);
                    }

                    if is_json {
                        let _ = writeln!(local_buffer,
                            "{{\"seq1\":\"{}\",\"seq2\":\"{}\",\"window_start\":{},\"window_end\":{},\"dN\":{},\"dS\":{},\"dN_dS\":{}}}",
                            json_escape(&ids[i]), json_escape(&ids[j]), start + 1, end,
                            format_json_f64(result.dn), format_json_f64(result.ds), format_json_f64(ratio));
                    } else {
                        let _ = writeln!(local_buffer, "{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                            name_field(&ids[i], sep), name_field(&ids[j], sep), start + 1, end,
                            norm_zero(result.dn), norm_zero(result.ds), norm_zero(ratio), s = sep);
                    }
                }
                pb.inc(num_windows as u64);
            }

            if let Some(stats) = summary {
                stats.flush_local(local_stats);
                local_stats.reset();
            }

            // One block per row i, in order, so the output is deterministic.
            sender.send((i, std::mem::take(local_buffer)))
                .expect("Writer thread channel closed unexpectedly");
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

