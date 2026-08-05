//! Pairwise result / neutrality / bootstrap table writers.

use super::*;

/// Writes pairwise results using a dedicated writer thread.
/// Computes pairs on-the-fly with lazy per-row caching (O(U) memory per thread).
/// Write a per-pair Nei-Gojobori neutrality-test table:
/// Seq1, Seq2, dN, dS, SE_dN, SE_dS, Z, P_value, where Z = (dN-dS)/√(V(dN)+V(dS))
/// and P_value is the two-sided normal p. `stats(u_i, u_j)` returns
/// `(dn, ds, var_dn, var_ds)`. SE/Z/P are NaN for the Li model (no NG variance).
pub fn write_pairwise_tests(
    ids: &[String],
    uidx_by_id: &[usize],
    stats: impl Fn(usize, usize) -> (f64, f64, f64, f64) + Sync,
    prefix: &str,
    sep: char,
    ext: &str,
) -> anyhow::Result<String> {
    use rayon::prelude::*;

    let n = ids.len();
    let pairs: Vec<(usize, usize)> =
        (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let is_json = ext == "json";
    let fmt = |v: f64| if v.is_finite() { format!("{:.6}", norm_zero(v)) } else { "NaN".to_string() };
    let jfmt = |v: f64| if v.is_finite() { format!("{:.6}", norm_zero(v)) } else { "null".to_string() };

    let rows: Vec<String> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (dn, ds, var_dn, var_ds) = stats(uidx_by_id[i], uidx_by_id[j]);
            let var_sum = var_dn + var_ds;
            let z = if var_sum > 0.0 { (dn - ds) / var_sum.sqrt() } else { f64::NAN };
            let p = crate::stats::normal_two_sided_p(z);
            if is_json {
                format!(
                    "  {{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dS\":{},\"SE_dN\":{},\"SE_dS\":{},\"Z\":{},\"P_value\":{}}}",
                    json_escape(&ids[i]), json_escape(&ids[j]), jfmt(dn), jfmt(ds), jfmt(var_dn.sqrt()), jfmt(var_ds.sqrt()), jfmt(z), jfmt(p)
                )
            } else {
                format!(
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    name_field(&ids[i], sep), name_field(&ids[j], sep),
                    fmt(dn), fmt(ds), fmt(var_dn.sqrt()), fmt(var_ds.sqrt()), fmt(z), fmt(p), s = sep
                )
            }
        })
        .collect();

    let output_path = format!("{}_pairwise_tests.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);
    if is_json {
        writeln!(file, "[")?;
        writeln!(file, "{}", rows.join(",\n"))?;
        writeln!(file, "]")?;
    } else {
        writeln!(file, "Seq1{s}Seq2{s}dN{s}dS{s}SE_dN{s}SE_dS{s}Z{s}P_value", s = sep)?;
        for row in &rows {
            writeln!(file, "{}", row)?;
        }
    }
    Ok(output_path)
}

/// Write a per-pair bootstrap-CI table:
/// Seq1, Seq2, dN, dN_CI_low, dN_CI_high, dS, dS_CI_low, dS_CI_high,
/// dN/dS, dNdS_CI_low, dNdS_CI_high. `stats(u_i, u_j)` returns
/// `(dn, ds, dn_lo, dn_hi, ds_lo, ds_hi, ratio_lo, ratio_hi)`.
pub fn write_pairwise_bootstrap(
    ids: &[String],
    uidx_by_id: &[usize],
    stats: impl Fn(usize, usize) -> (f64, f64, f64, f64, f64, f64, f64, f64) + Sync,
    prefix: &str,
    sep: char,
    ext: &str,
) -> anyhow::Result<String> {
    use rayon::prelude::*;

    let n = ids.len();
    let pairs: Vec<(usize, usize)> =
        (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let is_json = ext == "json";
    let fmt = |v: f64| if v.is_finite() { format!("{:.6}", norm_zero(v)) } else { "NaN".to_string() };
    let jfmt = |v: f64| if v.is_finite() { format!("{:.6}", norm_zero(v)) } else { "null".to_string() };

    let rows: Vec<String> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (dn, ds, dn_lo, dn_hi, ds_lo, ds_hi, r_lo, r_hi) = stats(uidx_by_id[i], uidx_by_id[j]);
            let ratio = if ds > 0.0 { dn / ds } else { f64::NAN };
            if is_json {
                format!(
                    "  {{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dN_CI_low\":{},\"dN_CI_high\":{},\"dS\":{},\"dS_CI_low\":{},\"dS_CI_high\":{},\"dN_dS\":{},\"dNdS_CI_low\":{},\"dNdS_CI_high\":{}}}",
                    json_escape(&ids[i]), json_escape(&ids[j]), jfmt(dn), jfmt(dn_lo), jfmt(dn_hi), jfmt(ds), jfmt(ds_lo), jfmt(ds_hi), jfmt(ratio), jfmt(r_lo), jfmt(r_hi)
                )
            } else {
                format!(
                    "{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
                    name_field(&ids[i], sep), name_field(&ids[j], sep),
                    fmt(dn), fmt(dn_lo), fmt(dn_hi), fmt(ds), fmt(ds_lo), fmt(ds_hi), fmt(ratio), fmt(r_lo), fmt(r_hi), s = sep
                )
            }
        })
        .collect();

    let output_path = format!("{}_pairwise_bootstrap.{}", prefix, ext);
    let mut file = BufWriter::new(File::create(&output_path)?);
    if is_json {
        writeln!(file, "[")?;
        writeln!(file, "{}", rows.join(",\n"))?;
        writeln!(file, "]")?;
    } else {
        writeln!(file, "Seq1{s}Seq2{s}dN{s}dN_CI_low{s}dN_CI_high{s}dS{s}dS_CI_low{s}dS_CI_high{s}dN/dS{s}dNdS_CI_low{s}dNdS_CI_high", s = sep)?;
        for row in &rows {
            writeln!(file, "{}", row)?;
        }
    }
    Ok(output_path)
}

pub fn write_pairwise(
    ids: &[String],
    uidx_by_id: &[usize],
    n_u: usize,
    compute_pair: impl Fn(usize, usize) -> DsDn + Sync,
    cfg: &OutputConfig,
) -> anyhow::Result<()> {
    let output_prefix = cfg.prefix;
    let model = cfg.model;
    let sep = cfg.sep;
    let ext = cfg.ext;
    let summary = cfg.summary;
    let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
    let pb_write = ProgressBar::new(total_pairs_to_write as u64);
    pb_write.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Computing & writing pairs: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let nan_count = AtomicUsize::new(0);
    // Each parallel task sends exactly one (row_index, block) tuple. The writer
    // reassembles blocks in ascending row order so the output is deterministic
    // regardless of thread scheduling (only the out-of-order window is buffered).
    let (tx, rx) = unbounded::<(usize, String)>();

    let is_json = ext == "json";
    // Create the output file in this (parent) thread so a create/permission error is
    // returned cleanly instead of panicking the workers on a closed channel later.
    let out_path = format!("{}_pairwise_results.{}", output_prefix, ext);
    let mut out_file = BufWriter::new(
        File::create(&out_path).map_err(|e| anyhow::anyhow!("Cannot create '{}': {}", out_path, e))?,
    );
    let writer_thread = thread::spawn({
        move || -> Result<(), std::io::Error> {
            let mut pending: std::collections::HashMap<usize, String> = std::collections::HashMap::new();
            let mut next = 0usize;
            if is_json {
                out_file.write_all(b"[\n")?;
                let mut first = true;
                let emit = |out_file: &mut BufWriter<File>, first: &mut bool, block: &str| -> Result<(), std::io::Error> {
                    for line in block.lines() {
                        if !*first { out_file.write_all(b",\n")?; }
                        out_file.write_all(line.as_bytes())?;
                        *first = false;
                    }
                    Ok(())
                };
                for (i, block) in rx {
                    pending.insert(i, block);
                    while let Some(b) = pending.remove(&next) { emit(&mut out_file, &mut first, &b)?; next += 1; }
                }
                while let Some(b) = pending.remove(&next) { emit(&mut out_file, &mut first, &b)?; next += 1; }
                out_file.write_all(b"\n]\n")?;
            } else {
                let header = match model {
                    Model::Li => format!("Seq1{s}Seq2{s}dN(Ka){s}dS(Ks){s}dN/dS\n", s = sep),
                    Model::Nei => format!("Seq1{s}Seq2{s}dN{s}dS{s}dN/dS\n", s = sep),
                };
                out_file.write_all(header.as_bytes())?;
                for (i, block) in rx {
                    pending.insert(i, block);
                    while let Some(b) = pending.remove(&next) { out_file.write_all(b.as_bytes())?; next += 1; }
                }
                while let Some(b) = pending.remove(&next) { out_file.write_all(b.as_bytes())?; next += 1; }
            }
            out_file.flush()?;
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
                FloatAccum::new(),
            )
        },
        |(sender, row_cache, gen_map, cur_gen, local_buffer, local_stats), i| {
            *cur_gen = cur_gen.wrapping_add(1);
            let gen = *cur_gen;
            let u_i = uidx_by_id[i];
            // Compute the self-comparison rather than hard-coding 0/0: an all-N /
            // all-gap sequence has no comparable codons, so its dN/dS is NaN, not a
            // spurious 0.0. compute_pair short-circuits the u==u case, so this is cheap.
            row_cache[u_i] = compute_pair(u_i, u_i);
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
                // An undefined dN or dS makes the ratio undefined too: a NaN numerator
                // must not slip through the ds==0 branch and print as +inf (which reads as
                // extreme positive selection). Saturated dN over zero dS is genuinely NaN.
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

                if is_json {
                    let _ = writeln!(local_buffer,
                        "{{\"seq1\":\"{}\",\"seq2\":\"{}\",\"dN\":{},\"dS\":{},\"dN_dS\":{}}}",
                        json_escape(&ids[i]), json_escape(&ids[j]),
                        format_json_f64(result.dn), format_json_f64(result.ds), format_json_f64(ratio));
                } else {
                    let _ = writeln!(local_buffer, "{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}",
                        name_field(&ids[i], sep), name_field(&ids[j], sep),
                        norm_zero(result.dn), norm_zero(result.ds), norm_zero(ratio), s = sep);
                }

            }

            // Flush thread-local stats to shared accumulator (once per row)
            if let Some(stats) = summary {
                stats.flush_local(local_stats);
                local_stats.reset();
            }

            // Send this row's block exactly once, tagged with its index, so the
            // writer can emit rows in order (empty for the last row, which has no pairs).
            sender.send((i, std::mem::take(local_buffer)))
                .expect("Writer thread channel closed unexpectedly");
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

