use super::*;

// ─── shared helpers ───────────────────────────────────────

/// A fresh temp dir; the returned handle removes its files when dropped.
fn cov_tmp() -> tempfile::TempDir {
    tempfile::tempdir().expect("create temp dir")
}

/// Owned "cov" output prefix inside a temp dir (writers append their own suffix).
fn cov_prefix(dir: &tempfile::TempDir) -> String {
    dir.path()
        .join("cov")
        .to_str()
        .expect("utf8 temp path")
        .to_string()
}

/// Read a file and return its non-empty lines.
fn cov_read_nonempty(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", path, e))
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.to_string())
        .collect()
}

/// Split a delimited line into fields.
fn cov_fields(line: &str, sep: char) -> Vec<String> {
    line.split(sep).map(|s| s.to_string()).collect()
}

/// Sequential ids g0..g{n-1}.
fn cov_ids(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("g{}", i)).collect()
}

/// Structurally validate a writer-emitted JSON array (no serde available):
/// opens with '[', closes with ']', has exactly `expected_objects` object lines
/// (each opening '{' and closing '}', ignoring a trailing comma), balanced braces.
fn cov_assert_json_array(content: &str, expected_objects: usize) {
    let lines: Vec<String> = content
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(lines.first().map(|s| s.as_str()), Some("["), "must open with '['");
    assert_eq!(lines.last().map(|s| s.as_str()), Some("]"), "must close with ']'");
    let objs: Vec<&String> = lines.iter().filter(|l| l.starts_with('{')).collect();
    assert_eq!(objs.len(), expected_objects, "unexpected object count");
    for o in &objs {
        assert!(
            o.trim_end_matches(',').ends_with('}'),
            "object must end with a brace: {}",
            o
        );
    }
    assert_eq!(
        content.matches('{').count(),
        content.matches('}').count(),
        "unbalanced braces"
    );
}

// Deterministic finite pair (ds > 0 so the ratio is finite).
fn cov_finite_pair(a: usize, b: usize) -> DsDn {
    DsDn { dn: 0.1 + a as f64 * 0.001, ds: 0.2 + b as f64 * 0.001 }
}

// Some pairs saturate (NaN), exercising the saturation counter/info path.
fn cov_nan_pair(a: usize, b: usize) -> DsDn {
    if (a + b).is_multiple_of(2) {
        DsDn { dn: f64::NAN, ds: 0.2 }
    } else {
        DsDn { dn: 0.1, ds: 0.2 }
    }
}

fn cov_tests_stats(_i: usize, _j: usize) -> (f64, f64, f64, f64) {
    (0.1, 0.2, 0.01, 0.02)
}

fn cov_boot_stats(_i: usize, _j: usize) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    (0.1, 0.2, 0.08, 0.12, 0.15, 0.25, 0.4, 0.6)
}

fn cov_win_finite(_a: &[u8], _b: &[u8]) -> DsDn {
    DsDn { dn: 0.1, ds: 0.2 }
}

fn cov_win_nan(a: &[u8], _b: &[u8]) -> DsDn {
    // Window-0 slices of a sequence that begins with 0 saturate.
    if a.first().copied() == Some(0) {
        DsDn { dn: f64::NAN, ds: 0.2 }
    } else {
        DsDn { dn: 0.1, ds: 0.2 }
    }
}

// ─── spawn_ordered_writer ───────────────────────────────────

#[test]
fn cov_spawn_ordered_writer_reassembles_in_index_order() {
    // Blocks sent 2,0,1 must be emitted 0,1,2 after the header — the core
    // determinism guarantee of the ordered writer.
    let dir = cov_tmp();
    let path = dir.path().join("ordered.txt").to_str().unwrap().to_string();
    let (tx, handle) =
        spawn_ordered_writer(path.clone(), "HEADER\n".to_string()).expect("spawn writer");
    tx.send((2, "c\n".to_string())).unwrap();
    tx.send((0, "a\n".to_string())).unwrap();
    tx.send((1, "b\n".to_string())).unwrap();
    drop(tx);
    handle.join().expect("writer join").expect("writer io");
    let content = std::fs::read_to_string(&path).expect("read ordered");
    assert_eq!(content, "HEADER\na\nb\nc\n");
}

// ─── write_pairwise_tests ───────────────────────────────────

#[test]
fn cov_pairwise_tests_tsv_header_rows_and_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();

    let path = write_pairwise_tests(&ids, &uidx, cov_tests_stats, prefix.as_str(), '\t', "tsv")
        .expect("write tests");
    let lines = cov_read_nonempty(&path);
    // Header + 4*(4-1)/2 = 6 data rows.
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Seq1", "Seq2", "dN", "dS", "SE_dN", "SE_dS", "Z", "P_value"]
    );
    assert_eq!(lines.len(), 1 + 6, "header + 6 pairs");
    let first = std::fs::read_to_string(&path).unwrap();

    // Determinism: identical inputs -> byte-identical output.
    let path2 = write_pairwise_tests(&ids, &uidx, cov_tests_stats, prefix.as_str(), '\t', "tsv")
        .expect("write tests 2");
    assert_eq!(path, path2);
    let second = std::fs::read_to_string(&path2).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

#[test]
fn cov_pairwise_tests_json_is_valid_array() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let path = write_pairwise_tests(&ids, &uidx, cov_tests_stats, prefix.as_str(), ',', "json")
        .expect("write tests json");
    let content = std::fs::read_to_string(&path).unwrap();
    cov_assert_json_array(&content, 6);
}

// ─── write_pairwise_bootstrap ───────────────────────────────

#[test]
fn cov_pairwise_bootstrap_tsv_header_rows_and_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let path = write_pairwise_bootstrap(&ids, &uidx, cov_boot_stats, prefix.as_str(), '\t', "tsv")
        .expect("write bootstrap");
    let lines = cov_read_nonempty(&path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec![
            "Seq1", "Seq2", "dN", "dN_CI_low", "dN_CI_high", "dS", "dS_CI_low",
            "dS_CI_high", "dN/dS", "dNdS_CI_low", "dNdS_CI_high",
        ]
    );
    assert_eq!(lines.len(), 1 + 6, "header + 6 pairs");
    let first = std::fs::read_to_string(&path).unwrap();

    let path2 = write_pairwise_bootstrap(&ids, &uidx, cov_boot_stats, prefix.as_str(), '\t', "tsv")
        .expect("write bootstrap 2");
    let second = std::fs::read_to_string(&path2).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

#[test]
fn cov_pairwise_bootstrap_json_is_valid_array() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    let path = write_pairwise_bootstrap(&ids, &uidx, cov_boot_stats, prefix.as_str(), ',', "json")
        .expect("write bootstrap json");
    let content = std::fs::read_to_string(&path).unwrap();
    cov_assert_json_array(&content, 3);
}

// ─── write_pairwise ───────────────────────────────────────

#[test]
fn cov_write_pairwise_nei_tsv_header_rows_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: None,
    };
    write_pairwise(&ids, &uidx, 4, cov_finite_pair, &cfg).expect("write pairwise");
    let out_path = format!("{}_pairwise_results.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Seq1", "Seq2", "dN", "dS", "dN/dS"]
    );
    assert_eq!(lines.len(), 1 + 6, "header + 6 pairs");
    let first = std::fs::read_to_string(&out_path).unwrap();
    write_pairwise(&ids, &uidx, 4, cov_finite_pair, &cfg).expect("write pairwise 2");
    let second = std::fs::read_to_string(&out_path).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

#[test]
fn cov_write_pairwise_li_header_differs() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Li,
        summary: None,
    };
    write_pairwise(&ids, &uidx, 3, cov_finite_pair, &cfg).expect("write pairwise li");
    let out_path = format!("{}_pairwise_results.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Seq1", "Seq2", "dN(Ka)", "dS(Ks)", "dN/dS"]
    );
    assert_eq!(lines.len(), 1 + 3, "header + 3 pairs");
}

#[test]
fn cov_write_pairwise_json_is_valid_array() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: ',',
        ext: "json",
        model: Model::Nei,
        summary: None,
    };
    write_pairwise(&ids, &uidx, 4, cov_finite_pair, &cfg).expect("write pairwise json");
    let out_path = format!("{}_pairwise_results.json", prefix);
    let content = std::fs::read_to_string(&out_path).unwrap();
    cov_assert_json_array(&content, 6);
}

#[test]
fn cov_write_pairwise_records_every_pair_with_nan() {
    // NaN results still produce one row per pair and one summary record per pair.
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let summary = SummaryStats::new();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: Some(&summary),
    };
    write_pairwise(&ids, &uidx, 4, cov_nan_pair, &cfg).expect("write pairwise nan");
    let out_path = format!("{}_pairwise_results.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(lines.len(), 1 + 6, "header + 6 pairs even with NaN results");
    assert_eq!(
        summary.total_count.load(Ordering::Relaxed),
        6,
        "summary records every written pair"
    );
}

#[test]
fn cov_write_pairwise_nan_dn_over_zero_ds_ratio_is_nan_not_inf() {
    // Regression: when dN is undefined (NaN from nonsynonymous saturation) and dS is a
    // finite zero, the ratio is undefined and must print as NaN, not +inf (which reads
    // as extreme positive selection).
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(2);
    let uidx: Vec<usize> = (0..2).collect();
    let cfg = OutputConfig {
        prefix: prefix.as_str(), sep: '\t', ext: "tsv", model: Model::Nei, summary: None,
    };
    let nan_over_zero = |_a: usize, _b: usize| DsDn { dn: f64::NAN, ds: 0.0 };
    write_pairwise(&ids, &uidx, 2, nan_over_zero, &cfg).expect("write pairwise");
    let out_path = format!("{}_pairwise_results.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    let f = cov_fields(&lines[1], '\t');
    assert_eq!(f[2], "NaN", "dN column");
    assert_eq!(f[4], "NaN", "dN/dS must be NaN, not inf");
}

// ─── write_lineage ────────────────────────────────────────

#[test]
fn cov_write_lineage_header_rows_plot_and_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(4);
    let uidx: Vec<usize> = (0..4).collect();
    let lineage_indices: Vec<usize> = vec![0, 0, 1, 1];
    let lineage_names: Vec<String> = vec!["A".to_string(), "B".to_string()];
    let summary = SummaryStats::new();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: Some(&summary),
    };
    let plot = write_lineage(
        &ids, &uidx, 4, cov_finite_pair, &lineage_indices, &lineage_names, &cfg,
    )
    .expect("write lineage");

    let out_path = format!("{}_lineage_summary.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"]
    );
    // Each of 4 genomes sees both lineages among the other 3 -> 2 rows each = 8.
    assert_eq!(lines.len(), 1 + 8, "header + 8 lineage rows");

    // One plot entry per emitted row, sorted by (genome, lineage).
    assert_eq!(plot.len(), 8);
    for w in plot.windows(2) {
        assert!(
            (w[0].0.as_str(), w[0].1.as_str()) <= (w[1].0.as_str(), w[1].1.as_str()),
            "plot data must be sorted by (genome, lineage)"
        );
    }
    for (_, _, r) in &plot {
        assert!(r.is_finite() && *r >= 0.0, "ratio finite and non-negative");
    }

    let first = std::fs::read_to_string(&out_path).unwrap();
    let _ = write_lineage(
        &ids, &uidx, 4, cov_finite_pair, &lineage_indices, &lineage_names, &cfg,
    )
    .expect("write lineage 2");
    let second = std::fs::read_to_string(&out_path).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

#[test]
fn cov_write_lineage_all_invalid_identical_genomes_excluded_not_zero() {
    // Regression: two identical all-N/all-gap genomes deduplicate to the same unique
    // index, so each genome's only cross-lineage neighbour is the self/diagonal. That
    // must be NaN (no comparable codons) and excluded — not a spurious 0.0 from a
    // hard-coded {0,0} diagonal reused across de-duplicated genomes.
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = vec!["L1_a".to_string(), "L2_b".to_string()];
    let uidx: Vec<usize> = vec![0, 0]; // identical -> same unique index
    let lineage_indices: Vec<usize> = vec![0, 1];
    let lineage_names: Vec<String> = vec!["L1".to_string(), "L2".to_string()];
    // compute_pair returns NaN for the all-invalid self/identity comparison (u == u).
    let compute = |a: usize, b: usize| {
        if a == b { DsDn { dn: f64::NAN, ds: f64::NAN } } else { DsDn { dn: 0.1, ds: 0.2 } }
    };
    let cfg = OutputConfig {
        prefix: prefix.as_str(), sep: '\t', ext: "tsv", model: Model::Nei, summary: None,
    };
    write_lineage(&ids, &uidx, 1, compute, &lineage_indices, &lineage_names, &cfg)
        .expect("write lineage");
    let out_path = format!("{}_lineage_summary.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(lines.len(), 1, "header only — all-invalid identical genomes are excluded, not 0.0");
}

#[test]
fn cov_write_lineage_no_summary_yields_empty_plot() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    let lineage_indices: Vec<usize> = vec![0, 1, 1];
    let lineage_names: Vec<String> = vec!["A".to_string(), "B".to_string()];
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: ',',
        ext: "csv",
        model: Model::Nei,
        summary: None,
    };
    let plot = write_lineage(
        &ids, &uidx, 3, cov_finite_pair, &lineage_indices, &lineage_names, &cfg,
    )
    .expect("write lineage no summary");
    assert!(plot.is_empty(), "no summary -> empty plot data");

    let out_path = format!("{}_lineage_summary.csv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], ','),
        vec!["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"]
    );
    // g0 (A) sees only lineage B -> 1 row (the count==0 skip fires for A);
    // g1, g2 each see A and B -> 2 rows. Total 5.
    assert_eq!(lines.len(), 1 + 5, "header + 5 lineage rows");
}

// ─── write_group_average ───────────────────────────────────

#[test]
fn cov_write_group_average_branches_header_and_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    // Group key = prefix before '_': G1 = {0,1}, G2 = {2}.
    let ids: Vec<String> =
        vec!["G1_a".to_string(), "G1_b".to_string(), "G2_x".to_string()];
    let uidx: Vec<usize> = (0..3).collect();
    let summary = SummaryStats::new();
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: Some(&summary),
    };
    let plot = write_group_average(&ids, &uidx, cov_finite_pair, false, &cfg)
        .expect("write group avg");

    let out_path = format!("{}_group_avg_dn_ds.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec![
            "Group1", "Group2", "NumSeqs1", "NumSeqs2", "NumComparisons",
            "Mean_dN/dS", "StdError", "95%CI",
        ]
    );
    // 2 groups -> 2*(2+1)/2 = 3 group pairs -> 3 data rows.
    assert_eq!(lines.len(), 1 + 3, "header + 3 group-pair rows");

    // Exactly one row per branch: single-comparison (N/A), multi-comparison
    // (numeric CI), and empty (0 comparisons, [NaN, NaN]).
    let data: Vec<Vec<String>> = lines[1..].iter().map(|l| cov_fields(l, '\t')).collect();
    for r in &data {
        assert_eq!(r.len(), 8, "each row has 8 columns");
    }
    let na_rows = data.iter().filter(|r| r[6] == "N/A").count();
    let empty_rows = data
        .iter()
        .filter(|r| r[4] == "0" && r[7] == "[NaN, NaN]")
        .count();
    let ci_rows = data
        .iter()
        .filter(|r| r[7].starts_with('[') && r[7] != "[NaN, NaN]")
        .count();
    assert_eq!(na_rows, 1, "one single-comparison (N/A) row");
    assert_eq!(empty_rows, 1, "one empty-group NaN row");
    assert_eq!(ci_rows, 1, "one multi-comparison CI row");

    assert_eq!(plot.len(), 3);
    for w in plot.windows(2) {
        assert!(w[0].label <= w[1].label, "group plot data sorted by label");
    }

    let first = std::fs::read_to_string(&out_path).unwrap();
    let _ = write_group_average(&ids, &uidx, cov_finite_pair, false, &cfg)
        .expect("write group avg 2");
    let second = std::fs::read_to_string(&out_path).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

// ─── write_pairwise_windows ────────────────────────────────

#[test]
fn cov_write_pairwise_windows_nei_header_rows_determinism() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    let seqs: Vec<Vec<u8>> = vec![
        vec![1u8, 2, 3, 4, 5, 6],
        vec![7u8, 8, 9, 10, 11, 12],
        vec![13u8, 14, 15, 16, 17, 18],
    ];
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: None,
    };
    // seq_len=6, window_size=3, step=3 -> (6-3)/3 + 1 = 2 windows.
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_finite, 3, 3, None, &cfg)
        .expect("write windows");
    let out_path = format!("{}_pairwise_windows.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Seq1", "Seq2", "Window_Start", "Window_End", "dN", "dS", "dN/dS"]
    );
    // 3 pairs * 2 windows = 6 rows.
    assert_eq!(lines.len(), 1 + 6, "header + 6 window rows");

    // 1-based inclusive window coordinates: (1,3) and (4,6).
    let data: Vec<Vec<String>> = lines[1..].iter().map(|l| cov_fields(l, '\t')).collect();
    for r in &data {
        assert_eq!(r.len(), 7, "each row has 7 columns");
        let start: usize = r[2].parse().unwrap();
        let end: usize = r[3].parse().unwrap();
        assert!(
            (start, end) == (1, 3) || (start, end) == (4, 6),
            "unexpected window {:?}",
            (start, end)
        );
    }

    let first = std::fs::read_to_string(&out_path).unwrap();
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_finite, 3, 3, None, &cfg)
        .expect("write windows 2");
    let second = std::fs::read_to_string(&out_path).unwrap();
    assert_eq!(first, second, "output must be deterministic");
}

#[test]
fn cov_write_pairwise_windows_all_invalid_identical_is_nan_not_zero() {
    // Regression: two identical all-N/all-gap sequences deduplicate to the same unique
    // index (u_i == u_j). Every window must be NaN (no comparable codons), not a
    // hard-coded 0.0 — even though the compute closure would return a finite value.
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = vec!["a".to_string(), "b".to_string()];
    let uidx: Vec<usize> = vec![0, 0]; // identical -> same unique index
    let inv = crate::codon::INVALID_CODON;
    let seqs: Vec<Vec<u8>> = vec![vec![inv, inv, inv, inv]]; // all-invalid unique sequence
    let cfg = OutputConfig {
        prefix: prefix.as_str(), sep: '\t', ext: "tsv", model: Model::Nei, summary: None,
    };
    // 2 windows of size 2; cov_win_finite would give 0.1/0.2 but the identity +
    // all-invalid guard must override every window to NaN.
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_finite, 2, 2, None, &cfg)
        .expect("write windows");
    let out_path = format!("{}_pairwise_windows.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(lines.len(), 1 + 2, "header + 2 window rows");
    for l in &lines[1..] {
        let f = cov_fields(l, '\t');
        assert_eq!(f[4], "NaN", "dN must be NaN for an all-invalid identical window");
        assert_eq!(f[5], "NaN", "dS must be NaN");
        assert_eq!(f[6], "NaN", "dN/dS must be NaN");
    }
}

#[test]
fn cov_write_pairwise_windows_li_header() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(2);
    let uidx: Vec<usize> = (0..2).collect();
    let seqs: Vec<Vec<u8>> = vec![vec![1u8, 2, 3, 4], vec![5u8, 6, 7, 8]];
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Li,
        summary: None,
    };
    // seq_len=4, window_size=2, step=2 -> 2 windows; 1 pair -> 2 rows.
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_finite, 2, 2, None, &cfg)
        .expect("write windows li");
    let out_path = format!("{}_pairwise_windows.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(
        cov_fields(&lines[0], '\t'),
        vec!["Seq1", "Seq2", "Window_Start", "Window_End", "dN(Ka)", "dS(Ks)", "dN/dS"]
    );
    assert_eq!(lines.len(), 1 + 2, "header + 2 window rows");
}

#[test]
fn cov_write_pairwise_windows_zero_windows_header_only() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    let seqs: Vec<Vec<u8>> = vec![
        vec![1u8, 2, 3, 4, 5, 6],
        vec![7u8, 8, 9, 10, 11, 12],
        vec![13u8, 14, 15, 16, 17, 18],
    ];
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: None,
    };
    // window_size (10) > seq_len (6) -> num_windows = 0 -> header only.
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_finite, 10, 3, None, &cfg)
        .expect("write windows zero");
    let out_path = format!("{}_pairwise_windows.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    assert_eq!(lines.len(), 1, "header only, no window rows");
}

#[test]
fn cov_write_pairwise_windows_nan_with_windowstats() {
    let dir = cov_tmp();
    let prefix = cov_prefix(&dir);
    let ids = cov_ids(3);
    let uidx: Vec<usize> = (0..3).collect();
    // Sequence 0 begins with 0 -> its window-0 slice saturates in cov_win_nan.
    let seqs: Vec<Vec<u8>> = vec![
        vec![0u8, 2, 3, 4, 5, 6],
        vec![7u8, 8, 9, 10, 11, 12],
        vec![13u8, 14, 15, 16, 17, 18],
    ];
    let ws = WindowStats::new(2);
    let cfg = OutputConfig {
        prefix: prefix.as_str(),
        sep: '\t',
        ext: "tsv",
        model: Model::Nei,
        summary: None,
    };
    write_pairwise_windows(&ids, &uidx, &seqs, cov_win_nan, 3, 3, Some(&ws), &cfg)
        .expect("write windows nan");
    let out_path = format!("{}_pairwise_windows.tsv", prefix);
    let lines = cov_read_nonempty(&out_path);
    // Every window still written despite NaN results.
    assert_eq!(lines.len(), 1 + 6, "header + 6 window rows even with NaN");
    let wd = ws.get_window_data();
    assert_eq!(wd.len(), 2, "one entry per window bin");
}


#[test]
fn json_escape_neutralizes_quotes_backslash_controls() {
    assert_eq!(json_escape("seq\"1"), "seq\\\"1");
    assert_eq!(json_escape("a\\b"), "a\\\\b");
    assert_eq!(json_escape("t\tn\n"), "t\\tn\\n");
    assert_eq!(json_escape("plain"), "plain");
    assert_eq!(json_escape("\u{1}"), "\\u0001");
}
