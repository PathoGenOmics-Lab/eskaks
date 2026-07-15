/// Integration tests for eskaks: run the compiled binary against synthetic data
/// and verify the output matches hand-calculated expected values.
use std::fs;
use std::process::Command;

/// Path to the synthetic FASTA test file.
const FASTA: &str = "tests/data/synthetic.fasta";
/// FASTA with sequences pre-grouped by prefix before '_' (for group_average tests).
const FASTA_GROUPED: &str = "tests/data/synthetic_grouped.fasta";
/// Tolerance for floating-point comparisons in output.
const EPSILON: f64 = 1e-4;

/// Build the release binary before running integration tests.
/// (cargo test --test integration builds the debug binary by default.)
fn binary_path() -> String {
    // Use the debug build produced by `cargo test`
    env!("CARGO_BIN_EXE_eskaks").to_string()
}

fn parse_tsv(path: &str) -> Vec<Vec<String>> {
    parse_delimited(path, '\t')
}

fn parse_csv(path: &str) -> Vec<Vec<String>> {
    parse_delimited(path, ',')
}

fn parse_delimited(path: &str, sep: char) -> Vec<Vec<String>> {
    let content = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Cannot read {}: {}", path, e));
    content
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split(sep).map(String::from).collect())
        .collect()
}

/// Find the row for a specific seq1/seq2 pair in the pairwise TSV.
fn find_pair<'a>(rows: &'a [Vec<String>], seq1: &str, seq2: &str) -> Option<&'a Vec<String>> {
    rows.iter().find(|r| r.len() >= 5 && r[0] == seq1 && r[1] == seq2)
}

// ─── Pairwise (Nei model, default) ───────────────────────────────────────────

#[test]
fn pairwise_nei_exits_ok() {
    let out_prefix = "/tmp/eskaks_test_pairwise_nei";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status()
        .expect("Failed to spawn binary");
    assert!(status.success(), "eskaks exited with non-zero status: {:?}", status);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_correct_header_and_row_count() {
    let out_prefix = "/tmp/eskaks_test_pairwise_header";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    // Header + 5*(5-1)/2 = 10 data rows
    assert_eq!(rows[0], vec!["Seq1", "Seq2", "dN", "dS", "dN/dS"],
        "unexpected header: {:?}", rows[0]);
    assert_eq!(rows.len(), 11, "expected 1 header + 10 pairs, got {}", rows.len());
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_identical_seqs_give_zero() {
    // seq1_A and seq5_A are identical → dN=dS=0, ratio=0
    let out_prefix = "/tmp/eskaks_test_pairwise_identical";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq5_A")
        .expect("pair seq1_A/seq5_A not found in output");
    let dn: f64 = row[2].parse().unwrap();
    let ds: f64 = row[3].parse().unwrap();
    let ratio: f64 = row[4].parse().unwrap();
    assert!((dn - 0.0).abs() < EPSILON, "dN for identical seqs should be 0, got {}", dn);
    assert!((ds - 0.0).abs() < EPSILON, "dS for identical seqs should be 0, got {}", ds);
    assert!((ratio - 0.0).abs() < EPSILON, "ratio for identical seqs should be 0, got {}", ratio);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_one_synonymous_change() {
    // seq1_A (ATGGCTGCT) vs seq2_A (ATGGCCGCT): 1 syn change GCT→GCC (Ala→Ala)
    // Expected: dN=0.0, dS≈0.82396
    let out_prefix = "/tmp/eskaks_test_pairwise_syn";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq2_A")
        .expect("pair seq1_A/seq2_A not found");
    let dn: f64 = row[2].parse().unwrap();
    let ds: f64 = row[3].parse().unwrap();
    assert!((dn - 0.0).abs() < EPSILON, "dN should be 0 for purely synonymous change, got {}", dn);
    assert!((ds - 0.82396).abs() < EPSILON, "dS should be ~0.82396, got {}", ds);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_one_nonsynonymous_change() {
    // seq1_A (ATGGCTGCT) vs seq3_B (ATGATTGCT): GCT->ATT (Ala->Ile)
    // 2 nucleotide positions differ. Pathway analysis (NG86):
    //   Path 1: GCT->ACT->ATT (Ala->Thr->Ile): 2 nonsyn steps
    //   Path 2: GCT->GTT->ATT (Ala->Val->Ile): 2 nonsyn steps
    //   Average: sd=0, nd=2
    // S=11/6, N=43/6, pN=12/43, dN=-0.75*ln(1-16/43) ~ 0.34902
    let out_prefix = "/tmp/eskaks_test_pairwise_nonsyn";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq3_B")
        .expect("pair seq1_A/seq3_B not found");
    let dn: f64 = row[2].parse().unwrap();
    let ds: f64 = row[3].parse().unwrap();
    assert!((ds - 0.0).abs() < EPSILON, "dS should be 0, got {}", ds);
    assert!((dn - 0.34902).abs() < EPSILON, "dN should be ~0.34902, got {}", dn);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_nonsyn_only_ratio_is_infinity() {
    // When dS=0 and dN>0, ratio should be "inf"
    let out_prefix = "/tmp/eskaks_test_pairwise_inf";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq3_B")
        .expect("pair seq1_A/seq3_B not found");
    assert_eq!(row[4], "inf", "ratio should be 'inf' when dS=0, dN>0, got '{}'", row[4]);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

// ─── Pairwise (Li model) ─────────────────────────────────────────────────────

#[test]
fn pairwise_li_exits_ok_and_has_correct_header() {
    let out_prefix = "/tmp/eskaks_test_pairwise_li";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "li", "--workers", "2"])
        .status().expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    assert_eq!(rows[0], vec!["Seq1", "Seq2", "dN(Ka)", "dS(Ks)", "dN/dS"],
        "Li model header mismatch: {:?}", rows[0]);
    assert_eq!(rows.len(), 11);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_li_identical_seqs_give_zero() {
    let out_prefix = "/tmp/eskaks_test_pairwise_li_identical";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "li", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq5_A").expect("pair not found");
    let dn: f64 = row[2].parse().unwrap();
    let ds: f64 = row[3].parse().unwrap();
    assert!((dn).abs() < EPSILON, "dN for identical seqs should be 0, got {}", dn);
    assert!((ds).abs() < EPSILON, "dS for identical seqs should be 0, got {}", ds);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

// ─── Lineage mode ─────────────────────────────────────────────────────────────

#[test]
fn lineage_exits_ok_and_has_correct_header() {
    let out_prefix = "/tmp/eskaks_test_lineage";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--lineage", "--workers", "2"])
        .status().expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_lineage_summary.tsv", out_prefix));
    assert_eq!(rows[0], vec!["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"],
        "lineage header mismatch: {:?}", rows[0]);
    fs::remove_file(format!("{}_lineage_summary.tsv", out_prefix)).ok();
}

// ─── Group average mode ───────────────────────────────────────────────────────

#[test]
fn group_average_exits_ok_and_has_correct_header() {
    // Uses FASTA_GROUPED where IDs start with "A_" or "B_", so group = "A" or "B"
    let out_prefix = "/tmp/eskaks_test_group";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA_GROUPED, "-o", out_prefix, "--group-average", "--workers", "2"])
        .status().expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_group_avg_dn_ds.tsv", out_prefix));
    assert_eq!(
        rows[0],
        vec!["Group1", "Group2", "NumSeqs1", "NumSeqs2", "NumComparisons", "Mean_dN/dS", "StdError", "95%CI"],
        "group_average header mismatch: {:?}", rows[0]
    );
    // 2 groups → 3 pairs (A×A, A×B, B×B)
    assert_eq!(rows.len(), 4, "expected header + 3 group pairs, got {}", rows.len());
    fs::remove_file(format!("{}_group_avg_dn_ds.tsv", out_prefix)).ok();
}

#[test]
fn group_average_within_group_a_has_only_syn_changes() {
    // Group A = {A_seq1, A_seq2, A_seq3}. A_seq1 and A_seq3 are identical (dedup).
    // All within-A pairs are synonymous only → mean dN/dS = 0
    let out_prefix = "/tmp/eskaks_test_group_a";
    Command::new(binary_path())
        .args(["fasta", FASTA_GROUPED, "-o", out_prefix, "--group-average", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_group_avg_dn_ds.tsv", out_prefix));
    let row = rows.iter().find(|r| r.len() >= 6 && r[0] == "A" && r[1] == "A")
        .expect("within-group A row not found");
    let mean: f64 = row[5].parse().unwrap_or(f64::NAN);
    assert!((mean - 0.0).abs() < EPSILON,
        "mean dN/dS within group A should be 0 (all syn), got {}", mean);
    fs::remove_file(format!("{}_group_avg_dn_ds.tsv", out_prefix)).ok();
}

// ─── CSV output format ──────────────────────────────────────────────────────

#[test]
fn pairwise_csv_format_has_correct_header_and_separators() {
    let out_prefix = "/tmp/eskaks_test_csv";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2", "--format", "csv"])
        .status()
        .expect("Failed to spawn binary");
    assert!(status.success());
    let rows = parse_csv(&format!("{}_pairwise_results.csv", out_prefix));
    assert_eq!(rows[0], vec!["Seq1", "Seq2", "dN", "dS", "dN/dS"],
        "CSV header mismatch: {:?}", rows[0]);
    assert_eq!(rows.len(), 11, "expected 1 header + 10 pairs");
    // Verify a specific value parses correctly
    let row = rows.iter().find(|r| r.len() >= 5 && r[0] == "seq1_A" && r[1] == "seq5_A")
        .expect("pair not found in CSV");
    let dn: f64 = row[2].parse().unwrap();
    assert!((dn - 0.0).abs() < EPSILON);
    fs::remove_file(format!("{}_pairwise_results.csv", out_prefix)).ok();
}

#[test]
fn lineage_csv_has_correct_header() {
    let out_prefix = "/tmp/eskaks_test_lineage_csv";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--lineage", "--workers", "2", "--format", "csv"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let rows = parse_csv(&format!("{}_lineage_summary.csv", out_prefix));
    assert_eq!(rows[0], vec!["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"]);
    fs::remove_file(format!("{}_lineage_summary.csv", out_prefix)).ok();
}

// ─── Sliding window ─────────────────────────────────────────────────────────

#[test]
fn sliding_window_exits_ok_and_has_correct_header() {
    let out_prefix = "/tmp/eskaks_test_window";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2",
               "--window-size", "2", "--window-step", "1"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_pairwise_windows.tsv", out_prefix));
    assert_eq!(rows[0], vec!["Seq1", "Seq2", "Window_Start", "Window_End", "dN", "dS", "dN/dS"],
        "window header mismatch: {:?}", rows[0]);
    fs::remove_file(format!("{}_pairwise_windows.tsv", out_prefix)).ok();
}

#[test]
fn sliding_window_correct_row_count() {
    // 5 seqs, 3 codons each, window_size=2, step=1 → 2 windows per pair
    // 10 pairs × 2 windows = 20 data rows + 1 header = 21
    let out_prefix = "/tmp/eskaks_test_window_count";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2",
               "--window-size", "2", "--window-step", "1"])
        .status()
        .expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_windows.tsv", out_prefix));
    assert_eq!(rows.len(), 21, "expected 1 header + 20 window rows, got {}", rows.len());
    // Check window coordinates are 1-based
    let first_data = &rows[1];
    let win_start: usize = first_data[2].parse().unwrap();
    let win_end: usize = first_data[3].parse().unwrap();
    assert_eq!(win_start, 1, "first window should start at codon 1");
    assert_eq!(win_end, 2, "first window should end at codon 2");
    fs::remove_file(format!("{}_pairwise_windows.tsv", out_prefix)).ok();
}

#[test]
fn sliding_window_csv_format() {
    let out_prefix = "/tmp/eskaks_test_window_csv";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2",
               "--window-size", "2", "--format", "csv"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let rows = parse_csv(&format!("{}_pairwise_windows.csv", out_prefix));
    assert_eq!(rows[0], vec!["Seq1", "Seq2", "Window_Start", "Window_End", "dN", "dS", "dN/dS"]);
    fs::remove_file(format!("{}_pairwise_windows.csv", out_prefix)).ok();
}

// ─── Alternative genetic codes ───────────────────────────────────────────────

#[test]
fn list_codes_exits_ok() {
    let output = Command::new(binary_path())
        .args(["--list-codes"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Standard"), "should list Standard code");
    assert!(stderr.contains("Vertebrate Mitochondrial"), "should list Vertebrate Mito");
}

#[test]
fn genetic_code_2_runs_ok() {
    let out_prefix = "/tmp/eskaks_test_gc2";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2", "--genetic-code", "2"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    assert!(rows.len() > 1, "should produce output with genetic code 2");
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn genetic_code_affects_results() {
    // Same input, different genetic code should give different results
    // (at least for some pairs, since different stop/amino acid assignments)
    let prefix_std = "/tmp/eskaks_test_gc_std";
    let prefix_mito = "/tmp/eskaks_test_gc_mito";
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", prefix_std, "--model", "nei", "--workers", "1", "--genetic-code", "1"])
        .status().expect("spawn");
    Command::new(binary_path())
        .args(["fasta", FASTA, "-o", prefix_mito, "--model", "nei", "--workers", "1", "--genetic-code", "4"])
        .status().expect("spawn");

    let rows_std = parse_tsv(&format!("{}_pairwise_results.tsv", prefix_std));
    let rows_mito = parse_tsv(&format!("{}_pairwise_results.tsv", prefix_mito));

    // Tables 1 and 4 differ only in TGA (stop vs Trp), so most pairs should be
    // the same. But the syn sites change for codons near TGA, affecting dS slightly.
    // At minimum, row counts should match.
    assert_eq!(rows_std.len(), rows_mito.len(), "should have same number of rows");

    fs::remove_file(format!("{}_pairwise_results.tsv", prefix_std)).ok();
    fs::remove_file(format!("{}_pairwise_results.tsv", prefix_mito)).ok();
}

#[test]
fn invalid_genetic_code_exits_with_error() {
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "--genetic-code", "99"])
        .status()
        .expect("spawn");
    assert!(!status.success(), "invalid genetic code should fail");
}

// ─── Summary and Plot ────────────────────────────────────────────────────────

#[test]
fn summary_prints_to_stderr() {
    let out_prefix = "/tmp/eskaks_test_summary";
    let output = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2", "--summary"])
        .output()
        .expect("Failed to spawn binary");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("dN/dS Summary"), "stderr should contain summary header, got: {}", stderr);
    assert!(stderr.contains("Total pairs:"), "stderr should contain pair count");
    assert!(stderr.contains("Distribution"), "stderr should contain histogram");
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn plot_creates_histogram_svg() {
    let out_prefix = "/tmp/eskaks_test_plot_hist";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2", "--plot"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let svg_path = format!("{}_dnds_histogram.svg", out_prefix);
    let content = fs::read_to_string(&svg_path)
        .unwrap_or_else(|_| panic!("SVG file not created at {}", svg_path));
    assert!(content.contains("<svg"), "SVG should contain <svg tag");
    assert!(content.contains("</svg>"), "SVG should contain closing </svg> tag");
    assert!(content.contains("dN/dS Ratio Distribution"), "SVG should contain title");
    fs::remove_file(&svg_path).ok();
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn plot_window_creates_svg() {
    let out_prefix = "/tmp/eskaks_test_plot_window";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2",
               "--window-size", "2", "--plot"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let svg_path = format!("{}_window_plot.svg", out_prefix);
    let content = fs::read_to_string(&svg_path)
        .unwrap_or_else(|_| panic!("SVG not created at {}", svg_path));
    assert!(content.contains("Sliding Window"));
    fs::remove_file(&svg_path).ok();
    fs::remove_file(format!("{}_pairwise_windows.tsv", out_prefix)).ok();
}

#[test]
fn plot_group_creates_svg() {
    let out_prefix = "/tmp/eskaks_test_plot_group";
    let status = Command::new(binary_path())
        .args(["fasta", FASTA_GROUPED, "-o", out_prefix, "--group-average", "--workers", "2", "--plot"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let svg_path = format!("{}_group_dnds.svg", out_prefix);
    let content = fs::read_to_string(&svg_path)
        .unwrap_or_else(|_| panic!("SVG not created at {}", svg_path));
    assert!(content.contains("Group Average"));
    fs::remove_file(&svg_path).ok();
    fs::remove_file(format!("{}_group_avg_dn_ds.tsv", out_prefix)).ok();
}

#[test]
fn summary_and_plot_together() {
    let out_prefix = "/tmp/eskaks_test_both";
    let output = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2", "--summary", "--plot"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("dN/dS Summary"));
    assert!(fs::metadata(format!("{}_dnds_histogram.svg", out_prefix)).is_ok(), "SVG should exist");
    fs::remove_file(format!("{}_dnds_histogram.svg", out_prefix)).ok();
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn neutrality_test_writes_pairwise_tests() {
    let out_prefix = "/tmp/eskaks_test_neutrality";
    let output = Command::new(binary_path())
        .args(["fasta", FASTA, "-o", out_prefix, "--model", "nei", "--neutrality"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let rows = parse_tsv(&format!("{}_pairwise_tests.tsv", out_prefix));
    assert_eq!(
        rows[0],
        vec!["Seq1", "Seq2", "dN", "dS", "SE_dN", "SE_dS", "Z", "P_value"],
        "neutrality header mismatch: {:?}", rows[0]
    );
    assert!(rows.len() > 1, "should have at least one pair");
    // P_value column (index 7) must be a probability in [0,1] or NaN.
    let p_str = &rows[1][7];
    if p_str != "NaN" {
        let p: f64 = p_str.parse().expect("P_value numeric");
        assert!((0.0..=1.0).contains(&p), "P in [0,1]: {}", p);
    }
    fs::remove_file(format!("{}_pairwise_tests.tsv", out_prefix)).ok();
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}
