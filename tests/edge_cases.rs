use std::fs;
use std::process::Command;

fn binary_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_eskaks"))
}

// ─── File handling ───────────────────────────────────────────────────────────

#[test]
fn nonexistent_file_gives_clear_error() {
    let output = Command::new(binary_path())
        .args(["nonexistent_file.fasta", "-o", "/tmp/eskaks_edge_nonexist"])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Failed to open input file"), "stderr: {}", stderr);
}

#[test]
fn empty_file_gives_clear_error() {
    let path = "/tmp/eskaks_test_empty.fasta";
    fs::write(path, "").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_empty"])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("No sequences found")
        || stderr.contains("at least 2 sequences")
        || stderr.contains("Failed to open")
        || stderr.contains("empty"),
        "stderr: {}", stderr);
    fs::remove_file(path).ok();
}

#[test]
fn single_sequence_gives_clear_error() {
    let path = "/tmp/eskaks_test_single.fasta";
    fs::write(path, ">seq1\nATGAAACCC\n").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_single"])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("at least 2 sequences"), "stderr: {}", stderr);
    fs::remove_file(path).ok();
}

// ─── Malformed input ─────────────────────────────────────────────────────────

#[test]
fn not_fasta_gives_error() {
    let path = "/tmp/eskaks_test_notfasta.txt";
    fs::write(path, "this is not a FASTA file\njust plain text\n").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_notfasta"])
        .output()
        .expect("spawn");
    assert!(!output.status.success(), "should fail on non-FASTA input");
    fs::remove_file(path).ok();
}

#[test]
fn sequence_not_divisible_by_3_still_works() {
    let path = "/tmp/eskaks_test_mod3.fasta";
    // 10 bases = 3 codons + 1 trailing base (should warn but not fail)
    fs::write(path, ">seq1\nATGAAACCCG\n>seq2\nATGAAACCCT\n").unwrap();
    let output = Command::new(binary_path())
        .env("RUST_LOG", "warn")
        .args([path, "-o", "/tmp/eskaks_edge_mod3"])
        .output()
        .expect("spawn");
    assert!(output.status.success(), "should succeed with warning");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("not divisible by 3") || stderr.contains("trailing"),
        "should warn about trailing bases. stderr: {}", stderr);
    fs::remove_file(path).ok();
    fs::remove_file("/tmp/eskaks_edge_mod3_pairwise_results.tsv").ok();
}

#[test]
fn all_gap_sequences_different_produce_nan() {
    let path = "/tmp/eskaks_test_allgap.fasta";
    // One all-gaps, one real → NaN because no valid codon pairs
    fs::write(path, ">seq1\n---------\n>seq2\nATGAAACCC\n").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_allgap"])
        .output()
        .expect("spawn");
    assert!(output.status.success(), "should succeed but produce NaN results");
    let results = fs::read_to_string("/tmp/eskaks_edge_allgap_pairwise_results.tsv").unwrap();
    // One seq all gaps → no valid codon pairs → NaN
    assert!(results.contains("NaN") || results.contains("nan"),
        "should contain NaN when one sequence is all gaps: {}", results);
    fs::remove_file(path).ok();
    fs::remove_file("/tmp/eskaks_edge_allgap_pairwise_results.tsv").ok();
}

#[test]
fn identical_sequences_produce_zero() {
    let path = "/tmp/eskaks_test_identical.fasta";
    fs::write(path, ">seq1\nATGAAACCCGGG\n>seq2\nATGAAACCCGGG\n").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_identical"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let results = fs::read_to_string("/tmp/eskaks_edge_identical_pairwise_results.tsv").unwrap();
    let lines: Vec<&str> = results.lines().collect();
    assert_eq!(lines.len(), 2, "should have header + 1 data row");
    // dN and dS should both be 0
    let data = lines[1];
    assert!(data.contains("0.0000") || data.contains("0\t") || data.contains("\t0."),
        "identical seqs should give 0 dN/dS: {}", data);
    fs::remove_file(path).ok();
    fs::remove_file("/tmp/eskaks_edge_identical_pairwise_results.tsv").ok();
}

// ─── Flag validation ─────────────────────────────────────────────────────────

#[test]
fn window_size_zero_fails() {
    let path = "tests/data/synthetic.fasta";
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_win0", "--window-size", "0"])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
}

#[test]
fn window_with_group_average_fails() {
    let path = "tests/data/synthetic.fasta";
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_wingrp", "--window-size", "10", "--group-average"])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("cannot be used with"), "stderr: {}", stderr);
}

// ─── Output format ──────────────────────────────────────────────────────────

#[test]
fn csv_format_produces_commas() {
    let path = "tests/data/synthetic.fasta";
    let status = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_csv", "--format", "csv"])
        .status()
        .expect("spawn");
    assert!(status.success());
    let results = fs::read_to_string("/tmp/eskaks_edge_csv_pairwise_results.csv").unwrap();
    assert!(results.contains(','), "CSV should use commas");
    assert!(!results.lines().next().unwrap().contains('\t'), "CSV should not use tabs");
    fs::remove_file("/tmp/eskaks_edge_csv_pairwise_results.csv").ok();
}

// ─── Large gap handling ─────────────────────────────────────────────────────

#[test]
fn mostly_gaps_still_computes() {
    let path = "/tmp/eskaks_test_gaps.fasta";
    // 4 codons: 3 gaps + 1 real
    fs::write(path, ">seq1\n---------ATG\n>seq2\n---------ATG\n").unwrap();
    let output = Command::new(binary_path())
        .args([path, "-o", "/tmp/eskaks_edge_gaps"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    fs::remove_file(path).ok();
    fs::remove_file("/tmp/eskaks_edge_gaps_pairwise_results.tsv").ok();
}
