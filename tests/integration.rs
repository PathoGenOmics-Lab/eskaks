/// Integration tests for eskaks: run the compiled binary against synthetic data
/// and verify the output matches hand-calculated expected values.
use std::fs;
use std::process::Command;

/// Path to the synthetic FASTA test file.
const FASTA: &str = "tests/data/synthetic.fasta";
/// FASTA with sequences pre-grouped by prefix before '_' (for group_average tests).
const FASTA_GROUPED: &str = "tests/data/synthetic_grouped.fasta";
/// Tolerance for floating-point comparisons in TSV output.
const EPSILON: f64 = 1e-4;

/// Build the release binary before running integration tests.
/// (cargo test --test integration builds the debug binary by default.)
fn binary_path() -> String {
    // Use the debug build produced by `cargo test`
    env!("CARGO_BIN_EXE_eskaks").to_string()
}

fn parse_tsv(path: &str) -> Vec<Vec<String>> {
    let content = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Cannot read {}: {}", path, e));
    content
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(String::from).collect())
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
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status()
        .expect("Failed to spawn binary");
    assert!(status.success(), "eskaks exited with non-zero status: {:?}", status);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_correct_header_and_row_count() {
    let out_prefix = "/tmp/eskaks_test_pairwise_header";
    Command::new(binary_path())
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
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
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
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
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
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
    // seq1_A (ATGGCTGCT) vs seq3_B (ATGATTGCT): 1 nonsyn change GCT→ATT (Ala→Ile)
    // Expected: dN≈0.15439, dS=0.0
    let out_prefix = "/tmp/eskaks_test_pairwise_nonsyn";
    Command::new(binary_path())
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    let row = find_pair(&rows, "seq1_A", "seq3_B")
        .expect("pair seq1_A/seq3_B not found");
    let dn: f64 = row[2].parse().unwrap();
    let ds: f64 = row[3].parse().unwrap();
    assert!((ds - 0.0).abs() < EPSILON, "dS should be 0 for purely nonsynonymous change, got {}", ds);
    assert!((dn - 0.15439).abs() < EPSILON, "dN should be ~0.15439, got {}", dn);
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn pairwise_nei_nonsyn_only_ratio_is_infinity() {
    // When dS=0 and dN>0, ratio should be "inf"
    let out_prefix = "/tmp/eskaks_test_pairwise_inf";
    Command::new(binary_path())
        .args([FASTA, "-o", out_prefix, "--model", "nei", "--workers", "2"])
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
        .args([FASTA, "-o", out_prefix, "--model", "li", "--workers", "2"])
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
        .args([FASTA, "-o", out_prefix, "--model", "li", "--workers", "2"])
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
        .args([FASTA, "-o", out_prefix, "--lineage", "--workers", "2"])
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
        .args([FASTA_GROUPED, "-o", out_prefix, "--group-average", "--workers", "2"])
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
        .args([FASTA_GROUPED, "-o", out_prefix, "--group-average", "--workers", "2"])
        .status().expect("spawn");
    let rows = parse_tsv(&format!("{}_group_avg_dn_ds.tsv", out_prefix));
    let row = rows.iter().find(|r| r.len() >= 6 && r[0] == "A" && r[1] == "A")
        .expect("within-group A row not found");
    let mean: f64 = row[5].parse().unwrap_or(f64::NAN);
    assert!((mean - 0.0).abs() < EPSILON,
        "mean dN/dS within group A should be 0 (all syn), got {}", mean);
    fs::remove_file(format!("{}_group_avg_dn_ds.tsv", out_prefix)).ok();
}
