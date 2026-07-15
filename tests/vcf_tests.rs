/// Integration tests for the VCF subcommand of eskaks.
use std::fs;
use std::process::Command;

const REF_FASTA: &str = "tests/data/test_ref.fasta";
const GFF3: &str = "tests/data/test_annotation.gff3";
const VCF: &str = "tests/data/test_variants.vcf";
const EPSILON: f64 = 1e-4;

fn binary_path() -> String {
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

fn parse_csv(path: &str) -> Vec<Vec<String>> {
    let content = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Cannot read {}: {}", path, e));
    content
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split(',').map(String::from).collect())
        .collect()
}

// ─── Basic VCF subcommand ────────────────────────────────────────────────────

#[test]
fn vcf_exits_ok() {
    let out_prefix = "/tmp/eskaks_test_vcf_basic";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .status()
        .expect("Failed to spawn binary");
    assert!(status.success(), "eskaks vcf exited with non-zero status: {:?}", status);
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

#[test]
fn vcf_produces_correct_header() {
    let out_prefix = "/tmp/eskaks_test_vcf_header";
    Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .status()
        .expect("spawn");
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    assert_eq!(
        rows[0],
        vec!["Gene", "Length_bp", "N_sites", "S_sites", "pN", "pS", "pN/pS",
             "Nonsyn_SNPs", "Syn_SNPs", "Total_SNPs",
             "Chrom", "Start", "End", "Strand", "Exp_N_frac",
             "P_value", "Q_value_BH", "P_Bonferroni"],
        "header mismatch: {:?}", rows[0]
    );
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

#[test]
fn vcf_produces_correct_gene_count() {
    let out_prefix = "/tmp/eskaks_test_vcf_genes";
    Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .status()
        .expect("spawn");
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    // 1 header + 2 gene rows
    assert_eq!(rows.len(), 3, "expected 1 header + 2 genes, got {}", rows.len());
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

#[test]
fn vcf_gene_a_has_correct_snp_classification() {
    // geneA: pos 4 G→A is nonsynonymous (GCT→ACT, Ala→Thr)
    //        pos 6 T→C is synonymous (GCT→GCC, Ala→Ala)
    let out_prefix = "/tmp/eskaks_test_vcf_genea";
    Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .status()
        .expect("spawn");
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    let gene_a = rows.iter().find(|r| r[0] == "geneA")
        .expect("geneA not found in output");
    let nonsyn: f64 = gene_a[7].parse().unwrap();
    let syn: f64 = gene_a[8].parse().unwrap();
    let total: f64 = gene_a[9].parse().unwrap();
    assert!((nonsyn - 1.0).abs() < 0.01, "geneA should have 1 nonsynonymous SNP, got {}", nonsyn);
    assert!((syn - 1.0).abs() < 0.01, "geneA should have 1 synonymous SNP, got {}", syn);
    assert!((total - 2.0).abs() < 0.01, "geneA should have 2 total SNPs, got {}", total);
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

#[test]
fn vcf_gene_b_has_one_syn_snp() {
    // geneB: pos 30 T→C is synonymous (GCT→GCC, Ala→Ala)
    let out_prefix = "/tmp/eskaks_test_vcf_geneb";
    Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .status()
        .expect("spawn");
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    let gene_b = rows.iter().find(|r| r[0] == "geneB")
        .expect("geneB not found in output");
    let nonsyn: f64 = gene_b[7].parse().unwrap();
    let syn: f64 = gene_b[8].parse().unwrap();
    assert!((nonsyn - 0.0).abs() < 0.01, "geneB should have 0 nonsynonymous SNPs, got {}", nonsyn);
    assert!((syn - 1.0).abs() < 0.01, "geneB should have 1 synonymous SNP, got {}", syn);
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

// ─── Output formats ─────────────────────────────────────────────────────────

#[test]
fn vcf_csv_format() {
    let out_prefix = "/tmp/eskaks_test_vcf_csv";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--format", "csv",
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    let rows = parse_csv(&format!("{}_pnps.csv", out_prefix));
    assert_eq!(rows[0][0], "Gene");
    assert!(rows.len() > 1);
    fs::remove_file(format!("{}_pnps.csv", out_prefix)).ok();
}

#[test]
fn vcf_json_format() {
    let out_prefix = "/tmp/eskaks_test_vcf_json";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--format", "json",
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    let content = fs::read_to_string(format!("{}_pnps.json", out_prefix)).unwrap();
    assert!(content.starts_with('['), "JSON should start with [");
    assert!(content.trim().ends_with(']'), "JSON should end with ]");
    assert!(content.contains("\"gene\""), "JSON should have gene key");
    fs::remove_file(format!("{}_pnps.json", out_prefix)).ok();
}

// ─── Genetic code ───────────────────────────────────────────────────────────

#[test]
fn vcf_with_genetic_code_11() {
    let out_prefix = "/tmp/eskaks_test_vcf_gc11";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--genetic-code", "11",
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

// ─── Filtering ──────────────────────────────────────────────────────────────

#[test]
fn vcf_pass_only_filter() {
    let out_prefix = "/tmp/eskaks_test_vcf_pass";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--pass-only",
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

#[test]
fn vcf_min_depth_filter() {
    let out_prefix = "/tmp/eskaks_test_vcf_mindp";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--min-depth", "100",
        ])
        .status()
        .expect("spawn");
    // Should succeed but filter out all SNPs (all have DP<100)
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    // All genes should have 0 SNPs
    for row in &rows[1..] {
        let total: f64 = row[9].parse().unwrap();
        assert!((total - 0.0).abs() < 0.01, "All SNPs should be filtered with min-depth 100");
    }
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

// ─── Plot ───────────────────────────────────────────────────────────────────

#[test]
fn vcf_plot_creates_svg() {
    let out_prefix = "/tmp/eskaks_test_vcf_plot";
    let status = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
            "--plot",
        ])
        .status()
        .expect("spawn");
    assert!(status.success());
    let svg_path = format!("{}_pnps_manhattan.svg", out_prefix);
    let content = fs::read_to_string(&svg_path)
        .unwrap_or_else(|_| panic!("SVG file not created at {}", svg_path));
    assert!(content.contains("<svg"), "SVG should contain <svg tag");
    assert!(content.contains("pN/pS"), "SVG should contain pN/pS in title");
    fs::remove_file(&svg_path).ok();
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

// ─── Error handling ─────────────────────────────────────────────────────────

#[test]
fn vcf_missing_ref_fails() {
    let output = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", "nonexistent.fasta",
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", "/tmp/eskaks_vcf_missing_ref",
        ])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
}

#[test]
fn vcf_missing_gff_fails() {
    let output = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", "nonexistent.gff3",
            "--vcf", VCF,
            "-o", "/tmp/eskaks_vcf_missing_gff",
        ])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
}

#[test]
fn vcf_missing_vcf_fails() {
    let output = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", "nonexistent.vcf",
            "-o", "/tmp/eskaks_vcf_missing_vcf",
        ])
        .output()
        .expect("spawn");
    assert!(!output.status.success());
}

// ─── Genome-wide (pooled) summary ────────────────────────────────────────────

/// Extract the float following a `label` on its line in captured stderr.
fn stderr_value_after(stderr: &str, label: &str) -> f64 {
    let line = stderr
        .lines()
        .find(|l| l.contains(label))
        .unwrap_or_else(|| panic!("no line containing {:?} in:\n{}", label, stderr));
    line.rsplit(':')
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap_or_else(|e| panic!("cannot parse value from {:?}: {}", line, e))
}

#[test]
fn vcf_summary_reports_genome_wide_pnps() {
    let out_prefix = "/tmp/eskaks_test_vcf_gw";
    let output = Command::new(binary_path())
        .args([
            "vcf",
            "--ref", REF_FASTA,
            "--gff", GFF3,
            "--vcf", VCF,
            "-o", out_prefix,
        ])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Genome-wide (pooled)"), "missing genome-wide block:\n{}", stderr);

    // The pooled ratio printed in the summary must equal the ratio recomputed
    // by pooling the per-gene TSV columns (N_sites, S_sites, Nonsyn, Syn).
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    let (mut n_sites, mut s_sites, mut nonsyn, mut syn) = (0.0, 0.0, 0.0, 0.0);
    for row in &rows[1..] {
        n_sites += row[2].parse::<f64>().unwrap();
        s_sites += row[3].parse::<f64>().unwrap();
        nonsyn += row[7].parse::<f64>().unwrap();
        syn += row[8].parse::<f64>().unwrap();
    }
    let expected = (nonsyn / n_sites) / (syn / s_sites);
    let reported = stderr_value_after(&stderr, "Overall pN/pS:");
    assert!(
        (reported - expected).abs() < EPSILON,
        "genome-wide pN/pS mismatch: summary={}, recomputed={}",
        reported, expected
    );
    // geneA (1N,1S) + geneB (0N,1S): purifying overall.
    assert!(stderr.contains("purifying"), "expected purifying call:\n{}", stderr);
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
}

// ─── Mutation-spectrum-weighted site counting (--kappa) ──────────────────────

/// Run the vcf subcommand with an optional --kappa and return (sum N_sites,
/// sum S_sites) from the per-gene TSV.
fn run_and_sum_sites(out_prefix: &str, kappa: Option<&str>) -> (f64, f64) {
    let mut args = vec![
        "vcf", "--ref", REF_FASTA, "--gff", GFF3, "--vcf", VCF, "-o", out_prefix,
    ];
    if let Some(k) = kappa {
        args.push("--kappa");
        args.push(k);
    }
    let status = Command::new(binary_path()).args(&args).status().expect("spawn");
    assert!(status.success());
    let rows = parse_tsv(&format!("{}_pnps.tsv", out_prefix));
    let (mut n, mut s) = (0.0, 0.0);
    for row in &rows[1..] {
        n += row[2].parse::<f64>().unwrap();
        s += row[3].parse::<f64>().unwrap();
    }
    fs::remove_file(format!("{}_pnps.tsv", out_prefix)).ok();
    (n, s)
}

#[test]
fn vcf_kappa_1_matches_default() {
    // --kappa 1 must reproduce the equal-rates default exactly.
    let (n_def, s_def) = run_and_sum_sites("/tmp/eskaks_kappa_def", None);
    let (n_k1, s_k1) = run_and_sum_sites("/tmp/eskaks_kappa_one", Some("1"));
    assert!((n_def - n_k1).abs() < EPSILON, "N: default {} vs kappa=1 {}", n_def, n_k1);
    assert!((s_def - s_k1).abs() < EPSILON, "S: default {} vs kappa=1 {}", s_def, s_k1);
}

#[test]
fn vcf_kappa_raises_s_lowers_n_conserves_total() {
    // Transition bias (kappa>1) must raise synonymous sites, lower
    // nonsynonymous sites, and conserve the total (each position is 1 site).
    let (n1, s1) = run_and_sum_sites("/tmp/eskaks_kappa_1", Some("1"));
    let (n5, s5) = run_and_sum_sites("/tmp/eskaks_kappa_5", Some("5"));
    assert!(s5 > s1, "S should rise with kappa: {} -> {}", s1, s5);
    assert!(n5 < n1, "N should fall with kappa: {} -> {}", n1, n5);
    assert!(((n1 + s1) - (n5 + s5)).abs() < EPSILON, "total sites must be conserved");
}

#[test]
fn vcf_invalid_kappa_fails() {
    for bad in ["0", "-2"] {
        let output = Command::new(binary_path())
            .args([
                "vcf", "--ref", REF_FASTA, "--gff", GFF3, "--vcf", VCF,
                "-o", "/tmp/eskaks_kappa_bad", "--kappa", bad,
            ])
            .output()
            .expect("spawn");
        assert!(!output.status.success(), "--kappa {} should be rejected", bad);
    }
}

// ─── FASTA subcommand still works ───────────────────────────────────────────

#[test]
fn fasta_subcommand_still_works() {
    let out_prefix = "/tmp/eskaks_test_vcf_fasta_compat";
    let status = Command::new(binary_path())
        .args([
            "fasta",
            "tests/data/synthetic.fasta",
            "-o", out_prefix,
            "--model", "nei",
            "--workers", "2",
        ])
        .status()
        .expect("spawn");
    assert!(status.success(), "fasta subcommand should still work");
    let rows = parse_tsv(&format!("{}_pairwise_results.tsv", out_prefix));
    assert!(rows.len() > 1, "should produce output");
    fs::remove_file(format!("{}_pairwise_results.tsv", out_prefix)).ok();
}

#[test]
fn list_codes_still_works_at_top_level() {
    let output = Command::new(binary_path())
        .args(["--list-codes"])
        .output()
        .expect("spawn");
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Standard"), "should list Standard code");
}
