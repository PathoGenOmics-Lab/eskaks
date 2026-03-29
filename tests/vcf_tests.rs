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
             "Nonsyn_SNPs", "Syn_SNPs", "Total_SNPs"],
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
    let nonsyn: u32 = gene_a[7].parse().unwrap();
    let syn: u32 = gene_a[8].parse().unwrap();
    let total: u32 = gene_a[9].parse().unwrap();
    assert_eq!(nonsyn, 1, "geneA should have 1 nonsynonymous SNP, got {}", nonsyn);
    assert_eq!(syn, 1, "geneA should have 1 synonymous SNP, got {}", syn);
    assert_eq!(total, 2, "geneA should have 2 total SNPs, got {}", total);
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
    let nonsyn: u32 = gene_b[7].parse().unwrap();
    let syn: u32 = gene_b[8].parse().unwrap();
    assert_eq!(nonsyn, 0, "geneB should have 0 nonsynonymous SNPs, got {}", nonsyn);
    assert_eq!(syn, 1, "geneB should have 1 synonymous SNP, got {}", syn);
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
        let total: u32 = row[9].parse().unwrap();
        assert_eq!(total, 0, "All SNPs should be filtered with min-depth 100");
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
