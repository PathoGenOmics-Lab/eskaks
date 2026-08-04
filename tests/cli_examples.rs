//! End-to-end CLI integration tests that run the compiled `eskaks` binary over the
//! bundled `examples/` datasets (the same commands the getting-started tutorial
//! documents). They lock in the tool's observable behaviour — output file names,
//! TSV/CSV/JSON/SVG/HTML structure, hand-checkable golden values, flag effects,
//! determinism, and error handling — so a regression in the wiring is caught, not
//! just in a unit.
//!
//! Every golden number here was captured from an actual run and independently
//! cross-checked (site counts via Nei-Gojobori, MK partition via
//! Pn+Dn==Nonsyn / Ps+Ds==Syn, 6 strains -> 15 pairs, 12 toy genes).

use std::path::PathBuf;
use std::process::{Command, Output};

// Example inputs, resolved relative to the crate root (we spawn with that cwd).
const FASTA_GENES: &str = "examples/genes.fasta";
const REF: &str = "examples/toy_genome/reference.fasta";
const GFF: &str = "examples/toy_genome/genes.gff3";
const VCF: &str = "examples/toy_genome/variants.vcf";
const DIVERGENCE: &str = "examples/toy_genome/divergence.tsv";

/// Tolerance for ratios/frequencies (output carries 6 decimals).
const EPS: f64 = 1e-4;
/// Tolerance for fractional site counts (output carries 4 decimals).
const EPS_SITES: f64 = 1e-3;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_eskaks"))
}

/// A throwaway output directory; files are written as `<prefix>_*`.
struct Run {
    _dir: tempfile::TempDir,
    prefix: String,
}

fn new_run() -> Run {
    let dir = tempfile::tempdir().expect("tempdir");
    let prefix = dir.path().join("run").to_str().unwrap().to_string();
    Run { _dir: dir, prefix }
}

/// Run the binary from the crate root (so `examples/...` resolves) and return the raw output.
fn run(args: &[&str]) -> Output {
    Command::new(bin())
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(args)
        .output()
        .expect("failed to spawn eskaks")
}

/// Run and assert success, returning the output (for stderr inspection).
fn run_ok(args: &[&str]) -> Output {
    let out = run(args);
    assert!(
        out.status.success(),
        "command failed: eskaks {}\nstderr:\n{}",
        args.join(" "),
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
}

fn split_rows(path: &str, sep: char) -> Vec<Vec<String>> {
    read(path)
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split(sep).map(|c| c.to_string()).collect())
        .collect()
}

fn tsv(path: &str) -> Vec<Vec<String>> {
    split_rows(path, '\t')
}

/// Find the data row whose first cell equals `name`.
fn row<'a>(rows: &'a [Vec<String>], name: &str) -> &'a Vec<String> {
    rows.iter()
        .find(|r| r.first().map(String::as_str) == Some(name))
        .unwrap_or_else(|| panic!("row {name:?} not found"))
}

/// Find a pairwise row by BOTH sequence ids (robust to output ordering).
fn find_pair<'a>(rows: &'a [Vec<String>], a: &str, b: &str) -> &'a Vec<String> {
    rows.iter()
        .find(|r| r.len() >= 2 && r[0] == a && r[1] == b)
        .unwrap_or_else(|| panic!("pair {a}/{b} not found"))
}

/// Find a group-average row by its two group labels, in either column order.
fn find_group<'a>(rows: &'a [Vec<String>], a: &str, b: &str) -> &'a Vec<String> {
    rows.iter()
        .find(|r| r.len() >= 2 && ((r[0] == a && r[1] == b) || (r[0] == b && r[1] == a)))
        .unwrap_or_else(|| panic!("group pair {a}/{b} not found"))
}

/// Distinct group labels appearing in a group-average table (columns 0 and 1).
fn group_labels(rows: &[Vec<String>]) -> std::collections::BTreeSet<String> {
    rows[1..].iter().flat_map(|r| [r[0].clone(), r[1].clone()]).collect()
}

fn f(cell: &str) -> f64 {
    cell.parse::<f64>()
        .unwrap_or_else(|_| panic!("not a float: {cell:?}"))
}

/// Grab the first numeric token appearing after `label` on its stderr line.
/// Robust to column padding and trailing annotations like `(100.0%)`.
fn num_after(err: &str, label: &str) -> f64 {
    let line = err
        .lines()
        .find(|l| l.contains(label))
        .unwrap_or_else(|| panic!("no stderr line containing {label:?}:\n{err}"));
    let after = &line[line.find(label).unwrap() + label.len()..];
    after
        .split_whitespace()
        .find_map(|t| t.trim_matches(|c: char| c == '%' || c == '(' || c == ')' || c == ',').parse::<f64>().ok())
        .unwrap_or_else(|| panic!("no number after {label:?} in {line:?}"))
}

fn assert_svg_well_formed(path: &str) {
    let s = read(path);
    assert!(s.starts_with("<?xml"), "{path}: missing XML prolog");
    assert!(s.contains("<svg "), "{path}: missing <svg>");
    assert!(s.trim_end().ends_with("</svg>"), "{path}: not closed");
    assert!(!s.contains("NaN"), "{path}: NaN leaked into SVG");
}

fn assert_report_html(path: &str) {
    let s = read(path);
    assert!(s.contains("<!DOCTYPE html>"), "{path}: missing doctype");
    assert!(s.contains("<title>eskaks report</title>"), "{path}: missing title");
    assert!(s.contains("<script"), "{path}: report has no embedded script");
    assert!(s.len() > 3000, "{path}: report suspiciously small ({} bytes)", s.len());
}

// ─────────────────────────── eskaks fasta (examples/genes.fasta) ───────────────────────────

#[test]
fn fasta_examples_pairwise_golden() {
    let r = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix]);
    let rows = tsv(&format!("{}_pairwise_results.tsv", r.prefix));
    assert_eq!(rows[0], ["Seq1", "Seq2", "dN", "dS", "dN/dS"]);
    // 6 strains -> C(6,2) = 15 unordered pairs.
    assert_eq!(rows.len() - 1, 15, "expected 15 pairs");

    // Golden values captured from a real run (Nei-Gojobori model, exclude-nonsense
    // site counting), looked up by BOTH ids so a future ordering change surfaces
    // as "pair not found".
    let ab = find_pair(&rows, "strain_A", "strain_B");
    assert!((f(&ab[2]) - 0.022990).abs() < EPS, "A/B dN {}", ab[2]);
    assert!((f(&ab[3]) - 0.164075).abs() < EPS, "A/B dS {}", ab[3]);
    assert!((f(&ab[4]) - 0.140120).abs() < EPS, "A/B dN/dS {}", ab[4]);
    let ac = find_pair(&rows, "strain_A", "strain_C");
    assert!((f(&ac[4]) - 0.077805).abs() < EPS, "A/C dN/dS {}", ac[4]);

    // Every pair is under purifying selection (dN/dS < 1).
    for p in &rows[1..] {
        assert!(f(&p[4]) < 1.0, "pair {}/{} not purifying: {}", p[0], p[1], p[4]);
    }
}

#[test]
fn fasta_examples_plot_and_report_files() {
    let r = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix, "--plot", "--report"]);
    assert_svg_well_formed(&format!("{}_dnds_histogram.svg", r.prefix));
    assert_report_html(&format!("{}_report.html", r.prefix));
}

#[test]
fn fasta_examples_deterministic() {
    let a = new_run();
    let b = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &a.prefix]);
    run_ok(&["fasta", FASTA_GENES, "-o", &b.prefix]);
    assert_eq!(
        read(&format!("{}_pairwise_results.tsv", a.prefix)),
        read(&format!("{}_pairwise_results.tsv", b.prefix)),
        "two identical fasta runs must produce byte-identical output"
    );
}

#[test]
fn fasta_examples_li_model() {
    let r = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix, "--model", "li"]);
    let rows = tsv(&format!("{}_pairwise_results.tsv", r.prefix));
    assert_eq!(rows[0], ["Seq1", "Seq2", "dN(Ka)", "dS(Ks)", "dN/dS"]);
    let ab = find_pair(&rows, "strain_A", "strain_B");
    assert!((f(&ab[2]) - 0.024933).abs() < EPS, "Ka {}", ab[2]);
    assert!((f(&ab[3]) - 0.164551).abs() < EPS, "Ks {}", ab[3]);
    // The purifying-selection headline must also hold on the Li (1993) code path.
    for p in &rows[1..] {
        assert!(f(&p[4]) < 1.0, "li pair {}/{} not purifying: {}", p[0], p[1], p[4]);
    }
}

#[test]
fn fasta_examples_summary_stderr() {
    let r = new_run();
    let out = run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix, "--summary"]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("dN/dS Summary"), "no summary block:\n{err}");
    // Pin the counts to their labelled lines (a bare "15" also appears in float values).
    assert_eq!(num_after(&err, "Total pairs:") as usize, 15);
    assert_eq!(num_after(&err, "Valid pairs:") as usize, 15);
}

#[test]
fn fasta_examples_csv_and_json_formats() {
    let r = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix, "--format", "csv"]);
    let csv = read(&format!("{}_pairwise_results.csv", r.prefix));
    assert!(csv.lines().next().unwrap().contains(','), "csv not comma-separated");

    let r2 = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r2.prefix, "--format", "json"]);
    let json = read(&format!("{}_pairwise_results.json", r2.prefix));
    let trimmed = json.trim();
    // Lightweight structural check (no JSON crate in the dependency set).
    assert!(trimmed.starts_with('['), "JSON should be a top-level array");
    assert!(trimmed.ends_with(']'), "JSON array not closed");
    assert!(json.contains("\"seq1\":\"strain_A\""), "missing expected record");
    assert!(json.contains("\"dN_dS\":"), "missing dN_dS key");
    assert_eq!(json.matches("\"seq1\":").count(), 15, "expected 15 pair records");
}

#[test]
fn fasta_examples_neutrality_and_bootstrap_reproducible() {
    let r = new_run();
    run_ok(&[
        "fasta", FASTA_GENES, "-o", &r.prefix,
        "--neutrality", "--bootstrap", "50", "--seed", "42",
    ]);
    let tests = read(&format!("{}_pairwise_tests.tsv", r.prefix));
    assert!(tests.lines().count() > 1, "neutrality output empty");
    let boot_a = read(&format!("{}_pairwise_bootstrap.tsv", r.prefix));

    // Same seed => byte-identical bootstrap CIs.
    let r2 = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r2.prefix, "--bootstrap", "50", "--seed", "42"]);
    let boot_b = read(&format!("{}_pairwise_bootstrap.tsv", r2.prefix));
    assert_eq!(boot_a, boot_b, "bootstrap must be reproducible for a fixed seed");
}

#[test]
fn fasta_examples_window_mode() {
    let r = new_run();
    run_ok(&[
        "fasta", FASTA_GENES, "-o", &r.prefix,
        "--window-size", "20", "--window-step", "10",
    ]);
    let rows = tsv(&format!("{}_pairwise_windows.tsv", r.prefix));
    assert_eq!(
        rows[0],
        ["Seq1", "Seq2", "Window_Start", "Window_End", "dN", "dS", "dN/dS"]
    );
    // 15 pairs x 5 windows (starts 1,11,21,31,41 over a 60-codon alignment).
    assert_eq!(rows.len() - 1, 75, "expected 75 window rows");
    let w = rows[1..]
        .iter()
        .find(|r| r[0] == "strain_A" && r[1] == "strain_B" && r[2] == "1")
        .expect("A/B window [1,20]");
    assert_eq!(w[3], "20");
    assert!((f(&w[4]) - 0.023104).abs() < EPS, "window dN {}", w[4]);
    assert!((f(&w[5]) - 0.215020).abs() < EPS, "window dS {}", w[5]);
    assert!((f(&w[6]) - 0.107451).abs() < EPS, "window dN/dS {}", w[6]);
}

#[test]
fn fasta_examples_lineage() {
    let r = new_run();
    run_ok(&["fasta", FASTA_GENES, "-o", &r.prefix, "--lineage"]);
    let rows = tsv(&format!("{}_lineage_summary.tsv", r.prefix));
    assert_eq!(
        rows[0],
        ["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"]
    );
    assert_eq!(rows.len() - 1, 6, "one row per genome");
    // All six IDs split on '_' to the single lineage "strain".
    let a = row(&rows, "strain_A");
    assert_eq!(a[1], "strain");
    assert!((f(&a[2]) - 0.026091).abs() < EPS, "Mean_dN {}", a[2]);
    assert!((f(&a[3]) - 0.155925).abs() < EPS, "Mean_dS {}", a[3]);
    assert!((f(&a[4]) - 0.167333).abs() < EPS, "dN/dS {}", a[4]);
}

// ─────────────────────────── eskaks vcf (examples/toy_genome) ───────────────────────────

/// The canonical toy-genome invocation, minus the report/plot side outputs.
fn vcf_core(prefix: &str, extra: &[&str]) -> Output {
    let mut args = vec![
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", VCF, "--genetic-code", "11", "-o", prefix,
    ];
    args.extend_from_slice(extra);
    run_ok(&args)
}

const PNPS_HEADER: [&str; 22] = [
    "Gene", "Length_bp", "N_sites", "S_sites", "pN", "pS", "pN/pS", "Nonsyn_SNPs",
    "Syn_SNPs", "Total_SNPs", "Chrom", "Start", "End", "Strand", "Exp_N_frac",
    "P_value", "Q_value_BH", "P_Bonferroni",
    "pN/pS_lo", "pN/pS_hi", "P_GC", "Q_GC_BH",
];

#[test]
fn vcf_examples_pnps_golden() {
    let r = new_run();
    vcf_core(&r.prefix, &[]);
    let rows = tsv(&format!("{}_pnps.tsv", r.prefix));
    assert_eq!(rows[0], PNPS_HEADER);
    assert_eq!(rows.len() - 1, 12, "12 toy genes");

    // gene01 golden — hand-checkable counts + Nei-Gojobori site counts + the full
    // neutrality-test column set (Exp_N_frac / P / BH-q / Bonferroni).
    let g = row(&rows, "gene01");
    assert_eq!(g[1], "399", "Length_bp");
    assert!((f(&g[2]) - 295.7798).abs() < EPS_SITES, "N_sites {}", g[2]);
    assert!((f(&g[3]) - 100.2202).abs() < EPS_SITES, "S_sites {}", g[3]);
    assert!((f(&g[6]) - 0.225889).abs() < EPS, "pN/pS {}", g[6]);
    assert!((f(&g[7]) - 10.0).abs() < EPS, "Nonsyn {}", g[7]);
    assert!((f(&g[8]) - 15.0).abs() < EPS, "Syn {}", g[8]);
    assert!((f(&g[9]) - 25.0).abs() < EPS, "Total {}", g[9]);
    assert_eq!(g[10], "chr1");
    assert_eq!(g[13], "+", "strand");
    // N_sites + S_sites reconstruct the CDS minus the excluded final codon (399-3=396).
    assert!((f(&g[2]) + f(&g[3]) - 396.0).abs() < 1e-2, "N+S should be 396");
    // Statistical columns (the recently added neutrality pipeline).
    assert!((f(&g[14]) - 0.746919).abs() < EPS, "Exp_N_frac {}", g[14]);
    assert!((f(&g[15]) - 4.967e-4).abs() < 1e-6, "P_value {}", g[15]);
    assert!((f(&g[16]) - 0.001987).abs() < 1e-5, "Q_value_BH {}", g[16]);
    assert!((f(&g[17]) - 0.005961).abs() < 1e-5, "P_Bonferroni {}", g[17]);

    // The two deliberately repetitive-named genes are present in the per-gene table.
    let names: Vec<&str> = rows[1..].iter().map(|r| r[0].as_str()).collect();
    assert!(names.contains(&"PPE_toy1"), "PPE_toy1 missing: {names:?}");
    assert!(names.contains(&"PE_PGRS_toy2"), "PE_PGRS_toy2 missing: {names:?}");
}

#[test]
fn vcf_examples_report_and_plots() {
    let r = new_run();
    vcf_core(&r.prefix, &["--report", "--plot", "--divergence", DIVERGENCE]);
    assert_svg_well_formed(&format!("{}_pnps_manhattan.svg", r.prefix));
    assert_svg_well_formed(&format!("{}_pvalue_manhattan.svg", r.prefix));
    assert_report_html(&format!("{}_report.html", r.prefix));
}

#[test]
fn vcf_examples_deterministic() {
    let a = new_run();
    let b = new_run();
    vcf_core(&a.prefix, &[]);
    vcf_core(&b.prefix, &[]);
    assert_eq!(
        read(&format!("{}_pnps.tsv", a.prefix)),
        read(&format!("{}_pnps.tsv", b.prefix)),
        "two identical vcf runs must produce byte-identical pN/pS output"
    );
}

#[test]
fn vcf_examples_multisample_merge_matches_single() {
    // Passing the same VCF twice exercises the multi-sample MERGE path (dedup
    // positions, recompute AF as fraction of samples) — the path that recently had
    // a determinism bug. Merging a sample with itself must equal the single-sample
    // run and be reproducible across runs.
    let single = new_run();
    vcf_core(&single.prefix, &[]);
    let merged = new_run();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", VCF, "--vcf", VCF,
        "--genetic-code", "11", "-o", &merged.prefix,
    ]);
    let merged2 = new_run();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", VCF, "--vcf", VCF,
        "--genetic-code", "11", "-o", &merged2.prefix,
    ]);
    let single_tsv = read(&format!("{}_pnps.tsv", single.prefix));
    let merged_tsv = read(&format!("{}_pnps.tsv", merged.prefix));
    assert_eq!(single_tsv, merged_tsv, "merge-with-self must equal the single-sample run");
    assert_eq!(
        merged_tsv,
        read(&format!("{}_pnps.tsv", merged2.prefix)),
        "the merge path must be deterministic across runs"
    );
}

#[test]
fn vcf_examples_vcf_list_matches_repeated_flag() {
    // A --vcf-list file (one path per line) must feed the same merge path as
    // repeated --vcf flags.
    let list = new_run();
    let list_file = format!("{}_samples.txt", list.prefix);
    std::fs::write(&list_file, format!("{VCF}\n{VCF}\n")).unwrap();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf-list", &list_file,
        "--genetic-code", "11", "-o", &list.prefix,
    ]);
    let repeated = new_run();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", VCF, "--vcf", VCF,
        "--genetic-code", "11", "-o", &repeated.prefix,
    ]);
    assert_eq!(
        read(&format!("{}_pnps.tsv", list.prefix)),
        read(&format!("{}_pnps.tsv", repeated.prefix)),
        "--vcf-list must match repeated --vcf"
    );
}

#[test]
fn vcf_examples_csv_and_json_formats() {
    let c = new_run();
    vcf_core(&c.prefix, &["--format", "csv"]);
    let rows = split_rows(&format!("{}_pnps.csv", c.prefix), ',');
    assert_eq!(rows[0][0], "Gene");
    assert_eq!(rows[0].len(), 22, "csv should have all 22 columns");
    let g = row(&rows, "gene01");
    assert_eq!(g[1], "399");
    assert!((f(&g[6]) - 0.225889).abs() < EPS, "csv gene01 pN/pS {}", g[6]);

    let j = new_run();
    vcf_core(&j.prefix, &["--format", "json"]);
    let json = read(&format!("{}_pnps.json", j.prefix));
    assert!(json.trim().starts_with('['), "vcf json should be an array");
    assert!(json.trim().ends_with(']'), "vcf json not closed");
    assert!(json.contains("\"gene\":\"gene01\""), "missing gene01 record");
    assert_eq!(json.matches("\"gene\":").count(), 12, "expected 12 gene records");
}

#[test]
fn vcf_examples_mk_table_with_invariant() {
    let r = new_run();
    vcf_core(&r.prefix, &["--mk"]);
    let rows = tsv(&format!("{}_mk.tsv", r.prefix));
    assert_eq!(
        rows[0],
        ["Gene", "Chrom", "Start", "End", "Strand", "Dn", "Ds", "Pn", "Ps", "NI", "alpha", "Fisher_p", "Fisher_q_BH"]
    );
    assert_eq!(rows.len() - 1, 12, "one MK row per gene");

    // gene01 golden partition and derived statistics.
    let g = row(&rows, "gene01");
    let (dn, ds, pn, ps) = (f(&g[5]), f(&g[6]), f(&g[7]), f(&g[8]));
    assert_eq!((dn, ds, pn, ps), (1.0, 1.0, 9.0, 14.0), "MK partition");
    assert!((f(&g[9]) - 0.642857).abs() < EPS, "NI {}", g[9]);
    assert!((f(&g[10]) - 0.357143).abs() < EPS, "alpha {}", g[10]);
    // First-principles cross-check: MK's nonsyn/syn totals must equal the pN/pS
    // SNP counts for the same gene (Pn+Dn == Nonsyn, Ps+Ds == Syn).
    assert_eq!(pn + dn, 10.0, "Pn+Dn must equal Nonsyn_SNPs (10)");
    assert_eq!(ps + ds, 15.0, "Ps+Ds must equal Syn_SNPs (15)");
}

#[test]
fn vcf_examples_min_snps_drops_low_count_genes() {
    let r = new_run();
    vcf_core(&r.prefix, &["--min-snps", "20"]);
    let rows = tsv(&format!("{}_pnps.tsv", r.prefix));
    // gene03 has 15 SNPs (< 20) and must be dropped from the per-gene table.
    assert_eq!(rows.len() - 1, 9, "min-snps 20 should leave 9 genes");
    let names: Vec<&str> = rows[1..].iter().map(|r| r[0].as_str()).collect();
    assert!(!names.contains(&"gene03"), "gene03 (15 SNPs) should be dropped");
}

#[test]
fn vcf_examples_kappa_raises_synonymous_sites() {
    let k1 = new_run();
    let k2 = new_run();
    vcf_core(&k1.prefix, &[]);
    vcf_core(&k2.prefix, &["--kappa", "2"]);
    let g1 = row(&tsv(&format!("{}_pnps.tsv", k1.prefix)), "gene01").clone();
    let g2 = row(&tsv(&format!("{}_pnps.tsv", k2.prefix)), "gene01").clone();
    // kappa=1 baseline vs kappa=2: transition up-weighting moves site mass from N to
    // S while conserving N+S=396.
    assert!((f(&g1[3]) - 100.2202).abs() < EPS_SITES, "k1 S_sites {}", g1[3]);
    assert!((f(&g2[2]) - 287.2636).abs() < EPS_SITES, "k2 N_sites {}", g2[2]);
    assert!((f(&g2[3]) - 108.7364).abs() < EPS_SITES, "k2 S_sites {}", g2[3]);
    assert!((f(&g2[2]) + f(&g2[3]) - 396.0).abs() < 1e-2, "k2 N+S should stay 396");
}

#[test]
fn vcf_examples_min_af_filters_snps() {
    let r = new_run();
    vcf_core(&r.prefix, &["--min-af", "0.5"]);
    let g = row(&tsv(&format!("{}_pnps.tsv", r.prefix)), "gene01").clone();
    // Dropping AF<0.5 variants nearly halves gene01's SNP load: 25 -> 12 (6 N, 6 S).
    assert!((f(&g[7]) - 6.0).abs() < EPS, "min-af Nonsyn {}", g[7]);
    assert!((f(&g[8]) - 6.0).abs() < EPS, "min-af Syn {}", g[8]);
    assert!((f(&g[9]) - 12.0).abs() < EPS, "min-af Total {}", g[9]);
}

#[test]
fn vcf_examples_min_depth_filters_snps() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &["--min-depth", "50"]);
    let err = String::from_utf8_lossy(&out.stderr);
    // Example DP spans ~20-90; requiring DP>=50 drops low-depth SNPs from the pool.
    assert_eq!(num_after(&err, "Total synonymous:") as usize, 71);
    assert_eq!(num_after(&err, "Total nonsynonymous:") as usize, 91);
}

#[test]
fn vcf_examples_max_af_drops_fixed_snps() {
    let r = new_run();
    vcf_core(&r.prefix, &["--max-af", "0.99"]);
    let rows = tsv(&format!("{}_pnps.tsv", r.prefix));
    let total: f64 = rows[1..].iter().map(|g| f(&g[9])).sum();
    // Excluding fixed variants (AF=1.0) removes 14 SNPs: 264 -> 250 across all genes.
    assert!((total - 250.0).abs() < EPS, "sum Total_SNPs {total}");
}

#[test]
fn vcf_examples_af_weighted_changes_ratio_and_label() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &["--af-weighted"]);
    let err = String::from_utf8_lossy(&out.stderr);
    // AF weighting computes πN/πS: gene01's ratio shifts from 0.225889 to 0.277661.
    let g = row(&tsv(&format!("{}_pnps.tsv", r.prefix)), "gene01").clone();
    assert!((f(&g[6]) - 0.277661).abs() < EPS, "af-weighted pN/pS {}", g[6]);
    assert!(err.contains("πN/πS"), "summary should switch to πN/πS:\n{err}");
}

#[test]
fn vcf_examples_exclude_repetitive_pools_core_only() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &["--exclude-repetitive"]);
    let err = String::from_utf8_lossy(&out.stderr);
    // The two repetitive genes leave the tested family and the pooled estimate...
    assert_eq!(num_after(&err, "Genes tested:") as usize, 10, "should test 10 core genes\n{err}");
    assert!((num_after(&err, "Overall pN/pS:") - 0.362316).abs() < EPS, "core-only pooled pN/pS\n{err}");
    // ...but are still LISTED in the per-gene table, flagged.
    let rows = tsv(&format!("{}_pnps.tsv", r.prefix));
    let names: Vec<&str> = rows[1..].iter().map(|g| g[0].as_str()).collect();
    assert!(names.contains(&"PPE_toy1") && names.contains(&"PE_PGRS_toy2"));
}

#[test]
fn vcf_examples_genomic_control_deflates_significance() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &["--genomic-control"]);
    let err = String::from_utf8_lossy(&out.stderr);
    // Dividing chi-square by the inflation factor collapses significance on this
    // clonal-like toy set: 7 significant genes -> 0 under BH-FDR·GC.
    assert!(err.contains("BH-FDR·GC"), "expected GC-corrected significance label:\n{err}");
    assert_eq!(num_after(&err, "Significant genes:") as usize, 0, "GC should deflate to 0\n{err}");
}

#[test]
fn vcf_examples_bootstrap_reproducible() {
    let a = new_run();
    let b = new_run();
    let out_a = vcf_core(&a.prefix, &["--bootstrap", "200", "--seed", "7"]);
    let out_b = vcf_core(&b.prefix, &["--bootstrap", "200", "--seed", "7"]);
    let ea = String::from_utf8_lossy(&out_a.stderr);
    let eb = String::from_utf8_lossy(&out_b.stderr);
    let ci = |e: &str| e.lines().find(|l| l.contains("CI")).unwrap_or("").to_string();
    assert_eq!(ci(&ea), ci(&eb), "bootstrap CI must be reproducible for a fixed seed");
    assert!(ci(&ea).contains("0.294088") && ci(&ea).contains("0.684979"), "CI value: {}", ci(&ea));
}

#[test]
fn vcf_examples_genome_wide_diversity_independent_of_min_snps() {
    // Regression: the genome-wide (pooled) diversity headline must pool over ALL
    // genes (like the pooled pN/pS), not the --min-snps-filtered subset, so
    // pi / theta_W / Tajima's D never silently depend on the --min-snps threshold.
    const MS: &str = "examples/toy_genome/variants_multisample.vcf";
    let base = [
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", MS, "--genetic-code", "11", "--diversity",
    ];

    let a = new_run();
    let mut args_a = base.to_vec();
    args_a.extend(["-o", a.prefix.as_str(), "--min-snps", "0"]);
    let ea = String::from_utf8_lossy(&run_ok(&args_a).stderr).into_owned();

    let b = new_run();
    let mut args_b = base.to_vec();
    args_b.extend(["-o", b.prefix.as_str(), "--min-snps", "20"]);
    let eb = String::from_utf8_lossy(&run_ok(&args_b).stderr).into_owned();

    for label in ["Segregating coding SNPs:", "piN/piS:", "Tajima's D:"] {
        assert_eq!(
            num_after(&ea, label),
            num_after(&eb, label),
            "genome-wide diversity '{label}' must not change with --min-snps\n\
             --min-snps 0:\n{ea}\n--min-snps 20:\n{eb}"
        );
    }

    // Sanity: --min-snps 20 genuinely filters the per-gene diversity table, so the
    // two runs differ; the headline just must not depend on that filtering.
    let rows0 = tsv(&format!("{}_diversity.tsv", a.prefix));
    let rows20 = tsv(&format!("{}_diversity.tsv", b.prefix));
    assert!(
        rows20.len() < rows0.len(),
        "--min-snps 20 should drop low-count genes from the per-gene diversity table \
         (rows0={}, rows20={})",
        rows0.len(),
        rows20.len()
    );
}

#[test]
fn vcf_report_never_leaks_raw_markup_from_a_gene_name() {
    // A gene name comes from an untrusted GFF, so it must never reach the report HTML
    // as live markup (a stored XSS the interactive tooltips previously allowed via an
    // attribute round-trip). Inject an HTML/script payload as a gene name and assert
    // the raw dangerous form never appears in the output.
    let gff = read("examples/toy_genome/genes.gff3")
        .replacen("gene=gene01", "gene=<img src=x onerror=alert(1)><script>alert(2)</script>", 1);
    let gff_path = std::env::temp_dir().join("eskaks_xss_test.gff3");
    std::fs::write(&gff_path, gff).unwrap();

    let r = new_run();
    run_ok(&[
        "vcf",
        "--ref", "examples/toy_genome/reference.fasta",
        "--gff", gff_path.to_str().unwrap(),
        "--vcf", "examples/toy_genome/variants_multisample.vcf",
        "--genetic-code", "11", "--report", "--mk",
        "-o", r.prefix.as_str(),
    ]);
    let html = read(&format!("{}_report.html", r.prefix));

    // The payload's markup must appear only in escaped form (Rust `esc()` turns the
    // name's `<`/`>` into `\uXXXX` in the JSON block). The bare tags below are unique
    // to the payload, so the report's own legitimate `<script>` tags don't false-match.
    for raw in ["<img src=x onerror", "<script>alert(2)"] {
        assert!(
            !html.contains(raw),
            "raw markup {raw:?} from a gene name leaked into the report HTML (XSS risk)"
        );
    }
}

#[test]
fn vcf_examples_summary_reports_pooled_and_significant() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &[]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("pN/pS Summary"), "no summary:\n{err}");
    assert_eq!(num_after(&err, "Genes analyzed:") as usize, 12);
    assert_eq!(num_after(&err, "Significant genes:") as usize, 7);
    assert!(err.contains("purifying selection"), "expected purifying pooled signal:\n{err}");
}

#[test]
fn vcf_examples_divergence_without_report_warns() {
    let r = new_run();
    let out = vcf_core(&r.prefix, &["--divergence", DIVERGENCE]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("only used by the interactive report"),
        "expected a warning that --divergence needs --report:\n{err}"
    );
}

// ─────────────── multi-lineage fasta (examples/lineages.fasta) ───────────────

const LINEAGES: &str = "examples/lineages.fasta";
const GROUP_AVG_HEADER: [&str; 8] = [
    "Group1", "Group2", "NumSeqs1", "NumSeqs2", "NumComparisons", "Mean_dN/dS", "StdError", "95%CI",
];

#[test]
fn fasta_lineages_group_average_default_split() {
    // IDs split on '_' into three lineages (Lineage2 / Lineage4 / Bovis), 2 isolates each.
    let r = new_run();
    run_ok(&["fasta", LINEAGES, "-o", &r.prefix, "--group-average"]);
    let rows = tsv(&format!("{}_group_avg_dn_ds.tsv", r.prefix));
    assert_eq!(rows[0], GROUP_AVG_HEADER);
    // 3 groups => 3 within-group + 3 between-group rows.
    assert_eq!(rows.len() - 1, 6, "3 groups -> 6 rows");
    assert_eq!(
        group_labels(&rows),
        ["Bovis", "Lineage2", "Lineage4"].iter().map(|s| s.to_string()).collect()
    );
    // Between-group comparison: 2 x 2 isolates = 4 pairwise comparisons.
    let l2_l4 = find_group(&rows, "Lineage2", "Lineage4");
    assert_eq!(l2_l4[4], "4", "NumComparisons for a between-lineage cell");
    assert!((f(&l2_l4[5]) - 0.621198).abs() < EPS, "Lineage2-Lineage4 mean {}", l2_l4[5]);
    // Within-lineage cell: C(2,2) = 1 comparison, SE reported as N/A.
    let bovis = find_group(&rows, "Bovis", "Bovis");
    assert_eq!(bovis[4], "1");
    assert_eq!(bovis[6], "N/A");
}

#[test]
fn fasta_lineages_first_letter_changes_grouping() {
    // --first-letter-lineage groups by the first character, so the two "Lineage*"
    // lineages merge into "L" while "Bovis" becomes "B": a genuinely different
    // partition than the default '_' split (3 groups -> 2 groups).
    let def = new_run();
    run_ok(&["fasta", LINEAGES, "-o", &def.prefix, "--group-average"]);
    let default_rows = tsv(&format!("{}_group_avg_dn_ds.tsv", def.prefix));

    let fl = new_run();
    run_ok(&["fasta", LINEAGES, "-o", &fl.prefix, "--group-average", "--first-letter-lineage"]);
    let rows = tsv(&format!("{}_group_avg_dn_ds.tsv", fl.prefix));

    assert_eq!(rows[0], GROUP_AVG_HEADER);
    // 2 groups => 2 within + 1 between = 3 rows (vs 6 for the default split).
    assert_eq!(rows.len() - 1, 3, "first-letter grouping should collapse to 3 rows");
    assert_ne!(rows.len(), default_rows.len(), "flag must change the grouping");
    assert_eq!(
        group_labels(&rows),
        ["B", "L"].iter().map(|s| s.to_string()).collect()
    );
    // "L" merges the 4 Lineage2+Lineage4 isolates: C(4,2)=6 within, 4x2=8 between.
    let ll = find_group(&rows, "L", "L");
    assert_eq!((ll[2].as_str(), ll[4].as_str()), ("4", "6"), "L within: 4 seqs, 6 comparisons");
    let lb = find_group(&rows, "L", "B");
    assert_eq!(lb[4], "8", "L-B: 4x2 = 8 comparisons");
}

#[test]
fn fasta_lineages_lineage_mode() {
    let r = new_run();
    run_ok(&["fasta", LINEAGES, "-o", &r.prefix, "--lineage"]);
    let rows = tsv(&format!("{}_lineage_summary.tsv", r.prefix));
    assert_eq!(
        rows[0],
        ["Genome", "Against_Lineage", "Mean_dN", "Mean_dS", "dN/dS_Ratio"]
    );
    // 6 genomes x 3 lineages each = 18 rows.
    assert_eq!(rows.len() - 1, 18, "6 genomes vs 3 lineages");
    let genomes: std::collections::BTreeSet<&str> = rows[1..].iter().map(|r| r[0].as_str()).collect();
    assert_eq!(genomes.len(), 6, "one block per genome");
}

// ─────────────── pass-only filter (examples/toy_genome/variants_mixed.vcf) ───────────────

#[test]
fn vcf_pass_only_drops_lowqual_variants() {
    // variants_mixed.vcf carries 12 gene01 SNPs: 8 PASS + 4 LowQual.
    const MIXED: &str = "examples/toy_genome/variants_mixed.vcf";
    let all = new_run();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", MIXED, "--genetic-code", "11", "-o", &all.prefix,
    ]);
    let pass = new_run();
    run_ok(&[
        "vcf", "--ref", REF, "--gff", GFF, "--vcf", MIXED, "--genetic-code", "11",
        "--pass-only", "-o", &pass.prefix,
    ]);
    let g_all = f(&row(&tsv(&format!("{}_pnps.tsv", all.prefix)), "gene01")[9]);
    let g_pass = f(&row(&tsv(&format!("{}_pnps.tsv", pass.prefix)), "gene01")[9]);
    assert_eq!(g_all, 12.0, "all 12 SNPs counted without --pass-only");
    assert_eq!(g_pass, 8.0, "--pass-only keeps only the 8 PASS SNPs");
}

// ─────────────────────────── error handling ───────────────────────────

#[test]
fn vcf_missing_required_ref_fails() {
    let r = new_run();
    let out = run(&["vcf", "--gff", GFF, "--vcf", VCF, "-o", &r.prefix]);
    assert!(!out.status.success(), "missing --ref must be a hard error");
}

#[test]
fn vcf_nonexistent_reference_fails_clearly() {
    let r = new_run();
    let out = run(&[
        "vcf", "--ref", "examples/toy_genome/does_not_exist.fasta",
        "--gff", GFF, "--vcf", VCF, "--genetic-code", "11", "-o", &r.prefix,
    ]);
    assert!(!out.status.success());
}


#[test]
fn vcf_examples_min_snps_counts_raw_not_af_weighted() {
    // --min-snps must filter on the real SNP count, so --af-weighted keeps the same
    // genes as plain mode (the AF-weighted fractional total is smaller but irrelevant).
    let plain = new_run();
    vcf_core(&plain.prefix, &["--min-snps", "20"]);
    let weighted = new_run();
    vcf_core(&weighted.prefix, &["--min-snps", "20", "--af-weighted"]);
    let n_plain = tsv(&format!("{}_pnps.tsv", plain.prefix)).len();
    let n_weighted = tsv(&format!("{}_pnps.tsv", weighted.prefix)).len();
    assert_eq!(n_plain, n_weighted, "--min-snps must not drop more genes under --af-weighted");
    assert!(n_plain > 1, "some genes should survive --min-snps 20");
}
