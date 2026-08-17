//! Golden-output snapshot test: run `eskaks vcf` on the bundled toy genome and assert
//! every table (pN/pS, variants, diversity, McDonald-Kreitman; TSV and JSON) matches the
//! committed reference in `tests/golden/`. The non-numeric skeleton (columns, keys,
//! ordering, formatting) must be byte-identical, so any structural drift (a shifted
//! column, a renamed key, a changed count, a stray token) still fails with a diff. Each
//! numeric token is compared with a tiny relative tolerance: several pN/pS columns
//! (`p_value`, `p_bonferroni`, `p_gc`, `q_gc_bh`, ...) are computed through transcendental
//! libm functions (`exp` via erfc, `ln` via ln_gamma) whose last bit differs between
//! macOS and glibc, so an exact byte compare fails on one OS of the CI matrix even when
//! nothing actually drifted. The tolerance (1e-9 relative) is far below any real change
//! yet well above cross-platform ULP noise.
//!
//! Regenerate the golden files after an *intended* output change, then review the diff
//! before committing:
//!     BLESS=1 cargo test --test golden

use std::path::PathBuf;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_eskaks").to_string()
}

fn manifest() -> &'static str {
    env!("CARGO_MANIFEST_DIR")
}

/// The deterministic command whose output is frozen. `--workers 1` removes any
/// parallel-ordering nondeterminism; `--genomic-control` populates the P_GC / Q_GC_BH
/// columns so they are covered too.
fn run(prefix: &str, format: &str) {
    let status = Command::new(bin())
        .current_dir(manifest())
        .args([
            "vcf",
            "--ref",
            "examples/toy_genome/reference.fasta",
            "--gff",
            "examples/toy_genome/genes.gff3",
            "--vcf",
            "examples/toy_genome/variants_multisample.vcf",
            "--genetic-code",
            "11",
            "--workers",
            "1",
            "--variants",
            "--diversity",
            "--mk",
            "--genomic-control",
            "--format",
            format,
            "-o",
            prefix,
            "--quiet",
        ])
        .status()
        .expect("failed to spawn eskaks");
    assert!(status.success(), "eskaks vcf ({format}) exited non-zero");
}

/// The per-codon recurrence scan is opt-in, so it gets its own invocation rather than a
/// new flag on `run` above: adding `--codon-scan` there would rewrite every existing
/// golden table and destroy the guarantee that a run without the flag is unchanged.
fn run_codon_scan(prefix: &str, format: &str) {
    let status = Command::new(bin())
        .current_dir(manifest())
        .args([
            "vcf",
            "--ref",
            "examples/toy_genome/reference.fasta",
            "--gff",
            "examples/toy_genome/genes.gff3",
            "--vcf",
            "examples/toy_genome/variants_multisample.vcf",
            "--genetic-code",
            "11",
            "--workers",
            "1",
            "--codon-scan",
            "--format",
            format,
            "-o",
            prefix,
            "--quiet",
        ])
        .status()
        .expect("failed to spawn eskaks");
    assert!(status.success(), "eskaks vcf --codon-scan ({format}) exited non-zero");
}

/// `--shared-codons` adds a column to the variants table, which would rewrite
/// `toy_variants.{tsv,json}` for every reader of the default output. It is therefore
/// opt-in, exactly like `--codon-scan`, and gets its own invocation and its own golden
/// pair so the new column is frozen without the existing table moving.
fn run_shared_codons(prefix: &str, format: &str) {
    let status = Command::new(bin())
        .current_dir(manifest())
        .args([
            "vcf",
            "--ref",
            "examples/toy_genome/reference.fasta",
            "--gff",
            "examples/toy_genome/genes.gff3",
            "--vcf",
            "examples/toy_genome/variants_multisample.vcf",
            "--genetic-code",
            "11",
            "--workers",
            "1",
            "--variants",
            "--shared-codons",
            "--format",
            format,
            "-o",
            prefix,
            "--quiet",
        ])
        .status()
        .expect("failed to spawn eskaks");
    assert!(status.success(), "eskaks vcf --shared-codons ({format}) exited non-zero");
}

/// `--tree` appends the independent-origin columns to BOTH the per-codon scan and the
/// per-variant table, so it gets its own invocation and its own golden pair for the same
/// reason `--shared-codons` does: the tables written without it must not move.
fn run_tree(prefix: &str, format: &str) {
    let status = Command::new(bin())
        .current_dir(manifest())
        .args([
            "vcf",
            "--ref",
            "examples/toy_genome/reference.fasta",
            "--gff",
            "examples/toy_genome/genes.gff3",
            "--vcf",
            "examples/toy_genome/variants_multisample.vcf",
            "--genetic-code",
            "11",
            "--workers",
            "1",
            "--codon-scan",
            "--variants",
            "--tree",
            "examples/toy_genome/samples.nwk",
            "--format",
            format,
            "-o",
            prefix,
            "--quiet",
        ])
        .status()
        .expect("failed to spawn eskaks");
    assert!(status.success(), "eskaks vcf --tree ({format}) exited non-zero");
}

/// Two numbers are "equal" within cross-platform libm noise. Real value changes are
/// orders of magnitude larger than this; last-bit differences from `exp`/`ln` are orders
/// of magnitude smaller.
fn approx_eq(a: f64, b: f64) -> bool {
    if a == b || (a.is_nan() && b.is_nan()) {
        return true;
    }
    (a - b).abs() <= 1e-12 + 1e-9 * b.abs()
}

/// Does a numeric literal start at `s[i]`? A number begins at a digit, or at a sign/dot
/// that immediately leads a digit. A lone `+`/`-`/`.` (e.g. the `+` strand field) is not
/// a number, so it is compared as an exact literal.
fn is_num_start(s: &[u8], i: usize) -> bool {
    let digit_at = |k: usize| k < s.len() && s[k].is_ascii_digit();
    match s[i] {
        c if c.is_ascii_digit() => true,
        b'-' | b'+' => digit_at(i + 1) || (i + 1 < s.len() && s[i + 1] == b'.' && digit_at(i + 2)),
        b'.' => digit_at(i + 1),
        _ => false,
    }
}

/// End index (exclusive) of the numeric literal beginning at `start`:
/// `[-+]? digits? (. digits?)? ([eE] [-+]? digits)?`.
fn num_end(s: &[u8], start: usize) -> usize {
    let mut i = start;
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        i += 1;
    }
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
    }
    if i < s.len() && s[i] == b'.' {
        i += 1;
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }
    }
    if i < s.len() && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < s.len() && (s[j] == b'-' || s[j] == b'+') {
            j += 1;
        }
        if j < s.len() && s[j].is_ascii_digit() {
            i = j;
            while i < s.len() && s[i].is_ascii_digit() {
                i += 1;
            }
        }
    }
    i
}

fn snippet(s: &str, at: usize) -> String {
    let lo = at.saturating_sub(24);
    let hi = (at + 24).min(s.len());
    format!("...{}...", s.get(lo..hi).unwrap_or("<?>"))
}

/// Compare `produced` against `golden`, allowing each numeric token a tiny relative
/// tolerance but requiring the surrounding text to be byte-identical. Returns a
/// description of the first real difference, so a genuine regression still reports the
/// offending value instead of hiding behind the tolerance.
fn compare_tolerant(produced: &str, golden: &str) -> Result<(), String> {
    let (p, g) = (produced.as_bytes(), golden.as_bytes());
    let (mut i, mut j) = (0usize, 0usize);
    while i < p.len() && j < g.len() {
        if is_num_start(p, i) && is_num_start(g, j) {
            let (pe, ge) = (num_end(p, i), num_end(g, j));
            match (produced[i..pe].parse::<f64>(), golden[j..ge].parse::<f64>()) {
                (Ok(a), Ok(b)) if approx_eq(a, b) => {
                    i = pe;
                    j = ge;
                    continue;
                }
                (Ok(a), Ok(b)) => {
                    let rel = if b != 0.0 { ((a - b) / b).abs() } else { f64::INFINITY };
                    return Err(format!(
                        "value {a} vs golden {b} (rel diff {rel:.2e}) near {}",
                        snippet(golden, j)
                    ));
                }
                _ => {} // not both parseable: fall through to a byte compare
            }
        }
        if p[i] != g[j] {
            return Err(format!(
                "text mismatch:\n    produced {}\n    golden   {}",
                snippet(produced, i),
                snippet(golden, j)
            ));
        }
        i += 1;
        j += 1;
    }
    if i != p.len() || j != g.len() {
        return Err(format!(
            "length mismatch: {} produced byte(s) and {} golden byte(s) left over",
            p.len() - i,
            g.len() - j
        ));
    }
    Ok(())
}

#[test]
fn golden_vcf_outputs_match() {
    let bless = std::env::var_os("BLESS").is_some();
    let prefix = std::env::temp_dir().join("eskaks_golden_run");
    let prefix = prefix.to_str().expect("temp path is valid UTF-8");
    let golden_dir = PathBuf::from(manifest()).join("tests/golden");

    let mut mismatches = Vec::new();
    for (format, ext) in [("tsv", "tsv"), ("json", "json")] {
        run(prefix, format);
        for table in ["pnps", "variants", "diversity", "mk"] {
            let produced = std::fs::read_to_string(format!("{prefix}_{table}.{ext}"))
                .unwrap_or_else(|e| panic!("cannot read produced {table}.{ext}: {e}"));
            let golden_path = golden_dir.join(format!("toy_{table}.{ext}"));
            if bless {
                std::fs::write(&golden_path, &produced).expect("write golden");
                continue;
            }
            let golden = std::fs::read_to_string(&golden_path).unwrap_or_else(|_| {
                panic!(
                    "missing golden {}; run `BLESS=1 cargo test --test golden` to create it",
                    golden_path.display()
                )
            });
            if let Err(reason) = compare_tolerant(&produced, &golden) {
                mismatches.push(format!("toy_{table}.{ext}: {reason}"));
            }
        }
    }

    assert!(
        mismatches.is_empty(),
        "output drifted from the golden snapshot:\n{}\n\
         If the change is intended, regenerate with `BLESS=1 cargo test --test golden` \
         and review the diff before committing.",
        mismatches.join("\n")
    );
}

#[test]
fn golden_codon_scan_matches() {
    let bless = std::env::var_os("BLESS").is_some();
    let prefix = std::env::temp_dir().join("eskaks_golden_codons");
    let prefix = prefix.to_str().expect("temp path is valid UTF-8");
    let golden_dir = PathBuf::from(manifest()).join("tests/golden");

    let mut mismatches = Vec::new();
    for (format, ext) in [("tsv", "tsv"), ("json", "json")] {
        run_codon_scan(prefix, format);
        let produced = std::fs::read_to_string(format!("{prefix}_codons.{ext}"))
            .unwrap_or_else(|e| panic!("cannot read produced codons.{ext}: {e}"));
        let golden_path = golden_dir.join(format!("toy_codons.{ext}"));
        if bless {
            std::fs::write(&golden_path, &produced).expect("write golden");
            continue;
        }
        let golden = std::fs::read_to_string(&golden_path).unwrap_or_else(|_| {
            panic!(
                "missing golden {}; run `BLESS=1 cargo test --test golden` to create it",
                golden_path.display()
            )
        });
        if let Err(reason) = compare_tolerant(&produced, &golden) {
            mismatches.push(format!("toy_codons.{ext}: {reason}"));
        }
    }

    assert!(
        mismatches.is_empty(),
        "the per-codon scan drifted from the golden snapshot:\n{}\n\
         If the change is intended, regenerate with `BLESS=1 cargo test --test golden` \
         and review the diff before committing.",
        mismatches.join("\n")
    );
}

#[test]
fn golden_shared_codons_matches() {
    let bless = std::env::var_os("BLESS").is_some();
    let prefix = std::env::temp_dir().join("eskaks_golden_shared");
    let prefix = prefix.to_str().expect("temp path is valid UTF-8");
    let golden_dir = PathBuf::from(manifest()).join("tests/golden");

    let mut mismatches = Vec::new();
    for (format, ext) in [("tsv", "tsv"), ("json", "json")] {
        run_shared_codons(prefix, format);
        let produced = std::fs::read_to_string(format!("{prefix}_variants.{ext}"))
            .unwrap_or_else(|e| panic!("cannot read produced variants.{ext}: {e}"));
        let golden_path = golden_dir.join(format!("toy_variants_shared.{ext}"));
        if bless {
            std::fs::write(&golden_path, &produced).expect("write golden");
            continue;
        }
        let golden = std::fs::read_to_string(&golden_path).unwrap_or_else(|_| {
            panic!(
                "missing golden {}; run `BLESS=1 cargo test --test golden` to create it",
                golden_path.display()
            )
        });
        if let Err(reason) = compare_tolerant(&produced, &golden) {
            mismatches.push(format!("toy_variants_shared.{ext}: {reason}"));
        }
    }

    assert!(
        mismatches.is_empty(),
        "the shared-codon marks drifted from the golden snapshot:\n{}\n\
         If the change is intended, regenerate with `BLESS=1 cargo test --test golden` \
         and review the diff before committing.",
        mismatches.join("\n")
    );
}

#[test]
fn golden_tree_origins_match() {
    let bless = std::env::var_os("BLESS").is_some();
    let prefix = std::env::temp_dir().join("eskaks_golden_tree");
    let prefix = prefix.to_str().expect("temp path is valid UTF-8");
    let golden_dir = PathBuf::from(manifest()).join("tests/golden");

    let mut mismatches = Vec::new();
    for (format, ext) in [("tsv", "tsv"), ("json", "json")] {
        run_tree(prefix, format);
        for (table, golden_name) in
            [("codons", "toy_codons_tree"), ("variants", "toy_variants_tree")]
        {
            let produced = std::fs::read_to_string(format!("{prefix}_{table}.{ext}"))
                .unwrap_or_else(|e| panic!("cannot read produced {table}.{ext}: {e}"));
            let golden_path = golden_dir.join(format!("{golden_name}.{ext}"));
            if bless {
                std::fs::write(&golden_path, &produced).expect("write golden");
                continue;
            }
            let golden = std::fs::read_to_string(&golden_path).unwrap_or_else(|_| {
                panic!(
                    "missing golden {}; run `BLESS=1 cargo test --test golden` to create it",
                    golden_path.display()
                )
            });
            if let Err(reason) = compare_tolerant(&produced, &golden) {
                mismatches.push(format!("{golden_name}.{ext}: {reason}"));
            }
        }
    }

    assert!(
        mismatches.is_empty(),
        "the independent-origin columns drifted from the golden snapshot:\n{}\n\
         If the change is intended, regenerate with `BLESS=1 cargo test --test golden` \
         and review the diff before committing.",
        mismatches.join("\n")
    );
}

/// The two appended columns must be the ONLY difference the flag makes: same rows, same
/// order, same values everywhere else. Without this, the opt-in promise is untested.
#[test]
fn shared_codons_only_appends_columns() {
    let plain_prefix = std::env::temp_dir().join("eskaks_golden_shared_plain");
    let plain_prefix = plain_prefix.to_str().expect("temp path is valid UTF-8");
    let flagged_prefix = std::env::temp_dir().join("eskaks_golden_shared_flagged");
    let flagged_prefix = flagged_prefix.to_str().expect("temp path is valid UTF-8");

    run(plain_prefix, "tsv");
    run_shared_codons(flagged_prefix, "tsv");
    let plain = std::fs::read_to_string(format!("{plain_prefix}_variants.tsv")).expect("plain");
    let flagged =
        std::fs::read_to_string(format!("{flagged_prefix}_variants.tsv")).expect("flagged");

    let plain_lines: Vec<&str> = plain.lines().collect();
    let flagged_lines: Vec<&str> = flagged.lines().collect();
    assert_eq!(plain_lines.len(), flagged_lines.len(), "same number of rows");
    let mut marked = 0usize;
    for (i, (p, f)) in plain_lines.iter().zip(&flagged_lines).enumerate() {
        let (rest, codon_change) = f.rsplit_once('\t').unwrap_or_else(|| panic!("row {i}: {f}"));
        let (prefix, mark) = rest.rsplit_once('\t').unwrap_or_else(|| panic!("row {i}: {f}"));
        assert_eq!(prefix, *p, "row {i} changed outside the two new columns");
        if i == 0 {
            assert_eq!((mark, codon_change), ("Codon_Shared", "Codon_Change"));
            continue;
        }
        assert!(
            matches!(mark, "true" | "false" | "NA"),
            "row {i} has an unexpected Codon_Shared value {mark:?}"
        );
        marked += usize::from(mark == "true");
        // `CTT>TTA`: two 3-base codons on the coding strand, and the change is real.
        let (from, to) = codon_change.split_once('>').expect("codon change is FROM>TO");
        assert!(
            from.len() == 3 && to.len() == 3 && from != to,
            "row {i} has a malformed Codon_Change {codon_change:?}"
        );
        // A codon change spanning more than one base can only come from a shared codon,
        // so the two columns cannot contradict each other. (The converse does not hold: an
        // allele most of whose carriers hold it ALONE is still marked, and still reports
        // its own single-base codon change.)
        let diffs = from.bytes().zip(to.bytes()).filter(|(a, b)| a != b).count();
        assert!(
            diffs == 1 || mark == "true",
            "row {i}: {codon_change} changes {diffs} bases but Codon_Shared={mark}"
        );
    }
    // The bundled toy VCF is genotyped, so nothing is unknowable and 13 codons are shared.
    assert!(
        !flagged_lines.iter().skip(1).any(|l| l.split('\t').nth(12) == Some("NA")),
        "a genotyped VCF must leave no row unknowable"
    );
    assert_eq!(marked, 27, "27 of the 264 rows sit on the 13 shared codons");
    // ...and at least one of them really is a multi-nucleotide codon change, or the joint
    // scoring never engaged and the assertions above are vacuous.
    assert!(
        flagged_lines.iter().skip(1).any(|l| {
            let cc = l.rsplit('\t').next().expect("last column");
            let (from, to) = cc.split_once('>').expect("codon change");
            from.bytes().zip(to.bytes()).filter(|(a, b)| a != b).count() > 1
        }),
        "no row reports a multi-nucleotide codon change"
    );
}

/// The opt-in scan must leave the default output byte-for-byte alone. `golden_vcf_outputs_match`
/// covers the tables it writes; this covers the one file that exists only under the flag.
#[test]
fn codon_scan_is_the_only_new_file() {
    let prefix = std::env::temp_dir().join("eskaks_golden_codons_optin");
    let prefix = prefix.to_str().expect("temp path is valid UTF-8");
    let codons = format!("{prefix}_codons.tsv");
    let _ = std::fs::remove_file(&codons);

    run(prefix, "tsv");
    assert!(
        !PathBuf::from(&codons).exists(),
        "a run WITHOUT --codon-scan must not write {codons}"
    );

    run_codon_scan(prefix, "tsv");
    assert!(PathBuf::from(&codons).exists(), "--codon-scan must write {codons}");
}

#[cfg(test)]
mod compare_tests {
    use super::*;
    use std::path::PathBuf;

    /// Scale every *float* token (those with a `.` or exponent, i.e. the computed values,
    /// not integer coordinates/counts) by `1 + rel`, to mimic a platform whose libm gives
    /// slightly different last bits.
    fn perturb_floats(s: &str, rel: f64) -> String {
        let b = s.as_bytes();
        let mut out = String::with_capacity(s.len());
        let mut i = 0;
        while i < b.len() {
            if is_num_start(b, i) {
                let e = num_end(b, i);
                let tok = &s[i..e];
                if tok.contains(['.', 'e', 'E']) {
                    if let Ok(v) = tok.parse::<f64>() {
                        out.push_str(&format!("{}", v * (1.0 + rel)));
                        i = e;
                        continue;
                    }
                }
                out.push_str(tok);
                i = e;
                continue;
            }
            out.push(b[i] as char);
            i += 1;
        }
        out
    }

    #[test]
    fn tolerance_absorbs_cross_platform_float_noise() {
        // The real pnps golden, with every computed value nudged by a bit more than the
        // observed macOS-vs-glibc difference, must still pass (an exact byte compare
        // would not) -- and a 1e-6 change must still fail. This is the exact scenario that
        // broke the byte-for-byte test on the ubuntu leg of the CI matrix.
        let golden = std::fs::read_to_string(
            PathBuf::from(manifest()).join("tests/golden/toy_pnps.json"),
        )
        .expect("read golden");
        let noisy = perturb_floats(&golden, 5e-14);
        assert_ne!(noisy, golden, "perturbation must change the bytes");
        assert!(
            compare_tolerant(&noisy, &golden).is_ok(),
            "cross-platform ULP noise must be tolerated"
        );
        let changed = perturb_floats(&golden, 1e-6);
        assert!(
            compare_tolerant(&changed, &golden).is_err(),
            "a real 1e-6 value change must still fail"
        );
    }

    #[test]
    fn accepts_last_bit_float_noise() {
        // The kind of difference macOS vs glibc produces in a p-value: identical to ~15
        // significant figures. Must pass.
        let a = r#"{"p":0.3118239208586769,"n":25.0000,"gene":"gene01"}"#;
        let b = r#"{"p":0.3118239208586771,"n":25.0000,"gene":"gene01"}"#;
        assert!(compare_tolerant(a, b).is_ok());
    }

    #[test]
    fn rejects_a_real_value_change() {
        let a = r#"{"p":0.31182}"#;
        let b = r#"{"p":0.30000}"#;
        assert!(compare_tolerant(a, b).is_err());
    }

    #[test]
    fn rejects_structural_drift() {
        // renamed key, extra column, changed count, reordered field: all caught.
        assert!(compare_tolerant(r#"{"p":0.5}"#, r#"{"q":0.5}"#).is_err());
        assert!(compare_tolerant("a\tb\t1.0", "a\t1.0").is_err());
        assert!(compare_tolerant(r#"{"n":25.0}"#, r#"{"n":24.0}"#).is_err());
    }

    #[test]
    fn digits_inside_identifiers_are_matched_exactly() {
        // "gene01"/"chr1"/"PE_PGRS_toy2" contain digits but are identical on both sides,
        // so a genuine p-value ULP diff alongside them still passes.
        let a = r#"{"gene":"PE_PGRS_toy2","chrom":"chr1","p":0.0004967449122598745}"#;
        let b = r#"{"gene":"PE_PGRS_toy2","chrom":"chr1","p":0.0004967449122598748}"#;
        assert!(compare_tolerant(a, b).is_ok());
        // but a changed identifier digit is a real mismatch
        let c = r#"{"gene":"PE_PGRS_toy2","chrom":"chr1","p":0.0004967449122598745}"#;
        let d = r#"{"gene":"PE_PGRS_toy3","chrom":"chr1","p":0.0004967449122598745}"#;
        assert!(compare_tolerant(c, d).is_err());
    }
}
