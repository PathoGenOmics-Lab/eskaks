//! Golden-output snapshot test: run `eskaks vcf` on the bundled toy genome and assert
//! every table (pN/pS, variants, diversity, McDonald-Kreitman; TSV and JSON) is
//! byte-identical to the committed reference in `tests/golden/`. Freezing the exact
//! bytes means any silent drift (a shifted column, a changed value, a stray `-0`, a
//! formatting regression) fails with an exact diff instead of slipping through.
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
            if produced != golden {
                mismatches.push(format!("toy_{table}.{ext}"));
            }
        }
    }

    assert!(
        mismatches.is_empty(),
        "output drifted from the golden snapshot: {}.\n\
         If the change is intended, regenerate with `BLESS=1 cargo test --test golden` \
         and review the diff before committing.",
        mismatches.join(", ")
    );
}
