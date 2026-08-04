//! Docs <-> code contract: keep the documentation in sync with the binary, so a new
//! output column or CLI flag cannot ship undocumented (the drift that once left four
//! pN/pS columns and several flags unmentioned).

use std::collections::BTreeSet;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_eskaks").to_string()
}

fn read(rel: &str) -> String {
    let path = format!("{}/{}", env!("CARGO_MANIFEST_DIR"), rel);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
}

/// The slice of a markdown doc from `heading` up to the next top-level `## ` heading.
fn section_of(doc: &str, heading: &str) -> String {
    let start = doc
        .find(heading)
        .unwrap_or_else(|| panic!("heading {heading:?} not found"));
    let rest = &doc[start + heading.len()..];
    let end = rest.find("\n## ").unwrap_or(rest.len());
    rest[..end].to_string()
}

/// `--long-flag` tokens appearing in help text (or any text).
fn long_flags(text: &str) -> BTreeSet<String> {
    let mut set = BTreeSet::new();
    for tok in text.split(|c: char| c.is_whitespace() || matches!(c, ',' | '<' | '=' | '[' | '|')) {
        if let Some(rest) = tok.strip_prefix("--") {
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '-')
                .collect();
            if name.len() >= 2 {
                set.insert(format!("--{name}"));
            }
        }
    }
    set
}

#[test]
fn every_pnps_column_is_documented() {
    // Source of truth: the exact header the binary writes (via the golden fixture).
    let header: Vec<String> = read("tests/golden/toy_pnps.tsv")
        .lines()
        .next()
        .expect("golden pnps has a header")
        .split('\t')
        .map(str::to_string)
        .collect();

    let section = section_of(&read("docs/vcf-analysis.md"), "## Output");
    let missing: Vec<&String> = header.iter().filter(|c| !section.contains(c.as_str())).collect();
    assert!(
        missing.is_empty(),
        "these pN/pS columns are written but not documented in vcf-analysis.md '## Output': {missing:?}"
    );
}

#[test]
fn every_vcf_help_flag_is_documented() {
    check_flags_documented("vcf", "docs/cli-reference.md");
}

#[test]
fn every_fasta_help_flag_is_documented() {
    check_flags_documented("fasta", "docs/cli-reference.md");
}

fn check_flags_documented(subcmd: &str, doc_rel: &str) {
    let out = Command::new(bin())
        .args([subcmd, "--help"])
        .output()
        .expect("spawn eskaks");
    let help = String::from_utf8_lossy(&out.stdout);
    let doc = read(doc_rel);

    // clap auto-generates --help/--version; they need no manual documentation.
    let auto: BTreeSet<String> = ["--help", "--version"].iter().map(|s| s.to_string()).collect();
    let missing: Vec<String> = long_flags(&help)
        .into_iter()
        .filter(|f| !auto.contains(f) && !doc.contains(f.as_str()))
        .collect();
    assert!(
        missing.is_empty(),
        "eskaks {subcmd} accepts these flags but {doc_rel} does not document them: {missing:?}"
    );
}
