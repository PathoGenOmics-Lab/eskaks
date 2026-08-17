//! Independent-origin counting (`--tree`), on a cohort whose truth is constructed.
//!
//! The gap this feature exists to close is precise: the per-codon scan counts DISTINCT
//! alleles, so it sees a residue that collected six different amino acids and is blind
//! to one allele that arose fifty times over. Carrier counts cannot rescue it, because
//! an ancestral mutation carried by an expanded clade has exactly the shape of many
//! independent origins.
//!
//! So this file builds a 64-sample cohort on a known tree containing, at the same
//! `A_c = 6` codon opportunity and in the same gene:
//!
//! * a **clonal** allele: one origin, carried by a whole 8-sample lineage;
//! * a **convergent** allele: the SAME single allele in eight separate lineages, eight
//!   origins, 16 carriers;
//! * a **calling-artefact** allele: five scattered singletons, which raw parsimony reads
//!   as five origins and `--min-origin-support 2` reads as none;
//! * a **multi-allelic** codon: four distinct alleles, one origin each, which is the
//!   shape the existing allele-count test already finds.
//!
//! Every expected number below is derived by hand from that construction, not read off
//! a run: the tree is explicit, the carriers are explicit, and Fitch parsimony on them
//! has one answer.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_eskaks").to_string()
}

/// 64 samples in 8 lineages of 8, `S01`..`S64`. Lineage `l` holds samples
/// `8l+1 ..= 8l+8`, and within it samples pair up as (1,2), (3,4), (5,6), (7,8).
const N_SAMPLES: usize = 64;
const LINEAGES: usize = 8;

fn sample_name(i: usize) -> String {
    format!("S{:02}", i + 1)
}

/// A balanced tree over the eight lineages, each lineage a polytomy of four sister
/// pairs. Written with branch lengths and a support value, because a real file has them.
///
/// The lineages are written in a SCRAMBLED order, so the first tip in the file is not
/// the first sample of the cohort. Every expected origin count below is invariant to
/// that, because a clade is a set of samples however the file orders it, and pinning the
/// invariance is the point: the answer must not depend on tip order. That the join is by
/// NAME and not by position is pinned separately, by `swapping_two_tip_labels_changes_the_answer`
/// below and by `tree::tests::the_join_maps_names_and_not_positions`.
const LINEAGE_ORDER: [usize; LINEAGES] = [3, 0, 5, 2, 7, 1, 6, 4];

fn newick() -> String {
    let pair = |a: usize, b: usize| format!("({}:0.001,{}:0.001):0.004", sample_name(a), sample_name(b));
    let lineage = |l: usize| {
        let base = l * 8;
        format!(
            "({},{},{},{})L{}:0.050",
            pair(base, base + 1),
            pair(base + 2, base + 3),
            pair(base + 4, base + 5),
            pair(base + 6, base + 7),
            l + 1
        )
    };
    let clade = |a: usize, b: usize| {
        format!("({},{})0.99:0.080", lineage(LINEAGE_ORDER[a]), lineage(LINEAGE_ORDER[b]))
    };
    format!(
        "(({},{}),({},{}));",
        clade(0, 1),
        clade(2, 3),
        clade(4, 5),
        clade(6, 7)
    )
}

/// 299 `GCT` (Ala) codons and a terminal stop. `GCT` has exactly 6 possible
/// nonsynonymous changes and 3 synonymous ones, none of them creating a stop, so every
/// codon of this gene has `A_c = 6` and the family is 299 codons / 1794 changes.
const CODONS: usize = 300;
const A_C: u64 = 6;
const FAMILY_M: u64 = (CODONS - 1) as u64;
const FAMILY_POSS: u64 = FAMILY_M * A_C;

/// First base (1-based) of codon `c` (1-based).
fn codon_pos(c: usize) -> usize {
    3 * c - 2
}

/// One VCF record: position, REF, the ALTs, and the samples carrying each ALT.
struct Record {
    pos: usize,
    ref_base: char,
    alts: Vec<(char, Vec<usize>)>,
}

/// The constructed cohort. Returns the records in position order.
fn records() -> Vec<Record> {
    let pair_of = |l: usize, p: usize| vec![l * 8 + 2 * p, l * 8 + 2 * p + 1];
    let mut recs = Vec::new();

    // 40 background alleles, one per codon 10..49, each carried by one sister pair:
    // one origin each, which is what makes the genome-wide rate an honest baseline.
    for (k, c) in (10..50).enumerate() {
        let (l, p) = (k % LINEAGES, (k / LINEAGES) % 4);
        recs.push(Record { pos: codon_pos(c), ref_base: 'G', alts: vec![('A', pair_of(l, p))] });
    }

    // CLONAL: one allele carried by all 8 samples of lineage 0. Eight carriers, and
    // exactly one origin, which is the case carrier counts get wrong.
    recs.push(Record {
        pos: codon_pos(100),
        ref_base: 'G',
        alts: vec![('C', (0..8).collect())],
    });

    // ARTEFACT: five scattered singletons, one per lineage. Raw parsimony calls that
    // five origins; requiring two carriers per origin calls it none.
    recs.push(Record {
        pos: codon_pos(150),
        ref_base: 'G',
        alts: vec![('T', (0..5).map(|l| l * 8).collect())],
    });

    // CONVERGENT: the SAME allele in the first sister pair of every lineage. Sixteen
    // carriers, eight independent origins. One allele, so the allele-count test cannot
    // see it at all.
    recs.push(Record {
        pos: codon_pos(200),
        ref_base: 'G',
        alts: vec![('T', (0..LINEAGES).flat_map(|l| pair_of(l, 0)).collect())],
    });

    // MULTI-ALLELIC: four distinct missense alleles at one codon, on disjoint carriers
    // (so they are not one multi-nucleotide event), one origin each. The shape the
    // allele-count test already finds, kept here so the two tests can be compared on a
    // codon where they agree.
    recs.push(Record {
        pos: codon_pos(250),
        ref_base: 'G',
        alts: vec![('A', pair_of(0, 1)), ('C', pair_of(0, 2)), ('T', pair_of(0, 3))],
    });
    recs.push(Record {
        pos: codon_pos(250) + 1,
        ref_base: 'C',
        alts: vec![('A', pair_of(1, 1))],
    });
    recs.sort_by_key(|r| r.pos);
    recs
}

fn write_inputs(dir: &Path) {
    let mut fasta = String::from(">chr1\n");
    let mut seq = "GCT".repeat(CODONS - 1);
    seq.push_str("TAA");
    for chunk in seq.as_bytes().chunks(60) {
        fasta.push_str(std::str::from_utf8(chunk).expect("ascii"));
        fasta.push('\n');
    }
    std::fs::write(dir.join("ref.fasta"), fasta).expect("write reference");

    let gff = format!(
        "##gff-version 3\nchr1\teskaks\tgene\t1\t{len}\t.\t+\t.\tID=g1;Name=geneA\n\
         chr1\teskaks\tCDS\t1\t{len}\t.\t+\t0\tParent=g1;gene=geneA\n",
        len = CODONS * 3
    );
    std::fs::write(dir.join("genes.gff3"), gff).expect("write gff");

    let mut vcf = std::io::BufWriter::new(
        std::fs::File::create(dir.join("cohort.vcf")).expect("create vcf"),
    );
    writeln!(vcf, "##fileformat=VCFv4.2").expect("write");
    writeln!(vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (haploid)\">")
        .expect("write");
    write!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT").expect("write");
    for i in 0..N_SAMPLES {
        write!(vcf, "\t{}", sample_name(i)).expect("write");
    }
    writeln!(vcf).expect("write");
    for r in records() {
        let alt_field: Vec<String> = r.alts.iter().map(|(a, _)| a.to_string()).collect();
        write!(
            vcf,
            "chr1\t{}\t.\t{}\t{}\t60\tPASS\tDP=50\tGT",
            r.pos,
            r.ref_base,
            alt_field.join(",")
        )
        .expect("write");
        for s in 0..N_SAMPLES {
            let gt = r
                .alts
                .iter()
                .position(|(_, carriers)| carriers.contains(&s))
                .map_or(0, |i| i + 1);
            write!(vcf, "\t{gt}").expect("write");
        }
        writeln!(vcf).expect("write");
    }
    drop(vcf);

    std::fs::write(dir.join("tree.nwk"), newick()).expect("write tree");
}

/// Run `eskaks vcf` on the constructed cohort with the given extra flags.
fn run(dir: &Path, prefix: &str, extra: &[&str]) -> std::process::Output {
    let mut cmd = Command::new(bin());
    cmd.args([
        "vcf",
        "--ref",
        dir.join("ref.fasta").to_str().expect("path"),
        "--gff",
        dir.join("genes.gff3").to_str().expect("path"),
        "--vcf",
        dir.join("cohort.vcf").to_str().expect("path"),
        "--genetic-code",
        "11",
        "--workers",
        "1",
        "--codon-scan",
        "--variants",
        "-o",
        dir.join(prefix).to_str().expect("path"),
    ]);
    cmd.args(extra);
    cmd.output().expect("spawn eskaks")
}

/// A delimited table as a list of (column name → value) maps.
fn table(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().expect("header").split('\t').map(str::to_string).collect();
    lines
        .map(|l| header.iter().cloned().zip(l.split('\t').map(str::to_string)).collect())
        .collect()
}

/// The codon-scan row for a residue.
fn row(rows: &[HashMap<String, String>], aa_pos: usize) -> HashMap<String, String> {
    rows.iter()
        .find(|r| r["AA_Pos"] == aa_pos.to_string())
        .unwrap_or_else(|| panic!("no codon-scan row for residue {aa_pos}"))
        .clone()
}

fn num(r: &HashMap<String, String>, col: &str) -> f64 {
    r[col].parse::<f64>().unwrap_or_else(|_| panic!("{col} = {:?} is not a number", r[col]))
}

#[test]
fn origins_separate_a_clonal_allele_from_a_convergent_one() {
    let dir = tempfile::tempdir().expect("tempdir");
    write_inputs(dir.path());
    let out = run(dir.path(), "tree", &["--tree", dir.path().join("tree.nwk").to_str().expect("path")]);
    assert!(out.status.success(), "run failed: {}", String::from_utf8_lossy(&out.stderr));

    let codons = table(&dir.path().join("tree_codons.tsv"));
    let (clonal, artefact, convergent, multi) =
        (row(&codons, 100), row(&codons, 150), row(&codons, 200), row(&codons, 250));

    // ── The construction, read back ─────────────────────────────────────────────
    // Every one of these codons has the same opportunity, so nothing below can be
    // explained by one of them being easier to mutate than another.
    for r in [&clonal, &artefact, &convergent, &multi] {
        assert_eq!(r["Poss_Nonsyn"], "6", "every GCT codon has A_c = 6");
        assert_eq!(r["Cooccurring"], "false", "no codon here is one multi-nucleotide event");
    }

    // CLONAL: 8 carriers, one origin. The carrier count is the second largest in the
    // cohort; the origin count is the smallest possible for an observed allele.
    assert_eq!(clonal["Nonsyn_Alleles"], "1");
    assert_eq!(clonal["Carriers_Max"], "8");
    assert_eq!(clonal["Nonsyn_Origins"], "1", "one clade is one origin, whatever its size");
    assert_eq!(clonal["Max_Allele_Origins"], "1");

    // CONVERGENT: the same single allele, 16 carriers, eight origins.
    assert_eq!(convergent["Nonsyn_Alleles"], "1", "ONE allele: the allele count cannot see this");
    assert_eq!(convergent["Carriers_Max"], "16");
    assert_eq!(convergent["Nonsyn_Origins"], "8", "one sister pair per lineage, eight lineages");
    assert_eq!(convergent["Max_Allele_Origins"], "8");

    // ARTEFACT: five scattered singleton calls. No origin subtends two carriers, so
    // under the default --min-origin-support 2 the codon scores zero origins.
    assert_eq!(artefact["Nonsyn_Alleles"], "1");
    assert_eq!(artefact["Carriers_Max"], "5");
    assert_eq!(artefact["Nonsyn_Origins"], "0", "singleton calls cannot support an origin");

    // MULTI-ALLELIC: four alleles, one origin each, so E_c == X_c exactly. This is the
    // reduction that makes the origin statistic a generalisation and not a rival.
    assert_eq!(multi["Nonsyn_Alleles"], "4");
    assert_eq!(multi["Nonsyn_Origins"], "4");
    assert_eq!(multi["Max_Allele_Origins"], "1");

    // ── The two nulls ───────────────────────────────────────────────────────────
    // theta = 47 distinct alleles (40 background + clonal + artefact + convergent + the
    // 4 of codon 250) over 1794 possible changes; lambda = 53 origins over the same
    // denominator (the artefact contributes 0). Both are printed in the run summary, so
    // both can be recomputed by hand from the table.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(&format!("{} codons, {} possible nonsyn changes", FAMILY_M, FAMILY_POSS)),
        "the family must be the whole coding gene: {stderr}"
    );
    assert!(stderr.contains("Distinct nonsyn alleles: 47"), "47 distinct alleles: {stderr}");
    assert!(stderr.contains("Origins (min support 2): 53"), "53 supported origins: {stderr}");

    let theta = 47.0 / FAMILY_POSS as f64;
    let lambda = 53.0 / FAMILY_POSS as f64;
    let exp_alleles = A_C as f64 * theta;
    let exp_origins = A_C as f64 * lambda;
    for r in [&clonal, &artefact, &convergent, &multi] {
        assert!((num(r, "Exp_Nonsyn_Alleles") - exp_alleles).abs() < 1e-6);
        assert!((num(r, "Exp_Nonsyn_Origins") - exp_origins).abs() < 1e-6);
    }

    // THE HEADLINE. The convergent codon is invisible to the allele-count test (it has
    // exactly one allele, like thousands of ordinary codons) and overwhelming on the
    // origin test. The clonal codon, with half the carriers and the same single allele,
    // is unremarkable on both. That is the whole claim, in four numbers.
    assert!(num(&convergent, "P_Recurrence") > 0.1, "one allele is never surprising by itself");
    assert_eq!(convergent["Q_Recurrence_BH"], "1.000000", "the allele-count test misses it");
    assert!(
        num(&convergent, "P_Origins") < 1e-9,
        "eight origins where 0.18 were expected: {}",
        convergent["P_Origins"]
    );
    assert!(
        num(&convergent, "Q_Origins") < 0.05,
        "and it must survive BH over the whole gene: {}",
        convergent["Q_Origins"]
    );

    assert!(num(&clonal, "P_Recurrence") > 0.1);
    assert!(num(&clonal, "P_Origins") > 0.1, "one origin is what the null expects");
    assert!(num(&clonal, "Q_Origins") > 0.05, "clonal expansion must NOT look like selection");
    assert!(
        num(&clonal, "Carriers_Max") >= num(&convergent, "Carriers_Max") / 2.0,
        "the two are comparable in carriers, so carriers cannot be doing the work"
    );

    // The artefact scores the largest p-value possible: P(E >= 0) is exactly 1.
    assert_eq!(artefact["P_Origins"], "1.000000");

    // Where the codon really did collect four alleles, both tests fire, and the origin
    // test says nothing new: E_c == X_c, so the two p-values differ only through their
    // rates. Neither replaces the other.
    assert!(num(&multi, "Q_Recurrence_BH") < 0.05);
    assert!(num(&multi, "Q_Origins") < 0.05);

    // ── Ranking ─────────────────────────────────────────────────────────────────
    // The table is scanned from the top, so the convergent codon has to be at the top.
    assert_eq!(codons[0]["AA_Pos"], "200", "the convergent codon must rank first");
    assert_eq!(codons[1]["AA_Pos"], "250", "then the four-allele codon");

    // ── The per-variant row ─────────────────────────────────────────────────────
    // This is where rpoB S450L would finally state its case: one row, one allele, eight
    // independent origins.
    let variants = table(&dir.path().join("tree_variants.tsv"));
    let by_pos = |pos: usize, alt: &str| {
        variants
            .iter()
            .find(|v| v["Pos"] == pos.to_string() && v["Alt"] == alt)
            .unwrap_or_else(|| panic!("no variant row at {pos} {alt}"))
            .clone()
    };
    assert_eq!(by_pos(codon_pos(200), "T")["Origins"], "8");
    assert_eq!(by_pos(codon_pos(100), "C")["Origins"], "1");
    assert_eq!(by_pos(codon_pos(150), "T")["Origins"], "0");
    // And every background allele arose once, which is what makes the rate a baseline.
    for c in 10..50 {
        assert_eq!(by_pos(codon_pos(c), "A")["Origins"], "1", "background codon {c}");
    }

    eprintln!(
        "CONSTRUCTED CASE  theta={theta:.6} lambda={lambda:.6}\n  \
         clonal      alleles=1 carriers=8  origins={} P_Recurrence={} P_Origins={} Q_Origins={}\n  \
         convergent  alleles=1 carriers=16 origins={} P_Recurrence={} P_Origins={} Q_Origins={}\n  \
         artefact    alleles=1 carriers=5  origins={} P_Origins={}\n  \
         four-allele alleles=4 origins={} Q_Recurrence={} Q_Origins={}",
        clonal["Nonsyn_Origins"], clonal["P_Recurrence"], clonal["P_Origins"], clonal["Q_Origins"],
        convergent["Nonsyn_Origins"], convergent["P_Recurrence"], convergent["P_Origins"],
        convergent["Q_Origins"],
        artefact["Nonsyn_Origins"], artefact["P_Origins"],
        multi["Nonsyn_Origins"], multi["Q_Recurrence_BH"], multi["Q_Origins"],
    );
}

#[test]
fn min_origin_support_is_what_stands_between_noise_and_a_genome_wide_hit() {
    // Five scattered false calls at one codon are five gains under raw parsimony. On
    // this family that is a BH-significant hit built entirely out of sequencing error,
    // and the ONLY thing that stops it is the support threshold.
    let dir = tempfile::tempdir().expect("tempdir");
    write_inputs(dir.path());
    let tree = dir.path().join("tree.nwk");
    let tree = tree.to_str().expect("path");

    let strict = run(dir.path(), "strict", &["--tree", tree]);
    assert!(strict.status.success(), "{}", String::from_utf8_lossy(&strict.stderr));
    let loose = run(dir.path(), "loose", &["--tree", tree, "--min-origin-support", "1"]);
    assert!(loose.status.success(), "{}", String::from_utf8_lossy(&loose.stderr));

    let strict_row = row(&table(&dir.path().join("strict_codons.tsv")), 150);
    let loose_row = row(&table(&dir.path().join("loose_codons.tsv")), 150);

    assert_eq!(strict_row["Nonsyn_Origins"], "0", "no origin subtends two carriers");
    assert_eq!(strict_row["Q_Origins"], "1.000000", "so the artefact is not a hit");
    assert_eq!(loose_row["Nonsyn_Origins"], "5", "raw parsimony counts every singleton");
    assert!(
        num(&loose_row, "Q_Origins") < 0.05,
        "and it becomes genome-wide significant: {}",
        loose_row["Q_Origins"]
    );

    // The real convergent signal survives the threshold either way: the filter removes
    // the failure mode without removing the feature.
    let strict_conv = row(&table(&dir.path().join("strict_codons.tsv")), 200);
    let loose_conv = row(&table(&dir.path().join("loose_codons.tsv")), 200);
    assert_eq!(strict_conv["Nonsyn_Origins"], "8");
    assert_eq!(loose_conv["Nonsyn_Origins"], "8");
    assert!(num(&strict_conv, "Q_Origins") < 0.05 && num(&loose_conv, "Q_Origins") < 0.05);

    eprintln!(
        "MIN-ORIGIN-SUPPORT  artefact codon 150: support 2 -> {} origins (Q {}), \
         support 1 -> {} origins (Q {})",
        strict_row["Nonsyn_Origins"], strict_row["Q_Origins"],
        loose_row["Nonsyn_Origins"], loose_row["Q_Origins"],
    );
}

#[test]
fn the_columns_exist_only_under_tree() {
    // The promise the goldens rest on: without --tree the tables are exactly what they
    // always were, column for column.
    let dir = tempfile::tempdir().expect("tempdir");
    write_inputs(dir.path());
    let plain = run(dir.path(), "plain", &[]);
    assert!(plain.status.success(), "{}", String::from_utf8_lossy(&plain.stderr));
    let tree = run(
        dir.path(),
        "withtree",
        &["--tree", dir.path().join("tree.nwk").to_str().expect("path")],
    );
    assert!(tree.status.success(), "{}", String::from_utf8_lossy(&tree.stderr));

    for (table_name, added) in [
        ("codons", vec![
            "Nonsyn_Origins",
            "Max_Allele_Origins",
            "Exp_Nonsyn_Origins",
            "P_Origins",
            "Q_Origins",
        ]),
        ("variants", vec!["Origins"]),
    ] {
        let plain_head = header(&dir.path().join(format!("plain_{table_name}.tsv")));
        let tree_head = header(&dir.path().join(format!("withtree_{table_name}.tsv")));
        for col in &added {
            assert!(!plain_head.contains(&col.to_string()), "{col} must be absent without --tree");
        }
        let mut want = plain_head.clone();
        want.extend(added.iter().map(|s| s.to_string()));
        assert_eq!(tree_head, want, "--tree may only APPEND columns to {table_name}");
    }

    // Row for row, the shared columns are identical: the flag adds information, it does
    // not change any answer already being given.
    let plain_rows = table(&dir.path().join("plain_codons.tsv"));
    let tree_rows = table(&dir.path().join("withtree_codons.tsv"));
    assert_eq!(plain_rows.len(), tree_rows.len());
    // Ranking may differ (the tree gives a second p-value to rank on), so compare rows
    // by residue rather than by position in the file.
    for p in &plain_rows {
        let t = row(&tree_rows, p["AA_Pos"].parse().expect("residue"));
        for (col, value) in p {
            assert_eq!(&t[col], value, "residue {} column {col} moved", p["AA_Pos"]);
        }
    }
}

#[test]
fn swapping_two_tip_labels_changes_the_answer() {
    // The tips are matched to samples by NAME. Swapping two labels leaves the topology
    // and the carriers untouched and changes only which sample sits where, so if the
    // origin counts did not move, the join would be reading positions instead.
    //
    // S01 is one of the convergent allele's carriers and sits in lineage 0's first
    // sister pair; S07 is in a different pair of the same lineage. After the swap the
    // two carriers in lineage 0 are no longer sisters, so that lineage's supported
    // origin disappears and E_c drops from 8 to 7.
    let dir = tempfile::tempdir().expect("tempdir");
    write_inputs(dir.path());
    let swapped =
        newick().replacen("S01", "@TMP@", 1).replacen("S07", "S01", 1).replace("@TMP@", "S07");
    assert_ne!(swapped, newick(), "the swap must actually change the file");
    std::fs::write(dir.path().join("swapped.nwk"), &swapped).expect("write tree");

    let base =
        run(dir.path(), "base", &["--tree", dir.path().join("tree.nwk").to_str().expect("path")]);
    assert!(base.status.success(), "{}", String::from_utf8_lossy(&base.stderr));
    let moved = run(
        dir.path(),
        "swapped",
        &["--tree", dir.path().join("swapped.nwk").to_str().expect("path")],
    );
    assert!(moved.status.success(), "{}", String::from_utf8_lossy(&moved.stderr));

    let base_conv = row(&table(&dir.path().join("base_codons.tsv")), 200);
    let moved_conv = row(&table(&dir.path().join("swapped_codons.tsv")), 200);
    assert_eq!(base_conv["Nonsyn_Origins"], "8");
    assert_eq!(moved_conv["Nonsyn_Origins"], "7", "one pair is no longer a pair");
}

fn header(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .lines()
        .next()
        .expect("header")
        .split('\t')
        .map(str::to_string)
        .collect()
}

#[test]
fn a_tree_that_does_not_match_the_cohort_is_a_hard_error() {
    // Silently dropping either side would change every origin count while changing
    // nothing a reader can see, so both directions have to be fatal, and the message has
    // to name the counts and some examples.
    let dir = tempfile::tempdir().expect("tempdir");
    write_inputs(dir.path());

    let cases: [(&str, &str, &str); 3] = [
        ("missing.nwk", "((S01,S02),(S03,S04));", "have no tip"),
        (
            "extra.nwk",
            &format!("({},GHOST_1,GHOST_2);", (0..N_SAMPLES).map(sample_name).collect::<Vec<_>>().join(",")),
            "match no sample",
        ),
        ("broken.nwk", "((S01,S02),(S03,S04);", "Newick"),
    ];
    for (name, text, needle) in cases {
        std::fs::write(dir.path().join(name), text).expect("write tree");
        let out = run(
            dir.path(),
            "bad",
            &["--tree", dir.path().join(name).to_str().expect("path")],
        );
        assert!(!out.status.success(), "{name} must fail the run");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(err.contains(needle), "{name}: expected {needle:?} in:\n{err}");
    }

    // A run whose tree is fine still succeeds, so the checks above are not just
    // "everything fails".
    let ok = run(dir.path(), "ok", &["--tree", dir.path().join("tree.nwk").to_str().expect("path")]);
    assert!(ok.status.success(), "{}", String::from_utf8_lossy(&ok.stderr));
}

/// The bundled toy genome, so the documented example is exercised end to end and the
/// answer it gives is pinned rather than described.
#[test]
fn the_bundled_toy_genome_runs_and_says_nothing_significant() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let prefix = std::env::temp_dir().join("eskaks_toy_tree_run");
    let out = Command::new(bin())
        .current_dir(&manifest)
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
            "-o",
            prefix.to_str().expect("path"),
        ])
        .output()
        .expect("spawn eskaks");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let err = String::from_utf8_lossy(&out.stderr);

    // 20 samples and 134 distinct nonsynonymous alleles come to 142 supported origins:
    // a handful of the toy cohort's alleles sit in two clades of this tree. That is the
    // whole of what 20 genomes can say, and the honest answer is that no residue clears
    // BH over a 1,421-codon family on either test.
    assert!(err.contains("Origins (min support 2): 142"), "{err}");
    let rows = table(&prefix.with_file_name(format!(
        "{}_codons.tsv",
        prefix.file_name().expect("prefix").to_string_lossy()
    )));
    let significant = rows
        .iter()
        .filter(|r| r["Q_Origins"].parse::<f64>().is_ok_and(|q| q < 0.05))
        .count();
    assert_eq!(significant, 0, "20 samples cannot establish convergence");
    // The ranking is still informative: the top row is a single allele that arose twice.
    assert_eq!(rows[0]["Nonsyn_Origins"], "2");
    assert_eq!(rows[0]["Nonsyn_Alleles"], "1", "an allele count of 1 that the tree splits in two");
}
