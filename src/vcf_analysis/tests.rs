use super::*;

fn make_gc() -> &'static GeneticCode {
    crate::genetic_code::get_table(1).unwrap()
}

#[test]
fn test_codon_to_aa() {
    let gc = make_gc();
    // ATG = Met
    assert_eq!(codon_to_aa(b"ATG", gc), Some(b'M'));
    // TAA = Stop
    assert_eq!(codon_to_aa(b"TAA", gc), Some(b'*'));
    // GCT = Ala
    assert_eq!(codon_to_aa(b"GCT", gc), Some(b'A'));
}

#[test]
fn test_complement() {
    assert_eq!(complement(b'A'), b'T');
    assert_eq!(complement(b'T'), b'A');
    assert_eq!(complement(b'C'), b'G');
    assert_eq!(complement(b'G'), b'C');
}

#[test]
fn test_reverse_complement() {
    assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
    assert_eq!(reverse_complement(b"AATTCC"), b"GGAATT");
}

#[test]
fn test_count_sites() {
    let gc = make_gc();
    // ATG GCT = Met Ala (2 codons)
    let cds = b"ATGGCT";
    let (n, s) = count_sites(cds, gc);
    // Each codon contributes exactly 3 sites total (S + N = 3 per codon, 6 total)
    let total = n + s;
    assert!((total - 6.0).abs() < 0.01, "Total sites should be 6.0, got {}", total);
    // ATG (Met): 0 syn, 8 nonsyn (1 change → stop TAA excluded) → N=3*8/8=3, S=0
    // GCT (Ala): 3 syn (GCA,GCC,GCG), 6 nonsyn, 0 stops → N=3*6/9=2, S=3*3/9=1
    assert!(n > 4.5 && n < 5.5, "N_sites: expected ~5.0, got {}", n);
    assert!(s > 0.5 && s < 1.5, "S_sites: expected ~1.0, got {}", s);
}

#[test]
fn test_is_transition() {
    assert!(is_transition(b'A', b'G'));
    assert!(is_transition(b'G', b'A'));
    assert!(is_transition(b'C', b'T'));
    assert!(is_transition(b'T', b'C'));
    assert!(!is_transition(b'A', b'C')); // transversion
    assert!(!is_transition(b'A', b'T')); // transversion
    assert!(!is_transition(b'G', b'T')); // transversion
}

#[test]
fn test_weighted_reduces_to_unweighted_at_kappa_1() {
    let gc = make_gc();
    // At kappa == 1 the weighted count must equal the unweighted count for
    // EVERY codon — including stop-adjacent ones (TAC, TCA), where changes
    // to stop codons are excluded — since both use codon-level pooling.
    for cds in [
        &b"ATGGCT"[..],       // Met Ala — no stop-adjacent changes
        &b"TACTGCTCACAA"[..], // Tyr Cys Ser Gln — TAC/TCA are stop-adjacent
        &b"TTTTATCATAAT"[..], // Phe Tyr His Asn — 2-fold codons
    ] {
        let (n0, s0) = count_sites(cds, gc);
        let (n1, s1) = count_sites_weighted(cds, gc, 1.0);
        assert!((n0 - n1).abs() < 1e-9, "N: unweighted {} vs weighted@1 {} for {:?}", n0, n1, cds);
        assert!((s0 - s1).abs() < 1e-9, "S: unweighted {} vs weighted@1 {} for {:?}", s0, s1, cds);
    }
}

#[test]
fn test_weighted_two_fold_site_follows_kappa_formula() {
    let gc = make_gc();
    // TTT = Phe. Only the 3rd position is (partly) synonymous: TTC (Phe) is a
    // transition (T->C); TTA/TTG (Leu) are transversions. Positions 1 and 2
    // are fully nonsynonymous with no stop changes.
    // So S = kappa/(kappa+2), N = 3 - S, for any kappa.
    for &k in &[0.5f64, 1.0, 2.0, 5.0] {
        let (n, s) = count_sites_weighted(b"TTT", gc, k);
        let expected_s = k / (k + 2.0);
        assert!((s - expected_s).abs() < 1e-9, "kappa={}: S expected {}, got {}", k, expected_s, s);
        assert!((n - (3.0 - expected_s)).abs() < 1e-9, "kappa={}: N expected {}, got {}", k, 3.0 - expected_s, n);
    }
}

#[test]
fn test_weighted_four_fold_site_is_kappa_invariant() {
    let gc = make_gc();
    // GCT = Ala, GCN all Ala: the 3rd position is 4-fold degenerate, so it is
    // fully synonymous regardless of ts/tv weighting. Its synonymous-site
    // contribution must be exactly 1.0 for every kappa.
    let s_at = |k: f64| count_sites_weighted(b"GCT", gc, k).1;
    let base = s_at(1.0);
    for &k in &[0.5f64, 2.0, 10.0] {
        assert!((s_at(k) - base).abs() < 1e-9, "4-fold S must be kappa-invariant (k={})", k);
    }
    // GCT: pos1 & pos2 fully nonsyn, pos3 fully syn → S == 1.0 exactly.
    assert!((base - 1.0).abs() < 1e-9, "GCT S should be 1.0, got {}", base);
}

#[test]
fn test_weighted_transition_bias_raises_s_lowers_n() {
    let gc = make_gc();
    // Across a mixed CDS, kappa>1 must raise total S and lower total N
    // relative to the equal-rates count (the whole point of the correction).
    let cds = b"TTTGCTAAACGT"; // Phe Ala Lys Arg
    let (n1, s1) = count_sites_weighted(cds, gc, 1.0);
    let (n5, s5) = count_sites_weighted(cds, gc, 5.0);
    assert!(s5 > s1, "S should rise with kappa: {} -> {}", s1, s5);
    assert!(n5 < n1, "N should fall with kappa: {} -> {}", n1, n5);
    // Sites are conserved: each position is exactly one site.
    assert!((n1 + s1 - (n5 + s5)).abs() < 1e-9, "total sites must be kappa-invariant");
}

#[test]
fn test_genomic_to_cds_offset_plus_strand() {
    let gene = Gene {
        name: "test".to_string(),
        seqid: "chr1".to_string(),
        strand: Strand::Plus,
        exons: vec![crate::gff::CdsExon {
            seqid: "chr1".to_string(),
            start: 100,
            end: 109,
            strand: Strand::Plus,
            phase: 0,
        }],
        length_bp: 10,
        start: 100,
    };

    assert_eq!(genomic_to_cds_offset(&gene, 100), Some(0));
    assert_eq!(genomic_to_cds_offset(&gene, 105), Some(5));
    assert_eq!(genomic_to_cds_offset(&gene, 109), Some(9));
    assert_eq!(genomic_to_cds_offset(&gene, 110), None);
    assert_eq!(genomic_to_cds_offset(&gene, 99), None);
}

/// Build a minimal `GenePnPs` for aggregation tests. Only the fields read
/// by `genome_wide_pn_ps` need to be meaningful.
fn gene_result(n_sites: f64, s_sites: f64, nonsyn: f64, syn: f64) -> GenePnPs {
    let pn = if n_sites > 0.0 { nonsyn / n_sites } else { 0.0 };
    let ps = if s_sites > 0.0 { syn / s_sites } else { 0.0 };
    let pn_ps = if ps > 0.0 {
        pn / ps
    } else if pn > 0.0 {
        f64::INFINITY
    } else {
        f64::NAN
    };
    GenePnPs {
        name: "g".to_string(),
        length_bp: 0,
        n_sites,
        s_sites,
        pn,
        ps,
        pn_ps,
        nonsyn_snps: nonsyn,
        syn_snps: syn,
        total_snps: nonsyn + syn,
        n_snps: (nonsyn + syn) as usize,
        genome_start: 0,
        genome_end: 0,
        strand: '+',
        chrom: "chr1".to_string(),
        p_value: f64::NAN,
        q_value: f64::NAN,
        p_bonferroni: f64::NAN,
        mk_dn: 0,
        mk_ds: 0,
        mk_pn: 0,
        mk_ps: 0,
        neglog10p: f64::NAN,
        pn_ps_lo: f64::NAN,
        pn_ps_hi: f64::NAN,
        p_gc: f64::NAN,
        q_gc: f64::NAN,
        repetitive: false,
        sfs_nonsyn: [0; SFS_NBINS],
        sfs_syn: [0; SFS_NBINS],
        variants: Vec::new(),
        scan_codons: 0,
        scan_poss_nonsyn: 0,
    }
}

#[test]
fn test_genome_wide_pools_counts_not_ratios() {
    // Gene 1: many sites, ratio 0.5. Gene 2: few sites, ratio 4.0.
    // Averaging ratios would give 2.25; pooling counts must not.
    let results = vec![
        gene_result(100.0, 100.0, 10.0, 20.0), // pN=0.10, pS=0.20
        gene_result(2.0, 2.0, 4.0, 1.0),       // pN=2.00, pS=0.50
    ];
    let gw = genome_wide_pn_ps(&results).expect("should aggregate");
    assert!((gw.n_sites - 102.0).abs() < 1e-9);
    assert!((gw.s_sites - 102.0).abs() < 1e-9);
    assert!((gw.nonsyn_snps - 14.0).abs() < 1e-9);
    assert!((gw.syn_snps - 21.0).abs() < 1e-9);
    // pooled pN = 14/102, pS = 21/102, ratio = 14/21
    assert!((gw.pn - 14.0 / 102.0).abs() < 1e-9);
    assert!((gw.ps - 21.0 / 102.0).abs() < 1e-9);
    assert!((gw.pn_ps - 14.0 / 21.0).abs() < 1e-9, "got {}", gw.pn_ps);
}

#[test]
fn test_genome_wide_empty_is_none() {
    assert!(genome_wide_pn_ps(&[]).is_none());
}

#[test]
fn test_genome_wide_no_variation_is_nan() {
    // Genes exist but carry no SNPs → ratio is NaN, not 0 or Inf.
    let results = vec![gene_result(50.0, 25.0, 0.0, 0.0)];
    let gw = genome_wide_pn_ps(&results).unwrap();
    assert_eq!(gw.pn, 0.0);
    assert_eq!(gw.ps, 0.0);
    assert!(gw.pn_ps.is_nan(), "expected NaN, got {}", gw.pn_ps);
}

#[test]
fn test_genome_wide_no_syn_is_infinite() {
    // Nonsynonymous variation but zero synonymous → +Infinity.
    let results = vec![gene_result(50.0, 25.0, 3.0, 0.0)];
    let gw = genome_wide_pn_ps(&results).unwrap();
    assert!(gw.pn_ps.is_infinite() && gw.pn_ps > 0.0);
}

#[test]
fn test_selection_label_bands() {
    assert!(selection_label(0.3).contains("purifying"));
    assert!(selection_label(1.0).contains("neutral"));
    assert!(selection_label(2.0).contains("diversifying"));
    assert!(selection_label(f64::NAN).contains("undetermined"));
    assert!(selection_label(f64::INFINITY).contains("synonymous"));
}

#[test]
fn test_format_ratio_special_values() {
    assert_eq!(format_ratio(f64::NAN), "NaN");
    assert_eq!(format_ratio(f64::INFINITY), "inf");
    assert_eq!(format_ratio(0.5), "0.500000");
}

#[test]
fn test_genomic_to_cds_offset_minus_strand() {
    let gene = Gene {
        name: "test".to_string(),
        seqid: "chr1".to_string(),
        strand: Strand::Minus,
        exons: vec![crate::gff::CdsExon {
            seqid: "chr1".to_string(),
            start: 100,
            end: 109,
            strand: Strand::Minus,
            phase: 0,
        }],
        length_bp: 10,
        start: 100,
    };

    // Minus strand: position 109 maps to CDS offset 0
    assert_eq!(genomic_to_cds_offset(&gene, 109), Some(0));
    assert_eq!(genomic_to_cds_offset(&gene, 100), Some(9));
}

#[test]
fn test_multi_exon_minus_strand_cds_consistency() {
    // Regression: for a multi-exon minus gene, extract_cds_sequence and
    // genomic_to_cds_offset must agree — the CDS base at a SNP's offset must be
    // the complement of the reference base there (the coding base). The old code
    // assembled the exons in swapped order, breaking this and mis-codon'ing SNPs.
    let ref_seq: Vec<u8> = (0..30).map(|i| b"ACGT"[i % 4]).collect();
    let ex = |start, end| crate::gff::CdsExon {
        seqid: "chr1".to_string(), start, end, strand: Strand::Minus, phase: 0,
    };
    // Stored descending for minus strand: high exon first.
    let gene = Gene {
        name: "m".to_string(), seqid: "chr1".to_string(), strand: Strand::Minus,
        exons: vec![ex(20, 25), ex(10, 15)], length_bp: 12, start: 10,
    };
    let cds = extract_cds_sequence(&gene, &ref_seq);
    assert_eq!(cds.len(), 12);
    for pos in [25usize, 24, 20, 15, 12, 10] {
        let off = genomic_to_cds_offset(&gene, pos).expect("position is in an exon");
        let expected = complement(ref_seq[pos - 1]).to_ascii_uppercase();
        assert_eq!(cds[off], expected, "pos {pos} -> off {off} mismatch");
    }
}

#[test]
fn test_is_repetitive_flags() {
    for n in ["PE_PGRS12", "PPE18", "Rv1234_PGRS", "IS6110", "some_transposase", "PE13", "maturase", "PE", "PPE"] {
        assert!(is_repetitive(n), "{n} should be repetitive");
    }
    // Real genes whose symbol merely starts with "PE"/"PPE" must NOT be flagged.
    for n in ["Rv0001", "katG", "rpoB", "gyrA", "pepN", "pepA", "pepD", "penA"] {
        assert!(!is_repetitive(n), "{n} should be core");
    }
}

#[test]
fn test_core_repetitive_split() {
    let mut a = gene_result(100.0, 100.0, 10.0, 20.0);
    a.repetitive = false;
    let mut b = gene_result(50.0, 50.0, 30.0, 5.0);
    b.repetitive = true;
    let (core, rep) = genome_wide_core_repetitive(&[a, b]);
    let core = core.expect("core present");
    let rep = rep.expect("rep present");
    // Core pools only the first gene, repetitive only the second.
    assert!((core.nonsyn_snps - 10.0).abs() < 1e-9 && (core.syn_snps - 20.0).abs() < 1e-9);
    assert!((rep.nonsyn_snps - 30.0).abs() < 1e-9 && (rep.syn_snps - 5.0).abs() < 1e-9);
}

#[test]
fn test_genomic_control_deflates_and_labels() {
    // Build a set of tested genes with a range of p-values.
    let mut genes: Vec<GenePnPs> = Vec::new();
    for (i, p) in [1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.7, 0.9].iter().enumerate() {
        let mut g = gene_result(100.0, 100.0, (10 + i) as f64, 20.0);
        g.p_value = *p;
        genes.push(g);
    }
    apply_multiple_testing(&mut genes, false);
    let lambda = genomic_inflation_lambda(&genes);
    assert!(lambda.is_finite() && lambda > 0.0);
    apply_genomic_control(&mut genes, lambda);
    // When λ > 1, every corrected p is >= the raw p (correction only deflates);
    // corrected values are valid probabilities.
    for g in &genes {
        assert!(g.p_gc.is_finite() && (0.0..=1.0).contains(&g.p_gc));
        if lambda > 1.0 {
            assert!(g.p_gc >= g.p_value - 1e-9, "GC should not make p smaller: {} < {}", g.p_gc, g.p_value);
        }
    }
    // λ < 2 tested genes → NaN.
    assert!(genomic_inflation_lambda(&genes[..1]).is_nan());
}

#[test]
fn test_genomic_control_survives_underflowed_pvalues() {
    // Genes so significant their exact p underflows to 0 but neglog10p is finite.
    // Older code took Φ⁻¹(1 − p/2) from the raw p and got +∞; λ and GC must stay finite.
    let mut genes: Vec<GenePnPs> = Vec::new();
    for i in 0..10 {
        let mut g = gene_result(100.0, 100.0, 90.0, 5.0);
        g.p_value = 0.0; // underflowed
        g.neglog10p = 40.0 + i as f64; // finite log-space statistic
        genes.push(g);
    }
    // add a couple of moderate genes so the set is not degenerate
    for p in [0.2_f64, 0.6_f64] {
        let mut g = gene_result(100.0, 100.0, 50.0, 50.0);
        g.p_value = p;
        g.neglog10p = -p.log10();
        genes.push(g);
    }
    apply_multiple_testing(&mut genes, false);
    let lambda = genomic_inflation_lambda(&genes);
    assert!(lambda.is_finite() && lambda > 0.0, "λ must be finite, got {lambda}");
    apply_genomic_control(&mut genes, lambda);
    for g in &genes {
        assert!(!g.p_gc.is_nan(), "p_gc must not be NaN");
        assert!((0.0..=1.0).contains(&g.p_gc));
    }
}

#[test]
fn cov_genomic_lambda_median_zero_is_nan_not_zero() {
    // Regression: when the median tested-gene χ² is 0 (the majority of genes have exact
    // two-sided p = 1 → χ² = 0), λ carries no inflation signal and must be NaN, not a
    // spurious 0.0 that renders as "λ 0.00" (extreme deflation). Mirrors the < 2-tested-
    // genes NaN sentinel.
    let mut genes: Vec<GenePnPs> = Vec::new();
    for _ in 0..3 {
        let mut g = gene_result(100.0, 100.0, 1.0, 1.0);
        g.p_value = 1.0;
        g.neglog10p = 0.0; // χ² = 0
        genes.push(g);
    }
    // One strong outlier so the set has a nonzero max but still a zero median.
    let mut strong = gene_result(100.0, 100.0, 90.0, 5.0);
    strong.p_value = 1e-8;
    strong.neglog10p = 8.0;
    genes.push(strong);
    apply_multiple_testing(&mut genes, false);
    assert!(
        genomic_inflation_lambda(&genes).is_nan(),
        "a median χ² of 0 must give NaN λ, not 0.0"
    );
}

#[test]
fn test_sfs_bins_edges() {
    assert_eq!(sfs_bin(0.05), 0);
    assert_eq!(sfs_bin(0.1), 0);
    assert_eq!(sfs_bin(0.15), 1);
    assert_eq!(sfs_bin(1.0), SFS_NBINS - 1);
    assert_eq!(sfs_bin(0.99), SFS_NBINS - 1);
}

// ---- output-writer / formatter coverage (prefix cov_out_) ----

/// Structural + numeric-safety invariants every non-empty SVG we emit must
/// hold. Mirrors the assertions used in plot.rs's own SVG tests, kept local
/// (a differently-named helper, so it never collides on merge).
fn cov_out_assert_svg(svg: &str) {
    assert!(svg.starts_with("<?xml"), "missing XML prolog");
    assert!(svg.contains("<svg "), "missing <svg> root");
    assert!(svg.trim_end().ends_with("</svg>"), "not closed with </svg>");
    assert!(!svg.contains("NaN"), "NaN leaked into SVG output");
    for bad in ["=\"inf\"", "=\"-inf\"", "=\"NaN\""] {
        assert!(!svg.contains(bad), "non-finite value in attribute ({bad})");
    }
}

#[test]
fn cov_out_parse_reference_fasta_roundtrip() {
    // A valid multi-record FASTA must round-trip to a map keyed by the FIRST
    // whitespace-delimited token after '>', with sequence bytes uppercased and
    // concatenated across wrapped lines.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("ref.fasta");
    std::fs::write(&path, ">seq1 some description\nACGT\nacgt\n>seq2\nTTTT\n").unwrap();

    let map = parse_reference_fasta(&path).unwrap();
    assert_eq!(map.len(), 2);
    // seq1: ACGT + acgt -> uppercased & concatenated = ACGTACGT
    assert_eq!(map.get("seq1").unwrap(), &b"ACGTACGT".to_vec());
    assert_eq!(map.get("seq2").unwrap(), &b"TTTT".to_vec());
    // Only the first token is used as the key, so "seq1 some description" is absent.
    assert!(!map.contains_key("seq1 some description"));
}

#[test]
fn cov_out_parse_reference_fasta_missing_file_is_err() {
    let missing = std::path::Path::new("/nonexistent_dir_cov_out_zzz/nofile.fasta");
    assert!(parse_reference_fasta(missing).is_err());
}

#[test]
fn cov_out_parse_reference_fasta_duplicate_contig_is_err() {
    // A repeated contig id would otherwise silently overwrite (last wins), scoring every
    // gene on that contig against the wrong sequence. It must be rejected.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("dup.fasta");
    std::fs::write(&path, ">chr1\nACGTACGT\n>chr1\nTTTTAAAA\n").unwrap();
    assert!(parse_reference_fasta(&path).is_err(), "duplicate contig id must be rejected");
}

#[test]
fn cov_core_rna_reference_matches_dna() {
    // A U (RNA) reference must yield the same N/S site counts as its DNA equivalent, since
    // the tool defines U == T. Regression for count_sites' self-skip comparing raw bytes.
    let gc = make_gc();
    let gene_d = cov_core_gene("g", "chr1", Strand::Plus, 1, 9, 0);
    let gene_r = cov_core_gene("g", "chr1", Strand::Plus, 1, 9, 0);
    let mut dna = HashMap::new();
    dna.insert("chr1".to_string(), b"TTTAAACCC".to_vec()); // Phe Lys Pro
    let mut rna = HashMap::new();
    rna.insert("chr1".to_string(), b"UUUAAACCC".to_vec()); // same, RNA
    let (rd, _) = compute_pn_ps(&dna, &[gene_d], &[], gc, false, 1.0, 0.95);
    let (rr, _) = compute_pn_ps(&rna, &[gene_r], &[], gc, false, 1.0, 0.95);
    assert_eq!(rd.len(), 1);
    assert_eq!(rr.len(), 1);
    assert!((rd[0].n_sites - rr[0].n_sites).abs() < 1e-9, "N_sites U {} vs T {}", rr[0].n_sites, rd[0].n_sites);
    assert!((rd[0].s_sites - rr[0].s_sites).abs() < 1e-9, "S_sites U {} vs T {}", rr[0].s_sites, rd[0].s_sites);
}

#[test]
fn cov_out_write_results_tsv_header_and_row() {
    // gene_result(100,100,10,20): pn=0.10, ps=0.20, pn_ps=0.5, total=30,
    // exp_n_frac = 100/200 = 0.5; p/q/bonferroni default to NaN -> "NA".
    let results = vec![gene_result(100.0, 100.0, 10.0, 20.0)];
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();

    let path = write_results(&results, prefix, &crate::models::OutputFormat::Tsv).unwrap();
    assert!(path.ends_with("_pnps.tsv"));
    let content = std::fs::read_to_string(&path).unwrap();
    let mut lines = content.lines();

    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "Gene\tLength_bp\tN_sites\tS_sites\tpN\tpS\tpN/pS\tNonsyn_SNPs\tSyn_SNPs\tTotal_SNPs\tChrom\tStart\tEnd\tStrand\tExp_N_frac\tP_value\tQ_value_BH\tP_Bonferroni\tpN/pS_lo\tpN/pS_hi\tP_GC\tQ_GC_BH"
    );

    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(row.len(), 22, "row = {row:?}");
    assert_eq!(row[0], "g");            // name
    assert_eq!(row[1], "0");            // length_bp
    assert_eq!(row[2], "100.0000");     // n_sites {:.4}
    assert_eq!(row[3], "100.0000");     // s_sites {:.4}
    assert_eq!(row[4], "0.100000");     // pn {:.6}
    assert_eq!(row[5], "0.200000");     // ps {:.6}
    assert_eq!(row[6], "0.500000");     // format_ratio(0.5)
    assert_eq!(row[7], "10.0000");      // nonsyn {:.4}
    assert_eq!(row[8], "20.0000");      // syn {:.4}
    assert_eq!(row[9], "30.0000");      // total {:.4}
    assert_eq!(row[10], "chr1");        // chrom
    assert_eq!(row[11], "0");           // start
    assert_eq!(row[12], "0");           // end
    assert_eq!(row[13], "+");           // strand
    assert_eq!(row[14], "0.500000");    // format_pval(exp_n_frac=0.5)
    assert_eq!(row[15], "NA");          // format_pval(NaN)
    assert_eq!(row[16], "NA");
    assert_eq!(row[17], "NA");
    assert_eq!(row[18], "NaN");         // pn_ps_lo: format_ratio(NaN)
    assert_eq!(row[19], "NaN");         // pn_ps_hi: format_ratio(NaN)
    assert_eq!(row[20], "NA");          // p_gc: format_pval(NaN)
    assert_eq!(row[21], "NA");          // q_gc
}

#[test]
fn cov_out_write_mk_results_tsv_header_row_and_filter() {
    // Only genes with mk_dn+mk_ds+mk_pn+mk_ps > 0 are written.
    // For (Dn,Ds,Pn,Ps) = (4,2,6,3):
    //   NI    = (Pn*Ds)/(Ps*Dn) = (6*2)/(3*4) = 12/12 = 1.0 -> "1.000000"
    //   alpha = 1 - (Ds*Pn)/(Dn*Ps) = 1 - (2*6)/(4*3) = 1 - 1 = 0.0 -> "0.000000"
    let mut g = gene_result(100.0, 100.0, 10.0, 20.0);
    g.mk_dn = 4;
    g.mk_ds = 2;
    g.mk_pn = 6;
    g.mk_ps = 3;
    // This gene has all mk counts == 0 (gene_result default) -> must be filtered out.
    let excluded = gene_result(50.0, 50.0, 0.0, 0.0);
    let results = vec![g, excluded];

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();

    let path = write_mk_results(&results, prefix, &crate::models::OutputFormat::Tsv).unwrap();
    assert!(path.ends_with("_mk.tsv"));
    let content = std::fs::read_to_string(&path).unwrap();
    let mut lines = content.lines();

    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "Gene\tChrom\tStart\tEnd\tStrand\tDn\tDs\tPn\tPs\tNI\talpha\tFisher_p\tFisher_q_BH"
    );

    let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
    // Only the informative gene is written.
    assert!(lines.next().is_none(), "filtered gene should not be written");
    assert_eq!(row.len(), 13, "row = {row:?}");
    assert_eq!(row[0], "g");
    assert_eq!(row[1], "chr1");
    assert_eq!(row[2], "0");
    assert_eq!(row[3], "0");
    assert_eq!(row[4], "+");
    assert_eq!(row[5], "4");   // Dn
    assert_eq!(row[6], "2");   // Ds
    assert_eq!(row[7], "6");   // Pn
    assert_eq!(row[8], "3");   // Ps
    assert_eq!(row[9], "1.000000");  // NI
    assert_eq!(row[10], "0.000000"); // alpha
    // Fisher p is a valid probability; with a single tested gene BH leaves q == p,
    // so the two rendered fields must be byte-identical.
    let p: f64 = row[11].parse().expect("fisher_p parses");
    assert!((0.0..=1.0).contains(&p), "fisher_p out of range: {p}");
    assert_eq!(row[11], row[12], "single-test BH: q must equal p");
}

#[test]
fn cov_out_write_pnps_plot_empty_skips_file() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();
    // No data points -> returns None and writes NO file (so the caller doesn't list a
    // phantom path in the output manifest).
    assert!(write_pnps_plot(&[], prefix, 0.05).unwrap().is_none());
    assert!(!std::path::Path::new(&format!("{}_pnps_manhattan.svg", prefix)).exists());
}

#[test]
fn cov_out_write_pnps_plot_wellformed_with_infinite_gene() {
    // A finite-ratio gene (plotted) plus an infinite-ratio gene (pS = 0, only
    // nonsynonymous variation): the infinite gene is excluded from the ratio
    // axis but the plot still renders cleanly with no NaN/inf attributes.
    let finite = gene_result(100.0, 100.0, 10.0, 20.0); // pn_ps = 0.5, total = 30
    let infinite = gene_result(100.0, 100.0, 30.0, 0.0); // pS = 0 -> pn_ps = +inf
    assert!(infinite.pn_ps.is_infinite());
    let results = vec![finite, infinite];

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();

    let path = write_pnps_plot(&results, prefix, 0.05).unwrap().expect("plot written");
    let svg = std::fs::read_to_string(&path).unwrap();
    cov_out_assert_svg(&svg);
}

#[test]
fn cov_out_write_pvalue_manhattan_empty_skips_file() {
    // gene_result leaves p_value = NaN, so no gene has a finite p-value: returns None
    // and writes NO file.
    let results = vec![gene_result(100.0, 100.0, 10.0, 20.0)];
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();
    assert!(write_pvalue_manhattan(&results, prefix, 0.05).unwrap().is_none());
    assert!(!std::path::Path::new(&format!("{}_pvalue_manhattan.svg", prefix)).exists());
}

#[test]
fn cov_out_write_pvalue_manhattan_wellformed_with_pvalues() {
    // Genes with finite p-values across several magnitudes exercise the
    // -log10(p) axis, the BH significance line, and per-point rendering.
    let mut results = Vec::new();
    for (i, p) in [1e-6_f64, 1e-3, 0.01, 0.05, 0.2].iter().enumerate() {
        let mut g = gene_result(100.0, 100.0, (10 + i) as f64, 20.0);
        g.p_value = *p;
        results.push(g);
    }
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("out");
    let prefix = prefix.to_str().unwrap();

    let path = write_pvalue_manhattan(&results, prefix, 0.05).unwrap().expect("plot written");
    let svg = std::fs::read_to_string(&path).unwrap();
    cov_out_assert_svg(&svg);
}

#[test]
fn cov_out_exp_n_frac_values() {
    // exp_n_frac = N/(N+S); undefined (NaN) when there are no sites.
    assert!((exp_n_frac(&gene_result(100.0, 100.0, 0.0, 0.0)) - 0.5).abs() < 1e-12);
    assert!((exp_n_frac(&gene_result(30.0, 10.0, 0.0, 0.0)) - 0.75).abs() < 1e-12);
    assert!(exp_n_frac(&gene_result(0.0, 0.0, 0.0, 0.0)).is_nan());
}

#[test]
fn cov_out_format_pval_strings() {
    // NaN -> "NA"; |v|<1e-3 (and v!=0) -> scientific {:.3e}; else {:.6}.
    assert_eq!(format_pval(f64::NAN), "NA");
    assert_eq!(format_pval(0.0), "0.000000");
    assert_eq!(format_pval(0.5), "0.500000");
    assert_eq!(format_pval(1e-6), "1.000e-6"); // Rust {:.3e}: no zero-padded exponent
    assert_eq!(format_pval(1e-3), "0.001000"); // 1e-3 is NOT < 1e-3, so fixed form
}

#[test]
fn cov_out_format_json_num_strings() {
    // Finite -> default Display literal (round-trips); non-finite -> "null".
    assert_eq!(format_json_num(0.0), "0");
    assert_eq!(format_json_num(0.5), "0.5");
    assert_eq!(format_json_num(1e-6), "0.000001");
    assert_eq!(format_json_num(f64::NAN), "null");
    assert_eq!(format_json_num(f64::INFINITY), "null");
    assert_eq!(format_json_num(f64::NEG_INFINITY), "null");
}

#[test]
fn cov_out_format_json_f64_strings() {
    // Non-finite -> "null"; exact zero -> fixed "0.000000"; else {:.6}.
    assert_eq!(format_json_f64(f64::NAN), "null");
    assert_eq!(format_json_f64(f64::INFINITY), "null");
    assert_eq!(format_json_f64(f64::NEG_INFINITY), "null");
    assert_eq!(format_json_f64(0.0), "0.000000");
    assert_eq!(format_json_f64(0.5), "0.500000");
    assert_eq!(format_json_f64(2.0), "2.000000");
}

#[test]
fn cov_out_format_ratio_more() {
    // Extends the existing special-values test (distinct name, no duplication).
    assert_eq!(format_ratio(f64::NEG_INFINITY), "inf"); // sign is dropped by design
    assert_eq!(format_ratio(2.5), "2.500000");
    assert_eq!(format_ratio(0.0), "0.000000");
}

// ---- cov_core_: numeric-core coverage (compute_pn_ps, bootstrap, site counting) ----

/// Single-exon gene builder for the cov_core_ end-to-end tests. Uses only the
/// public `Gene`/`CdsExon` fields, mirroring the inline construction in
/// `test_genomic_to_cds_offset_plus_strand`.
fn cov_core_gene(name: &str, seqid: &str, strand: Strand, start: usize, end: usize, phase: u8) -> Gene {
    Gene {
        name: name.to_string(),
        seqid: seqid.to_string(),
        strand,
        exons: vec![crate::gff::CdsExon {
            seqid: seqid.to_string(),
            start,
            end,
            strand,
            phase,
        }],
        length_bp: end - start + 1,
        start,
    }
}

/// Minimal single-ALT `VcfSnp` builder (all public fields populated).
fn cov_core_snp(chrom: &str, pos: usize, ref_allele: u8, alt: u8, af: f64) -> VcfSnp {
    VcfSnp {
        chrom: chrom.to_string(),
        pos,
        ref_allele,
        alt_alleles: vec![alt],
        alt_freqs: vec![af],
        gt_counts: None,
        filter: "PASS".to_string(),
        depth: None,
    }
}

#[test]
fn cov_core_compute_plus_strand_hand_checked() {
    let gc = make_gc();
    // Reference CDS = codons TTT(Phe) GCT(Ala) on the plus strand.
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"TTTGCT".to_vec());
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 6, 0);
    // pos 3: REF T -> ALT C  => TTT->TTC  Phe->Phe   SYNONYMOUS
    // pos 4: REF G -> ALT A  => GCT->ACT  Ala->Thr   NONSYNONYMOUS
    let snps = vec![
        cov_core_snp("chr1", 3, b'T', b'C', 0.8),
        cov_core_snp("chr1", 4, b'G', b'A', 0.4),
    ];
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.5);
    assert_eq!(res.len(), 1);
    let g = &res[0];

    // Empirical numerators: exactly one syn and one nonsyn ALT.
    assert_eq!(g.nonsyn_snps, 1.0);
    assert_eq!(g.syn_snps, 1.0);
    assert_eq!(g.total_snps, 2.0);
    assert_eq!(g.strand, '+');

    // Nei-Gojobori site counts (single-nt changes, stop changes excluded):
    //   TTT(Phe): pos3 T->C is the only syn change; syn=1, nonsyn=8
    //             => S=3*1/9=1/3, N=3*8/9=8/3
    //   GCT(Ala): 3rd position 4-fold degenerate; syn=3, nonsyn=6
    //             => S=3*3/9=1, N=3*6/9=2
    //   totals: N = 8/3 + 2 = 14/3,  S = 1/3 + 1 = 4/3
    assert!((g.n_sites - 14.0 / 3.0).abs() < 1e-9, "n_sites {}", g.n_sites);
    assert!((g.s_sites - 4.0 / 3.0).abs() < 1e-9, "s_sites {}", g.s_sites);
    // pn = 1/(14/3) = 3/14 ; ps = 1/(4/3) = 3/4 ; pn/ps = (3/14)/(3/4) = 2/7
    assert!((g.pn - 3.0 / 14.0).abs() < 1e-9, "pn {}", g.pn);
    assert!((g.ps - 3.0 / 4.0).abs() < 1e-9, "ps {}", g.ps);
    assert!((g.pn_ps - 2.0 / 7.0).abs() < 1e-9, "pn_ps {}", g.pn_ps);

    // MK split at mk_fixed_af=0.5: syn af 0.8 >= 0.5 -> Ds; nonsyn af 0.4 < 0.5 -> Pn.
    assert_eq!(g.mk_ds, 1);
    assert_eq!(g.mk_pn, 1);
    assert_eq!(g.mk_dn, 0);
    assert_eq!(g.mk_ps, 0);

    // SFS_EDGES = [0.1,0.2,0.4,0.6,0.8,1.0]; af 0.8 -> bin 4 (syn), af 0.4 -> bin 2 (nonsyn).
    assert_eq!(g.sfs_syn[4], 1);
    assert_eq!(g.sfs_nonsyn[2], 1);

    // Non-AF-weighted with SNPs present: neutrality test is defined.
    assert!(g.p_value.is_finite());
}

#[test]
fn cov_core_zero_synonymous_sites_gives_nan_ps_not_zero() {
    // A translatable gene of only Met (ATG) and Trp (TGG) codons has zero synonymous
    // sites, so pS is 0/0 (undefined). With a nonsynonymous SNP present it must report
    // pS = NaN, not a spurious 0.0 that reads as an observed synonymous rate of zero.
    let gc = make_gc();
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"ATGTGG".to_vec()); // Met Trp — no synonymous sites
    let gene = cov_core_gene("metg", "chr1", Strand::Plus, 1, 6, 0);
    // pos 4: REF T -> ALT C  => TGG(Trp) -> CGG(Arg)  NONSYNONYMOUS
    let snps = vec![cov_core_snp("chr1", 4, b'T', b'C', 0.5)];
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.5);
    assert_eq!(res.len(), 1, "a translatable Met/Trp gene is retained");
    let g = &res[0];
    assert_eq!(g.s_sites, 0.0, "Met/Trp CDS has zero synonymous sites");
    assert!(g.n_sites > 0.0, "nonsynonymous sites still exist");
    assert!(g.ps.is_nan(), "pS must be NaN for zero synonymous sites, got {}", g.ps);
    assert!(g.pn.is_finite() && g.pn > 0.0, "pN stays defined, got {}", g.pn);
}

#[test]
fn cov_core_untranslatable_all_n_cds_is_skipped() {
    // A CDS that is entirely N in the reference has no translatable codons: pN and pS
    // would both be 0/0. The gene must be skipped, not emitted as an all-zero row that
    // reads as a real gene under total constraint.
    let gc = make_gc();
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"NNNNNNNNN".to_vec());
    let gene = cov_core_gene("nng", "chr1", Strand::Plus, 1, 9, 0);
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &[], gc, false, 1.0, 0.95);
    assert!(res.is_empty(), "all-N CDS must be skipped, got {} gene(s)", res.len());
}

#[test]
fn cov_core_pooled_zero_s_sites_gives_nan_ps_not_zero() {
    // The genome-wide pooled pS is undefined (NaN), not 0.0, when the retained set has
    // no synonymous sites at all — otherwise the summary prints a defined synonymous
    // rate of exactly 0 where none exists.
    let g = gene_result(30.0, 0.0, 3.0, 0.0); // n_sites>0, s_sites==0, some nonsyn SNPs
    let pooled = genome_wide_pn_ps(&[g]).expect("one gene pools");
    assert!(pooled.ps.is_nan(), "pooled pS must be NaN when S_sites == 0, got {}", pooled.ps);
    assert!(pooled.pn_ps.is_infinite(), "pN>0 over pS=0 pools to +inf, got {}", pooled.pn_ps);
}

#[test]
fn cov_core_compute_minus_strand_hand_checked() {
    let gc = make_gc();
    // Same coding sequence TTT GCT, but on the minus strand the reference holds
    // its reverse complement: reverse_complement("TTTGCT") == "AGCAAA".
    let mut reference = HashMap::new();
    reference.insert("chrM".to_string(), b"AGCAAA".to_vec());
    let gene = cov_core_gene("gm", "chrM", Strand::Minus, 1, 6, 0);
    // Minus-strand mapping over exon 1..6: cds_offset = 6 - pos.
    //   syn:    coding TTT->TTC at cds offset 2 => genomic pos 4.
    //           expected VCF REF = complement(cds[2]='T') = 'A';
    //           coding ALT 'C' = complement('G'), so VCF ALT = 'G'.
    //   nonsyn: coding GCT->ACT at cds offset 3 => genomic pos 3.
    //           expected VCF REF = complement(cds[3]='G') = 'C';
    //           coding ALT 'A' = complement('T'), so VCF ALT = 'T'.
    let snps = vec![
        cov_core_snp("chrM", 4, b'A', b'G', 0.5),
        cov_core_snp("chrM", 3, b'C', b'T', 0.5),
    ];
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);
    assert_eq!(res.len(), 1);
    let g = &res[0];
    assert_eq!(g.strand, '-');
    assert_eq!(g.nonsyn_snps, 1.0);
    assert_eq!(g.syn_snps, 1.0);
    // Same coding sequence => same Nei-Gojobori site counts and ratio as the plus case.
    assert!((g.n_sites - 14.0 / 3.0).abs() < 1e-9, "n_sites {}", g.n_sites);
    assert!((g.s_sites - 4.0 / 3.0).abs() < 1e-9, "s_sites {}", g.s_sites);
    assert!((g.pn_ps - 2.0 / 7.0).abs() < 1e-9, "pn_ps {}", g.pn_ps);
}

#[test]
fn cov_core_compute_af_weighted_uses_frequencies() {
    let gc = make_gc();
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"TTTGCT".to_vec());
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 6, 0);
    let snps = vec![
        cov_core_snp("chr1", 3, b'T', b'C', 0.8), // synonymous, weight 0.8
        cov_core_snp("chr1", 4, b'G', b'A', 0.4), // nonsynonymous, weight 0.4
    ];
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &snps, gc, true, 1.0, 0.5);
    assert_eq!(res.len(), 1);
    let g = &res[0];
    // Under --af-weighted each ALT contributes its allele frequency.
    assert!((g.syn_snps - 0.8).abs() < 1e-9, "syn {}", g.syn_snps);
    assert!((g.nonsyn_snps - 0.4).abs() < 1e-9, "nonsyn {}", g.nonsyn_snps);
    assert!((g.total_snps - 1.2).abs() < 1e-9, "total {}", g.total_snps);
    // pn = 0.4/(14/3), ps = 0.8/(4/3);
    // pn/ps = (0.4/(14/3)) / (0.8/(4/3)) = (0.4*4)/(0.8*14) = 1.6/11.2 = 1/7
    assert!((g.pn_ps - 1.0 / 7.0).abs() < 1e-9, "pn_ps {}", g.pn_ps);
    // Neutrality test and Wilson CI are undefined under AF weighting.
    assert!(g.p_value.is_nan());
    assert!(g.pn_ps_lo.is_nan() && g.pn_ps_hi.is_nan());
}

#[test]
fn cov_core_compute_skips_gene_with_missing_reference() {
    let gc = make_gc();
    // Reference has no entry for the gene's seqid -> gene is filtered out.
    let reference: HashMap<String, Vec<u8>> = HashMap::new();
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 6, 0);
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &[], gc, false, 1.0, 0.95);
    assert!(res.is_empty());
}

#[test]
fn cov_core_compute_skips_short_cds() {
    let gc = make_gc();
    // CDS shorter than one codon (2 bp) -> gene is skipped.
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"AT".to_vec());
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 2, 0);
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &[], gc, false, 1.0, 0.95);
    assert!(res.is_empty());
}

#[test]
fn cov_core_compute_ref_mismatch_snp_skipped() {
    let gc = make_gc();
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"TTTGCT".to_vec());
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 6, 0);
    // VCF claims REF=A at pos 3, but the reference has T there -> SNP is skipped.
    let snps = vec![cov_core_snp("chr1", 3, b'A', b'C', 0.5)];
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &snps, gc, false, 1.0, 0.95);
    assert_eq!(res.len(), 1);
    let g = &res[0];
    assert_eq!(g.total_snps, 0.0);
    assert_eq!(g.nonsyn_snps, 0.0);
    assert_eq!(g.syn_snps, 0.0);
}

#[test]
fn cov_core_bootstrap_determinism_and_bounds() {
    // Every gene has syn>0 and nonsyn>0, so every resample pools ps>0 and pn>0
    // => a finite, positive pooled pN/pS ratio.
    let results = vec![
        gene_result(100.0, 100.0, 10.0, 20.0),
        gene_result(80.0, 90.0, 12.0, 8.0),
        gene_result(120.0, 60.0, 5.0, 15.0),
        gene_result(40.0, 55.0, 9.0, 6.0),
    ];
    let a = bootstrap_genome_wide_ci(&results, 200, 12345, 0.95).expect("finite CI");
    let b = bootstrap_genome_wide_ci(&results, 200, 12345, 0.95).expect("finite CI");
    // (a) Determinism: same seed -> bit-identical bounds.
    assert_eq!(a.0, b.0);
    assert_eq!(a.1, b.1);
    // (b) Ordering and (c) finiteness.
    assert!(a.0 <= a.1, "lo {} > hi {}", a.0, a.1);
    assert!(a.0.is_finite() && a.1.is_finite());
    assert!(a.0 > 0.0, "pooled pN/pS is strictly positive here, got {}", a.0);
}

#[test]
fn cov_core_bootstrap_none_for_empty_and_zero_boot() {
    let results = vec![gene_result(100.0, 100.0, 10.0, 20.0)];
    // (d) None for empty input / n_boot == 0.
    assert!(bootstrap_genome_wide_ci(&[], 100, 1, 0.95).is_none());
    assert!(bootstrap_genome_wide_ci(&results, 0, 1, 0.95).is_none());
}

#[test]
fn cov_core_extract_cds_applies_phase() {
    // phase=1 on the first exon drops the single leading base of the assembled CDS.
    let gene = cov_core_gene("g", "chr1", Strand::Plus, 1, 6, 1);
    let ref_seq = b"ATTTGC".to_vec();
    let cds = extract_cds_sequence(&gene, &ref_seq);
    // Assembled "ATTTGC", skip phase=1 -> "TTTGC".
    assert_eq!(cds, b"TTTGC");
}

#[test]
fn cov_core_extract_cds_skips_out_of_bounds_exon() {
    // Second exon starts past the end of the reference -> skipped, not indexed.
    let ex = |start, end| crate::gff::CdsExon {
        seqid: "chr1".to_string(), start, end, strand: Strand::Plus, phase: 0,
    };
    let gene = Gene {
        name: "g".to_string(), seqid: "chr1".to_string(), strand: Strand::Plus,
        exons: vec![ex(1, 6), ex(100, 105)], length_bp: 12, start: 1,
    };
    let ref_seq = b"TTTGCT".to_vec(); // only 6 bp available
    let cds = extract_cds_sequence(&gene, &ref_seq);
    assert_eq!(cds, b"TTTGCT");
}

#[test]
fn cov_core_genomic_to_cds_offset_in_phase_region_is_none() {
    // phase=2: the first two coding positions lie within the phase region and
    // map to None; the third position maps to CDS offset 0.
    let gene = cov_core_gene("g", "chr1", Strand::Plus, 1, 9, 2);
    assert_eq!(genomic_to_cds_offset(&gene, 1), None);
    assert_eq!(genomic_to_cds_offset(&gene, 2), None);
    assert_eq!(genomic_to_cds_offset(&gene, 3), Some(0));
    assert_eq!(genomic_to_cds_offset(&gene, 4), Some(1));
}

#[test]
fn cov_core_genomic_to_cds_offset_multi_exon_plus() {
    // Exons 1..3 (CDS offsets 0..2) and 10..12 (CDS offsets 3..5), phase 0.
    let ex = |start, end| crate::gff::CdsExon {
        seqid: "chr1".to_string(), start, end, strand: Strand::Plus, phase: 0,
    };
    let gene = Gene {
        name: "g".to_string(), seqid: "chr1".to_string(), strand: Strand::Plus,
        exons: vec![ex(1, 3), ex(10, 12)], length_bp: 6, start: 1,
    };
    assert_eq!(genomic_to_cds_offset(&gene, 3), Some(2));
    assert_eq!(genomic_to_cds_offset(&gene, 10), Some(3));
    assert_eq!(genomic_to_cds_offset(&gene, 12), Some(5));
    // Intergenic gap between the two exons -> None.
    assert_eq!(genomic_to_cds_offset(&gene, 5), None);
}

#[test]
fn cov_core_count_sites_skips_stop_and_ambiguous_codons() {
    let gc = make_gc();
    let just_phe = count_sites(b"TTT", gc);
    // TAA is a stop codon and contributes no sites.
    let with_stop = count_sites(b"TTTTAA", gc);
    assert!((with_stop.0 - just_phe.0).abs() < 1e-9);
    assert!((with_stop.1 - just_phe.1).abs() < 1e-9);
    // A codon carrying an ambiguous base (N) is skipped entirely.
    let with_ambig = count_sites(b"TTTNGC", gc);
    assert!((with_ambig.0 - just_phe.0).abs() < 1e-9);
    assert!((with_ambig.1 - just_phe.1).abs() < 1e-9);
}

#[test]
fn cov_core_count_sites_excludes_changes_to_stop() {
    // The convention that determines agreement with KaKs_Calculator / SNPGenie:
    // a single-nt change that turns a *coding* codon into a stop is excluded from
    // the site counts (not counted as nonsynonymous, as classic Nei-Gojobori 1986
    // would). TGG (Trp) is the sharp case — of its 9 changes, 2nd→TAG and 3rd→TGA
    // are stop, so only 7 valid changes remain, all nonsynonymous:
    //   S = 3·0/7 = 0.0,  N = 3·7/7 = 3.0.
    // (Under the count-nonsense-as-nonsyn convention this would still be 0/3, but
    // 2-fold codons adjacent to a stop, e.g. TAC/TGC, would differ — see below.)
    let gc = make_gc();
    let (n, s) = count_sites(b"TGG", gc);
    assert!((n - 3.0).abs() < 1e-9, "TGG N sites: got {}, want 3.0", n);
    assert!((s - 0.0).abs() < 1e-9, "TGG S sites: got {}, want 0.0", s);

    // TGC (Cys): 3rd pos → TGA is stop (excluded); TGT is syn (Cys). Of the 8 valid
    // changes, 1 syn and 7 nonsyn ⇒ S = 3·1/8 = 0.375. If changes-to-stop were kept
    // (9 changes, 1 syn) S would be 3·1/9 = 0.333 — so this value confirms exclusion.
    let (_n2, s2) = count_sites(b"TGC", gc);
    assert!((s2 - 0.375).abs() < 1e-9, "TGC S sites: got {}, want 0.375 (stop excluded)", s2);
}

#[test]
fn cov_core_count_sites_weighted_skips_stop_and_ambiguous_codons() {
    let gc = make_gc();
    let just_phe = count_sites_weighted(b"TTT", gc, 2.0);
    // Whole stop / ambiguous codons contribute nothing (guard at the codon level).
    let with_stop = count_sites_weighted(b"TTTTAA", gc, 2.0);
    let with_ambig = count_sites_weighted(b"TTTNGC", gc, 2.0);
    assert!((with_stop.0 - just_phe.0).abs() < 1e-9 && (with_stop.1 - just_phe.1).abs() < 1e-9);
    assert!((with_ambig.0 - just_phe.0).abs() < 1e-9 && (with_ambig.1 - just_phe.1).abs() < 1e-9);
}

#[test]
fn cov_core_codon_to_aa_handles_ambiguous_and_rna() {
    let gc = make_gc();
    // Any ambiguous base collapses the whole codon to None.
    assert_eq!(codon_to_aa(b"ANG", gc), None);
    assert_eq!(codon_to_aa(b"NNN", gc), None);
    // Lowercase and RNA (U) bases still decode (U indexes like T).
    assert_eq!(codon_to_aa(b"atg", gc), Some(b'M'));
    assert_eq!(codon_to_aa(b"AUG", gc), Some(b'M'));
}

#[test]
fn cov_core_base_to_li_and_complement_edges() {
    // base_to_li: U/u map to index 3 (like T); non-bases -> None.
    assert_eq!(base_to_li(b'U'), Some(3));
    assert_eq!(base_to_li(b'u'), Some(3));
    assert_eq!(base_to_li(b'N'), None);
    assert_eq!(base_to_li(b'-'), None);
    // complement: U/u -> A; unknown bases collapse to N; lowercase complements too.
    assert_eq!(complement(b'U'), b'A');
    assert_eq!(complement(b'u'), b'A');
    assert_eq!(complement(b'N'), b'N');
    assert_eq!(complement(b'Z'), b'N');
    assert_eq!(complement(b'a'), b'T');
    assert_eq!(complement(b'g'), b'C');
}

#[test]
fn cov_core_sfs_bin_above_one_saturates_last_bin() {
    // An out-of-range AF (> 1.0) falls through the edge loop to the last bin.
    assert_eq!(sfs_bin(1.5), SFS_NBINS - 1);
}


#[test]
fn cov_bootstrap_keeps_diversifying_extreme_replicates() {
    // Genes with syn=0 but nonsyn>0 pool to pN/pS = +inf (diversifying extreme).
    // These replicates must be kept toward the upper CI, not discarded as undefined:
    // the CI must exist (not None) and its upper bound must reach +inf.
    let genes: Vec<GenePnPs> = (0..5).map(|_| gene_result(100.0, 50.0, 5.0, 0.0)).collect();
    let ci = bootstrap_genome_wide_ci(&genes, 100, 42, 0.95);
    assert!(ci.is_some(), "CI must exist for diversifying-extreme replicates, not None");
    assert!(ci.unwrap().1.is_infinite(), "upper CI bound must reach +inf");
}

#[test]
fn cov_core_gene_exon_past_reference_is_skipped() {
    // A CDS exon overrunning the reference used to desync the offset mapping and drop
    // SNPs silently; such a gene must now be skipped outright.
    let gc = make_gc();
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"TTTGCT".to_vec()); // 6 bp
    let gene = cov_core_gene("g", "chr1", Strand::Plus, 1, 12, 0); // exon 1..12 > 6 bp
    let (res, _diag) = compute_pn_ps(&reference, &[gene], &[], gc, false, 1.0, 0.95);
    assert!(res.is_empty(), "a gene whose CDS overruns the reference must be skipped");
}


#[test]
fn cov_core_diagnostics_report_in_cds_snps_and_internal_stops() {
    let gc = make_gc();
    // Clean CDS TTT GCT, one SNP inside it.
    let mut reference = HashMap::new();
    reference.insert("chr1".to_string(), b"TTTGCT".to_vec());
    let gene = cov_core_gene("g1", "chr1", Strand::Plus, 1, 6, 0);
    let (_r, diag) = compute_pn_ps(&reference, &[gene], &[cov_core_snp("chr1", 3, b'T', b'C', 0.8)], gc, false, 1.0, 0.5);
    assert_eq!(diag.snps_in_cds, 1);
    assert_eq!(diag.genes_with_internal_stops, 0);

    // CDS starting with a stop codon (TAA) -> flagged as an internal stop.
    let mut ref2 = HashMap::new();
    ref2.insert("chr1".to_string(), b"TAAGCT".to_vec());
    let gene2 = cov_core_gene("g2", "chr1", Strand::Plus, 1, 6, 0);
    let (_r2, diag2) = compute_pn_ps(&ref2, &[gene2], &[], gc, false, 1.0, 0.5);
    assert_eq!(diag2.genes_with_internal_stops, 1);
}
