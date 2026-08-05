//! NCBI genetic code tables for alternative translation tables.
//!
//! Supports 20 NCBI translation tables covering all known genetic codes.
//! Default is table 1 (Standard). Common alternatives include:
//! - Table 2: Vertebrate Mitochondrial
//! - Table 4: Mycoplasma/Spiroplasma
//! - Table 11: Bacterial, Archaeal and Plant Plastid
//!
//! Tables generated programmatically from the NCBI codon table definitions:
//! <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
//!
//! All amino acid strings use Li indexing: `16*b1 + 4*b2 + b3` where A=0, C=1, G=2, T=3.
//! Stop codons are encoded as `*`.

/// A genetic code translation table.
#[derive(Clone)]
pub struct GeneticCode {
    /// NCBI table number.
    pub id: u8,
    /// Human-readable name.
    pub name: &'static str,
    /// 64-entry amino acid array indexed by Li codon index.
    /// Stop codons are `b'*'`.
    pub aa_table: [u8; 64],
}

/// All supported NCBI genetic code tables.
///
/// Codon order (Li indexing, A=0 C=1 G=2 T=3):
/// ```text
/// AAA AAC AAG AAT  ACA ACC ACG ACT  AGA AGC AGG AGT  ATA ATC ATG ATT
/// CAA CAC CAG CAT  CCA CCC CCG CCT  CGA CGC CGG CGT  CTA CTC CTG CTT
/// GAA GAC GAG GAT  GCA GCC GCG GCT  GGA GGC GGG GGT  GTA GTC GTG GTT
/// TAA TAC TAG TAT  TCA TCC TCG TCT  TGA TGC TGG TGT  TTA TTC TTG TTT
/// ```
static TABLES: &[GeneticCode] = &[
    GeneticCode {
        id: 1,
        name: "Standard",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 2,
        name: "Vertebrate Mitochondrial",
        aa_table: *b"KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 3,
        name: "Yeast Mitochondrial",
        aa_table: *b"KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 4,
        name: "Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma/Spiroplasma",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 5,
        name: "Invertebrate Mitochondrial",
        aa_table: *b"KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 6,
        name: "Ciliate, Dasycladacean and Hexamita Nuclear",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 9,
        name: "Echinoderm and Flatworm Mitochondrial",
        aa_table: *b"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 10,
        name: "Euplotid Nuclear",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
    },
    GeneticCode {
        id: 11,
        name: "Bacterial, Archaeal and Plant Plastid",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 12,
        name: "Alternative Yeast Nuclear",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 13,
        name: "Ascidian Mitochondrial",
        aa_table: *b"KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        aa_table: *b"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 21,
        name: "Trematode Mitochondrial",
        aa_table: *b"NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF",
    },
    GeneticCode {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF",
    },
    GeneticCode {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        aa_table: *b"KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
    },
    GeneticCode {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF",
    },
    GeneticCode {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        aa_table: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLALEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
    },
    GeneticCode {
        id: 33,
        name: "Cephalodiscidae Mitochondrial UAA-Tyr",
        aa_table: *b"KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
    },
];

/// Look up a genetic code table by NCBI ID.
pub fn get_table(id: u8) -> Option<&'static GeneticCode> {
    TABLES.iter().find(|t| t.id == id)
}

/// List all available genetic code table IDs and names.
pub fn list_tables() -> Vec<(u8, &'static str)> {
    TABLES.iter().map(|t| (t.id, t.name)).collect()
}

/// Get all codon indices that encode stop codons for a given genetic code and model.
/// Returns indices in the model's native indexing (Li or Nei).
pub fn stop_codon_indices(gc: &GeneticCode, model: crate::models::Model) -> Vec<u8> {
    match model {
        crate::models::Model::Li => {
            (0u8..64)
                .filter(|&i| gc.aa_table[i as usize] == b'*')
                .collect()
        }
        crate::models::Model::Nei => {
            let nei_aa = li_to_nei_aa(&gc.aa_table);
            (0u8..64)
                .filter(|&i| nei_aa[i as usize] == '*')
                .collect()
        }
    }
}

/// Convert a Li-indexed amino acid table to Nei indexing.
///
/// Li index: `16*b1 + 4*b2 + b3` where A=0, C=1, G=2, T=3.
/// Nei index: `16*remap(b2) + 4*remap(b1) + remap(b3)` where remap: A→2, C→1, G→3, T→0.
pub fn li_to_nei_aa(li_table: &[u8; 64]) -> [char; 64] {
    let mut nei_table = ['?'; 64];
    let remap = |b: usize| -> usize {
        match b { 0 => 2, 1 => 1, 2 => 3, 3 => 0, _ => unreachable!() }
    };

    for (li_idx, &aa) in li_table.iter().enumerate() {
        let b1 = li_idx / 16;
        let b2 = (li_idx / 4) % 4;
        let b3 = li_idx % 4;
        let nei_idx = 16 * remap(b2) + 4 * remap(b1) + remap(b3);
        nei_table[nei_idx] = aa as char;
    }
    nei_table
}

/// Per-codon synonymous SITES under the Nei-Gojobori "exclude changes-to-stop"
/// convention — the same rule this crate's VCF pN/pS counter
/// ([`crate::vcf_analysis::sites::count_sites`]) and the reference tools
/// (KaKs_Calculator, SNPGenie) use, so both eskaks paths share one convention.
///
/// For each codon, of its nine single-nucleotide changes we count how many are
/// synonymous vs nonsynonymous, EXCLUDING any change that creates a stop codon,
/// then split the codon's three sites in that ratio: `S = 3·syn/(syn+nonsyn)`
/// (and the caller takes `N = 3 − S`). Classic Nei-Gojobori 1986 instead keeps a
/// change-to-stop as nonsynonymous (fixed denominator 9); excluding it makes S
/// slightly larger at stop-adjacent codons. Stop codons themselves return 0.
pub fn compute_syn_sites(nei_aa: &[char; 64]) -> [f64; 64] {
    let mut syn_sites = [0.0f64; 64];

    for idx in 0..64usize {
        let aa = nei_aa[idx];
        if aa == '*' {
            continue; // stop codons are never summed; keep 0 to mirror count_sites
        }
        let slots = [idx / 16, (idx / 4) % 4, idx % 4];
        let (mut syn, mut nonsyn) = (0usize, 0usize);
        for pos in 0..3 {
            for alt in 0..4usize {
                if alt == slots[pos] {
                    continue;
                }
                let mut new_slots = slots;
                new_slots[pos] = alt;
                let new_idx = new_slots[0] * 16 + new_slots[1] * 4 + new_slots[2];
                match nei_aa[new_idx] {
                    '*' => {}                 // exclude changes to a stop codon
                    a if a == aa => syn += 1, // synonymous
                    _ => nonsyn += 1,         // nonsynonymous (non-stop)
                }
            }
        }
        let valid = (syn + nonsyn) as f64;
        if valid > 0.0 {
            syn_sites[idx] = 3.0 * syn as f64 / valid;
        }
    }
    syn_sites
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_table_exists() {
        let t = get_table(1).expect("Table 1 should exist");
        assert_eq!(t.name, "Standard");
    }

    #[test]
    fn all_tables_have_64_entries() {
        for t in TABLES {
            assert_eq!(t.aa_table.len(), 64, "Table {} has wrong length", t.id);
        }
    }

    #[test]
    fn standard_table_matches_hardcoded_nei_aa() {
        let t = get_table(1).unwrap();
        let nei_aa = li_to_nei_aa(&t.aa_table);
        // This is the hardcoded AA_ARRAY in nei.rs
        let expected: [char; 64] = [
            'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
            'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
            'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
            'C', 'C', '*', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G',
        ];
        assert_eq!(nei_aa, expected, "Standard code Nei conversion mismatch");
    }

    #[test]
    fn compute_syn_sites_matches_hand_derived_exclude_nonsense() {
        let t = get_table(1).unwrap();
        let syn = compute_syn_sites(&li_to_nei_aa(&t.aa_table));
        // Nei index of a codon (see the module header): 16·remap(b2)+4·remap(b1)+remap(b3),
        // remap A→2 C→1 G→3 T→0. This lets us look up specific codons by name.
        let remap = |b: u8| match b { b'A' => 2, b'C' => 1, b'G' => 3, b'T' => 0, _ => unreachable!() };
        let ix = |c: &[u8; 3]| 16 * remap(c[1]) + 4 * remap(c[0]) + remap(c[2]);
        let s = |c: &[u8; 3]| syn[ix(c)];

        // Hand-derived synonymous SITES = 3·syn/(syn+nonsyn), stop changes excluded:
        //   TTT (Phe): 1 syn (TTC), 8 nonsyn, 0 stop           → 3·1/9 = 1/3
        //   GCT (Ala): 3 syn (4-fold 3rd pos), 6 nonsyn, 0 stop → 3·3/9 = 1.0
        //   ATG (Met): unique codon, 0 syn                       → 0.0
        //   TGG (Trp): 2nd→TAG, 3rd→TGA stop (excluded); 7 nonsyn → 0.0
        //   TAT (Tyr): 3rd→TAA & →TAG stop (excluded); 1 syn (TAC), 6 nonsyn → 3·1/7 = 3/7
        //   TGC (Cys): 3rd→TGA stop (excluded); 1 syn (TGT), 7 nonsyn        → 3·1/8 = 0.375
        // The last two are the whole point: under classic Nei-Gojobori (keep nonsense)
        // TAT would be 3·1/9 = 1/3 and TGC 3·1/9 = 1/3, so these values prove exclusion.
        assert!((s(b"TTT") - 1.0 / 3.0).abs() < 1e-12, "TTT {}", s(b"TTT"));
        assert!((s(b"GCT") - 1.0).abs() < 1e-12, "GCT {}", s(b"GCT"));
        assert!((s(b"ATG") - 0.0).abs() < 1e-12, "ATG {}", s(b"ATG"));
        assert!((s(b"TGG") - 0.0).abs() < 1e-12, "TGG {}", s(b"TGG"));
        assert!((s(b"TAT") - 3.0 / 7.0).abs() < 1e-12, "TAT {}", s(b"TAT"));
        assert!((s(b"TGC") - 0.375).abs() < 1e-12, "TGC {}", s(b"TGC"));
        // Every non-stop codon splits into exactly 3 sites (S + N = 3); a stop returns 0.
        let nei = li_to_nei_aa(&t.aa_table);
        for (idx, &site) in syn.iter().enumerate() {
            assert!((0.0..=3.0).contains(&site), "codon {} S out of range: {}", idx, site);
            if nei[idx] == '*' {
                assert_eq!(site, 0.0, "stop codon {} must contribute 0 sites", idx);
            }
        }
    }

    #[test]
    fn compute_syn_sites_agrees_with_vcf_count_sites() {
        // Unification guard: the fasta (Nei) path's precomputed per-codon S-site table
        // must equal, codon for codon, the vcf path's `count_sites` — one convention,
        // two code paths. Both are Nei-Gojobori with changes-to-stop excluded.
        let t = get_table(1).unwrap();
        let syn = compute_syn_sites(&li_to_nei_aa(&t.aa_table));
        // Nei index → bases: inverse of remap (A→2 C→1 G→3 T→0) ⇒ 0→T 1→C 2→A 3→G.
        let inv = b"TCAG";
        for (idx, &s_fasta) in syn.iter().enumerate() {
            let (s0, s1, s2) = (idx / 16, (idx / 4) % 4, idx % 4);
            let codon = [inv[s1], inv[s0], inv[s2]]; // (b1, b2, b3)
            let (_n, s_cs) = crate::vcf_analysis::count_sites(&codon, t);
            assert!(
                (s_fasta - s_cs).abs() < 1e-12,
                "codon {} (nei idx {}): fasta S={} vs count_sites S={}",
                std::str::from_utf8(&codon).unwrap(),
                idx,
                s_fasta,
                s_cs
            );
        }
    }

    #[test]
    fn vertebrate_mito_differs_from_standard() {
        let std_t = get_table(1).unwrap();
        let mito = get_table(2).unwrap();
        assert_ne!(std_t.aa_table, mito.aa_table);
        // Vertebrate mito: AGA=* (stop), AGG=* (stop), ATA=M, TGA=W
        // Li indices: AGA=8, AGG=10, ATA=12, TGA=56
        assert_eq!(mito.aa_table[8], b'*');  // AGA = stop
        assert_eq!(std_t.aa_table[8], b'R'); // AGA = Arg
        assert_eq!(mito.aa_table[12], b'M'); // ATA = Met
        assert_eq!(std_t.aa_table[12], b'I'); // ATA = Ile
    }

    #[test]
    fn bacterial_table_same_as_standard() {
        let std_t = get_table(1).unwrap();
        let bac = get_table(11).unwrap();
        assert_eq!(std_t.aa_table, bac.aa_table);
    }

    #[test]
    fn list_tables_returns_all() {
        let tables = list_tables();
        assert_eq!(tables.len(), 20, "Should have 20 tables");
        assert!(tables.iter().any(|(id, _)| *id == 1));
        assert!(tables.iter().any(|(id, _)| *id == 2));
        assert!(tables.iter().any(|(id, _)| *id == 11));
        assert!(tables.iter().any(|(id, _)| *id == 33));
    }

    #[test]
    fn all_tables_have_valid_amino_acids() {
        let valid = b"ACDEFGHIKLMNPQRSTVWY*";
        for t in TABLES {
            for (i, &aa) in t.aa_table.iter().enumerate() {
                assert!(valid.contains(&aa),
                    "Table {} index {}: invalid amino acid '{}'", t.id, i, aa as char);
            }
        }
    }

    #[test]
    fn mito_table_triggers_nei_fallback() {
        // Vertebrate mitochondrial: AGA and AGG are stop codons.
        // This creates pairs where ALL pathways go through stop intermediates,
        // unlike the standard genetic code where fallbacks are unreachable.
        let mito = get_table(2).unwrap();
        let nei_aa = li_to_nei_aa(&mito.aa_table);

        // Check that some 2-diff pairs now have all-stop intermediates
        let mut fallback_count = 0;
        for idx_a in 0u8..64 {
            for idx_b in (idx_a + 1)..64 {
                if nei_aa[idx_a as usize] == '*' || nei_aa[idx_b as usize] == '*' { continue; }
                let slots_a = [idx_a / 16, (idx_a / 4) % 4, idx_a % 4];
                let slots_b = [idx_b / 16, (idx_b / 4) % 4, idx_b % 4];
                let mut diffs = Vec::new();
                for pos in 0..3 {
                    if slots_a[pos as usize] != slots_b[pos as usize] { diffs.push(pos); }
                }
                if diffs.len() == 2 {
                    let mut all_stop = true;
                    for &first in &diffs {
                        let mut inter = slots_a;
                        inter[first as usize] = slots_b[first as usize];
                        let inter_idx = inter[0] as usize * 16 + inter[1] as usize * 4 + inter[2] as usize;
                        if nei_aa[inter_idx] != '*' { all_stop = false; break; }
                    }
                    if all_stop { fallback_count += 1; }
                }
            }
        }
        assert!(fallback_count > 0,
            "Vertebrate mito should have at least one 2-diff pair triggering fallback, got 0");
    }

    #[test]
    fn tables_24_and_33_aga_is_ser() {
        // Regression: NCBI transl_table 24 and 33 define AGA=Ser, AGG=Lys. The Li
        // index for a codon is 16*b1+4*b2+b3 (A=0,C=1,G=2,T=3), so AGA=8 and AGG=10.
        for id in [24u8, 33] {
            let t = get_table(id).unwrap();
            assert_eq!(t.aa_table[8], b'S', "code {id}: AGA must be Ser");
            assert_eq!(t.aa_table[10], b'K', "code {id}: AGG must be Lys");
            assert_eq!(t.aa_table[9], b'S', "code {id}: AGC must be Ser");
            assert_eq!(t.aa_table[11], b'S', "code {id}: AGT must be Ser");
        }
    }

    #[test]
    fn cov_get_table_unknown_id_is_none() {
        // This crate defines no NCBI tables 0, 7, 8, or 99 => lookup returns None.
        assert!(get_table(0).is_none());
        assert!(get_table(7).is_none());
        assert!(get_table(8).is_none());
        assert!(get_table(99).is_none());
    }

    #[test]
    fn cov_stop_codon_indices_standard_li_and_nei() {
        // Standard code (table 1) stops are TAA, TAG, TGA.
        // Li indexing 16*b1+4*b2+b3 (A=0,C=1,G=2,T=3):
        //   TAA = 16*3+4*0+0 = 48, TAG = 16*3+4*0+2 = 50, TGA = 16*3+4*2+0 = 56.
        // stop_codon_indices iterates 0..64 ascending, so the order is fixed.
        let gc = get_table(1).unwrap();
        assert_eq!(stop_codon_indices(gc, crate::models::Model::Li), vec![48u8, 50, 56]);
        // Nei indexing 16*remap(b2)+4*remap(b1)+remap(b3), remap A->2 C->1 G->3 T->0:
        //   TAA -> 16*2+4*0+2 = 34, TAG -> 16*2+4*0+3 = 35, TGA -> 16*3+4*0+2 = 50.
        assert_eq!(stop_codon_indices(gc, crate::models::Model::Nei), vec![34u8, 35, 50]);
    }

    #[test]
    fn cov_stop_codon_indices_count_matches_star_entries() {
        // Both model indexings must return exactly one index per '*' entry of the
        // Li-indexed aa_table, and every Li stop index must map back to b'*'.
        for &id in &[1u8, 2, 6, 11] {
            let gc = get_table(id).unwrap();
            let star = gc.aa_table.iter().filter(|&&a| a == b'*').count();
            let li = stop_codon_indices(gc, crate::models::Model::Li);
            let nei = stop_codon_indices(gc, crate::models::Model::Nei);
            assert_eq!(li.len(), star, "code {id}: Li stop count");
            assert_eq!(nei.len(), star, "code {id}: Nei stop count");
            for i in &li {
                assert_eq!(gc.aa_table[*i as usize], b'*', "code {id}: Li idx {i}");
            }
        }
    }
}
