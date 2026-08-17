//! How many times each allele arose, given a phylogeny over the cohort.
//!
//! This is the plumbing between three things that already exist: the per-sample carrier
//! sets ([`crate::vcf::CarrierSet`]), the tree ([`crate::tree`]), and the per-codon
//! recurrence scan. It answers one question per ALT allele, how many independent
//! origins does parsimony require, and hands the answer to the variants table and to
//! the codon scan, where it generalises the tested statistic from "distinct alleles" to
//! "independent origins" (see [`crate::vcf_analysis::codons`]).
//!
//! # Sample identity
//!
//! eskaks otherwise knows samples only as column positions. A tree knows them only as
//! names. The join is therefore the fragile part, and it is deliberately strict: see
//! [`crate::tree::SampleTree::join`], which refuses any disagreement rather than
//! dropping the unmatched side, because a silently dropped sample changes every origin
//! count without changing anything a reader can see.
//!
//! The naming rule is fixed and documented so a run is reproducible:
//!
//! * **one VCF**: the sample column names of the `#CHROM` header, in column order;
//! * **many VCFs** (`--vcf-list` or repeated `--vcf`): each file is one sample, named by
//!   its single inner sample column when it has exactly one, else by the file's
//!   basename with `.vcf` / `.vcf.gz` stripped.

use super::*;
use crate::tree::{OriginScratch, SampleTree};

/// An ALT allele within one contig: `(POS, ALT)`. An allele is identified by
/// `(CHROM, POS, ALT)` everywhere in eskaks; the contig is the outer key here so a
/// lookup can borrow the contig name instead of allocating a copy of it per variant row.
type AlleleKey = (usize, u8);

/// Independent-origin counts, per contig, keyed by [`AlleleKey`].
#[derive(Debug, Default)]
pub struct AlleleOrigins {
    counts: HashMap<String, HashMap<AlleleKey, u32>>,
    /// ALT alleles that had a carrier set, and so received a count.
    pub counted: usize,
    /// ALT alleles with no carrier set (a site where every genotype was a no-call, say).
    /// They stay `NA` rather than being scored as zero origins.
    pub uncounted: usize,
}

impl AlleleOrigins {
    /// Count independent origins for every ALT allele that carries per-sample identity.
    ///
    /// Cost per allele is its number of carriers plus one walk of the tree, so this is
    /// linear in (alleles x tree size) with a tiny constant and no pairwise anything.
    pub fn compute(snps: &[VcfSnp], tree: &SampleTree, min_support: u32) -> AlleleOrigins {
        // One scratch per rayon worker, reused across that worker's alleles: the
        // buffers are the size of the tree, so allocating per allele would cost more
        // than the parsimony passes themselves.
        let per_snp: Vec<(&str, Vec<(AlleleKey, u32)>)> = snps
            .par_iter()
            .map_init(
                || OriginScratch::for_tree(tree.tree()),
                |scratch, snp| {
                    let Some(carriers) = snp.carriers.as_ref() else {
                        return (snp.chrom.as_str(), Vec::new());
                    };
                    let row = snp
                        .alt_alleles
                        .iter()
                        .enumerate()
                        .filter_map(|(i, alt)| {
                            let cs = carriers.get(i)?;
                            let n = tree.origins(cs.samples(), min_support, scratch);
                            Some(((snp.pos, *alt), n))
                        })
                        .collect();
                    (snp.chrom.as_str(), row)
                },
            )
            .collect();

        let mut counts: HashMap<String, HashMap<AlleleKey, u32>> = HashMap::new();
        for (chrom, row) in per_snp {
            if row.is_empty() {
                continue;
            }
            let per_chrom = counts.entry(chrom.to_string()).or_default();
            for (key, n) in row {
                // A duplicated (CHROM, POS, ALT) record describes ONE allele; its
                // carriers are the same set, so keeping the larger count is a no-op in
                // practice and is the safe rule if two records ever disagree.
                let slot = per_chrom.entry(key).or_insert(0);
                *slot = (*slot).max(n);
            }
        }
        let counted: usize = counts.values().map(HashMap::len).sum();
        let total_alts: usize = snps.iter().map(|s| s.alt_alleles.len()).sum();
        AlleleOrigins { counts, counted, uncounted: total_alts.saturating_sub(counted) }
    }

    /// Origins of one ALT allele, or `None` when it had no per-sample carriers.
    ///
    /// The single-allele lookup this module's own tests assert on. `attach` below walks
    /// the per-contig map directly, so without the attribute this reads as dead code in
    /// the binary target, which compiles the module a second time.
    #[allow(dead_code)]
    pub fn get(&self, chrom: &str, pos: usize, alt: u8) -> Option<u32> {
        self.counts.get(chrom)?.get(&(pos, alt)).copied()
    }

    /// Write each variant's origin count onto it, keyed back to the SNP records.
    ///
    /// Keying rather than plumbing a `CarrierSet` into every [`Variant`] is deliberate:
    /// a codon in two overlapping genes produces two rows for one allele, and a bitset
    /// per row would duplicate the cohort per row for nothing.
    pub fn attach(&self, results: &mut [GenePnPs]) {
        for g in results.iter_mut() {
            // One contig lookup per gene, then one (POS, ALT) lookup per row.
            let Some(per_chrom) = self.counts.get(g.chrom.as_str()) else {
                for v in g.variants.iter_mut() {
                    v.origins = None;
                }
                continue;
            };
            for v in g.variants.iter_mut() {
                v.origins = per_chrom.get(&(v.pos, v.alt_allele)).copied();
            }
        }
    }
}

/// The cohort's sample names, in sample-index order, under the documented rule (see the
/// module docs). `vcf_paths` is the run's VCF list in the order the merge used, so the
/// position in the returned vector IS the sample index the carrier sets were built with.
pub fn cohort_sample_names(vcf_paths: &[String]) -> anyhow::Result<Vec<String>> {
    if vcf_paths.len() == 1 {
        let names = crate::vcf::sample_names(std::path::Path::new(&vcf_paths[0]))?;
        if names.is_empty() {
            anyhow::bail!(
                "{}: the VCF has no per-sample genotype columns, so there are no samples to \
                 match a tree's tips to. --tree needs genotypes.",
                vcf_paths[0]
            );
        }
        return Ok(names);
    }
    let mut out = Vec::with_capacity(vcf_paths.len());
    for p in vcf_paths {
        let path = std::path::Path::new(p);
        let inner = crate::vcf::sample_names(path)?;
        out.push(if inner.len() == 1 {
            inner.into_iter().next().expect("length checked")
        } else {
            basename_sample(p)
        });
    }
    Ok(out)
}

/// A per-sample VCF's name when its header does not supply exactly one: the file's
/// basename with a `.gz`/`.bgz` wrapper and then a `.vcf`/`.bcf` extension removed, so
/// `data/ERR1234.vcf.gz` is the sample `ERR1234`.
fn basename_sample(path: &str) -> String {
    let base = std::path::Path::new(path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string());
    let base = base
        .strip_suffix(".gz")
        .or_else(|| base.strip_suffix(".bgz"))
        .unwrap_or(&base)
        .to_string();
    base.strip_suffix(".vcf")
        .or_else(|| base.strip_suffix(".bcf"))
        .unwrap_or(&base)
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::Tree;
    use crate::vcf::CarrierSet;

    fn snp(pos: usize, alt: u8, carriers: &[usize], n: usize) -> VcfSnp {
        let mut cs = CarrierSet::new(n);
        for &s in carriers {
            cs.insert(s);
        }
        VcfSnp {
            chrom: "chr1".to_string(),
            pos,
            ref_allele: b'A',
            alt_alleles: vec![alt],
            alt_freqs: vec![carriers.len() as f64 / n as f64],
            gt_counts: None,
            carriers: Some(vec![cs]),
            filter: "PASS".to_string(),
            depth: None,
        }
    }

    fn eight_tip_tree() -> SampleTree {
        let names: Vec<String> = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let tree = Tree::parse_newick("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)));").expect("newick");
        SampleTree::join(tree, &names).expect("join")
    }

    #[test]
    fn a_clonal_allele_and_a_convergent_one_are_told_apart() {
        let tree = eight_tip_tree();
        let snps = vec![
            // One clade of four: one origin, four carriers.
            snp(10, b'T', &[0, 1, 2, 3], 8),
            // The same four carriers, spread over four clades: four origins.
            snp(20, b'G', &[0, 2, 4, 6], 8),
        ];
        let o = AlleleOrigins::compute(&snps, &tree, 1);
        assert_eq!(o.get("chr1", 10, b'T'), Some(1));
        assert_eq!(o.get("chr1", 20, b'G'), Some(4));
        assert_eq!((o.counted, o.uncounted), (2, 0));
        // Carrier COUNT cannot tell these apart; that is the whole point.
        assert_eq!(snps[0].carriers.as_ref().unwrap()[0].len(), 4);
        assert_eq!(snps[1].carriers.as_ref().unwrap()[0].len(), 4);
    }

    #[test]
    fn an_allele_without_carriers_is_na_and_not_zero() {
        let tree = eight_tip_tree();
        let mut s = snp(10, b'T', &[0, 1], 8);
        s.carriers = None; // an AF-only record, or a site of pure no-calls
        let o = AlleleOrigins::compute(&[s], &tree, 1);
        assert_eq!(o.get("chr1", 10, b'T'), None, "unknown must not read as zero origins");
        assert_eq!((o.counted, o.uncounted), (0, 1));
    }

    #[test]
    fn origins_land_on_the_right_variant_rows() {
        // The attach step keys on (CHROM, POS, ALT). A gene on another contig with the
        // same coordinates must not collect the wrong count.
        let tree = eight_tip_tree();
        let snps = vec![snp(10, b'T', &[0, 1, 2, 3], 8), snp(10, b'G', &[0, 2, 4, 6], 8)];
        let o = AlleleOrigins::compute(&snps, &tree, 1);
        assert_eq!(o.get("chr1", 10, b'T'), Some(1));
        assert_eq!(o.get("chr1", 10, b'G'), Some(4), "same position, different ALT");
        assert_eq!(o.get("chr2", 10, b'T'), None, "same position, different contig");
        assert_eq!(o.get("chr1", 11, b'T'), None);
    }

    #[test]
    fn the_naming_rule_for_a_multi_vcf_run_is_the_documented_one() {
        // Basename, with the wrappers stripped, is what a per-sample VCF with no usable
        // header name falls back to.
        assert_eq!(basename_sample("/data/runs/ERR1234.vcf.gz"), "ERR1234");
        assert_eq!(basename_sample("ERR1234.vcf"), "ERR1234");
        assert_eq!(basename_sample("./x/ERR1234.bcf"), "ERR1234");
        assert_eq!(basename_sample("ERR1234"), "ERR1234");
        // A name with dots in it keeps them: only the recognised wrappers come off.
        assert_eq!(basename_sample("/d/S1.filtered.vcf.gz"), "S1.filtered");
    }
}
