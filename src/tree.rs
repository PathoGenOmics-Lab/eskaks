//! A phylogeny over the cohort, and how many times an allele arose on it.
//!
//! # Why a tree at all
//!
//! The per-codon recurrence scan ([`crate::vcf_analysis::compute_codon_scan`]) counts
//! **distinct alleles**, because each distinct allele is at minimum one independent
//! mutational event and that is the only recurrence claim a phylogeny-free tool can
//! make. It therefore sees *gyrA* D94 (six distinct alleles at one residue) and misses
//! *rpoB* S450L, whose whole evidence is that ONE allele arose many times over. One
//! allele is one allele whether it arose once and spread to two thousand isolates or
//! arose fifty times independently; nothing in a VCF separates those two.
//!
//! A tree separates them. Given which samples carry an allele, Fitch parsimony returns
//! the minimum number of times the allele must have arisen on that tree, and summing
//! that over a codon's missense alleles generalises the existing statistic exactly:
//! when every allele arose once, the sum **is** the count of distinct alleles.
//!
//! # What is counted, and the one knob that makes it usable
//!
//! Raw parsimony is brutally fragile here. A genuine single-origin clade allele plus
//! three random false calls scores four gains, which at genome scale is a significant
//! hit built out of three sequencing errors. Counting only origins whose clade carries
//! the allele in at least `min_support` samples removes that failure mode outright
//! (singleton false calls can never form a supported clade) and, in simulation, also
//! demotes a pure calling artefact below the real signals it otherwise outranks. It
//! also absorbs the origin inflation an unresolved (polytomous) tree produces.
//!
//! The filter biases **conservative**: a mutation that genuinely arose once in one
//! sampled genome goes uncounted. That is the right direction for a genome-wide scan.
//!
//! # Rooting
//!
//! The root is taken as a non-carrier, which is what makes "how many times did this
//! arise" well defined without an outgroup: gains are counted downward from a
//! non-carrier ancestor. Where the root's own parsimony state is nevertheless
//! carrier-only (an allele in every sampled genome), that counts as one origin above
//! the root rather than none.

use std::collections::HashMap;

/// One node of a parsed tree. Tips are the nodes with no children.
#[derive(Debug, Clone)]
struct Node {
    /// Index of the parent, or `usize::MAX` for the root.
    parent: usize,
    children: Vec<usize>,
    /// The label as written, unquoted. `None` for an unlabelled node.
    label: Option<String>,
}

/// The sentinel [`Node::parent`] of the root.
const NO_PARENT: usize = usize::MAX;

/// A rooted phylogeny over named tips, with the traversal orders the parsimony passes
/// need precomputed so counting an allele's origins costs one walk of the nodes.
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    root: usize,
    /// Node indices of the tips, in the order they appear in the Newick string.
    tips: Vec<usize>,
    /// Tip slot of each node (its index in `tips`), or `u32::MAX` for an internal node.
    tip_slot: Vec<u32>,
    /// Every node, children before parents. Reversed, it is a valid parent-first order.
    postorder: Vec<usize>,
}

/// Reusable per-thread scratch for [`Tree::count_origins`], so a genome-scale run does
/// not allocate three vectors per allele.
#[derive(Debug, Clone)]
pub struct OriginScratch {
    /// Fitch state set per node: bit 0 = non-carrier, bit 1 = carrier.
    state: Vec<u8>,
    /// Carrier tips in each node's subtree.
    subtree: Vec<u32>,
    /// The state finally assigned to each node by the parent-first pass.
    assigned: Vec<u8>,
    /// Carrier flag per tip slot, cleared after every allele.
    tip_state: Vec<bool>,
    /// The slots set for the current allele, so clearing costs the carriers rather
    /// than the cohort.
    touched: Vec<usize>,
}

impl OriginScratch {
    /// Scratch sized for `tree`.
    pub fn for_tree(tree: &Tree) -> Self {
        OriginScratch {
            state: vec![0; tree.nodes.len()],
            subtree: vec![0; tree.nodes.len()],
            assigned: vec![0; tree.nodes.len()],
            tip_state: vec![false; tree.tips.len()],
            touched: Vec::new(),
        }
    }
}

impl Tree {
    /// Number of tips.
    pub fn n_tips(&self) -> usize {
        self.tips.len()
    }

    /// Number of nodes, tips included.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Every tip's label, by slot. `None` for an unlabelled tip, which is an error for
    /// any caller that has to join tips to samples.
    pub fn tip_labels(&self) -> impl Iterator<Item = Option<&str>> + '_ {
        self.tips.iter().map(|&n| self.nodes[n].label.as_deref())
    }

    /// How many times the allele marked in `scratch` must have arisen on this tree.
    ///
    /// Fitch parsimony over the two-state carrier character, with the root taken as a
    /// non-carrier, counting the edges that go non-carrier → carrier. An origin is
    /// counted only when the clade below it carries the allele in at least
    /// `min_support` samples, which is what makes the count robust to isolated false
    /// calls and to unresolved nodes (see the module docs).
    ///
    /// The carrier marks are read, never cleared: the caller owns them (and clears the
    /// few it set) so the cost is the number of carriers rather than the cohort size.
    pub fn count_origins(&self, min_support: u32, scratch: &mut OriginScratch) -> u32 {
        let (state, subtree, assigned) = (
            &mut scratch.state,
            &mut scratch.subtree,
            &mut scratch.assigned,
        );
        // Pass 1, children before parents: the Fitch state set and the number of
        // carrier tips below each node. The multifurcating rule (keep the states held
        // by the most children) is what makes a polytomy cost one change rather than
        // one per branch, and it reduces to the usual intersection-else-union rule at
        // a bifurcation.
        for &v in &self.postorder {
            let node = &self.nodes[v];
            if node.children.is_empty() {
                let carries = scratch.tip_state[self.tip_slot[v] as usize];
                state[v] = if carries { 0b10 } else { 0b01 };
                subtree[v] = u32::from(carries);
                continue;
            }
            let (mut absent, mut present, mut below) = (0u32, 0u32, 0u32);
            for &c in &node.children {
                absent += u32::from(state[c] & 0b01 != 0);
                present += u32::from(state[c] & 0b10 != 0);
                below += subtree[c];
            }
            let best = absent.max(present);
            state[v] = u8::from(absent == best) | (u8::from(present == best) << 1);
            subtree[v] = below;
        }

        // Pass 2, parents before children: resolve each set to one state. The root
        // prefers non-carrier wherever its set allows, which is what "the root is a
        // non-carrier" means when the data leave the choice open; every other node
        // inherits its parent's state whenever that is in its own set, so a change is
        // placed only where parsimony forces one.
        let mut origins = 0u32;
        for &v in self.postorder.iter().rev() {
            let node = &self.nodes[v];
            if v == self.root {
                assigned[v] = u8::from(state[v] & 0b01 == 0);
                // An allele whose parsimony state is "carrier" at the root arose above
                // the root: one origin, not none.
                if assigned[v] == 1 && subtree[v] >= min_support {
                    origins += 1;
                }
                continue;
            }
            let up = assigned[node.parent];
            assigned[v] = if state[v] & (1 << up) != 0 {
                up
            } else {
                u8::from(state[v] & 0b01 == 0)
            };
            if up == 0 && assigned[v] == 1 && subtree[v] >= min_support {
                origins += 1;
            }
        }
        origins
    }

    /// Parse a Newick string into a rooted tree.
    ///
    /// Hand-rolled and iterative on purpose: it adds no dependency, and a recursive
    /// descent would overflow the stack on the ladder-shaped trees a clonal cohort
    /// produces. Branch lengths and internal labels (including bootstrap support) are
    /// parsed and discarded (parsimony uses neither), but a malformed length is an
    /// error rather than something to skip past, since it means the file is not what
    /// the caller thinks it is.
    ///
    /// Square-bracket comments (including NHX annotations) are skipped anywhere.
    /// Labels may be single- or double-quoted, with a doubled quote as the escape.
    /// Underscores in unquoted labels are kept **literal**: the Newick convention of
    /// reading them as spaces would silently rename samples such as `ERR_1234`.
    pub fn parse_newick(text: &str) -> anyhow::Result<Tree> {
        let mut p = Parser { src: text.as_bytes(), i: 0 };
        let mut nodes: Vec<Node> = Vec::new();
        let mut stack: Vec<usize> = Vec::new();
        // The node completed at the current level: what a label, a branch length or a
        // ')' should attach to. `None` means the next token opens a fresh node.
        let mut cur: Option<usize> = None;
        let mut finished = false;

        /// Create a node, attached to whichever group is currently open.
        fn push(nodes: &mut Vec<Node>, stack: &[usize]) -> usize {
            let idx = nodes.len();
            let parent = stack.last().copied().unwrap_or(NO_PARENT);
            nodes.push(Node { parent, children: Vec::new(), label: None });
            if parent != NO_PARENT {
                nodes[parent].children.push(idx);
            }
            idx
        }

        while let Some(tok) = p.next_token()? {
            match tok {
                Token::Open => {
                    if cur.is_some() {
                        anyhow::bail!("Newick: '(' directly after a completed node at byte {}", p.i);
                    }
                    let n = push(&mut nodes, &stack);
                    stack.push(n);
                }
                Token::Label(name) => match cur {
                    // A label straight after ')' names the group that just closed.
                    Some(n) if !nodes[n].children.is_empty() && nodes[n].label.is_none() => {
                        nodes[n].label = Some(name);
                    }
                    Some(_) => anyhow::bail!("Newick: two labels in a row near byte {}", p.i),
                    None => {
                        let n = push(&mut nodes, &stack);
                        nodes[n].label = Some(name);
                        cur = Some(n);
                    }
                },
                Token::Length => {
                    // `(:0.1,:0.2)` is legal Newick: an unnamed tip with a length.
                    if cur.is_none() {
                        cur = Some(push(&mut nodes, &stack));
                    }
                }
                Token::Comma => {
                    if stack.is_empty() {
                        anyhow::bail!(
                            "Newick: ',' outside any group at byte {}: a file must hold ONE tree",
                            p.i
                        );
                    }
                    if cur.is_none() {
                        push(&mut nodes, &stack);
                    }
                    cur = None;
                }
                Token::Close => {
                    if cur.is_none() && stack.last().is_some_and(|&s| !nodes[s].children.is_empty()) {
                        push(&mut nodes, &stack);
                    }
                    let n = stack
                        .pop()
                        .ok_or_else(|| anyhow::anyhow!("Newick: unbalanced ')' at byte {}", p.i))?;
                    cur = Some(n);
                }
                Token::Semicolon => {
                    finished = true;
                    break;
                }
            }
        }

        if !stack.is_empty() {
            anyhow::bail!("Newick: {} unclosed '(': the tree ends mid-group", stack.len());
        }
        if nodes.is_empty() {
            anyhow::bail!("Newick: the file contains no tree");
        }
        if finished {
            // Anything after the terminating ';' is a second tree (a bootstrap set, say).
            // Counting origins on the wrong tree is exactly the silent error to refuse.
            if let Some(tok) = p.next_token()? {
                let _ = tok;
                anyhow::bail!(
                    "Newick: content after the terminating ';': supply a single tree, not a \
                     multi-tree file"
                );
            }
        }

        let root = cur.unwrap_or(0);
        if nodes[root].parent != NO_PARENT {
            anyhow::bail!("Newick: the tree has no single root");
        }

        // Tips, and the postorder walk, computed once. The walk is iterative for the
        // same reason the parser is: a clonal ladder is as deep as it is wide.
        let mut tips = Vec::new();
        let mut tip_slot = vec![u32::MAX; nodes.len()];
        let mut postorder = Vec::with_capacity(nodes.len());
        let mut work: Vec<(usize, usize)> = vec![(root, 0)];
        while let Some((v, child)) = work.pop() {
            if child < nodes[v].children.len() {
                work.push((v, child + 1));
                work.push((nodes[v].children[child], 0));
            } else {
                if nodes[v].children.is_empty() {
                    tip_slot[v] = tips.len() as u32;
                    tips.push(v);
                }
                postorder.push(v);
            }
        }
        if postorder.len() != nodes.len() {
            anyhow::bail!("Newick: the tree is disconnected (parsed a cycle or a stray node)");
        }
        if tips.is_empty() {
            anyhow::bail!("Newick: the tree has no tips");
        }
        Ok(Tree { nodes, root, tips, tip_slot, postorder })
    }
}

/// One Newick token. Branch lengths carry no information for parsimony, so `Length`
/// records only that one was consumed (and validated).
enum Token {
    Open,
    Close,
    Comma,
    Semicolon,
    Length,
    Label(String),
}

struct Parser<'a> {
    src: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    /// Skip whitespace and `[...]` comments, which Newick allows anywhere.
    fn skip_filler(&mut self) -> anyhow::Result<()> {
        loop {
            while self.i < self.src.len() && self.src[self.i].is_ascii_whitespace() {
                self.i += 1;
            }
            if self.i < self.src.len() && self.src[self.i] == b'[' {
                let start = self.i;
                self.i += 1;
                while self.i < self.src.len() && self.src[self.i] != b']' {
                    self.i += 1;
                }
                if self.i >= self.src.len() {
                    anyhow::bail!("Newick: unterminated '[' comment opened at byte {start}");
                }
                self.i += 1;
                continue;
            }
            return Ok(());
        }
    }

    fn next_token(&mut self) -> anyhow::Result<Option<Token>> {
        self.skip_filler()?;
        let Some(&c) = self.src.get(self.i) else {
            return Ok(None);
        };
        match c {
            b'(' => {
                self.i += 1;
                Ok(Some(Token::Open))
            }
            b')' => {
                self.i += 1;
                Ok(Some(Token::Close))
            }
            b',' => {
                self.i += 1;
                Ok(Some(Token::Comma))
            }
            b';' => {
                self.i += 1;
                Ok(Some(Token::Semicolon))
            }
            b':' => {
                self.i += 1;
                self.skip_filler()?;
                let start = self.i;
                while self.i < self.src.len()
                    && !matches!(self.src[self.i], b'(' | b')' | b',' | b';' | b'[')
                    && !self.src[self.i].is_ascii_whitespace()
                {
                    self.i += 1;
                }
                let text = String::from_utf8_lossy(&self.src[start..self.i]).into_owned();
                if text.parse::<f64>().is_err() {
                    anyhow::bail!(
                        "Newick: branch length {:?} at byte {} is not a number",
                        text,
                        start
                    );
                }
                Ok(Some(Token::Length))
            }
            b'\'' | b'"' => {
                let quote = c;
                let start = self.i;
                self.i += 1;
                let mut out = String::new();
                loop {
                    let Some(&b) = self.src.get(self.i) else {
                        anyhow::bail!("Newick: unterminated quoted label opened at byte {start}");
                    };
                    self.i += 1;
                    if b == quote {
                        // A doubled quote is a literal one, the Newick escape.
                        if self.src.get(self.i) == Some(&quote) {
                            out.push(quote as char);
                            self.i += 1;
                            continue;
                        }
                        break;
                    }
                    out.push(b as char);
                }
                Ok(Some(Token::Label(out)))
            }
            _ => {
                let start = self.i;
                while self.i < self.src.len()
                    && !matches!(self.src[self.i], b'(' | b')' | b',' | b';' | b':' | b'[')
                    && !self.src[self.i].is_ascii_whitespace()
                {
                    self.i += 1;
                }
                if self.i == start {
                    anyhow::bail!("Newick: unexpected byte {:?} at {}", c as char, start);
                }
                Ok(Some(Token::Label(
                    String::from_utf8_lossy(&self.src[start..self.i]).into_owned(),
                )))
            }
        }
    }
}

/// A tree whose tips have been matched to the cohort's sample indices.
///
/// The join is where a silent corruption would otherwise live: a tip that matches no
/// sample would be scored as a non-carrier and a sample that matches no tip would
/// vanish, and either one changes the origin count without changing anything visible.
/// [`SampleTree::join`] therefore refuses both.
#[derive(Debug, Clone)]
pub struct SampleTree {
    tree: Tree,
    /// Sample index carried by each tip slot.
    tip_sample: Vec<usize>,
    /// Tip slot of each sample index (the inverse of `tip_sample`).
    sample_slot: Vec<usize>,
}

impl SampleTree {
    /// Match the tree's tips to `samples` (sample index = position in the slice).
    ///
    /// Both sides must name exactly the same set. Anything else is an error naming the
    /// counts and a few examples, because the alternative is an origin count computed
    /// over a cohort the user did not supply.
    pub fn join(tree: Tree, samples: &[String]) -> anyhow::Result<SampleTree> {
        let unnamed = tree.tip_labels().filter(|l| l.is_none()).count();
        if unnamed > 0 {
            anyhow::bail!(
                "the tree has {unnamed} unnamed tip(s), which cannot be matched to a sample. \
                 Every tip must carry the sample's name."
            );
        }
        let mut by_name: HashMap<&str, usize> = HashMap::new();
        for (slot, label) in tree.tip_labels().enumerate() {
            let name = label.expect("checked above");
            if by_name.insert(name, slot).is_some() {
                anyhow::bail!("the tree has more than one tip named {name:?}");
            }
        }
        let mut seen: HashMap<&str, usize> = HashMap::new();
        for (idx, s) in samples.iter().enumerate() {
            if seen.insert(s.as_str(), idx).is_some() {
                anyhow::bail!(
                    "two samples are both named {s:?}; sample names must be unique to match \
                     tree tips"
                );
            }
        }

        let missing: Vec<&str> = samples
            .iter()
            .map(String::as_str)
            .filter(|s| !by_name.contains_key(s))
            .collect();
        let extra: Vec<&str> = by_name.keys().copied().filter(|t| !seen.contains_key(t)).collect();
        if !missing.is_empty() || !extra.is_empty() {
            let mut extra = extra;
            extra.sort_unstable();
            anyhow::bail!(
                "--tree does not match the cohort: {} of {} sample(s) have no tip ({}), and {} \
                 of {} tip(s) match no sample ({}). Names must agree exactly; prune the tree to \
                 the samples analysed, or rename them. Dropping either side silently would \
                 change every origin count.",
                missing.len(),
                samples.len(),
                examples(&missing),
                extra.len(),
                by_name.len(),
                examples(&extra),
            );
        }

        let mut tip_sample = vec![0usize; tree.n_tips()];
        let mut sample_slot = vec![0usize; samples.len()];
        for (name, slot) in &by_name {
            let idx = seen[name];
            tip_sample[*slot] = idx;
            sample_slot[idx] = *slot;
        }
        Ok(SampleTree { tree, tip_sample, sample_slot })
    }

    /// The underlying tree.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Sample index of a tip slot (used by the tests that check the join direction).
    #[allow(dead_code)]
    pub fn sample_of_slot(&self, slot: usize) -> usize {
        self.tip_sample[slot]
    }

    /// How many times an allele carried by `carriers` (sample indices) arose.
    ///
    /// Cost is the number of carriers plus one walk of the tree, so a cohort of
    /// thousands where almost every allele is rare stays affordable. The scratch is
    /// left clear, so one buffer serves every allele of a run.
    pub fn origins(
        &self,
        carriers: impl IntoIterator<Item = usize>,
        min_support: u32,
        scratch: &mut OriginScratch,
    ) -> u32 {
        scratch.touched.clear();
        for s in carriers {
            let Some(&slot) = self.sample_slot.get(s) else { continue };
            if !scratch.tip_state[slot] {
                scratch.tip_state[slot] = true;
                scratch.touched.push(slot);
            }
        }
        let n = self.tree.count_origins(min_support, scratch);
        for i in 0..scratch.touched.len() {
            let slot = scratch.touched[i];
            scratch.tip_state[slot] = false;
        }
        scratch.touched.clear();
        n
    }
}

/// Up to three names from a list, for a diagnostic message.
fn examples(names: &[&str]) -> String {
    if names.is_empty() {
        return "none".to_string();
    }
    let shown = names.len().min(3);
    let mut out = names[..shown].join(", ");
    if names.len() > shown {
        out.push_str(&format!(", … (+{} more)", names.len() - shown));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tree(s: &str) -> Tree {
        Tree::parse_newick(s).unwrap_or_else(|e| panic!("parse {s:?}: {e}"))
    }

    /// Origins of the named carriers on a tree whose tips are the given names.
    fn origins(newick: &str, samples: &[&str], carriers: &[&str], min_support: u32) -> u32 {
        let names: Vec<String> = samples.iter().map(|s| s.to_string()).collect();
        let st = SampleTree::join(tree(newick), &names).expect("join");
        let mut scratch = OriginScratch::for_tree(st.tree());
        let idx: Vec<usize> = carriers
            .iter()
            .map(|c| samples.iter().position(|s| s == c).expect("carrier is a sample"))
            .collect();
        st.origins(idx, min_support, &mut scratch)
    }

    #[test]
    fn newick_shapes_that_must_parse() {
        // Lengths, internal labels, bootstrap values, quoted names, comments, a
        // unary node, a polytomy, and no trailing semicolon.
        let t = tree("((A:0.1,B:0.2)0.95:0.3,(C,D,E));");
        assert_eq!(t.n_tips(), 5);
        assert_eq!(t.tip_labels().flatten().collect::<Vec<_>>(), vec!["A", "B", "C", "D", "E"]);
        assert_eq!(tree("('a b','c,d')").n_tips(), 2);
        assert_eq!(
            tree("('it''s',B);").tip_labels().flatten().collect::<Vec<_>>(),
            vec!["it's", "B"]
        );
        assert_eq!(tree("(A[a comment],B)[root];").n_tips(), 2);
        assert_eq!(tree("((A));").n_tips(), 1, "a unary node is still one tip");
        assert_eq!(tree("A;").n_tips(), 1, "a one-tip tree");
        // An underscore is part of the name, NOT a space: renaming ERR_1234 to
        // "ERR 1234" would break the join against the VCF header.
        assert_eq!(tree("(ERR_1234,B);").tip_labels().flatten().next(), Some("ERR_1234"));
        // Unnamed tips with lengths are legal Newick and must parse (the join then
        // rejects them, which is a much clearer error than a parse failure).
        assert_eq!(tree("(:0.1,:0.2);").n_tips(), 2);
        assert_eq!(t.n_nodes(), 5 + 3, "5 tips, 2 internal groups, 1 root");
    }

    #[test]
    fn malformed_newick_is_refused_not_guessed() {
        for bad in [
            "((A,B);",              // unclosed group
            "(A,B));",              // unbalanced close
            "A,B;",                 // two roots
            "(A,B)(C,D);",          // '(' after a completed node
            "(A B,C);",             // two labels in a row
            "(A:x,B);",             // a branch length that is not a number
            "(A,B);(C,D);",         // two trees in one file
            "(A[unterminated,B);",  // unterminated comment
            "('A,B);",              // unterminated quote
            "",                     // nothing at all
        ] {
            assert!(
                Tree::parse_newick(bad).is_err(),
                "{bad:?} must be refused, not silently reinterpreted"
            );
        }
    }

    #[test]
    fn a_clonal_allele_has_one_origin_and_a_convergent_one_has_many() {
        // Two lineages of four. The whole point of the feature, in eight tips.
        let nwk = "(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)));";
        let all = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"];
        // One clade: one origin, however many samples carry it.
        assert_eq!(origins(nwk, &all, &["A1", "A2", "A3", "A4"], 1), 1);
        assert_eq!(origins(nwk, &all, &["A1", "A2"], 1), 1);
        // The same allele in two separate clades: two origins, and the count does not
        // care that the carrier COUNT is identical to the clonal case above.
        assert_eq!(origins(nwk, &all, &["A1", "A2", "B1", "B2"], 1), 2);
        // Scattered one per pair: four origins from four carriers.
        assert_eq!(origins(nwk, &all, &["A1", "A3", "B1", "B3"], 1), 4);
        // Every sample carries it: it arose once, above the root, not zero times.
        assert_eq!(origins(nwk, &all, &all, 1), 1);
        // Nobody carries it: no origins.
        assert_eq!(origins(nwk, &all, &[], 1), 0);
    }

    /// Sixteen tips in four lineages of four, which is enough structure for a handful
    /// of false calls to stay a minority.
    const FOUR_LINEAGES: &str = "(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),\
                                 ((C1,C2),(C3,C4)),((D1,D2),(D3,D4)));";
    const SIXTEEN: [&str; 16] = [
        "A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3",
        "D4",
    ];

    #[test]
    fn min_support_kills_the_false_call_failure_mode() {
        // The measured failure: a genuine single-origin clade allele plus three random
        // false calls scores four origins, which at genome scale is a significant hit
        // manufactured out of sequencing error.
        let real_plus_noise = ["A1", "A2", "A3", "A4", "B1", "C1", "D1"];
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &real_plus_noise, 1), 4, "parsimony is fooled");
        assert_eq!(
            origins(FOUR_LINEAGES, &SIXTEEN, &real_plus_noise, 2),
            1,
            "requiring two carriers per origin leaves only the real clade"
        );
        // Twenty-five singleton false calls would be no different in kind: every one of
        // them is its own tip, and no tip can subtend two carriers.
        let noise: Vec<&str> = SIXTEEN.iter().copied().filter(|s| s.ends_with('1')).collect();
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &noise, 2), 0, "pure noise scores nothing");
        // And the filter is conservative in the other direction too: a genuine
        // singleton origin goes uncounted, which is the right way to be wrong.
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &["A1"], 2), 0);
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &["A1"], 1), 1);
        // Supported origins survive the filter untouched: this is the shape of a real
        // convergent site, four independent pairs.
        let convergent = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2"];
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &convergent, 2), 4);
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &convergent, 1), 4);
        // A higher bar demands bigger clades, and drops those four two-carrier origins.
        assert_eq!(origins(FOUR_LINEAGES, &SIXTEEN, &convergent, 3), 0);
    }

    #[test]
    fn a_hole_inside_a_carrier_clade_does_not_manufacture_origins() {
        // Missing data or a reversion inside a clade: parsimony can call this one gain
        // plus one loss, which must stay ONE origin rather than becoming two.
        let nwk = "(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)));";
        let all = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"];
        assert_eq!(origins(nwk, &all, &["A1", "A2", "A4"], 1), 1);
        assert_eq!(origins(nwk, &all, &["A1", "A2", "A4"], 2), 1);
    }

    #[test]
    fn a_polytomy_costs_one_origin_per_clade_not_one_per_branch() {
        // An unresolved node is exactly where naive parsimony inflates origins. The
        // multifurcating Fitch rule keeps a majority-carrier polytomy at one gain.
        let nwk = "(((A1,A2,A3,A4,A5),(B1,B2,B3,B4,B5)));";
        let all = ["A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"];
        assert_eq!(origins(nwk, &all, &["A1", "A2", "A3", "A4", "A5"], 1), 1);
        // Three of the five tips of one polytomy: still one gain there (the node's
        // state is carrier by majority) rather than three.
        assert_eq!(origins(nwk, &all, &["A1", "A2", "A3"], 1), 1);
    }

    #[test]
    fn the_join_refuses_every_way_the_names_can_disagree() {
        let nwk = "((A,B),(C,D));";
        let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        // Exact agreement is fine, in any order.
        assert!(SampleTree::join(tree(nwk), &names(&["D", "C", "B", "A"])).is_ok());
        // A sample with no tip, a tip with no sample, a duplicated tip, a duplicated
        // sample, and an unnamed tip: all fatal.
        for (t, s) in [
            (nwk, vec!["A", "B", "C", "D", "E"]),
            (nwk, vec!["A", "B", "C"]),
            ("((A,A),(C,D));", vec!["A", "C", "D"]),
            (nwk, vec!["A", "B", "C", "C"]),
            ("((A,B),(C,:0.1));", vec!["A", "B", "C", "D"]),
        ] {
            let err = SampleTree::join(tree(t), &names(&s))
                .expect_err("mismatched names must be fatal");
            let msg = err.to_string();
            assert!(!msg.is_empty(), "the error must say what went wrong");
        }
    }

    #[test]
    fn the_join_maps_names_and_not_positions() {
        // Tip order and sample order are independent, so a positional join would be
        // wrong here and silently so.
        let st = SampleTree::join(
            tree("((A,B),(C,D));"),
            &["D".into(), "C".into(), "B".into(), "A".into()],
        )
        .expect("join");
        // Tip slots are A,B,C,D; samples are D,C,B,A.
        assert_eq!(st.sample_of_slot(0), 3, "tip A is sample index 3");
        assert_eq!(st.sample_of_slot(3), 0, "tip D is sample index 0");
        // And the origin count follows the names: samples 3 and 2 are tips A and B,
        // one clade, one origin.
        let mut scratch = OriginScratch::for_tree(st.tree());
        assert_eq!(st.origins([3usize, 2], 1, &mut scratch), 1);
        // Samples 3 and 0 are tips A and D, in different clades: two origins.
        assert_eq!(st.origins([3usize, 0], 1, &mut scratch), 2);
    }

    #[test]
    fn scratch_is_reusable_across_alleles() {
        // The carrier marks are cleared by the caller, so a second allele must not see
        // the first one's carriers. A leak here would silently merge alleles.
        let st = SampleTree::join(
            tree("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)));"),
            &(1..=8)
                .map(|i| format!("{}{}", if i <= 4 { "A" } else { "B" }, (i - 1) % 4 + 1))
                .collect::<Vec<_>>(),
        )
        .expect("join");
        let mut scratch = OriginScratch::for_tree(st.tree());
        assert_eq!(st.origins([0usize, 1], 1, &mut scratch), 1);
        assert_eq!(st.origins([4usize, 5], 1, &mut scratch), 1);
        assert_eq!(st.origins([0usize, 1, 4, 5], 1, &mut scratch), 2);
        assert_eq!(st.origins([], 1, &mut scratch), 0);
    }

    #[test]
    fn a_deep_ladder_tree_does_not_overflow_the_stack() {
        // A clonal cohort produces ladder-shaped trees: 5,000 nested groups would blow
        // a recursive parser or a recursive traversal.
        let n = 5_000;
        let mut nwk = String::new();
        for i in 0..n {
            nwk.push('(');
            nwk.push_str(&format!("t{i},"));
        }
        nwk.push_str("tlast");
        for _ in 0..n {
            nwk.push(')');
        }
        nwk.push(';');
        let t = tree(&nwk);
        assert_eq!(t.n_tips(), n + 1);
        let names: Vec<String> =
            (0..n).map(|i| format!("t{i}")).chain(std::iter::once("tlast".to_string())).collect();
        let st = SampleTree::join(t, &names).expect("join");
        let mut scratch = OriginScratch::for_tree(st.tree());
        // The two deepest tips are sisters: one origin.
        assert_eq!(st.origins([n - 1, n], 1, &mut scratch), 1);
    }
}
