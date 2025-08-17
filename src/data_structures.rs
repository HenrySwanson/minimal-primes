use std::collections::VecDeque;

use crate::digits::DigitSeq;

/// items lower than `min_allowed_weight` must be empty
pub struct WeightedVec<T> {
    elements: Vec<VecDeque<T>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

pub struct CandidateSequences {
    // items in here never change their index. removing an
    // item is just replacing it with None. indices also can't
    // be re-used.
    inner: Vec<Option<DigitSeq>>,
}

#[derive(Debug, Clone)]
pub struct CandidateIndices {
    /// Indices of candidates we're intentionally including
    idxs: Vec<usize>,
    /// Since the set of candidates might grow in the meantime, we
    /// need to track the start of where new candidates are.
    start_unknown: usize,
}

#[derive(Debug, Clone)]
pub struct AppendTree<T> {
    nodes: Vec<AppendTreeNode<T>>,
}

#[derive(Debug, Clone, Copy)]
pub struct AppendTreeNodeID(usize);

#[derive(Debug, Clone)]
struct AppendTreeNode<T> {
    contents: Vec<Content<T>>,
    // TODO: parent?
}

#[derive(Debug, Clone)]
enum Content<T> {
    Item(T),
    Child { tag: T, idx: usize },
}

impl<T> WeightedVec<T> {
    pub fn new() -> Self {
        Self {
            elements: vec![],
            min_allowed_weight: 0,
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.elements.iter().flatten()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.iter().all(|layer| layer.is_empty())
    }

    pub fn len(&self) -> usize {
        self.elements.iter().map(|layer| layer.len()).sum()
    }

    /// Returns the minimum weight present in this collection. This can be
    /// different from [self.min_allowed_weight], because that layer might
    /// be empty (or even the layers above).
    pub fn min_weight(&self) -> Option<usize> {
        self.elements.iter().position(|layer| !layer.is_empty())
    }

    pub fn put(&mut self, item: T, weight: usize) {
        debug_assert!(weight >= self.min_allowed_weight);

        loop {
            // Keep appending empty deques until we reach the right weight
            match self.elements.get_mut(weight) {
                Some(layer) => {
                    layer.push_back(item);
                    break;
                }
                None => self.elements.push(VecDeque::new()),
            }
        }
    }

    // TODO: should this be public?
    pub fn find_first_non_empty_layer_mut(&mut self) -> Option<&mut VecDeque<T>> {
        for (i, layer) in self.elements.iter_mut().enumerate() {
            if i < self.min_allowed_weight {
                assert!(
                    layer.is_empty(),
                    "all layers below min_allowed_weight must be empty"
                );
            }

            if !layer.is_empty() {
                return Some(layer);
            }
        }

        None
    }
}

impl CandidateSequences {
    pub fn new() -> Self {
        Self { inner: vec![] }
    }

    pub fn len(&self) -> usize {
        self.iter().count()
    }

    pub fn iter(&self) -> impl Iterator<Item = &DigitSeq> {
        self.inner.iter().flatten()
    }

    pub fn insert(&mut self, seq: DigitSeq) {
        // Does this contain an existing candidate? Reject it.
        for other in self.iter() {
            // we shouldn't be inserting duplicates ever
            assert_ne!(&seq, other);
            if is_proper_substring(other, &seq) {
                return;
            }
        }

        // Okay, we're definitely going to insert this. Remove any candidates that
        // contain this.
        // We can't merge this loop with the one above! We need to make a complete
        // decision first before we start modifying things.
        for slot in self.inner.iter_mut() {
            if let Some(other) = slot {
                if is_proper_substring(&seq, other) {
                    *slot = None;
                }
            }
        }

        // Insert
        self.inner.push(Some(seq));
    }

    // TODO: ugly as sin but i can deal with it later after giving this thing some stable indices
    pub fn clone_and_sort_and_iter(&self) -> impl Iterator<Item = &DigitSeq> {
        let mut primes: Vec<_> = self.iter().collect();
        primes.sort();
        primes.into_iter()
    }

    pub fn new_indices(&self) -> CandidateIndices {
        CandidateIndices {
            idxs: vec![],
            start_unknown: self.inner.len(),
        }
    }

    pub fn get_many<'a, 'b>(
        &'a self,
        indices: &'b CandidateIndices,
    ) -> impl Iterator<Item = (usize, &'a DigitSeq)> + 'b
    where
        'a: 'b,
    {
        indices
            .idxs
            .iter()
            .copied()
            .chain(indices.start_unknown..self.inner.len())
            .flat_map(|idx| self.inner[idx].as_ref().map(|val| (idx, val)))
    }
}

impl Default for CandidateSequences {
    fn default() -> Self {
        Self::new()
    }
}

impl CandidateIndices {
    pub fn empty() -> Self {
        Self {
            idxs: vec![],
            start_unknown: 0,
        }
    }

    pub fn add(&mut self, idx: usize) {
        self.idxs.push(idx);
    }
}

impl<T> AppendTree<T> {
    pub fn new() -> Self {
        let root = AppendTreeNode::new();
        Self { nodes: vec![root] }
    }

    pub fn root(&self) -> AppendTreeNodeID {
        AppendTreeNodeID(0)
    }

    pub fn append(&mut self, node_id: AppendTreeNodeID, item: T) -> Result<(), T> {
        match self.nodes.get_mut(node_id.0) {
            Some(node) => {
                node.contents.push(Content::Item(item));
                Ok(())
            }
            None => Err(item),
        }
    }

    pub fn make_child(&mut self, node_id: AppendTreeNodeID, tag: T) -> Result<AppendTreeNodeID, T> {
        // gotta get this before the mutable borrow begins
        let num_nodes = self.nodes.len();

        match self.nodes.get_mut(node_id.0) {
            Some(node) => {
                // Add a new node to the whole tree, and push it into
                // this node's children.
                let child_idx = num_nodes;
                node.contents.push(Content::Child {
                    tag,
                    idx: child_idx,
                });
                self.nodes.push(AppendTreeNode::new());
                Ok(AppendTreeNodeID(child_idx))
            }
            None => Err(tag),
        }
    }
}

// is this principled? no. does it work well? yeah
impl<T: std::fmt::Display> AppendTree<T> {
    pub fn pretty_print_to_stdout(&self) {
        self.pretty_print_helper(0, 0);
    }

    fn pretty_print_helper(&self, node_idx: usize, indent: usize) {
        let node = &self.nodes[node_idx];

        for content in &node.contents {
            match content {
                Content::Item(t) => println!("{:indent$}{}", "", t, indent = indent * 2),
                Content::Child {
                    tag,
                    idx: child_idx,
                } => {
                    println!("{:indent$}{}", "", tag, indent = indent * 2);
                    self.pretty_print_helper(*child_idx, indent + 1);
                }
            }
        }
    }
}

impl<T> AppendTreeNode<T> {
    pub fn new() -> Self {
        Self { contents: vec![] }
    }
}

pub fn is_proper_substring(needle: &DigitSeq, haystack: &DigitSeq) -> bool {
    // Save some time when the needle is too large, and also, rule out identical
    // strings.
    if needle.0.len() >= haystack.0.len() {
        return false;
    }

    let mut iter = haystack.0.iter().copied();
    for d in needle.0.iter().copied() {
        // Chomp iter until we find that digit
        loop {
            match iter.next() {
                Some(d2) if d == d2 => break,
                Some(_) => {}
                None => return false,
            }
        }
    }
    // If we got here, then hooray, this is a match!
    true
}
