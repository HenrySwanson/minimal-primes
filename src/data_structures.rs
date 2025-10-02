use std::collections::VecDeque;

use crate::digits::DigitSeq;

/// TODO: this is an implementation detail that should be removed
pub struct WeightedVec<T> {
    elements: Vec<VecDeque<T>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

/// This struct contains all of the (potentially) minimal primes we discovered
/// so far.
///
/// It is possible for this struct to contain a prime P that turns out not to
/// be minimal, i.e., contains another prime Q. However, when Q is added to this
/// container, P will be removed.
#[derive(Debug)]
pub struct CandidateSequences {
    // items in here never change their index. removing an
    // item is just replacing it with None. indices also can't
    // be re-used.
    inner: Vec<Option<DigitSeq>>,
}

/// A collection of indices for [CandidateSequences] that automatically extends
/// to include new elements added to the container. This makes it useful for
/// tracking elements that may have a particular property.
///
/// As an example, say we have a `CandidateSequences` with five elements and
/// a `CandidateIndices` with indexes 0, 2, and 3. If we add two more elements
/// to the container, then the next time we iterate through it with the indices
/// object, we'll get indexes 0, 2, 3, 5, and 6.
#[derive(Debug, Clone)]
pub struct CandidateIndices {
    /// Indices of candidates we're intentionally including
    idxs: Vec<usize>,
    /// Since the set of candidates might grow in the meantime, we
    /// need to track the start of where new candidates are.
    start_unknown: usize,
}

/// A n-ary tree containing `T`s, where each node contains a `T` and an ordered
/// list of children, which can be more `T`s or more nodes (potentially
/// interleaved).
///
/// This allows structures like:
/// ```
/// Node(T1)   (mix of Ts and child nodes)
/// - T2
/// - T3
/// - Node(T4) (no children)
/// - Node(T5) (only Ts)
///   - T6
///   - T7
///   - T8
/// - Node(T9) (lots of T-less nodes)
///   - Node(T10)
///   - Node(T11)
/// - T12
/// ```
///
/// Used for tracking the evolution of the search tree as we look for minimal
/// primes.
#[derive(Debug, Clone)]
pub struct AppendTree<T> {
    nodes: Vec<AppendTreeNode<T>>,
}

/// Identifies nodes in a [AppendTree].
#[derive(Debug, Clone, Copy)]
pub struct AppendTreeNodeID(usize);

/// Internal nodes in an [AppendTree].
#[derive(Debug, Clone)]
struct AppendTreeNode<T> {
    contents: Vec<Content<T>>,
    // TODO: parent?
}

/// Children of an [AppendTreeNode]
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

    pub fn len(&self) -> usize {
        self.elements.iter().map(|layer| layer.len()).sum()
    }

    /// Returns the minimum weight present in this collection. This can be
    /// different from [self.min_allowed_weight], because that layer might
    /// be empty (or even the layers above).
    pub fn min_weight(&self) -> Option<usize> {
        self.elements.iter().position(|layer| !layer.is_empty())
    }

    /// Inserts a new item into the vector. The weight must be at least the
    /// minimum weight (no backtracking!).
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
            if seq.properly_contains(other) {
                return;
            }
        }

        // Okay, we're definitely going to insert this. Remove any candidates that
        // contain this.
        // We can't merge this loop with the one above! We need to make a complete
        // decision first before we start modifying things.
        for slot in self.inner.iter_mut() {
            if let Some(other) = slot {
                if other.properly_contains(&seq) {
                    *slot = None;
                }
            }
        }

        // Insert
        self.inner.push(Some(seq));
    }

    /// Return a sorted list of the primes contained in this struct.
    ///
    /// TODO: ugly as sin but i can deal with it later after giving this
    /// thing some stable indices
    pub fn clone_and_sort_and_iter(&self) -> impl Iterator<Item = &DigitSeq> {
        let mut primes: Vec<_> = self.iter().collect();
        primes.sort();
        primes.into_iter()
    }

    /// Returns a [CandidateIndices] containing none of the current elements.
    pub fn empty_indices(&self) -> CandidateIndices {
        CandidateIndices {
            idxs: vec![],
            start_unknown: self.inner.len(),
        }
    }

    /// Returns an iterator over the elements represented by the given
    /// [CandidateIndices].
    pub fn get_many<'slf, 'idx>(
        &'slf self,
        indices: &'idx CandidateIndices,
    ) -> impl Iterator<Item = (usize, &'slf DigitSeq)> + 'idx
    where
        'slf: 'idx,
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

// TODO: this should work negatively, i.e., remove instead of add!
impl CandidateIndices {
    /// Returns a set of indices that's completely empty; nothing has
    /// been ruled out.
    pub fn zero() -> Self {
        Self {
            idxs: vec![],
            start_unknown: 0,
        }
    }

    /// Inserts an index into this struct.
    pub fn add(&mut self, idx: usize) {
        self.idxs.push(idx);
    }
}

impl<T> AppendTree<T> {
    /// Creates a tree with a single root node.
    pub fn new() -> Self {
        let root = AppendTreeNode::new();
        Self { nodes: vec![root] }
    }

    /// Returns the root of this tree.
    pub fn root(&self) -> AppendTreeNodeID {
        AppendTreeNodeID(0)
    }

    /// Appends a new item to the given subtree. If the node ID doesn't exist,
    /// returns the item in the `Err` variant.
    pub fn append(&mut self, node_id: AppendTreeNodeID, item: T) -> Result<(), T> {
        match self.nodes.get_mut(node_id.0) {
            Some(node) => {
                node.contents.push(Content::Item(item));
                Ok(())
            }
            None => Err(item),
        }
    }

    /// Appends a new child node to the given subtree, with `tag` as the first
    /// item. If the node ID doesn't exist, returns the item in the `Err` variant.
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
