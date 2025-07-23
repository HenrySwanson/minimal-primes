use std::collections::VecDeque;

use crate::digits::DigitSeq;

pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    // usize is the node ID, might remove this later
    by_weight: WeightedVec<(T, usize)>,
    /// optional tree for tracing
    /// TODO: make optional
    tree_tracer: TreeTracer,
}

struct TreeTracer {
    nodes: Vec<TreeTracerNode>,
}

struct TreeTracerNode {
    tag: String, // TODO: change that?
    children: Vec<usize>,
}

pub trait Weight {
    fn weight(&self) -> usize;
}

/// items lower than `min_allowed_weight` must be empty
struct WeightedVec<T> {
    elements: Vec<VecDeque<T>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

pub struct CandidateSequences {
    inner: Vec<DigitSeq>,
}

// TODO: remove ToString bound
impl<T: Weight + ToString> Frontier<T> {
    pub fn new(initial: T) -> Self {
        let mut ret = Self {
            by_weight: WeightedVec::new(),
            tree_tracer: TreeTracer { nodes: vec![] },
        };
        ret.tree_tracer.nodes.push(TreeTracerNode {
            tag: initial.to_string(),
            children: vec![],
        });
        ret.put(initial, 0);
        ret
    }

    pub fn len(&self) -> usize {
        self.by_weight.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_weight.is_empty()
    }

    /// Returns the minimum weight present in this collection. This can be
    /// different from [self.min_allowed_weight], because that layer might
    /// be empty (or even the layers above).
    pub fn min_weight(&self) -> Option<usize> {
        self.by_weight.min_weight()
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.by_weight.iter().map(|(value, _)| value)
    }

    fn put(&mut self, node: T, idx: usize) {
        let weight = node.weight();
        self.by_weight.put((node, idx), weight);
    }

    // TODO: return some richer type from the closure?
    /// Pops the next node, and applies the given function to it,
    /// inserting the returned nodes into itself.
    ///
    /// Returned nodes must have weight at least as large as the
    /// popped node.
    ///
    /// Returns false if this structure is empty and there are no
    /// nodes to explore.
    pub fn explore_next(&mut self, f: impl FnOnce(T) -> Vec<T>) -> bool {
        // Pop out an element of least weight
        let layer = match self.by_weight.find_first_non_empty_layer_mut() {
            Some(layer) => layer,
            None => return false,
        };

        let (node, idx) = layer.pop_front().expect("non-empty layer");
        for child in f(node) {
            let child_tag = child.to_string();
            let child_idx = self.tree_tracer.nodes.len();

            self.tree_tracer.nodes[idx].children.push(child_idx);
            self.tree_tracer.nodes.push(TreeTracerNode {
                tag: child_tag,
                children: vec![],
            });

            self.put(child, child_idx);
        }
        true
    }

    pub fn explore_one_level(&mut self, mut f: impl FnMut(T) -> Vec<T>) -> bool {
        let layer = match self.by_weight.find_first_non_empty_layer_mut() {
            Some(layer) => std::mem::take(layer),
            None => return false,
        };

        for (node, idx) in layer {
            for child in f(node) {
                let child_tag = child.to_string();
                let child_idx = self.tree_tracer.nodes.len();

                self.tree_tracer.nodes[idx].children.push(child_idx);
                self.tree_tracer.nodes.push(TreeTracerNode {
                    tag: child_tag,
                    children: vec![],
                });

                self.put(child, child_idx);
            }
        }

        true
    }

    pub fn print_tree_to_stdout(&self) {
        self.tree_tracer.print_tree_to_stdout();
    }
}

impl TreeTracer {
    fn print_tree_to_stdout(&self) {
        self.print_to_stdout_helper(0, 0);
    }

    fn print_to_stdout_helper(&self, node_idx: usize, indent: usize) {
        let node = &self.nodes[node_idx];

        println!("{:indent$}{}", "", node.tag, indent = indent * 2);
        for child_idx in &node.children {
            self.print_to_stdout_helper(*child_idx, indent + 1);
        }
    }
}

impl<T> WeightedVec<T> {
    fn new() -> Self {
        Self {
            elements: vec![],
            min_allowed_weight: 0,
        }
    }

    fn iter(&self) -> impl Iterator<Item = &T> {
        self.elements.iter().flatten()
    }

    fn is_empty(&self) -> bool {
        self.elements.iter().all(|layer| layer.is_empty())
    }

    fn len(&self) -> usize {
        self.elements.iter().map(|layer| layer.len()).sum()
    }

    fn min_weight(&self) -> Option<usize> {
        self.elements.iter().position(|layer| !layer.is_empty())
    }

    fn put(&mut self, item: T, weight: usize) {
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

    fn find_first_non_empty_layer_mut(&mut self) -> Option<&mut VecDeque<T>> {
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
        self.inner.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &DigitSeq> {
        self.inner.iter()
    }

    pub fn insert(&mut self, seq: DigitSeq) {
        // Does this contain an existing candidate? Reject it.
        for other in self.inner.iter_mut() {
            if is_proper_substring(other, &seq) {
                return;
            }
        }

        // Remove any candidates that contain this (if any)
        self.inner.retain(|other| !is_proper_substring(&seq, other));

        // Insert
        self.inner.push(seq);
    }

    pub fn sort(&mut self) {
        self.inner.sort();
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
