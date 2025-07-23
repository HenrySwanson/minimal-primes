use std::collections::VecDeque;

use crate::digits::DigitSeq;

pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    /// items lower than `min_allowed_weight` must be empty
    by_weight: Vec<VecDeque<T>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

pub trait Weight {
    fn weight(&self) -> usize;
}

pub struct CandidateSequences {
    inner: Vec<DigitSeq>,
}

impl<T: Weight> Frontier<T> {
    pub fn new(initial: T) -> Self {
        let mut ret = Self {
            by_weight: vec![],
            min_allowed_weight: initial.weight(),
        };
        ret.put(initial);
        ret
    }

    pub fn len(&self) -> usize {
        self.by_weight.iter().map(|layer| layer.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.by_weight.iter().all(|layer| layer.is_empty())
    }

    /// Returns the minimum weight present in this collection. This can be
    /// different from [self.min_allowed_weight], because that layer might
    /// be empty (or even the layers above).
    pub fn min_weight(&self) -> Option<usize> {
        self.by_weight.iter().position(|layer| !layer.is_empty())
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.by_weight.iter().flatten()
    }

    fn put(&mut self, node: T) {
        let weight = node.weight();
        debug_assert!(weight >= self.min_allowed_weight);

        loop {
            match self.by_weight.get_mut(weight) {
                Some(layer) => {
                    layer.push_back(node);
                    break;
                }
                None => self.by_weight.push(VecDeque::new()),
            }
        }
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
        let layer = match self.find_first_non_empty_layer_mut() {
            Some(layer) => layer,
            None => return false,
        };

        let node = layer.pop_front().expect("non-empty layer");
        for child in f(node) {
            self.put(child);
        }
        true
    }

    pub fn explore_one_level(&mut self, mut f: impl FnMut(T) -> Vec<T>) -> bool {
        let layer = match self.find_first_non_empty_layer_mut() {
            Some(layer) => std::mem::take(layer),
            None => return false,
        };

        for node in layer {
            for child in f(node) {
                self.put(child);
            }
        }

        true
    }

    fn find_first_non_empty_layer_mut(&mut self) -> Option<&mut VecDeque<T>> {
        for (i, layer) in self.by_weight.iter_mut().enumerate() {
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
