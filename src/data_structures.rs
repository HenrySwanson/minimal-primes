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
    inner: Vec<DigitSeq>,
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
