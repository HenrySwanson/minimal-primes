use crate::sequences::DigitSeq;

pub struct Frontier<T> {
    /// maps weight to nodes; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    /// items lower than `min_allowed_weight` must be empty
    by_weight: Vec<Vec<T>>,
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
    pub fn new() -> Self {
        Self {
            by_weight: vec![],
            min_allowed_weight: 0,
        }
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

    /// Removes the families of the least weight from the structure.
    /// Once this method is called, elements of weight < self.min_weight
    /// should not be inserted! Elements with exactly the minimum weight
    /// are allowed though (lateral exploration).
    pub fn pop(&mut self) -> Option<Vec<T>> {
        for (i, layer) in self.by_weight.iter_mut().enumerate() {
            if i < self.min_allowed_weight {
                assert!(layer.is_empty())
            }

            // Take the first non-empty layer we see
            if !layer.is_empty() {
                return Some(std::mem::take(layer));
            }
        }
        None
    }

    pub fn extend(&mut self, iter: impl IntoIterator<Item = T>) {
        for node in iter {
            self.put(node);
        }
    }

    pub fn put(&mut self, node: T) {
        let weight = node.weight();
        debug_assert!(weight >= self.min_allowed_weight);

        loop {
            match self.by_weight.get_mut(weight) {
                Some(layer) => {
                    layer.push(node);
                    break;
                }
                None => self.by_weight.push(vec![]),
            }
        }
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
        self.inner
            .retain(|other| !is_proper_substring(&seq, other));

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
