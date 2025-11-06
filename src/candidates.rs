use crate::digits::DigitSeq;

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
    pub fn indices_none(&self) -> CandidateIndices {
        CandidateIndices {
            idxs: vec![],
            start_unknown: self.inner.len(),
        }
    }

    /// Returns a [CandidateIndices] containing all of the current elements.
    pub fn indices_all(&self) -> CandidateIndices {
        CandidateIndices {
            idxs: vec![],
            start_unknown: 0,
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

impl CandidateIndices {
    /// Inserts an index into this struct.
    pub fn add(&mut self, idx: usize) {
        self.idxs.push(idx);
    }
}
