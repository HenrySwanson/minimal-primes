use crate::digits::{DigitSeq, DigitSet};

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
    /// tracks the digits presents in `inner`, so that we can very quickly
    /// eliminate candidates during containment testing. removed entries
    /// (those that are `None`) get a mask of [REMOVED].
    masks: Vec<DigitSet>,
    /// How many of `inner` are still present. Tracked as we go, because
    /// counting them is O(n) and [CandidateSequences::len] is called often.
    num_present: usize,
}

/// A candidate minimal prime, together with its index and precomputed digit set.
#[derive(Debug, Clone, Copy)]
pub struct Candidate<'a> {
    pub idx: usize,
    pub seq: &'a DigitSeq,
    pub mask: DigitSet,
}

/// The mask stored for a vacated slot. Every bit is set, so it never looks
/// like a subset of anything and scans pass it by.
const REMOVED: DigitSet = DigitSet::from_mask(u64::MAX);

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
        Self {
            inner: vec![],
            masks: vec![],
            num_present: 0,
        }
    }

    /// The number of primes in this container. Will differ from [CandidateSequences::upper_bound]
    /// if elements have been removed.
    pub fn len(&self) -> usize {
        self.num_present
    }

    /// One past the largest index any current element has.
    ///
    /// This is what you want when you have checked all known primes, but want to
    /// be aware of any new ones that arise. It will differ from [CandidateSequences::len]
    /// if elements have been removed.
    pub fn upper_bound(&self) -> usize {
        self.inner.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &DigitSeq> {
        self.inner.iter().flatten()
    }

    // TODO: do we need to test this here? given how dense this vector is, it might
    // be worth just keeping the extra candidates and eliminate them all at the end.
    // also we don't really need the "does contain" check; callers all check against
    // "possible contained primes"... might get to throw away some work here :o
    pub fn insert(&mut self, seq: DigitSeq) {
        let our_mask = seq.digit_set();

        // Does the new sequence contain any of our existing sequences?
        for (idx, &mask) in self.masks.iter().enumerate() {
            // `properly_contains` is slow, if our digits aren't a subset of
            // their digits, we can skip this test. matches are very rare, so
            // this prefilter gets rid of like 99% of the work :)
            if !mask.is_subset_of(our_mask) {
                continue;
            }

            // skip removed elements
            let Some(other) = &self.inner[idx] else {
                continue;
            };

            // we shouldn't be inserting duplicates ever
            assert_ne!(&seq, other);

            // do the actual expensive containment check
            if seq.properly_contains(other) {
                return;
            }
        }

        // Okay, we're definitely going to insert this. But due to some peculiarities
        // in how we explore the space, we might be inserted a sequence that's contained
        // in some of our existing entries, which we need to remove.
        //
        // We can't merge this loop with the one above! We need to make a complete
        // decision first before we start modifying the list.
        for (idx, mask) in self.masks.iter_mut().enumerate() {
            // This time we're checking whether we are a subset of the other element,
            // so the arguments are flipped.
            if !mask.is_superset_of(our_mask) {
                continue;
            }
            if let Some(other) = &self.inner[idx]
                && other.properly_contains(&seq)
            {
                self.inner[idx] = None;
                *mask = REMOVED;
                self.num_present -= 1;
            }
        }

        // Insert
        self.inner.push(Some(seq));
        self.masks.push(our_mask);
        self.num_present += 1;
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
    ) -> impl Iterator<Item = Candidate<'slf>> + 'idx
    where
        'slf: 'idx,
    {
        indices
            .idxs
            .iter()
            .copied()
            .chain(indices.start_unknown..self.inner.len())
            .flat_map(|idx| self.get(idx))
    }

    /// Returns an iterator over the elements from `start` onwards.
    pub fn get_tail<'slf, 'idx>(
        &'slf self,
        start: usize,
    ) -> impl Iterator<Item = Candidate<'slf>> + 'idx
    where
        'slf: 'idx,
    {
        (start..self.inner.len()).flat_map(|idx| self.get(idx))
    }

    fn get(&self, idx: usize) -> Option<Candidate<'_>> {
        self.inner[idx].as_ref().map(|seq| Candidate {
            idx,
            seq,
            mask: self.masks[idx],
        })
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
