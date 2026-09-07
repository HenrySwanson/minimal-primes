use crate::digits::{DigitSeq, DigitSet};
use stable_vec::StableVec;

/// This struct contains all of the candidates for minimal primes we've
/// discovered so far.
///
/// There are two useful properties it has:
/// - Stable indices: after an element is inserted, its index never changes,
///   even when other elements are removed. This is needed to make
///   [CandidateIndices] work correctly.
/// - Antichain: no two sequences in this chain contain the other. This helps
///   is necessary, but not sufficient (see below) to ensure we're collecting
///   the minimal primes. See the documentation for [CandidateSequences::insert]
///   for more details.
///
/// The reason that we use the term "candidate", instead of "minimal", is
/// that, depending on how we explore the search space, we may not discover
/// primes in ascending order by length. If we did, callers could easily
/// guarantee that primes are minimal before adding them to this list. However,
/// that is not the case, so we simply maintain a set of candidates, and it's
/// up to the caller to ensure that all smaller primes are discovered.
#[derive(Debug)]
pub struct CandidateSequences {
    // i was handrolling this before, but let's just use an existing crate...
    inner: StableVec<(DigitSeq, DigitSet)>,
}

/// A candidate minimal prime, together with its index and precomputed digit set.
#[derive(Debug, Clone, Copy)]
pub struct Candidate<'a> {
    pub idx: usize,
    pub seq: &'a DigitSeq,
    pub mask: DigitSet,
}

/// A collection of indices for [CandidateSequences] that automatically extends
/// to include new elements added to the container. This makes it useful for
/// tracking a filtered list of elements; if new elements are added, subsequent
/// calls to [CandidateSequences::get_many] will surface them.
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
            inner: StableVec::new(),
        }
    }

    /// The number of primes in this container. Will differ from [CandidateSequences::upper_bound]
    /// if elements have been removed.
    pub fn num_elements(&self) -> usize {
        self.inner.num_elements()
    }

    /// One past the largest index any current element has.
    ///
    /// This is what you want when you have checked all known primes, but want to
    /// be aware of any new ones that arise. It will differ from [CandidateSequences::num_elements]
    /// if elements have been removed.
    pub fn upper_bound(&self) -> usize {
        self.inner.next_push_index()
    }

    pub fn iter(&self) -> impl Iterator<Item = &DigitSeq> {
        self.inner.values().map(|(seq, _)| seq)
    }

    /// Inserts the given sequence into the list, but preserves the antichain property.
    ///
    /// If this sequence contains any existing sequence in the list, it is discarded,
    /// and conversely, if it is a subsequence of any sequence(s) in the list, those
    /// are discarded, and the given sequence is inserted instead.
    pub fn insert(&mut self, seq: DigitSeq) {
        // TODO: do we need to test this here? given how dense this vector is, it might
        // be worth just keeping the extra candidates and eliminate them all at the end.
        // also we don't really need the "does contain" check; callers all check against
        // "possible contained primes"... might get to throw away some work here :o

        let our_mask = seq.digit_set();

        // Does the new sequence contain any of our existing sequences? If so, exit now.
        for (other_seq, other_mask) in self.inner.values() {
            // we shouldn't be inserting duplicates ever
            debug_assert_ne!(
                seq, *other_seq,
                "Attempted to insert {seq} twice; something is wrong with our search"
            );

            // `properly_contains` is slow, so let's guard it by comparing the digit
            // sets. this filters out like 99% of all candidates, so it's definitely
            // worth doing!
            if our_mask.is_superset_of(*other_mask) && seq.properly_contains(other_seq) {
                return;
            }
        }

        // Okay, we're definitely going to insert this. But due to some peculiarities
        // in how we explore the space, we might be inserted a sequence that's contained
        // in some of our existing entries, which we need to remove.
        //
        // Note that we can't merge this loop with the one above! We need to know for
        // sure we're going to insert this entry before we modify the list at all.
        self.inner.retain(|(other_seq, other_mask)| {
            // This time we're checking whether we are a subset of the other element,
            // so the arguments are flipped.
            let contains_seq =
                other_mask.is_superset_of(our_mask) && other_seq.properly_contains(&seq);
            // `retain` *keeps* the element if this returns `true`, so flip it. makes
            // more sense to me to calculate it like this
            !contains_seq
        });

        // Insert the new element!
        self.inner.push((seq, our_mask));
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
            start_unknown: self.inner.next_push_index(),
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
            .chain(indices.start_unknown..self.inner.next_push_index())
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
        (start..self.inner.next_push_index()).flat_map(|idx| self.get(idx))
    }

    fn get(&self, idx: usize) -> Option<Candidate<'_>> {
        let (seq, set) = self.inner.get(idx)?;
        Some(Candidate {
            idx,
            seq,
            mask: *set,
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
