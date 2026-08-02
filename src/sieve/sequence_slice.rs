use std::ops::Range;

use bitvec::prelude::BitVec;

use crate::sequence::Sequence;

/// Tracks which n in a range still have a candidate that hasn't been ruled
/// out yet, for a single sequence.
#[derive(Debug)]
pub struct SequenceSlice {
    pub seq: Sequence,
    n_lo: usize,
    n_bitvec: BitVec,
}

impl SequenceSlice {
    pub fn new(seq: Sequence, range: Range<usize>) -> Self {
        Self {
            seq,
            n_lo: range.start,
            n_bitvec: BitVec::repeat(true, range.len()),
        }
    }

    pub fn n_lo(&self) -> usize {
        self.n_lo
    }

    /// The total width of n this slice covers, regardless of how many
    /// candidates remain within it.
    pub fn range_len(&self) -> usize {
        self.n_bitvec.len()
    }

    #[cfg(test)]
    pub fn check_n(&self, n: usize) -> bool {
        self.n_bitvec[n - self.n_lo]
    }

    pub fn num_remaining(&self) -> usize {
        self.n_bitvec.count_ones()
    }

    pub fn iter_remaining(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.n_bitvec.iter_ones().map(|i| self.n_lo + i)
    }

    pub fn eliminate_multiple(&mut self, p: u64, base: u64, start: usize, spacing: usize) {
        let mut idx = start - self.n_lo;

        // It's possible the term we're about to eliminate is actually p itself.
        // Let's avoid that, if so.
        // TODO: actually, this should report a prime immediately! not that that'll
        // happen in the interesting cases, but still!
        if self.seq.check_term_equal(base, p, start) {
            idx += spacing;
        }
        // Insane edge case: it could also be zero! In that case, bump it up twice.
        else if self.seq.check_term_equal(base, 0, start) {
            idx += 2 * spacing;
        }
        // Can it be negative? No, that would not be meaningful for the kinds
        // of sequences we're considering.
        if self.seq.c < 0 {
            // 0th term is (k*1+c) / d, which can only go negative if c is large and negative
            debug_assert!(self.seq.k > self.seq.c.unsigned_abs())
        }

        while let Some(mut slot) = self.n_bitvec.get_mut(idx) {
            slot.set(false);
            idx += spacing;
        }
    }
}
