use itertools::Itertools;
use num_bigint::BigUint;

use crate::digits::{Digit, DigitSeq};

/// A *family* is a subset of digit sequences, specified by concatenating fixed
/// digit sequences and cores, for example, `1[78]*23[9]*4`. Some sequences in
/// this family are:
/// - `1234`
/// - `188882394`
/// - `1787239994`
#[derive(Debug, Clone, PartialEq)]
pub struct Family {
    // invariant: digitseqs.len() = cores.len() + 1
    pub digitseqs: Vec<DigitSeq>,
    pub cores: Vec<Core>,
}

/// A *core* is an unordered set of digits, representing any sequence, of any
/// length made from those digits. We notate cores, and families, with a
/// regex-like syntax, so the core `[134]*` represents any strings made out
/// of only 1s, 3s, and 4s (including the empty string).
#[derive(Debug, Clone, PartialEq)]
pub struct Core {
    digits: Vec<Digit>,
}

/// A *simple* family is a family with exactly one core, which contains only
/// one digit.
///
/// For example, `4[6]*7` is a simple family.
#[derive(Debug, Clone)]
pub struct SimpleFamily {
    pub before: DigitSeq,
    pub center: Digit,
    pub min_repeats: usize,
    pub after: DigitSeq,
}

impl Family {
    /// Creates the family `[0..B]*`, where B is the given base. This family
    /// contains all digit sequences.
    pub fn any(base: u8) -> Self {
        Self {
            digitseqs: vec![DigitSeq::new(), DigitSeq::new()],
            cores: vec![Core::full(base)],
        }
    }

    /// Returns the weight of the sequence, i.e., the sum of the lengths of the
    /// fixed digit sequences. Equivalently, the length of the smallest string in
    /// this family.
    pub fn weight(&self) -> usize {
        self.digitseqs.iter().map(|seq| seq.0.len()).sum()
    }

    /// Reduces the family to an equivalent but simpler form. For example,
    /// it removes empty cores (which can only expand to the empty string).
    pub fn simplify(&mut self) {
        debug_assert_eq!(self.digitseqs.len(), self.cores.len() + 1);

        // Contract out any empty cores
        // Normal for loop won't work because we're mutating the thing
        // we're iterating over.
        let mut i = 0;
        while let Some(core) = self.cores.get(i) {
            if core.is_empty() {
                // join x[]z into xz, and don't increment i, since
                // we've skooched everything one to the left
                let rhs = self.digitseqs.remove(i + 1);
                self.digitseqs[i] += rhs;
                self.cores.remove(i);
            } else {
                i += 1;
            }
        }

        // Delete [0]* at the beginning of a family
        if let Some(first_core) = self.cores.first() {
            let first_seq = &self.digitseqs[0];
            if first_seq.0.is_empty()
                && first_core.len() == 1
                && first_core.iter().next().unwrap() == Digit(0)
            {
                self.digitseqs.remove(0);
                self.cores.remove(0);
            }
        }
    }

    /// Returns the sequence gotten by removing all the cores (equivalently,
    /// replacing them with empty strings).
    pub fn contract(&self) -> DigitSeq {
        DigitSeq(
            self.digitseqs
                .iter()
                .flat_map(|seq| &seq.0)
                .copied()
                .collect(),
        )
    }

    /// Returns the sequence gotten by substituting the given digit for the
    /// specified core, and deleting all the others.
    ///
    /// Does not check that the digit is in that core.
    pub fn substitute(&self, slot: usize, digit: Digit) -> DigitSeq {
        self.substitute_multiple(slot, [digit])
    }

    /// Returns the sequence gotten by substituting the given digits for the
    /// specified core, and deleting all the others.
    ///
    /// Does not check that the digits are in that core.
    pub fn substitute_multiple(
        &self,
        slot: usize,
        digits: impl IntoIterator<Item = Digit>,
    ) -> DigitSeq {
        let mut output = DigitSeq::new();
        for i in 0..=slot {
            output += &self.digitseqs[i];
        }
        for d in digits {
            output += d;
        }
        for i in (slot + 1)..self.digitseqs.len() {
            output += &self.digitseqs[i];
        }
        output
    }

    /// Works like [Self::substitute], but with two slots and two digits. If
    /// the slots are the same, `digit_i` precedes `digit_j`.
    pub fn substitute_two(
        &self,
        slot_i: usize,
        digit_i: Digit,
        slot_j: usize,
        digit_j: Digit,
    ) -> DigitSeq {
        let mut output = DigitSeq::new();
        for (k, fixed) in self.digitseqs.iter().enumerate() {
            output += fixed;
            if k == slot_i {
                output += digit_i;
            }
            if k == slot_j {
                output += digit_j;
            }
        }
        output
    }

    /// Expands a family on the ith core "to the left", meaning, if the family
    /// is xLz, with L = {y1, y2, ...}, this returns the families xz, xy1Lz,
    /// xy2Lz, ....
    ///
    /// This is Lemma 19 in Bright, and in his code, it's called "exploring".
    pub fn expand_left(&self, slot: usize) -> Vec<Self> {
        self.cores[slot]
            .iter()
            // skip 0 if it'd be the first digit
            .filter(|digit| !(digit.0 == 0 && slot == 0 && self.digitseqs[0].0.is_empty()))
            .map(|digit| {
                // Insert the new digit into the fixed part of this segment
                let mut new = self.clone();
                new.digitseqs[slot] += digit;
                new
            })
            .collect()
    }

    /// Expands a family on the ith core "to the right", meaning, if the family
    /// is xLz, with L = {y1, y2, ...}, this returns the families xz, xLy1z,
    /// xLy2Lz, ....
    ///
    /// This is Lemma 19 in Bright, and in his code, it's called "exploring".
    pub fn expand_right(&self, slot: usize) -> Vec<Self> {
        self.cores[slot]
            .iter()
            .map(|digit| {
                // If the invariant is true, we'll definitely have an entry
                // at slot + 1.
                let mut new = self.clone();
                new.digitseqs[slot + 1].0.insert(0, digit);
                new
            })
            .collect()
    }

    /// Returns true if this family contains `needle`.
    pub fn could_contain(&self, needle: &DigitSeq) -> bool {
        let mut needle_iter = needle.0.iter().copied().peekable();

        for (seq, core) in self.digitseqs.iter().zip(&self.cores) {
            // Consume the known digits first...
            for d in &seq.0 {
                match needle_iter.peek() {
                    Some(d2) if d == d2 => {
                        needle_iter.next();
                    }
                    Some(_) => {}
                    None => return true,
                }
            }

            // ...then use the core to consume as much of the needle
            // as possible.
            loop {
                match needle_iter.peek() {
                    Some(d2) if core.digits.contains(d2) => {
                        needle_iter.next();
                    }
                    Some(_) => break,
                    None => return true,
                }
            }
        }

        // Get the last core too
        for d in &self.digitseqs.last().unwrap().0 {
            match needle_iter.peek() {
                Some(d2) if d == d2 => {
                    needle_iter.next();
                }
                Some(_) => {}
                None => return true,
            }
        }

        // could have exhausted it on the very last round, still gotta check
        needle_iter.peek().is_none()
    }
}

impl Core {
    pub fn new(digits: Vec<Digit>) -> Self {
        Self { digits }
    }

    pub fn full(base: u8) -> Self {
        Self {
            digits: (0..base).map(Digit).collect(),
        }
    }

    pub fn remove(&mut self, d: Digit) {
        self.digits.retain(|d2| d != *d2);
    }

    pub fn without(mut self, d: Digit) -> Self {
        self.remove(d);
        self
    }

    pub fn clear(&mut self) {
        self.digits.clear();
    }

    pub fn iter(&self) -> impl Iterator<Item = Digit> + Clone + '_ {
        self.digits.iter().copied()
    }

    pub fn is_empty(&self) -> bool {
        self.digits.is_empty()
    }

    pub fn len(&self) -> usize {
        self.digits.len()
    }
}

impl SimpleFamily {
    /// Returns the smallest member of this family, i.e., `before + center *
    /// min_repeats + after`.
    pub fn contract(&self) -> DigitSeq {
        let mut seq = self.before.clone();
        for _ in 0..self.min_repeats {
            seq += self.center;
        }
        seq += &self.after;
        seq
    }

    /// Returns the value of [Self::contract], interpreted in the given base.
    pub fn value(&self, base: u8) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.before.0 {
            value = value * base + d.0;
        }
        for _ in 0..self.min_repeats {
            value = value * base + self.center.0;
        }
        for d in &self.after.0 {
            value = value * base + d.0;
        }
        value
    }

    /// Returns the smallest n for which this family will contain the given
    /// digit sequence as a substring, or None if no such n exists.
    pub fn will_contain_at(&self, needle: &DigitSeq) -> Option<usize> {
        let mut needle_iter = needle.0.iter().copied().peekable();
        let mut repeats_required = 0;

        // Three stages: go through before, then center, then after.
        // Try to consume the whole needle.
        for d in self.before.0.iter().copied() {
            match needle_iter.peek() {
                Some(d2) if d == *d2 => {
                    needle_iter.next();
                }
                Some(_) => {}
                None => break,
            }
        }

        // For the center, consume as many digits as we can, even if it's
        // more than we currently have.
        loop {
            match needle_iter.peek() {
                Some(d2) if self.center == *d2 => {
                    repeats_required += 1;
                    needle_iter.next();
                }
                // different digit, time to leave
                Some(_) => break,
                // done with the needle!
                None => break,
            }
        }

        for d in self.after.0.iter().copied() {
            match needle_iter.peek() {
                Some(d2) if d == *d2 => {
                    needle_iter.next();
                }
                Some(_) => {}
                None => break,
            }
        }

        if needle_iter.peek().is_some() {
            None
        } else {
            Some(repeats_required)
        }
    }

    #[cfg(test)]
    pub fn pattern(&self) -> String {
        format!("{}{}*{}", self.before, self.center, self.after)
    }
}

impl TryFrom<Family> for SimpleFamily {
    type Error = Family;

    fn try_from(mut family: Family) -> Result<Self, Self::Error> {
        // We need to have exactly one core, and only one digit in it.
        if family.cores.len() != 1 || family.cores[0].len() != 1 {
            return Err(family);
        }

        // Great! Pull out the bits we want.
        let center = family.cores[0].iter().next().unwrap();
        let mut after = family.digitseqs.pop().unwrap();
        let mut before = family.digitseqs.pop().unwrap();

        // It's quite likely we've got some repeated digits next to
        // the center. Let's merge those in.
        let mut num_repeats = 0;
        while before.0.last() == Some(&center) {
            num_repeats += 1;
            before.0.pop();
        }
        while after.0.first() == Some(&center) {
            num_repeats += 1;
            after.0.remove(0);
        }

        Ok(SimpleFamily {
            before,
            center,
            min_repeats: num_repeats,
            after,
        })
    }
}

impl std::fmt::Display for Family {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        debug_assert_eq!(self.digitseqs.len(), self.cores.len() + 1);
        for i in 0..self.cores.len() {
            write!(
                f,
                "{}[{}]*",
                self.digitseqs[i],
                self.cores[i].iter().format("")
            )?
        }
        write!(f, "{}", self.digitseqs.last().expect("digitseqs nonempty"))
    }
}

impl std::fmt::Display for SimpleFamily {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}*{} -- x{}",
            self.before, self.center, self.after, self.min_repeats
        )
    }
}
