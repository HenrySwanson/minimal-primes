use std::fmt::Write;

use itertools::Itertools;
use num_bigint::BigUint;
use serde::{Deserialize, Serialize};

/// A digit. The base is not specified and is provided as another parameter
/// in the necessary methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Digit(pub u8);

/// A set of digits. The base is not specified, and is provided as another
/// parameter in the necessary methods (uncommon).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DigitSet {
    mask: u64,
}

/// A sequence of digits. The base is not specified and is provided as
/// another parameter in the necessary methods.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DigitSeq(pub Vec<Digit>);

impl DigitSet {
    /// Returns a core containing no digits.
    pub fn new(digits: impl IntoIterator<Item = Digit>) -> Self {
        let mut core = Self { mask: 0 };
        for d in digits {
            core.insert(d);
        }
        core
    }

    /// Returns a core containing all the digits in the given base, i.e., 0, 1,
    /// ..., `base` - 1.
    ///
    /// Panics if `base >= 64`.
    pub fn full(base: u8) -> Self {
        let Some(mask) = 1_u64.checked_shl(u32::from(base)) else {
            panic!("Base {base} is too large to be represented as a u64 bitmask!")
        };
        // 100...00 -> 011...11
        Self { mask: mask - 1 }
    }

    fn bit(d: Digit) -> u64 {
        debug_assert!(
            u32::from(d.0) < u64::BITS,
            "digit {d} too big for a u64 bitmask"
        );
        1 << d.0
    }

    pub fn insert(&mut self, d: Digit) {
        self.mask |= Self::bit(d);
    }

    pub fn contains(&self, d: Digit) -> bool {
        self.mask & Self::bit(d) != 0
    }

    pub fn remove(&mut self, d: Digit) {
        self.mask &= !Self::bit(d);
    }

    pub fn without(mut self, d: Digit) -> Self {
        self.remove(d);
        self
    }

    pub fn clear(&mut self) {
        self.mask = 0;
    }

    /// Iterates the digits in this core, smallest first.
    pub fn iter(&self) -> impl Iterator<Item = Digit> + Clone {
        let mut mask = self.mask;
        std::iter::from_fn(move || {
            if mask == 0 {
                return None;
            }

            // We could iteratively shift to find the 1, but nowadays, CPUs
            // have nice instructions for this kind of thing, so let's just
            // use one of those instructions, then erase that 1.
            let d = Digit(mask.trailing_zeros() as u8);
            // This trick is still pretty good though, AFAIK. Clear the lowest 1.
            mask &= mask - 1;
            Some(d)
        })
    }

    pub fn is_empty(&self) -> bool {
        self.mask == 0
    }

    pub fn len(&self) -> usize {
        self.mask.count_ones() as usize
    }
}

impl DigitSeq {
    /// Creates an empty sequence of digits.
    pub fn new() -> Self {
        Self(vec![])
    }

    /// Returns the value of this sequence, interpreted in the given base.
    pub fn value(&self, base: u8) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.0 {
            value *= base;
            value += d.0;
        }
        value
    }

    /// Returns the value of the concatenation of these sequences, interpreted
    /// in the given base.
    pub fn concat_value<'a>(seqs: impl IntoIterator<Item = &'a DigitSeq>, base: u8) -> BigUint {
        let mut value = BigUint::ZERO;
        for seq in seqs.into_iter() {
            for d in &seq.0 {
                value *= base;
                value += d.0;
            }
        }
        value
    }

    /// Returns true if `needle` appears as a proper subsequence of this
    /// sequence.
    pub fn properly_contains(&self, needle: &DigitSeq) -> bool {
        // Save some time when the needle is too large, and also, rule out identical
        // strings.
        if needle.0.len() >= self.0.len() {
            return false;
        }

        let mut iter = self.0.iter().copied();
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
}

impl Default for DigitSeq {
    fn default() -> Self {
        Self::new()
    }
}

impl std::ops::Add for DigitSeq {
    type Output = DigitSeq;

    fn add(mut self, mut rhs: Self) -> Self::Output {
        self.0.append(&mut rhs.0);
        self
    }
}

impl std::ops::Add<Digit> for DigitSeq {
    type Output = DigitSeq;

    fn add(mut self, rhs: Digit) -> Self::Output {
        self.0.push(rhs);
        self
    }
}

impl std::ops::Add<DigitSeq> for Digit {
    type Output = DigitSeq;

    fn add(self, mut rhs: DigitSeq) -> Self::Output {
        rhs.0.insert(0, self);
        rhs
    }
}

impl std::ops::AddAssign for DigitSeq {
    fn add_assign(&mut self, mut rhs: Self) {
        self.0.append(&mut rhs.0)
    }
}

impl std::ops::AddAssign<&DigitSeq> for DigitSeq {
    fn add_assign(&mut self, rhs: &Self) {
        self.0.extend(&rhs.0)
    }
}

impl std::ops::AddAssign<Digit> for DigitSeq {
    fn add_assign(&mut self, rhs: Digit) {
        self.0.push(rhs)
    }
}

impl From<Digit> for DigitSeq {
    fn from(d: Digit) -> Self {
        Self(vec![d])
    }
}

impl From<Vec<Digit>> for DigitSeq {
    fn from(digits: Vec<Digit>) -> Self {
        Self(digits)
    }
}

impl Ord for DigitSeq {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First compare length; the longer one is larger
        self.0.len().cmp(&other.0.len()).then_with(|| {
            // If they're equal length, sort lexicographically
            self.0.cmp(&other.0)
        })
    }
}

impl PartialOrd for DigitSeq {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::fmt::Display for Digit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let d = self.0;
        let ch = match d {
            0..=9 => d + b'0',
            10..=35 => d - 10 + b'A',
            _ => {
                return write!(f, "({d})");
            }
        };
        f.write_char(ch as char)
    }
}

impl std::fmt::Display for DigitSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.iter().format(""))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[track_caller]
    fn check_digitset(core: DigitSet, digits: &[u8]) {
        let digits: Vec<_> = digits.iter().copied().map(Digit).collect();
        assert_eq!(core.iter().collect::<Vec<_>>(), digits);
        assert_eq!(core.len(), digits.len());
    }

    #[test]
    fn empty_digitset_has_no_digits() {
        let core = DigitSet::new(vec![]);

        assert_eq!(core.iter().count(), 0);
        assert!(!core.contains(Digit(0)));
        assert!(!core.contains(Digit(63)));
    }

    #[test]
    fn full_digitset_has_every_digit() {
        for base in [2, 10, 31, 36] {
            let core = DigitSet::full(base);
            assert_eq!(core.len(), usize::from(base));
            assert!(core.contains(Digit(0)));
            assert!(core.contains(Digit(base - 1)));
            assert!(!core.contains(Digit(base)));
        }
    }

    #[test]
    fn digitset_round_trips_digits() {
        let digits = vec![Digit(0), Digit(3), Digit(9)];
        let core = DigitSet::new(digits.clone());
        check_digitset(core, &[0, 3, 9]);
        assert!(core.contains(Digit(3)));
        assert!(!core.contains(Digit(4)));
    }

    #[test]
    fn core_insert_and_remove() {
        let digits = vec![Digit(0), Digit(3), Digit(9)];
        let mut core = DigitSet::new(digits.clone());
        check_digitset(core, &[0, 3, 9]);

        // Insert a new digit
        core.insert(Digit(5));
        assert!(core.contains(Digit(5)));
        check_digitset(core, &[0, 3, 5, 9]);

        // Remove a digit
        core.remove(Digit(0));
        assert!(!core.contains(Digit(0)));
        check_digitset(core, &[3, 5, 9]);

        // Removing a digit that's not present is a no-op
        core.remove(Digit(2));
        assert!(!core.contains(Digit(2)));
        check_digitset(core, &[3, 5, 9]);

        // Inserting a digit that's already there is also a no-op
        core.insert(Digit(9));
        assert!(core.contains(Digit(9)));
        check_digitset(core, &[3, 5, 9]);
    }
}
