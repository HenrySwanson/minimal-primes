
use itertools::Itertools;
use num_bigint::BigUint;

use crate::digits::{Digit, DigitSeq};

#[derive(Debug, Clone, PartialEq)]
pub struct Family {
    // invariant: digitseqs.len() = cores.len() + 1
    pub digitseqs: Vec<DigitSeq>,
    pub cores: Vec<Vec<Digit>>,
}

#[derive(Debug, Clone)]
pub struct SimpleFamily {
    pub before: DigitSeq,
    pub center: Digit,
    pub num_repeats: usize,
    pub after: DigitSeq,
}

impl Family {
    pub fn any(base: u8) -> Self {
        Self {
            digitseqs: vec![DigitSeq::new(), DigitSeq::new()],
            cores: vec![(0..base).map(Digit).collect()],
        }
    }

    pub fn weight(&self) -> usize {
        self.digitseqs.iter().map(|seq| seq.0.len()).sum()
    }

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
            if first_seq.0.is_empty() && first_core.len() == 1 && first_core[0] == Digit(0) {
                self.digitseqs.remove(0);
                self.cores.remove(0);
            }
        }
    }

    pub fn contract(&self) -> DigitSeq {
        DigitSeq(
            self.digitseqs
                .iter()
                .flat_map(|seq| &seq.0)
                .copied()
                .collect(),
        )
    }

    pub fn substitute(&self, slot: usize, digit: Digit) -> DigitSeq {
        self.substitute_multiple(slot, [digit])
    }

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

    pub fn split_left(&self, slot: usize) -> Vec<Self> {
        self.cores[slot]
            .iter()
            .copied()
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

    pub fn split_right(&self, slot: usize) -> Vec<Self> {
        self.cores[slot]
            .iter()
            .copied()
            .map(|digit| {
                // If the invariant is true, we'll definitely have an entry
                // at slot + 1.
                let mut new = self.clone();
                new.digitseqs[slot + 1].0.insert(0, digit);
                new
            })
            .collect()
    }
}

impl SimpleFamily {
    pub fn sequence_with(&self, n: usize) -> DigitSeq {
        let mut seq = self.before.clone();
        for _ in 0..n {
            seq += self.center;
        }
        seq += &self.after;
        seq
    }

    pub fn value(&self, base: u8) -> BigUint {
        self.value_with(base, self.num_repeats)
    }

    pub fn value_with(&self, base: u8, n: usize) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.before.0 {
            value = value * base + d.0;
        }
        for _ in 0..n {
            value = value * base + self.center.0;
        }
        for d in &self.after.0 {
            value = value * base + d.0;
        }
        value
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
        let center = family.cores[0][0];
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
            num_repeats,
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
            self.before, self.center, self.after, self.num_repeats
        )
    }
}
