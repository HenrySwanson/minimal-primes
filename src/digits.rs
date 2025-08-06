use std::fmt::Write;

use itertools::Itertools;
use num_bigint::BigUint;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Digit(pub u8);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DigitSeq(pub Vec<Digit>);

impl DigitSeq {
    pub fn new() -> Self {
        Self(vec![])
    }

    pub fn value(&self, base: u8) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.0 {
            value *= base;
            value += d.0;
        }
        value
    }

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
