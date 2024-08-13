use std::fmt::Write;

use itertools::Itertools;
use num_bigint::BigUint;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Digit(pub u8);

#[derive(Debug, Clone, PartialEq)]
pub struct DigitSeq(pub Vec<Digit>);

#[derive(Debug, Clone, PartialEq)]
pub struct Pattern {
    pub before: DigitSeq,
    pub center: Vec<Digit>,
    pub after: DigitSeq,
}

impl DigitSeq {
    pub fn new() -> Self {
        Self(vec![])
    }

    pub fn single(d: Digit) -> Self {
        Self(vec![d])
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

impl From<Digit> for DigitSeq {
    fn from(d: Digit) -> Self {
        Self(vec![d])
    }
}

impl Pattern {
    pub fn any(base: u8) -> Self {
        Self {
            before: DigitSeq::new(),
            center: (0..base).map(Digit).collect(),
            after: DigitSeq::new(),
        }
    }

    pub fn weight(&self) -> usize {
        self.before.0.len() + self.after.0.len()
    }

    pub fn substitute(self, digit: Digit) -> DigitSeq {
        self.before + digit + self.after
    }

    pub fn split_left(&self) -> Vec<Self> {
        self.center
            .iter()
            .copied()
            // skip 0 if it'd be the first digit
            .filter(|digit| digit.0 != 0 || !self.before.0.is_empty())
            .map(|digit| {
                let new_before = self.before.clone() + (digit);

                Self {
                    before: new_before,
                    center: self.center.clone(),
                    after: self.after.clone(),
                }
            })
            .collect()
    }

    pub fn split_right(&self) -> Vec<Self> {
        self.center
            .iter()
            .copied()
            .map(|digit| {
                let new_after = (digit) + self.after.clone();
                Self {
                    before: self.before.clone(),
                    center: self.center.clone(),
                    after: new_after,
                }
            })
            .collect()
    }
}

impl std::fmt::Display for Digit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let d = self.0;
        let ch = match d {
            0..=9 => d + b'0',
            10..=35 => d - 10 + b'A',
            _ => {
                return write!(f, "({})", d);
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

impl std::fmt::Display for Pattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}[{}]+{}",
            self.before,
            self.center.iter().format(""),
            self.after
        )
    }
}
