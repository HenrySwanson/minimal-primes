use std::fmt::Write;

use itertools::Itertools;
use num_bigint::BigUint;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Digit(pub u8);

#[derive(Debug, Clone, PartialEq)]
pub struct DigitSeq(pub Vec<Digit>);

#[derive(Debug, Clone, PartialEq)]
pub struct Pattern {
    pub segments: Vec<Segment>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub fixed: DigitSeq,
    pub core: Vec<Digit>,
}

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

impl Pattern {
    pub fn any(base: u8) -> Self {
        Self {
            segments: vec![Segment {
                fixed: DigitSeq::new(),
                core: (0..base).map(Digit).collect(),
            }],
        }
    }

    pub fn weight(&self) -> usize {
        self.segments.iter().map(|seg| seg.fixed.0.len()).sum()
    }

    pub fn substitute(self, slot: usize, digit: Digit) -> DigitSeq {
        let mut output = DigitSeq::new();
        for (i, seg) in self.segments.into_iter().enumerate() {
            output += seg.fixed;
            if i == slot {
                debug_assert!(seg.core.contains(&digit));
                output += digit;
            }
        }
        output
    }

    pub fn split_left(&self, slot: usize) -> Vec<Self> {
        let segment = &self.segments[slot];
        segment
            .core
            .iter()
            .copied()
            // skip 0 if it'd be the first digit
            .filter(|digit| !(digit.0 == 0 && slot == 0 && segment.fixed.0.is_empty()))
            .map(|digit| {
                // Insert the new digit into the fixed part of this segment
                let mut new = self.clone();
                new.segments[slot].fixed += digit;
                new
            })
            .collect()
    }

    pub fn split_right(&self, slot: usize) -> Vec<Self> {
        let segment = &self.segments[slot];
        segment
            .core
            .iter()
            .copied()
            .map(|digit| {
                // Copy and insert another digit to the right of this pattern.
                // This may involve creating a new segment.
                let mut new = self.clone();
                match new.segments.get_mut(slot + 1) {
                    Some(seg) => seg.fixed.0.insert(0, digit),
                    None => new.segments.push(Segment {
                        fixed: digit.into(),
                        core: vec![],
                    }),
                }
                new
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
        for seg in &self.segments {
            write!(f, "{}[{}]*", seg.fixed, seg.core.iter().format(""))?
        }
        Ok(())
    }
}
