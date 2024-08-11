use itertools::Itertools;
use num_bigint::BigUint;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Digit(pub u32);

#[derive(Debug, Clone)]
pub struct DigitSeq(pub Vec<Digit>);

#[derive(Debug, Clone)]
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

    pub fn value(&self, base: u32) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.0 {
            value *= base;
            value += d.0;
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

impl Pattern {
    pub fn any(base: u32) -> Self {
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
        self.before + DigitSeq::single(digit) + self.after
    }

    pub fn split_left(&self) -> Vec<Self> {
        self.center
            .iter()
            .copied()
            // skip 0 if it'd be the first digit
            .filter(|digit| digit.0 != 0 || !self.before.0.is_empty())
            .map(|digit| {
                let new_before = self.before.clone() + DigitSeq::single(digit);

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
                let new_after = DigitSeq::single(digit) + self.after.clone();
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
        write!(f, "{:x}", self.0)
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
