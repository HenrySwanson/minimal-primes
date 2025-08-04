use itertools::Itertools;
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;

use crate::digits::{Digit, DigitSeq};

#[derive(Debug, Clone, PartialEq)]
pub struct Family {
    // invariant: digitseqs.len() = cores.len() + 1
    pub digitseqs: Vec<DigitSeq>,
    pub cores: Vec<Core>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Core {
    digits: Vec<Digit>,
}

#[derive(Debug, Clone)]
pub struct SimpleFamily {
    pub before: DigitSeq,
    pub center: Digit,
    pub min_repeats: usize,
    pub after: DigitSeq,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sequence {
    pub k: u64,
    pub c: i64,
    pub d: u64,
}

impl Family {
    pub fn any(base: u8) -> Self {
        Self {
            digitseqs: vec![DigitSeq::new(), DigitSeq::new()],
            cores: vec![Core::full(base)],
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
            if first_seq.0.is_empty()
                && first_core.len() == 1
                && first_core.iter().next().unwrap() == Digit(0)
            {
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

    pub fn split_left(&self, slot: usize) -> Vec<Self> {
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

    pub fn split_right(&self, slot: usize) -> Vec<Self> {
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
    pub fn sequence(&self) -> DigitSeq {
        let mut seq = self.before.clone();
        for _ in 0..self.min_repeats {
            seq += self.center;
        }
        seq += &self.after;
        seq
    }

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

impl Sequence {
    pub fn new(k: u64, c: i64, d: u64) -> Self {
        assert_ne!(k, 0);
        assert_ne!(c, 0);
        assert_ne!(d, 0);
        // k and c have to be opposites mod d
        assert_eq!(k.checked_add_signed(c).unwrap() % d, 0);

        // Do some quick reduction to put it in lowest terms
        let gcd = k.gcd(&c.unsigned_abs()).gcd(&d);

        Self {
            k: k / gcd,
            // casting is okay because 0 < gcd <= |c|
            c: c / (gcd as i64),
            d: d / gcd,
        }
    }

    // TODO: better error type
    pub fn try_from_family(simple: &SimpleFamily, base: u8) -> Result<Self, String> {
        // Compute the sequence for this family: xy*z
        let x = simple.before.value(base);
        let y = simple.center.0;
        let z = simple.after.value(base);

        let b_z = BigUint::from(base).pow(simple.after.0.len() as u32);
        let d = u64::from(base) - 1;
        let k = (x * d + y) * &b_z;
        let c = BigInt::from(d * z) - BigInt::from(y * b_z);

        // Try to fit it into the appropriate ranges
        let k = u64::try_from(k)
            .map_err(|e| format!("Can't convert {} to u64 for {}", e.into_original(), simple))?;
        let c = i64::try_from(c)
            .map_err(|e| format!("Can't convert {} to i64 for {}", e.into_original(), simple))?;

        Ok(Sequence::new(k, c, d))
    }

    pub fn compute_term(&self, n: u32, base: u64) -> BigUint {
        let bn = BigUint::from(base).pow(n);
        let kbnc = if self.c > 0 {
            self.k * bn + self.c.unsigned_abs()
        } else {
            self.k * bn - self.c.unsigned_abs()
        };
        let (q, r) = kbnc.div_rem(&self.d.into());
        debug_assert_eq!(r, BigUint::ZERO);
        q
    }

    pub fn check_term_equal(&self, base: u64, p: u64, n: usize) -> bool {
        let mut x = u128::from(p);
        x *= u128::from(self.d);
        if self.c > 0 {
            let c = u128::from(self.c.unsigned_abs());

            x = match x.checked_sub(c) {
                Some(x) => x,
                None => return false,
            };
        } else {
            x += u128::from(self.c.unsigned_abs());
        }

        if x % u128::from(self.k) != 0 {
            return false;
        }
        x /= u128::from(self.k);

        for _ in 0..n {
            if x % u128::from(base) != 0 {
                return false;
            }
            x /= u128::from(base);
        }

        x == 1
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

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}*b^n+{})/{}", self.k, self.c, self.d)
    }
}
