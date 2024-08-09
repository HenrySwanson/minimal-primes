use itertools::Itertools;
use num_bigint::BigUint;
use std::collections::VecDeque;

const BASE: u32 = 10;

fn main() {
    println!("Hello, world!");

    let init = Pattern::any();
    let mut branches = VecDeque::new();
    branches.push_back(init);

    let mut minimal_primes = vec![];

    let mut counter = 0;
    while let Some(mut pattern) = branches.pop_front() {
        counter += 1;

        if counter > 80 {
            branches.push_back(pattern);
            break;
        }

        // Say our pattern is ABC[PQR]+XYZ.
        // We want to explore all possible children that have weight one more than this one.
        // There's many ways we can do that.
        println!("Investigating {}", pattern);

        // Hack! Fix this
        if hacky_test(&pattern) {
            println!("Delet: {}", pattern);
            continue;
        }

        // Can any of PQR be eliminated? Also, take care of ABC[PQR]XYZ (no center).
        let mut allowed_digits = vec![];
        for digit in pattern.center.iter().copied() {
            let seq = pattern.clone().substitute(digit);

            match find_contained_prime(&seq) {
                Some(p) => {
                    assert_ne!(seq.0, p.0);
                    println!("Discarding {}, contains primes {}", seq, p);
                }
                None => {
                    let value = seq.value();
                    if is_prime(&value) {
                        println!("Found a minimal prime {}", value);
                        minimal_primes.push(value);
                    } else {
                        allowed_digits.push(digit);
                    }
                }
            }
        }

        // Here's the new pattern
        pattern.center = allowed_digits;
        println!("Reduced pattern: {}", pattern);

        // Now that we've reduced the center, let's split it, left or right.
        // On vibes, let's split left when it's the second digit.
        if pattern.weight() == 1 {
            branches.extend(pattern.split_left());
        } else {
            branches.extend(pattern.split_right());
        }
    }

    println!("Stopping now!");
    println!("Remaining branches:");
    for branch in branches {
        println!("  {}", branch);
    }
    println!("Minimal primes: ({} total)", minimal_primes.len());
    minimal_primes.sort();
    for p in minimal_primes {
        println!("  {}", p);
    }
}

fn hacky_test(pattern: &Pattern) -> bool {
    // Divisible by 2
    if pattern.after.0.last().is_some_and(|d| d.0 % 2 == 0) {
        return true;
    }
    // Nothing but 0s left
    if pattern.before.0.is_empty() && pattern.center == vec![Digit(0)] {
        return true;
    }
    {
        // Divisible by 3
        let mut three_sum_const = 0;
        for digit in &pattern.before.0 {
            three_sum_const += digit.0;
        }
        for digit in &pattern.after.0 {
            three_sum_const += digit.0;
        }
        if three_sum_const % 3 != 0 {
            return false;
        }
        for digit in &pattern.center {
            if digit.0 % 3 != 0 {
                return false;
            }
        }
        return true;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Digit(u32);

#[derive(Debug, Clone)]
struct DigitSeq(Vec<Digit>);

#[derive(Debug, Clone)]
struct Pattern {
    before: DigitSeq,
    center: Vec<Digit>,
    after: DigitSeq,
}

impl DigitSeq {
    fn new() -> Self {
        Self(vec![])
    }

    fn single(d: Digit) -> Self {
        Self(vec![d])
    }

    fn value(&self) -> BigUint {
        let mut value = BigUint::ZERO;
        for d in &self.0 {
            value *= BASE;
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
    fn any() -> Self {
        Self {
            before: DigitSeq::new(),
            center: all_digits(),
            after: DigitSeq::new(),
        }
    }

    fn weight(&self) -> usize {
        self.before.0.len() + self.after.0.len()
    }

    fn substitute(self, digit: Digit) -> DigitSeq {
        self.before + DigitSeq::single(digit) + self.after
    }

    fn split_left(&self) -> Vec<Self> {
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

    fn split_right(&self) -> Vec<Self> {
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

fn find_contained_prime(seq: &DigitSeq) -> Option<DigitSeq> {
    for subseq in seq.0.iter().copied().powerset() {
        // Skip the trivial cases
        if subseq.is_empty() || subseq.len() == seq.0.len() {
            continue;
        }

        // Check if it's prime
        let subseq = DigitSeq(subseq);
        if is_prime(&subseq.value()) {
            return Some(subseq);
        }
    }

    None
}

fn is_prime(n: &BigUint) -> bool {
    let two: BigUint = BigUint::from(2_u32);

    // 0 and 1 are not prime, but 2 is
    match n.cmp(&two) {
        std::cmp::Ordering::Less => return false,
        std::cmp::Ordering::Equal => return true,
        std::cmp::Ordering::Greater => {}
    }

    // Check divisibility by 2
    if !n.bit(0) {
        return false;
    }

    // Now go up by 3, 5, 7, ...
    let mut factor = BigUint::from(3_u32);
    loop {
        if &factor * &factor > *n {
            return true;
        }

        if n.clone() % &factor == BigUint::ZERO {
            return false;
        }

        factor += 2_u32;
    }
}

fn all_digits() -> Vec<Digit> {
    (0..BASE).map(Digit).collect()
}

impl std::fmt::Display for Digit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
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
