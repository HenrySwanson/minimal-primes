use data::{DigitSeq, Pattern, BASE};
use itertools::Itertools;
use num_bigint::BigUint;
use std::collections::VecDeque;

mod data;

fn main() {
    depth_first();
}

struct Indenter {
    stack: Vec<usize>,
}

impl Indenter {
    const BAR: &'static str = " │  ";
    const TEE: &'static str = " ├─ ";
    const EMPTY: &'static str = "    ";

    fn new() -> Self {
        Self { stack: vec![] }
    }

    fn leader(&self) -> String {
        let mut s = String::new();
        for (idx, n) in self.stack.iter().copied().enumerate() {
            for _ in 1..n {
                s += Self::EMPTY;
            }
            if idx == self.stack.len() - 1 {
                s += Self::TEE;
            } else {
                s += Self::BAR;
            }
        }
        s
    }

    fn indent(&mut self) {
        self.stack.push(1);
    }

    fn increase_indent(&mut self) {
        if let Some(last) = self.stack.last_mut() {
            *last += 1;
        }
    }

    fn dedent(&mut self) {
        self.stack.pop();
    }
}

fn depth_first() {
    let mut minimal_primes = vec![];
    let mut cutoff_patterns = vec![];

    let mut top = VecDeque::new();
    top.push_back(Pattern::any());

    let mut indenter = Indenter::new();

    let mut stack = vec![top];
    while let Some(top) = stack.last_mut() {
        // Take the front off the vector, if present
        let mut pattern = match top.pop_front() {
            Some(pattern) => pattern,
            None => {
                // Nothing left to investigate; drop down to the next one
                stack.pop();
                indenter.dedent();
                continue;
            }
        };

        println!("{}{}", indenter.leader(), pattern);
        indenter.indent();

        // Cutoff point
        if pattern.weight() > 10 {
            println!("{}...cutoff!", indenter.leader());
            cutoff_patterns.push(pattern);
            indenter.dedent();
            continue;
        }

        // Say our pattern is ABC[PQR]+XYZ.
        // We want to explore all possible children that have weight one more than this one.
        // There's many ways we can do that.

        // First, we try to reduce the center, by seeing whether any of P,Q,R can be
        // substituted into the center.
        // If the resulting sequences contains or is a prime, we can eliminate that
        // branch.
        let mut allowed_digits = vec![];
        for digit in pattern.center.iter().copied() {
            let seq = pattern.clone().substitute(digit);

            match find_contained_prime(&seq) {
                Some(p) => {
                    assert_ne!(seq.0, p.0);
                    println!("{}{} contains prime: {}", indenter.leader(), seq, p);
                }
                None => {
                    let value = seq.value();
                    if is_prime(&value) {
                        println!("{}{} minimal prime!", indenter.leader(), seq);
                        minimal_primes.push((seq, value));
                    } else {
                        allowed_digits.push(digit);
                    }
                }
            }
        }

        // Now we've reduced the center, and have a new pattern.
        if allowed_digits.is_empty() {
            println!("{}eliminated all digits", indenter.leader());
            indenter.dedent();
            continue;
        } else {
            pattern.center = allowed_digits;
        }

        // Let's check whether this pattern is guaranteed to never be prime,
        // i.e., if there's a perpetual factor.
        if let Some(f) = find_perpetual_factor(&pattern) {
            println!(
                "{}{} always divisible by {}!",
                indenter.leader(),
                pattern,
                f
            );
            indenter.dedent();
            continue;
        }

        // If we couldn't eliminate it, let's split it, left or right.
        // On vibes, let's split left when it's the second digit.
        if pattern.weight() == 1 {
            println!("{}{} reduced, split left", indenter.leader(), pattern);
            stack.push(VecDeque::from(pattern.split_left()));
        } else {
            println!("{}{} reduced, split right", indenter.leader(), pattern);
            stack.push(VecDeque::from(pattern.split_right()));
        }
        indenter.increase_indent();
    }

    println!("Stopping now!");
    println!("Remaining branches:");
    for pattern in cutoff_patterns {
        println!("  {}", pattern);
    }
    println!("Minimal primes: ({} total)", minimal_primes.len());
    minimal_primes.sort_by_key(|(_, val)| val.clone());
    for (seq, p) in minimal_primes {
        println!("  {} = ({} in base 10)", seq, p);
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

fn find_perpetual_factor(pattern: &Pattern) -> Option<BigUint> {
    // This function is used to eliminate patterns that will always result
    // in composite numbers, letting us cut off infinite branches of the
    // search space.
    // There are a few possible ways this can happen. We'll use base 10
    // in the comments for familiarity/

    // p divides BASE (e.g., 2, 5)
    // ---------------------------
    // This happens when the last digit is not coprime with the base.
    // If so, the GCD of the last digit and BASE divides the pattern.
    if let Some(d) = pattern.after.0.last() {
        let gcd = gcd(d.0.into(), BASE.into());
        assert!(gcd != BigUint::ZERO);
        if gcd != BigUint::from(1_u32) {
            return Some(gcd);
        }
    }

    // p does not divide BASE (e.g. 7)
    // -------------------------------
    // This is how we detect patterns like 4[6]9 being divisible by 7.
    //
    // If 49, 469, 4669, etc are all divisible by 7, it means that:
    // - 49 is divisible by 7
    // - 4666...6669 is equivalent to 49 mod 7 (they're both 0)
    //
    // This means that 46666...60 is equivalent to 40 mod 7, but since 7
    // and 10 are coprime, that means 466...666 and 4 are equivalent.
    //
    // It suffices to show that 4 and 46 are equivalent! If they are, then
    // (10*x+6) maps them to 46 and 466, making 4, 46, and 466 equivalent,
    // etc.
    //
    // So we need to prove two facts:
    // - 49 is divisible by 7
    // - 4 and 46 are equivalent mod 7, i.e., 46 - 4 is divisible by 7
    //
    // More generally, for a pattern a[b]c:
    // - ac is divisible by p
    // - ab - a is divisible by p
    //
    // Even better, we don't have to try a bunch of different p; what
    // we're actually looking for is whether ac and (ab - a) have some
    // non-trivial common factor, i.e., we can just use the GCD!
    //
    // Lastly, for patterns with more than one digit in the center,
    // we can try all of them individually and lump them into the same GCD.
    // This is what will let us eliminate patterns like 2[369]1.
    let a = pattern.before.value();
    let ac = (pattern.before.clone() + pattern.after.clone()).value();
    let mut gcd_accumulated = ac;
    for b in &pattern.center {
        let ab = (pattern.before.clone() + DigitSeq::single(*b)).value();
        gcd_accumulated = gcd(gcd_accumulated, ab - &a);
    }
    assert!(gcd_accumulated != BigUint::ZERO);
    if gcd_accumulated != BigUint::from(1_u32) {
        Some(gcd_accumulated)
    } else {
        None
    }
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

        if n % &factor == BigUint::ZERO {
            return false;
        }

        factor += 2_u32;
    }
}

fn gcd(mut a: BigUint, mut b: BigUint) -> BigUint {
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }

    while b != BigUint::ZERO {
        (b, a) = (a % &b, b);
    }

    a
}
