use clap::Parser;
use data::{DigitSeq, Pattern};
use itertools::Itertools;
use num_bigint::BigUint;
use std::collections::VecDeque;

mod data;
mod tree_format;

#[derive(Parser)]
struct Args {
    /// Base, e.g., decimal, binary, etc.
    #[arg(default_value_t = 10)]
    base: u8,

    /// How far to explore into a branch before we give up
    #[arg(short, long, default_value_t = 10)]
    cutoff: usize,
}

fn main() {
    let args = Args::parse();

    let mut ctx = SearchContext::new(args.base);
    for i in 0..args.cutoff {
        println!("Level {} - {} branches", i, ctx.queue.len());
        ctx.search_one_level();
    }

    ctx.minimal_primes.sort_by_key(|(_, p)| p.clone());
    println!("---- BRANCHES REMAINING ----");
    for pat in ctx.queue.iter() {
        println!("{}", pat);
    }
    println!("---- MINIMAL PRIMES ({}) ----", ctx.minimal_primes.len());
    println!(
        "{}",
        ctx.minimal_primes.iter().map(|(seq, _)| seq).format(", ")
    );
    println!("------------");
}

struct SearchContext {
    base: u8,

    minimal_primes: Vec<(DigitSeq, BigUint)>,
    queue: VecDeque<Pattern>,
}

impl SearchContext {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            minimal_primes: vec![],
            queue: vec![Pattern::any(base)].into(),
        }
    }

    pub fn search_one_level(&mut self) {
        let old_queue = std::mem::take(&mut self.queue);
        for mut pattern in old_queue {
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

                match self.test_for_contained_prime(&seq) {
                    Some(p) => {
                        assert_ne!(seq.0, p.0);
                    }
                    None => {
                        let value = seq.value(self.base);
                        if is_prime(&value) {
                            self.minimal_primes.push((seq, value));
                        } else {
                            allowed_digits.push(digit);
                        }
                    }
                }
            }

            // Now we've reduced the center, and have a new pattern.
            if allowed_digits.is_empty() {
                continue;
            } else {
                pattern.center = allowed_digits;
            }

            // Let's check whether this pattern is guaranteed to never be prime,
            // i.e., if there's a perpetual factor.
            if let Some(f) = self.find_perpetual_factor(&pattern) {
                continue;
            }

            // If we couldn't eliminate it, let's split it, left or right.
            if pattern.weight() == 1 {
                self.queue.extend(VecDeque::from(pattern.split_left()));
            } else {
                self.queue.extend(VecDeque::from(pattern.split_right()));
            }
        }
    }

    fn test_for_contained_prime(&self, seq: &DigitSeq) -> Option<DigitSeq> {
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        'outer: for (subseq, _) in self.minimal_primes.iter() {
            let mut iter = seq.0.iter().copied();
            for d in subseq.0.iter().copied() {
                // Chomp iter until we find that digit
                match iter.any(|d2| d == d2) {
                    true => {}
                    false => continue 'outer,
                }
            }

            // If we got here, then hooray, this is a match!
            return Some(subseq.clone());
        }
        None
    }

    fn find_perpetual_factor(&self, pattern: &Pattern) -> Option<BigUint> {
        // This function is used to eliminate patterns that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity.

        // p divides BASE (e.g., 2, 5)
        // ---------------------------
        // This happens when the last digit is not coprime with the base.
        // If so, the GCD of the last digit and BASE divides the pattern.
        if let Some(d) = pattern.after.0.last() {
            let gcd = gcd(d.0.into(), self.base.into());
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
        let a = pattern.before.value(self.base);
        let ac = (pattern.before.clone() + pattern.after.clone()).value(self.base);
        let mut gcd_accumulated = ac.clone();
        for b in &pattern.center {
            let ab = (pattern.before.clone() + *b).value(self.base);
            gcd_accumulated = gcd(gcd_accumulated, ab - &a);
        }
        assert!(gcd_accumulated != BigUint::ZERO);
        if gcd_accumulated != BigUint::from(1_u32) {
            return Some(gcd_accumulated);
        }

        // Okay, let's try it again but checking two possible divisors; one for
        // ac, abbc, abbbbc, and one for abc, abbbc, etc.
        let mut gcd_accumulated_1 = ac;
        let mut gcd_accumulated_2 = BigUint::ZERO;
        for b in &pattern.center {
            let abb = (pattern.before.clone() + *b + *b).value(self.base);
            let abc = pattern.clone().substitute(*b).value(self.base);
            let abb_minus_a = abb - &a;
            gcd_accumulated_1 = gcd(gcd_accumulated_1, abb_minus_a.clone());
            gcd_accumulated_2 = gcd(gcd_accumulated_2, abb_minus_a);
            gcd_accumulated_2 = gcd(gcd_accumulated_2, abc);
        }
        assert!(gcd_accumulated_1 != BigUint::ZERO);
        assert!(gcd_accumulated_2 != BigUint::ZERO);
        if gcd_accumulated_1 != BigUint::from(1_u32) && gcd_accumulated_2 != BigUint::from(1_u32) {
            return Some(gcd_accumulated_1);
        }

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
    // Since we're working with a binary representation, it's more efficient
    // to use this algorithm than the standard Euclidean one.

    // First, pull out all the factors of 2.
    // gcd(2^i a, 2^j b) = 2^k gcd(a, b) when a, b odd and k = min(i, j);
    // But since x.trailing_zeros() is None when x is 0, take this opportunity
    // to check for the edge case of a = 0 or b = 0.
    let i = match a.trailing_zeros() {
        Some(i) => i,
        None => return b,
    };
    let j = match b.trailing_zeros() {
        Some(j) => j,
        None => return a,
    };

    // How many of those 2s do we want to keep?
    let k = i.min(j);
    a >>= i;
    b >>= j;

    loop {
        // Now they're both odd
        debug_assert!(a.bit(0));
        debug_assert!(b.bit(0));

        // Swap so a is larger
        if a < b {
            std::mem::swap(&mut a, &mut b);
        }

        // Subtract, just like with Euclid
        // gcd(u, v) = gcd(u, v-u)
        a -= &b;

        // Now a is even; remove its 2s (again checking for the case of a = 0).
        // gcd(2^i a, b) = gcd(a, b) when b is odd
        match a.trailing_zeros() {
            Some(i) => {
                a >>= i;
            }
            None => {
                // gcd(0, b) = b, then add in the 2^k we removed earlier
                return b << k;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KNOWN_MINIMAL_SETS: [(u8, &str); 9]= [
        (2, "10, 11"),
        (3, "2, 10, 111"),
        (4, "2, 3, 11"),
        (5, "2, 3, 10, 111, 401, 414, 14444, 44441"),
        (6, "2, 3, 5, 11, 4401, 4441, 40041"),
        (7, "2, 3, 5, 10, 14, 16, 41, 61, 11111"),
        (8, "2, 3, 5, 7, 111, 141, 161, 401, 661, 4611, 6101, 6441, 60411, 444641, 444444441"),
        (9, "2, 3, 5, 7, 14, 18, 41, 81, 601, 661, 1011, 1101"),
        (10, "2, 3, 5, 7, 11, 19, 41, 61, 89, 409, 449, 499, 881, 991, 6469, 6949, 9001, 9049, 9649, 9949, 60649, 666649, 946669, 60000049, 66000049, 66600049"),
    ];

    #[test]
    fn test_known_bases() {
        for (base, output) in KNOWN_MINIMAL_SETS {
            let mut ctx = SearchContext::new(base);

            for _ in 0..10 {
                ctx.search_one_level();
            }

            ctx.minimal_primes.sort_by_key(|(_, p)| p.clone());
            let primes = ctx.minimal_primes.iter().map(|(seq, _)| seq).join(", ");
            assert_eq!(primes, output);

            // Except for bases 8 and 9, everything should be solvable
            if base != 8 && base != 9 {
                assert!(
                    ctx.queue.is_empty(),
                    "Found unexpected branches for base {}: {}",
                    base,
                    ctx.queue.iter().join("\n")
                );
            }
        }
    }
}
