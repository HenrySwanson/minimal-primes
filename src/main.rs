use clap::Parser;
use data::{DigitSeq, Pattern};
use itertools::Itertools;
use num_bigint::BigUint;
use std::collections::VecDeque;
use tree_format::Indenter;

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
    // depth_first(args);

    let mut ctx = SearchContext::new(args.base);
    for i in 0..12 {
        println!("------------");
        println!("Level {}", i);
        println!("------------");
        println!(
            "{}",
            ctx.minimal_primes.iter().map(|(seq, _)| seq).format(", ")
        );
        println!("------------");
        for pat in ctx.queue.iter() {
            println!("{}", pat);
        }
        println!("------------");
        ctx.search_one_level();
    }
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

                match find_contained_prime(&seq, self.base) {
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
            if let Some(f) = find_perpetual_factor(&pattern, self.base) {
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
}

fn depth_first(args: Args) {
    let mut minimal_primes = vec![];
    let mut cutoff_patterns = vec![];

    let mut indenter = Indenter::new();

    let mut top = VecDeque::new();
    top.push_back(Pattern::any(args.base));
    let mut stack = vec![top];
    while let Some(top) = stack.last_mut() {
        // Take the front off the vector, if present
        let mut pattern = match top.pop_front() {
            Some(pattern) => pattern,
            None => {
                // Nothing left to investigate; drop down to the next one
                stack.pop();
                indenter.pop();
                continue;
            }
        };

        if top.is_empty() {
            indenter.close();
        }

        indenter.write_line(&format!("{}", pattern));
        indenter.push_new();

        // Cutoff point
        if pattern.weight() > args.cutoff {
            indenter.prepare_pop();
            indenter.write_line("...cutoff!");
            cutoff_patterns.push(pattern);
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

            match find_contained_prime(&seq, args.base) {
                Some(p) => {
                    assert_ne!(seq.0, p.0);
                    indenter.write_line(&format!("{} contains prime: {}", seq, p));
                }
                None => {
                    let value = seq.value(args.base);
                    if is_prime(&value) {
                        indenter.write_line(&format!("{} minimal prime", seq));
                        minimal_primes.push((seq, value));
                    } else {
                        allowed_digits.push(digit);
                    }
                }
            }
        }

        // Now we've reduced the center, and have a new pattern.
        if allowed_digits.is_empty() {
            indenter.prepare_pop();
            indenter.write_line("eliminated all digits");
            continue;
        } else {
            pattern.center = allowed_digits;
        }

        // Let's check whether this pattern is guaranteed to never be prime,
        // i.e., if there's a perpetual factor.
        if let Some(f) = find_perpetual_factor(&pattern, args.base) {
            indenter.prepare_pop();
            indenter.write_line(&format!("{} always divisible by {}!", pattern, f));
            continue;
        }

        // If we couldn't eliminate it, let's split it, left or right.
        // On vibes, let's split left when it's the second digit.
        indenter.close();
        if pattern.weight() == 1 {
            indenter.write_line(&format!("{} reduced, split left", pattern));
            stack.push(VecDeque::from(pattern.split_left()));
        } else {
            indenter.write_line(&format!("{} reduced, split right", pattern));
            stack.push(VecDeque::from(pattern.split_right()));
        }
        indenter.push_new();
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

fn find_contained_prime(seq: &DigitSeq, base: u8) -> Option<DigitSeq> {
    for subseq in seq.0.iter().copied().powerset() {
        // Skip the trivial cases
        if subseq.is_empty() || subseq.len() == seq.0.len() {
            continue;
        }

        // Check if it's prime
        let subseq = DigitSeq(subseq);
        if is_prime(&subseq.value(base)) {
            return Some(subseq);
        }
    }

    None
}

fn find_perpetual_factor(pattern: &Pattern, base: u8) -> Option<BigUint> {
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
        let gcd = gcd(d.0.into(), base.into());
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
    let a = pattern.before.value(base);
    let ac = (pattern.before.clone() + pattern.after.clone()).value(base);
    let mut gcd_accumulated = ac;
    for b in &pattern.center {
        let ab = (pattern.before.clone() + *b).value(base);
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