use itertools::Itertools;
use num_bigint::BigUint;
use std::collections::VecDeque;

const BASE: u32 = 10;

fn main() {
    depth_first();
}

fn depth_first() {
    let mut minimal_primes = vec![];
    let mut cutoff_patterns = vec![];

    let mut top = VecDeque::new();
    top.push_back(Pattern::any());

    let mut stack = vec![top];
    while let Some(top) = stack.last_mut() {
        // Take the front off the vector, if present
        let mut pattern = match top.pop_front() {
            Some(pattern) => pattern,
            None => {
                // Nothing left to investigate; drop down to the next one
                stack.pop();
                continue;
            }
        };

        // For formatting
        let leader = " ".repeat(pattern.weight());
        println!("{}{}", leader, pattern);

        // Cutoff point
        if pattern.weight() > 10 {
            println!("{} ...cutoff!", leader);
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

            match find_contained_prime(&seq) {
                Some(p) => {
                    assert_ne!(seq.0, p.0);
                    println!(" {}{} contains prime: {}", leader, seq, p);
                }
                None => {
                    let value = seq.value();
                    if is_prime(&value) {
                        println!(" {}{} minimal prime!", leader, seq);
                        minimal_primes.push(value);
                    } else {
                        allowed_digits.push(digit);
                    }
                }
            }
        }

        // Now we've reduced the center, and have a new pattern.
        if allowed_digits.is_empty() {
            println!(" {}eliminated all digits", leader);
            continue;
        } else {
            pattern.center = allowed_digits;
        }

        // Let's check whether this pattern is guaranteed to never be prime,
        // i.e., if there's a perpetual factor.
        match find_perpetual_factor(&pattern) {
            Some(f) => {
                println!(" {}{} always divisible by {}!", leader, pattern, f);
                continue;
            }
            None => println!(" {}{}", leader, pattern),
        }

        // If we couldn't eliminate it, let's split it, left or right.
        // On vibes, let's split left when it's the second digit.
        if pattern.weight() == 1 {
            stack.push(VecDeque::from(pattern.split_left()));
        } else {
            stack.push(VecDeque::from(pattern.split_right()));
        }
    }

    println!("Stopping now!");
    println!("Remaining branches:");
    for pattern in cutoff_patterns {
        println!("  {}", pattern);
    }
    println!("Minimal primes: ({} total)", minimal_primes.len());
    minimal_primes.sort();
    for p in minimal_primes {
        println!("  {}", p);
    }
}

fn breadth_first() {
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

    // TODO: turn this into a real "neverprime" machine.
    // there are a few kinds of divisibility tests, all of which are necessary
    // * last-digit: 2,5
    //   * only the last digit matters
    //   * kills branches like [468]+4 that'll never contain primes
    //   * only works for primes dividing BASE
    // * mod-sum: 3
    //   * hyper-specific, only affects factors of BASE-1
    //   * can deal with things like 8[0]+1
    //   * can work on large centers, e.g., [069]+6669, or 1[0369]+2.
    // * tricky: lots of primes?
    //   * this is the evil case that is necessary but idk exactly how to do it yet
    //   * only way to kill 4[6]+9
    //   * mostly works on singleton centers, but i bet there's a rare case where
    //     you get things equivalent mod p in the center
    //   * approach: 4 is 4 mod 7, and when you tack on a 6, you get 4*10+6=5+6=4 again.
    //     this tells you you can collapse the 6s. go check 49.
    //   * the repdigit depends on the beginning; 3, 36, 366, don't have the same remainder
    //     mod 7, but 3, 31, 311, do. (so 3[1]+5 would be never prime)
    //   * i think this works with arbitrarily large primes. consider 4[3]+42 mod 13.
    //     4*10+3=1+3=4, so check 442, which is 2*13*17
    //   * maybe i can tackle it without knowing p?
    //     4*10+6 = 4 for what p? it's factors of 4*9+6 i.e. 42 (hey there's 7!)
    //     4*10+3 = 4 for factors of 4*9+3 = 39 (hey, there's 13!)
    //   * maybe i can roll mod-3 into that then...
    //     for 5[0369]28, 5*9+x, where x=0369, gives 45,48,51,54, all of which are
    //     multiples of 3. so we can collapse them and check 528 mod 3, which is 0.
    //   * i think that clinches it! i don't have to try arbitrarily high primes, i just
    //     have to track factors of (before)*(BASE-1)+(center) :)
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
