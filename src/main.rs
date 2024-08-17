use clap::Parser;
use data::{DigitSeq, Pattern};
use itertools::Itertools;
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use std::sync::atomic::{AtomicBool, Ordering};

mod data;
mod tree_format;

static LOGGING_ENABLED: AtomicBool = AtomicBool::new(false);

macro_rules! debug_println {
    ($($arg:tt)*) => {
        if LOGGING_ENABLED.load(Ordering::Relaxed) {
            println!($($arg)*);
        }
    };
}

#[derive(Parser)]
struct Args {
    /// Base, e.g., decimal, binary, etc.
    #[arg(default_value_t = 10)]
    base: u8,

    /// Stop exploring when patterns get above this weight.
    /// If neither --max-weight or --max-iter is specified,
    /// --max-weight defaults to 10.
    #[arg(long)]
    max_weight: Option<usize>,

    /// Stop exploring after a specific number of iterations.
    /// If neither --max-weight or --max-iter is specified,
    /// --max-weight defaults to 10.
    #[arg(long)]
    max_iter: Option<usize>,

    /// Enable logging
    #[arg(long)]
    log: bool,
}

fn main() {
    let mut args = Args::parse();
    if args.max_iter.is_none() && args.max_weight.is_none() {
        args.max_weight = Some(10);
    }

    LOGGING_ENABLED.store(args.log, Ordering::Relaxed);

    let mut ctx = SearchContext::new(args.base);
    while let Some(weight) = ctx.frontier.min_weight() {
        if let Some(max) = args.max_weight {
            if weight > max {
                println!("Reached weight cutoff; stopping...");
                break;
            }
        }

        if let Some(max) = args.max_iter {
            if ctx.iter >= max {
                println!("Reached iteration cutoff; stopping...");
                break;
            }
        }

        println!(
            "Iteration {} - Weight {} - {} branches",
            ctx.iter,
            weight,
            ctx.frontier.len()
        );

        ctx.search_one_level();
    }

    ctx.minimal_primes.sort_by_key(|(_, p)| p.clone());
    println!("---- BRANCHES REMAINING ----");
    for pat in ctx.frontier.iter() {
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

    /// iteration counter; corresponds to the weight of the patterns
    /// we're looking at
    iter: usize,
    /// patterns we haven't explored yet
    frontier: Frontier,
    /// primes we've discovered so far, in two different formats
    minimal_primes: Vec<(DigitSeq, BigUint)>,
}

struct Frontier {
    /// maps weight to pattern; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    /// items lower than `min_allowed_weight` must be empty
    by_weight: Vec<Vec<Pattern>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

impl Frontier {
    pub fn new(base: u8) -> Self {
        let initial_pattern = Pattern::any(base);
        debug_assert_eq!(initial_pattern.weight(), 0);
        Self {
            by_weight: vec![vec![initial_pattern]],
            min_allowed_weight: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.by_weight.iter().map(|layer| layer.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.by_weight.iter().all(|layer| layer.is_empty())
    }

    /// Returns the minimum weight present in this collection. This can be
    /// different from [self.min_allowed_weight], because that layer might
    /// be empty (or even the layers above).
    pub fn min_weight(&self) -> Option<usize> {
        self.by_weight.iter().position(|layer| !layer.is_empty())
    }

    pub fn iter(&self) -> impl Iterator<Item = &Pattern> {
        self.by_weight.iter().flat_map(|layer| layer)
    }

    /// Removes the patterns of the least weight from the structure.
    /// Once this method is called, elements of weight < self.min_weight
    /// should not be inserted! Elements with exactly the minimum weight
    /// are allowed though (lateral exploration).
    pub fn pop(&mut self) -> Option<Vec<Pattern>> {
        for (i, layer) in self.by_weight.iter_mut().enumerate() {
            if i < self.min_allowed_weight {
                assert!(layer.is_empty())
            }

            // Take the first non-empty layer we see
            if !layer.is_empty() {
                return Some(std::mem::take(layer));
            }
        }
        None
    }

    pub fn extend(&mut self, iter: impl IntoIterator<Item = Pattern>) {
        for pat in iter {
            let weight = pat.weight();
            debug_assert!(weight >= self.min_allowed_weight);
            loop {
                match self.by_weight.get_mut(weight) {
                    Some(layer) => {
                        layer.push(pat);
                        break;
                    }
                    None => self.by_weight.push(vec![]),
                }
            }
        }
    }
}

impl SearchContext {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            iter: 0,
            frontier: Frontier::new(base),
            minimal_primes: vec![],
        }
    }

    pub fn search_one_level(&mut self) {
        let old_queue = self.frontier.pop().unwrap_or_default();
        for pattern in old_queue {
            // Say our pattern is xL*z.
            // We want to explore all possible children with weight one more than this one.
            debug_println!(" Exploring {}", pattern);

            // First, we try to reduce the cores.
            let mut pattern = self.reduce_cores(pattern);
            pattern.simplify();
            if pattern.cores.is_empty() {
                debug_println!("  {} was reduced to trivial string", pattern);
                continue;
            }

            // Now, run some tests to see whether this pattern is guaranteed to
            // be composite.
            if self.test_for_perpetual_composite(&pattern) {
                debug_println!("  Discarding {}, is always composite", pattern);
                continue;
            }

            // Let's see if we can split it in an interesting way
            if let Some(children) = self.split_on_repeat(&pattern) {
                self.frontier.extend(children);
                continue;
            }

            // If we couldn't eliminate the pattern, let's split it, left or right.
            // We can't split on a non-empty core, but after we simplify, we shouldn't
            // have to worry about that.
            let slot = self.iter % pattern.cores.len();
            debug_assert!(!pattern.cores[slot].is_empty());
            if pattern.weight() == 1 {
                debug_println!("  Splitting {} left", pattern);
                self.frontier.extend(pattern.split_left(slot));
            } else {
                debug_println!("  Splitting {} right", pattern);
                self.frontier.extend(pattern.split_right(slot));
            }
        }

        self.iter += 1;
    }

    fn reduce_cores(&mut self, mut pattern: Pattern) -> Pattern {
        let old_pat = pattern.clone();
        for (i, core) in pattern.cores.iter_mut().enumerate() {
            // Substitute elements from the core into the string to see if any
            // of them contain or are a prime.
            // NOTE: this is where we generate minimal primes of (weight + 1), so
            // next loop, those should all be available.
            let mut allowed_digits = vec![];
            for digit in core.iter().copied() {
                let seq = old_pat.substitute(i, digit);

                match self.test_for_contained_prime(&seq) {
                    Some(p) => {
                        assert_ne!(&seq, p);
                        debug_println!("  Discarding {}, contains prime {}", seq, p);
                    }
                    None => {
                        debug_println!("  Testing for primality {}", seq);
                        let value = seq.value(self.base);
                        if is_prime(&value, None).probably() {
                            debug_println!("  Discarding {}, is minimal prime", seq);
                            self.minimal_primes.push((seq, value));
                        } else {
                            allowed_digits.push(digit);
                        }
                    }
                }
            }

            *core = allowed_digits;
        }
        // Now we've reduced the core, and have a new pattern.
        debug_println!("  Reducing {} to {}", old_pat, pattern);
        pattern
    }

    fn test_for_contained_prime(&self, seq: &DigitSeq) -> Option<&DigitSeq> {
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        self.minimal_primes
            .iter()
            .map(|(subseq, _)| subseq)
            .find(|subseq| is_substring(subseq, seq))
    }

    fn test_for_perpetual_composite(&self, pattern: &Pattern) -> bool {
        // This function is used to eliminate patterns that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity, unless specified otherwise.

        // p divides BASE (e.g., 2, 5)
        if let Some(_) = self.shares_factor_with_base(pattern) {
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect patterns like 4[6]9 being divisible by 7.
        for stride in 1..=2 {
            if let Some(_) = self.find_perpetual_factor(pattern, stride) {
                return true;
            }
        }
        false
    }

    fn shares_factor_with_base(&self, pattern: &Pattern) -> Option<BigUint> {
        // Get the last digit of the pattern
        let last_seq = pattern.digitseqs.last().expect("digitseqs nonempty");
        let d = last_seq.0.last()?;

        let gcd = gcd(d.0.into(), self.base.into());
        debug_assert!(gcd != BigUint::ZERO);
        if gcd != BigUint::from(1_u32) {
            debug_println!("  {} has divisor {}", pattern, gcd);
            Some(gcd)
        } else {
            None
        }
    }

    fn find_perpetual_factor(&self, pattern: &Pattern, stride: usize) -> Option<Vec<BigUint>> {
        // If 49, 469, 4669, etc are all divisible by 7, it means that:
        // - 49 is divisible by 7
        // - 4666...6669 is equivalent to 49 mod 7 (they're both 0)
        //
        // For the second fact, it suffices to show that 49 and 469 are equivalent!
        // If they are, then 10*(x-9)+69 maps them to 469 and 4669, making 49,
        // 469, and 4669 equivalent.
        //
        // (It would also suffice to show that 4 and 46 are equivalent mod 7.)
        //
        // So we need to prove two facts:
        // - 49 is divisible by 7
        // = 469 is equivalent to 49 mod 7, i.e., it's divisible by 7
        // - 4 and 46 are equivalent mod 7, i.e., 46 - 4 is divisible by 7
        //
        // More generally, for a pattern a[b]*c:
        // - ac and abc are divisible by p
        //
        // Even better, we don't have to try a bunch of different p; what
        // we're actually looking for is whether ac and abc have some
        // non-trivial common factor, i.e., we can just use the GCD!
        //
        // For patterns with more than one digit in the center,
        // we can try all of them individually and lump them into the same GCD.
        // This is what will let us eliminate patterns like 2[369]1.
        //
        // Lastly, we can extend this idea to situations where ac, abbc, abbbbc,
        // have one divisor, and abc, abbbc, abbbbbc, have a different one.
        let one = BigUint::from(1_u32);

        // TODO: generalize
        if pattern.cores.len() != 1 {
            return None;
        }

        let mut gcds = vec![BigUint::ZERO; stride];
        for i in 0..stride {
            // The smaller of the two sets: xL^iz
            for center in pattern.cores[0]
                .iter()
                .copied()
                .combinations_with_replacement(i)
            {
                let value = DigitSeq::concat_value(
                    [&pattern.digitseqs[0], &center.into(), &pattern.digitseqs[1]],
                    self.base,
                );
                // Update the GCD. If we ever see a 1, it's always going to
                // be that way, so bail out instantly.
                gcds[i] = gcd(gcds[i].clone(), value);
                if gcds[i] == one {
                    return None;
                }
            }
            // The larger of the two sets: xL^(i+stride)z
            for center in pattern.cores[0]
                .iter()
                .copied()
                .combinations_with_replacement(i + stride)
            {
                let value = DigitSeq::concat_value(
                    [&pattern.digitseqs[0], &center.into(), &pattern.digitseqs[1]],
                    self.base,
                );
                gcds[i] = gcd(gcds[i].clone(), value);
                if gcds[i] == one {
                    return None;
                }
            }
        }
        for g in &gcds {
            debug_assert_ne!(*g, one);
        }
        debug_println!("  {} is divisible by {}", pattern, gcds.iter().format(", "));
        Some(gcds)
    }

    fn split_on_repeat(&self, pattern: &Pattern) -> Option<Vec<Pattern>> {
        debug_println!(" Trying to split {}", pattern);
        for (i, core) in pattern.cores.iter().enumerate() {
            for d in core.iter().copied() {
                // Check whether x yy z contains a prime subword
                let seq = pattern.substitute_multiple(i, &[d, d]);
                if let Some(p) = self.test_for_contained_prime(&seq) {
                    assert_ne!(&seq, p);
                    debug_println!("  {} contains a prime {}", seq, p);

                    // Split into two patterns, x(L-y)z and x(L-y)y(L-y)z
                    let yless_core: Vec<_> = core.iter().copied().filter(|x| *x != d).collect();
                    let mut pattern_1 = pattern.clone();
                    pattern_1.cores[i] = yless_core.clone();
                    let mut pattern_2 = pattern_1.clone();
                    pattern_2.digitseqs.insert(i + 1, d.into());
                    pattern_2.cores.insert(i + 1, yless_core);
                    pattern_1.simplify();
                    pattern_2.simplify();

                    debug_println!("  {} split into {} and {}", pattern, pattern_1, pattern_2);
                    return Some(vec![pattern_1, pattern_2]);
                }
            }
        }
        None
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

fn is_substring(needle: &DigitSeq, haystack: &DigitSeq) -> bool {
    let mut iter = haystack.0.iter().copied();
    for d in needle.0.iter().copied() {
        // Chomp iter until we find that digit
        match iter.any(|d2| d == d2) {
            true => {}
            false => return false,
        }
    }
    // If we got here, then hooray, this is a match!
    true
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

            for _ in 0..12 {
                ctx.search_one_level();
            }

            ctx.minimal_primes.sort_by_key(|(_, p)| p.clone());
            let primes = ctx.minimal_primes.iter().map(|(seq, _)| seq).join(", ");
            assert_eq!(primes, output);

            // Except for bases 8 and 9, everything should be solvable
            if base != 8 && base != 9 {
                assert!(
                    ctx.frontier.is_empty(),
                    "Found unexpected branches for base {}: {}",
                    base,
                    ctx.frontier.iter().join("\n")
                );
            }
        }
    }
}
