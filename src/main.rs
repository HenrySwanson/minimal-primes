use clap::Parser;
use composite::{
    big_one, find_even_odd_factor, find_perpetual_factor, gcd, shares_factor_with_base,
};
use data::{DigitSeq, Pattern};
use itertools::Itertools;
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use std::sync::atomic::AtomicBool;

mod composite;
mod data;

static LOGGING_ENABLED: AtomicBool = AtomicBool::new(false);

const DEFAULT_MAX_WEIGHT: usize = 10;

macro_rules! debug_println {
    ($($arg:tt)*) => {
        if crate::LOGGING_ENABLED.load(std::sync::atomic::Ordering::Relaxed) {
            println!($($arg)*);
        }
    };
}

pub(crate) use debug_println;

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
        args.max_weight = Some(DEFAULT_MAX_WEIGHT);
    }

    LOGGING_ENABLED.store(args.log, std::sync::atomic::Ordering::Relaxed);

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
    println!("---- MINIMAL PRIMES ----");
    println!(
        "{}",
        ctx.minimal_primes.iter().map(|(seq, _)| seq).format(", ")
    );
    println!("------------");
    println!(
        "{} primes found, {} branches unresolved",
        ctx.minimal_primes.len(),
        ctx.frontier.len()
    );
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
        self.by_weight.iter().flatten()
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

            // Test this for primality
            // TODO: normally we've tested this already, in reduce_cores,
            // but split_on_repeat can produce strings we've never tested :/
            // What's a better way to avoid this redundancy?
            let seq = pattern.contract();
            match self.test_for_contained_prime(&seq) {
                Some(p) => {
                    assert_ne!(&seq, p);
                    debug_println!("  Discarding {}, contains prime {}", pattern, p);
                }
                None => {
                    debug_println!("  Testing for primality {}", seq);
                    let value = seq.value(self.base);
                    if is_prime(&value, None).probably() {
                        debug_println!("  Discarding {}, contracts to minimal prime", pattern);
                        self.minimal_primes.push((seq, value));
                        continue;
                    }
                }
            }

            // Then, we try to reduce the cores.
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
            if let Some(children) = self.split_on_repeat(&pattern, 3) {
                self.frontier.extend(children);
                continue;
            }

            if pattern.weight() > 6 {
                if let Some(child) = self.split_on_necessary_digit(&pattern) {
                    self.frontier.extend([child]);
                    continue;
                }
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
            // We also need to consider the case where the chosen core expands to
            // the empty string. However, in the case where there's one core, this
            // is pretty redundant with the work we're doing in reduce_core().
            // For example: if we reduce a[xyz]c, we test the primality of axc, ayc
            // and azc. So after we split, and get ax[xyz]c, there's no need to
            // test ax[]c again.
            if pattern.cores.len() > 1 {
                pattern.cores[slot].clear();
                pattern.simplify();
                self.frontier.extend([pattern]);
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
        if let Some(factor) = shares_factor_with_base(self.base, pattern) {
            debug_println!("  {} has divisor {}", pattern, factor);
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect patterns like 4[6]9 being divisible by 7.
        for stride in 1..=2 {
            if let Some(factors) = find_perpetual_factor(self.base, pattern, stride) {
                debug_println!(
                    "  {} is divisible by {}",
                    pattern,
                    factors.iter().format(", ")
                );
                return true;
            }
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(self.base, pattern) {
            debug_println!(
                "  {} is divisible by either {} or {}",
                pattern,
                even_factor,
                odd_factor
            );
            return true;
        }

        false
    }

    fn split_on_repeat(&self, pattern: &Pattern, max_repeats: usize) -> Option<Vec<Pattern>> {
        debug_println!(" Trying to split {}", pattern);
        for (i, core) in pattern.cores.iter().enumerate() {
            for d in core.iter().copied() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = pattern.substitute_multiple(i, std::iter::repeat(d).take(n));
                    if let Some(p) = self.test_for_contained_prime(&seq) {
                        assert_ne!(&seq, p);
                        debug_println!("  {} contains a prime {}", seq, p);

                        // Split into n patterns, x (L-y) (y (L-y))^i z for i in 0..n
                        let yless_core: Vec<_> = core.iter().copied().filter(|x| *x != d).collect();
                        let mut first_child = pattern.clone();
                        first_child.cores[i] = yless_core.clone();

                        let mut children = vec![first_child];

                        while children.len() < n {
                            let mut new = children.last().unwrap().clone();
                            new.digitseqs.insert(i + 1, d.into());
                            new.cores.insert(i + 1, yless_core.clone());
                            children.push(new);
                        }

                        // Simplify everything (don't do it while we're generating patterns),
                        // since that'd mess with indices).
                        for child in children.iter_mut() {
                            child.simplify();
                        }

                        debug_println!(
                            "  {} split into {}",
                            pattern,
                            children.iter().format(" and ")
                        );
                        return Some(children);
                    }
                }
            }
        }
        None
    }

    fn split_on_necessary_digit(&self, pattern: &Pattern) -> Option<Pattern> {
        // There's a case in base 11 (and probably others) where we have
        // just one core, where all the digits except one are even, and so
        // is the rest of the number.
        // This tells me that we are required to have at least one of that digit,
        // or else we'll forever be even.
        // This function detects that situation and splits the pattern accordingly.
        // TODO: generalize to multiple cores!
        // TODO: generalize to multiple digits?
        // TODO: does this belong in composite? not quite i think

        if pattern.cores.len() != 1 {
            return None;
        }

        if pattern.cores[0].len() <= 1 {
            return None;
        }

        let contracted = pattern.contract().value(self.base);

        for d in pattern.cores[0].iter().copied() {
            let mut g = contracted.clone();
            // Try everything except this digit
            for d2 in pattern.cores[0].iter().copied() {
                if d == d2 {
                    continue;
                }
                g = gcd(g, pattern.substitute(0, d2).value(self.base));
            }

            if g > big_one() {
                // Got a match! Return xLyLz
                let mut new = pattern.clone();
                let d_less_core = pattern.cores[0]
                    .iter()
                    .copied()
                    .filter(|d2| *d2 != d)
                    .collect();

                new.digitseqs.insert(1, d.into());
                new.cores.insert(1, d_less_core);
                debug_println!("  {} must have a {}, transforming into {}", pattern, d, new);
                return Some(new);
            }
        }
        None
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
    use std::io;

    use super::*;

    enum Status {
        /// Completely solved; all branches eliminated.
        Complete,
        /// Can get all the minimal primes, but there's some branches
        /// we can't realize are composite.
        StrayBranches,
        /// Not solved yet. Limit the number of iterations.
        Unsolved(usize),
    }

    macro_rules! declare_test_for_base {
        ($name:ident, $base:literal, $status:expr) => {
            #[test]
            fn $name() {
                test_for_base($base, $status);
            }
        };
    }

    declare_test_for_base!(test_base_2, 2, Status::Complete);
    declare_test_for_base!(test_base_3, 3, Status::Complete);
    declare_test_for_base!(test_base_4, 4, Status::Complete);
    declare_test_for_base!(test_base_5, 5, Status::Complete);
    declare_test_for_base!(test_base_6, 6, Status::Complete);
    declare_test_for_base!(test_base_7, 7, Status::Complete);
    declare_test_for_base!(test_base_8, 8, Status::StrayBranches);
    declare_test_for_base!(test_base_9, 9, Status::StrayBranches);
    declare_test_for_base!(test_base_10, 10, Status::Complete);
    declare_test_for_base!(test_base_11, 11, Status::Unsolved(10));
    declare_test_for_base!(test_base_12, 12, Status::Complete);
    declare_test_for_base!(test_base_13, 13, Status::Unsolved(7));
    declare_test_for_base!(test_base_14, 14, Status::Complete);
    declare_test_for_base!(test_base_15, 15, Status::Complete);
    declare_test_for_base!(test_base_16, 16, Status::Unsolved(10));
    declare_test_for_base!(test_base_17, 17, Status::Unsolved(4));
    declare_test_for_base!(test_base_18, 18, Status::Complete);
    declare_test_for_base!(test_base_19, 19, Status::Unsolved(3));
    // 20 is solvable, but requires --release to be practical
    declare_test_for_base!(test_base_20, 20, Status::Unsolved(50));
    declare_test_for_base!(test_base_21, 21, Status::Unsolved(4));
    declare_test_for_base!(test_base_22, 22, Status::Unsolved(5));
    declare_test_for_base!(test_base_23, 23, Status::Unsolved(3));
    declare_test_for_base!(test_base_24, 24, Status::StrayBranches);
    declare_test_for_base!(test_base_25, 25, Status::Unsolved(3));
    declare_test_for_base!(test_base_26, 26, Status::Unsolved(3));
    declare_test_for_base!(test_base_27, 27, Status::Unsolved(3));
    declare_test_for_base!(test_base_28, 28, Status::Unsolved(3));
    declare_test_for_base!(test_base_29, 29, Status::Unsolved(2));
    // 30 is solvable, but takes too long, even with --release
    declare_test_for_base!(test_base_30, 30, Status::Unsolved(50));

    fn test_for_base(base: u8, status: Status) {
        let (max_weight, all_primes_found, no_stray_branches) = match status {
            Status::Complete => {
                // Simulate it for the full duration
                let max_weight = get_max_weight(base);
                (max_weight, true, true)
            }
            Status::StrayBranches => {
                // Simulate it for the full duration
                let max_weight = get_max_weight(base);
                (max_weight, true, false)
            }
            Status::Unsolved(max_weight) => (max_weight, false, false),
        };

        // Calculate as many primes as we ask
        let final_ctx = calculate(base, max_weight);

        // Check that we have the right primes
        compare_primes(&final_ctx, all_primes_found);

        // Check that we've eliminated all branches
        if no_stray_branches {
            assert!(
                final_ctx.frontier.is_empty(),
                "Some branches were not eliminated!\n{}",
                final_ctx.frontier.iter().format("\n")
            );
        } else {
            assert!(
                !final_ctx.frontier.is_empty(),
                "All branches were eliminated, this test should be marked Complete!"
            )
        }
    }

    fn calculate(base: u8, max_weight: usize) -> SearchContext {
        let mut ctx = SearchContext::new(base);
        while let Some(weight) = ctx.frontier.min_weight() {
            if weight > max_weight {
                break;
            }
            ctx.search_one_level();
        }
        ctx.minimal_primes.sort_by_key(|(_, p)| p.clone());
        ctx
    }

    fn iter_ground_truth(base: u8) -> impl Iterator<Item = String> {
        use io::BufRead;
        let file_path = format!("mepn-data/minimal.{}.txt", base);
        let file = std::fs::File::open(file_path).expect("open ground truth file");
        io::BufReader::new(file)
            .lines()
            .map(|line| line.expect("read line"))
    }

    fn get_max_weight(base: u8) -> usize {
        iter_ground_truth(base)
            .map(|s| s.len())
            .max()
            .expect("nonempty ground truth")
    }

    fn compare_primes(ctx: &SearchContext, require_all: bool) {
        let mut truth_iter = iter_ground_truth(ctx.base).peekable();
        let mut iter = ctx
            .minimal_primes
            .iter()
            .map(|(seq, _)| seq.to_string())
            .peekable();

        let mut fail = false;
        loop {
            // We want to consume both iterators until they're both gone.
            //
            // If there's a mismatch, we need to be intelligent about which
            // iterator to increase (the one with the smaller prime.)
            // It turns out to be simpler to combine th
            let cmp = match (truth_iter.peek(), iter.peek()) {
                // If we get two different primes, we only want to advance
                // the iterator with the lesser prime, so that we can re-sync.
                (Some(p_truth), Some(p_got)) => {
                    // First compare them by length, then lexicographically;
                    // otherwise 9 will compare larger than 10.
                    p_truth
                        .len()
                        .cmp(&p_got.len())
                        .then_with(|| p_truth.cmp(&p_got))
                }
                // It turns out to be easier to combine this case
                // with the one above by treating None as the biggest
                // prime.
                (None, Some(_)) => std::cmp::Ordering::Greater,
                // same
                (Some(_), None) => std::cmp::Ordering::Less,
                // both exhausted! break
                (None, None) => break,
            };
            match cmp {
                // truth < actual, so we skipped something we were
                // supposed to see. this might be okay.
                std::cmp::Ordering::Less => {
                    let p = truth_iter.next().unwrap();
                    if require_all {
                        println!("Didn't see expected prime: {}", p);
                        fail = true;
                    }
                }
                // truth > actual, so there's something in actual that
                // shouldn't be there. this is always a failure.
                std::cmp::Ordering::Greater => {
                    let p = iter.next().unwrap();
                    println!("Got extra unexpected prime: {}", p);
                    fail = true;
                }
                // all is well, just increment both iterators
                std::cmp::Ordering::Equal => {
                    iter.next();
                    truth_iter.next();
                }
            }
        }

        assert!(!fail, "Some mismatches between primes");
    }
}
