use std::cell::RefCell;
use std::time::{Duration, Instant};

use itertools::Itertools;
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;

use crate::composite::{find_even_odd_factor, find_perpetual_factor, shares_factor_with_base};
use crate::data::{DigitSeq, Pattern, SimplePattern};
use crate::debug_println;
use crate::math::{big_one, gcd_reduce};

#[derive(Debug, Default)]
pub struct Stats {
    pub num_primality_checks: usize,
    pub duration_primality_checks: Duration,
    pub num_substring_checks: usize,
    pub duration_substring_checks: Duration,
    pub num_simple_substring_checks: usize,
    pub duration_simple_substring_checks: Duration,
    pub num_branches_explored: usize,
}

pub struct SearchContext {
    pub base: u8,

    /// iteration counter; corresponds to the weight of the patterns
    /// we're looking at
    pub iter: usize,
    /// patterns we haven't explored yet
    pub frontier: Frontier,
    /// primes we've discovered so far, in two different formats
    /// depending on the order we discovered these, they may not be minimal!
    pub primes: Vec<(DigitSeq, BigUint)>,

    /// For potentially getting insight into what's going on
    pub stats: RefCell<Stats>,
}

pub struct Frontier {
    /// maps weight to pattern; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    /// items lower than `min_allowed_weight` must be empty
    by_weight: Vec<Vec<SearchNode>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

#[derive(Debug)]
pub enum SearchNode {
    Arbitrary(Pattern),
    Simple(SimplePattern),
}

impl Frontier {
    pub fn new(base: u8) -> Self {
        let initial_pattern = Pattern::any(base);
        debug_assert_eq!(initial_pattern.weight(), 0);
        Self {
            by_weight: vec![vec![SearchNode::Arbitrary(initial_pattern)]],
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

    pub fn iter(&self) -> impl Iterator<Item = &SearchNode> {
        self.by_weight.iter().flatten()
    }

    /// Removes the patterns of the least weight from the structure.
    /// Once this method is called, elements of weight < self.min_weight
    /// should not be inserted! Elements with exactly the minimum weight
    /// are allowed though (lateral exploration).
    pub fn pop(&mut self) -> Option<Vec<SearchNode>> {
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
            let node = SearchNode::Arbitrary(pat);
            self.put(node);
        }
    }

    pub fn put(&mut self, node: SearchNode) {
        let weight = match &node {
            SearchNode::Arbitrary(pattern) => pattern.weight(),
            SearchNode::Simple(simple) => {
                simple.before.0.len() + simple.num_repeats + simple.after.0.len()
            }
        };
        debug_assert!(weight >= self.min_allowed_weight);

        loop {
            match self.by_weight.get_mut(weight) {
                Some(layer) => {
                    layer.push(node);
                    break;
                }
                None => self.by_weight.push(vec![]),
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
            primes: vec![],
            stats: RefCell::new(Stats::default()),
        }
    }

    pub fn minimize_primes(&mut self) {
        // TODO: this seems kind of ad-hoc. is there a better approach?
        let minimized: Vec<_> = self
            .primes
            .iter()
            .filter(|(seq, _)| self.test_for_contained_prime(seq).is_none())
            .cloned()
            .collect();
        self.primes = minimized;
    }

    pub fn search_one_level(&mut self) {
        let old_queue = self.frontier.pop().unwrap_or_default();
        for pattern in old_queue {
            // Say our pattern is xL*z.
            // We want to explore all possible children with weight one more than this one.
            match pattern {
                SearchNode::Arbitrary(pattern) => {
                    debug_println!(" Exploring {}", pattern);
                    self.explore_pattern(pattern)
                }
                SearchNode::Simple(pattern) => {
                    debug_println!(" Exploring simple {}", pattern,);
                    self.explore_simple_pattern(pattern)
                }
            }
            self.stats.borrow_mut().num_branches_explored += 1;
        }
        self.iter += 1;
    }

    fn explore_pattern(&mut self, pattern: Pattern) {
        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = pattern.contract();
        if let Some(p) = self.test_for_contained_prime(&seq) {
            assert_ne!(&seq, p);
            debug_println!("  Discarding {}, contains prime {}", pattern, p);
            return;
        }

        debug_println!("  Testing for primality {}", seq);
        let value = seq.value(self.base);
        if self.test_for_prime(&value) {
            debug_println!("  Saving {}, contracts to prime", pattern);
            self.primes.push((seq, value));
            return;
        }

        // Then, we try to reduce the cores.
        let mut pattern = self.reduce_cores(pattern);
        pattern.simplify();
        if pattern.cores.is_empty() {
            debug_println!("  {} was reduced to trivial string", pattern);
            return;
        }

        // Now, run some tests to see whether this pattern is guaranteed to
        // be composite.
        if self.test_for_perpetual_composite(&pattern) {
            debug_println!("  Discarding {}, is always composite", pattern);
            return;
        }

        // TODO: is this right?
        // Check if this pattern is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        if pattern.cores.len() == 1 && pattern.cores[0].len() == 1 {
            let after = pattern.digitseqs.pop().unwrap();
            let before = pattern.digitseqs.pop().unwrap();
            let node = SearchNode::Simple(SimplePattern {
                before,
                center: pattern.cores[0][0],
                num_repeats: 0,
                after,
            });
            self.frontier.put(node);
            return;
        }

        // Let's see if we can split it in an interesting way
        if let Some(children) = self.split_on_repeat(&pattern, 3) {
            self.frontier.extend(children);
            return;
        }

        if pattern.weight() >= 2 {
            if let Some(child) = self.split_on_necessary_digit(&pattern) {
                self.frontier.extend([child]);
                return;
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

    fn explore_simple_pattern(&mut self, mut pattern: SimplePattern) {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.
        // Even composite testing isn't very useful here, since we
        // can't get here without passing through a composite test.
        // All we can do is add another digit and see if it becomes
        // prime or not. Or contains another prime.

        // Test if it contains a prime
        let start = Instant::now();
        for (prime, _) in &self.primes {
            let result = is_substring_of_simple(&prime, &pattern);
            self.stats.borrow_mut().num_simple_substring_checks += 1;
            if result {
                debug_println!("  Discarding {}, contains prime {}", pattern, prime);
                self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();
                return;
            }
        }
        self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();

        // Test if it is a prime
        let mut value = BigUint::ZERO;
        for d in &pattern.before.0 {
            value = value * self.base + d.0;
        }
        for _ in 0..pattern.num_repeats {
            value = value * self.base + pattern.center.0;
        }
        for d in &pattern.after.0 {
            value = value * self.base + d.0;
        }

        if self.test_for_prime(&value) {
            debug_println!("  Saving {}, is prime", pattern);

            let mut seq = pattern.before.clone();
            for _ in 0..pattern.num_repeats {
                seq += pattern.center;
            }
            seq += pattern.after;

            self.primes.push((seq, value));
            return;
        }

        pattern.num_repeats += 1;
        self.frontier.put(SearchNode::Simple(pattern));
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

                if let Some(p) = self.test_for_contained_prime(&seq) {
                    assert_ne!(&seq, p);
                    debug_println!("  Discarding {}, contains prime {}", seq, p);
                    continue;
                }

                debug_println!("  Testing for primality {}", seq);
                let value = seq.value(self.base);
                if self.test_for_prime(&value) {
                    debug_println!("  Saving {}, is prime", seq);
                    self.primes.push((seq, value));
                } else {
                    allowed_digits.push(digit);
                }
            }

            *core = allowed_digits;
        }
        // Now we've reduced the core, and have a new pattern.
        debug_println!("  Reducing {} to {}", old_pat, pattern);
        pattern
    }

    fn test_for_contained_prime(&self, seq: &DigitSeq) -> Option<&DigitSeq> {
        let start = Instant::now();
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        let result = self.primes.iter().map(|(subseq, _)| subseq).find(|subseq| {
            self.stats.borrow_mut().num_substring_checks += 1;
            is_proper_substring(subseq, seq)
        });

        self.stats.borrow_mut().duration_substring_checks += start.elapsed();
        result
    }

    fn test_for_prime(&mut self, value: &BigUint) -> bool {
        let start = Instant::now();
        let result = is_prime(value, None).probably();
        self.stats.borrow_mut().duration_primality_checks += start.elapsed();
        self.stats.borrow_mut().num_primality_checks += 1;
        result
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
            let g = gcd_reduce(
                // We want to try "no digits" and "all digits except d"
                std::iter::once(contracted.clone()).chain(
                    pattern.cores[0]
                        .iter()
                        .copied()
                        .filter(|d2| *d2 != d)
                        .map(|d2| pattern.substitute(0, d2).value(self.base)),
                ),
            );

            if g != big_one() {
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

fn is_proper_substring(needle: &DigitSeq, haystack: &DigitSeq) -> bool {
    // Save some time when the needle is too large, and also, rule out identical
    // strings.
    if needle.0.len() >= haystack.0.len() {
        return false;
    }

    let mut iter = haystack.0.iter().copied();
    for d in needle.0.iter().copied() {
        // Chomp iter until we find that digit
        loop {
            match iter.next() {
                Some(d2) if d == d2 => break,
                Some(_) => {}
                None => return false,
            }
        }
    }
    // If we got here, then hooray, this is a match!
    true
}

fn is_substring_of_simple(needle: &DigitSeq, haystack: &SimplePattern) -> bool {
    let mut needle_iter = needle.0.iter().copied().peekable();

    // Three stages: go through before, then center, then after.
    // Try to consume the whole needle.
    for d in haystack.before.0.iter().copied() {
        match needle_iter.peek() {
            Some(d2) if d == *d2 => {
                needle_iter.next();
            }
            Some(_) => {}
            None => return true,
        }
    }

    for _ in 0..haystack.num_repeats {
        match needle_iter.peek() {
            Some(d2) if haystack.center == *d2 => {
                needle_iter.next();
            }
            // it's the same digit repeated, if it doesn't match
            // right now, just leave
            Some(_) => break,
            None => return true,
        }
    }

    for d in haystack.after.0.iter().copied() {
        match needle_iter.peek() {
            Some(d2) if d == *d2 => {
                needle_iter.next();
            }
            Some(_) => {}
            None => return true,
        }
    }

    needle_iter.peek().is_none()
}

impl std::fmt::Display for SearchNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SearchNode::Arbitrary(p) => p.fmt(f),
            SearchNode::Simple(p) => p.fmt(f),
        }
    }
}
