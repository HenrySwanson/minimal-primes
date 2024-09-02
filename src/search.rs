use std::cell::RefCell;
use std::time::{Duration, Instant};

use itertools::Itertools;
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use num_traits::One;

use crate::composite::{find_even_odd_factor, find_perpetual_factor, shares_factor_with_base};
use crate::data::{DigitSeq, Family, SimpleFamily};
use crate::debug_println;
use crate::math::gcd_reduce;

pub fn search_for_simple_families(
    base: u8,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    stop_when_simple: bool,
) -> SearchContext {
    let mut ctx = SearchContext::new(base);
    while let Some(weight) = ctx.frontier.min_weight() {
        if let Some(max) = max_weight {
            if weight > max {
                println!("Reached weight cutoff; stopping...");
                break;
            }
        }

        if let Some(max) = max_iter {
            if ctx.iter >= max {
                println!("Reached iteration cutoff; stopping...");
                break;
            }
        }

        if stop_when_simple {
            if ctx.frontier.all_simple() {
                println!("All remaining families are simple; stopping...");
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

    ctx.minimize_primes();
    ctx.primes.sort_by_key(|(_, p)| p.clone());
    ctx
}

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

    /// iteration counter; corresponds to the weight of the families
    /// we're looking at
    pub iter: usize,
    /// families we haven't explored yet
    pub frontier: Frontier,
    /// primes we've discovered so far, in two different formats
    /// depending on the order we discovered these, they may not be minimal!
    pub primes: Vec<(DigitSeq, BigUint)>,

    /// For potentially getting insight into what's going on
    pub stats: RefCell<Stats>,
}

pub struct Frontier {
    /// maps weight to family; used to ensure we're exploring the
    /// search space in (non-strictly) increasing order.
    /// items lower than `min_allowed_weight` must be empty
    by_weight: Vec<Vec<SearchNode>>,
    /// the 'ratchet' that enforces that we can't backtrack to an
    /// element of lower weight.
    min_allowed_weight: usize,
}

#[derive(Debug, Clone)]
pub enum SearchNode {
    Arbitrary(Family),
    Simple(SimpleFamily),
}

impl Frontier {
    pub fn new(base: u8) -> Self {
        let initial_family = Family::any(base);
        debug_assert_eq!(initial_family.weight(), 0);
        Self {
            by_weight: vec![vec![SearchNode::Arbitrary(initial_family)]],
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

    /// Removes the families of the least weight from the structure.
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

    pub fn extend(&mut self, iter: impl IntoIterator<Item = Family>) {
        for family in iter {
            let node = SearchNode::Arbitrary(family);
            self.put(node);
        }
    }

    pub fn put(&mut self, node: SearchNode) {
        let weight = match &node {
            SearchNode::Arbitrary(family) => family.weight(),
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

    pub fn all_simple(&self) -> bool {
        self.iter().all(|f| matches!(f, SearchNode::Simple(_)))
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
        for family in old_queue {
            // Say our family is xL*z.
            // We want to explore all possible children with weight one more than this one.
            match family {
                SearchNode::Arbitrary(family) => {
                    debug_println!(" Exploring {}", family);
                    self.explore_family(family)
                }
                SearchNode::Simple(family) => {
                    debug_println!(" Exploring simple {}", family,);
                    self.explore_simple_family(family)
                }
            }
            self.stats.borrow_mut().num_branches_explored += 1;
        }
        self.iter += 1;
    }

    fn explore_family(&mut self, family: Family) {
        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = family.contract();
        if let Some(p) = self.test_for_contained_prime(&seq) {
            assert_ne!(&seq, p);
            debug_println!("  Discarding {}, contains prime {}", family, p);
            return;
        }

        debug_println!("  Testing for primality {}", seq);
        let value = seq.value(self.base);
        if self.test_for_prime(&value) {
            debug_println!("  Saving {}, contracts to prime", family);
            self.primes.push((seq, value));
            return;
        }

        // Then, we try to reduce the cores.
        let mut family = self.reduce_cores(family);
        family.simplify();
        if family.cores.is_empty() {
            debug_println!("  {} was reduced to trivial string", family);
            return;
        }

        // Now, run some tests to see whether this family is guaranteed to
        // be composite.
        if self.test_for_perpetual_composite(&family) {
            debug_println!("  Discarding {}, is always composite", family);
            return;
        }

        // TODO: is this right?
        // Check if this family is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        let mut family = match SimpleFamily::try_from(family) {
            Ok(simple) => {
                self.frontier.put(SearchNode::Simple(simple));
                return;
            }
            // can't convert, put it back to normal
            Err(f) => f,
        };

        // Let's see if we can split it in an interesting way
        if let Some(children) = self.split_on_repeat(&family, 3) {
            self.frontier.extend(children);
            return;
        }

        if family.weight() >= 2 {
            if let Some(child) = self.split_on_necessary_digit(&family) {
                self.frontier.extend([child]);
                return;
            }
        }

        // If we couldn't eliminate the family, let's split it, left or right.
        // We can't split on a non-empty core, but after we simplify, we shouldn't
        // have to worry about that.
        let slot = self.iter % family.cores.len();
        debug_assert!(!family.cores[slot].is_empty());
        if family.weight() == 1 {
            debug_println!("  Splitting {} left", family);
            self.frontier.extend(family.split_left(slot));
        } else {
            debug_println!("  Splitting {} right", family);
            self.frontier.extend(family.split_right(slot));
        }
        // We also need to consider the case where the chosen core expands to
        // the empty string. However, in the case where there's one core, this
        // is pretty redundant with the work we're doing in reduce_core().
        // For example: if we reduce a[xyz]c, we test the primality of axc, ayc
        // and azc. So after we split, and get ax[xyz]c, there's no need to
        // test ax[]c again.
        if family.cores.len() > 1 {
            family.cores[slot].clear();
            family.simplify();
            self.frontier.extend([family]);
        }
    }

    fn explore_simple_family(&mut self, mut family: SimpleFamily) {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.
        // Even composite testing isn't very useful here, since we
        // can't get here without passing through a composite test.
        // All we can do is add another digit and see if it becomes
        // prime or not. Or contains another prime.

        // Test if it contains a prime
        let start = Instant::now();
        for (prime, _) in &self.primes {
            let result = is_substring_of_simple(prime, &family);
            self.stats.borrow_mut().num_simple_substring_checks += 1;
            if result {
                debug_println!("  Discarding {}, contains prime {}", family, prime);
                self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();
                return;
            }
        }
        self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();

        // Test if it is a prime
        let mut value = BigUint::ZERO;
        for d in &family.before.0 {
            value = value * self.base + d.0;
        }
        for _ in 0..family.num_repeats {
            value = value * self.base + family.center.0;
        }
        for d in &family.after.0 {
            value = value * self.base + d.0;
        }

        if self.test_for_prime(&value) {
            debug_println!("  Saving {}, is prime", family);

            let mut seq = family.before.clone();
            for _ in 0..family.num_repeats {
                seq += family.center;
            }
            seq += family.after;

            self.primes.push((seq, value));
            return;
        }

        family.num_repeats += 1;
        self.frontier.put(SearchNode::Simple(family));
    }

    fn reduce_cores(&mut self, mut family: Family) -> Family {
        let old_family = family.clone();
        for (i, core) in family.cores.iter_mut().enumerate() {
            // Substitute elements from the core into the string to see if any
            // of them contain or are a prime.
            // NOTE: this is where we generate minimal primes of (weight + 1), so
            // next loop, those should all be available.
            let mut allowed_digits = vec![];
            for digit in core.iter().copied() {
                let seq = old_family.substitute(i, digit);

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
        // Now we've reduced the core, and have a new family.
        debug_println!("  Reducing {} to {}", old_family, family);
        family
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

    fn test_for_perpetual_composite(&self, family: &Family) -> bool {
        // This function is used to eliminate families that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity, unless specified otherwise.

        // p divides BASE (e.g., 2, 5)
        if let Some(factor) = shares_factor_with_base(self.base, family) {
            debug_println!("  {} has divisor {}", family, factor);
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect families like 4[6]9 being divisible by 7.
        for stride in 1..=2 {
            if let Some(factors) = find_perpetual_factor(self.base, family, stride) {
                debug_println!(
                    "  {} is divisible by {}",
                    family,
                    factors.iter().format(", ")
                );
                return true;
            }
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(self.base, family) {
            debug_println!(
                "  {} is divisible by either {} or {}",
                family,
                even_factor,
                odd_factor
            );
            return true;
        }

        false
    }

    fn split_on_repeat(&self, family: &Family, max_repeats: usize) -> Option<Vec<Family>> {
        debug_println!(" Trying to split {}", family);
        for (i, core) in family.cores.iter().enumerate() {
            for d in core.iter().copied() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = family.substitute_multiple(i, std::iter::repeat(d).take(n));
                    if let Some(p) = self.test_for_contained_prime(&seq) {
                        assert_ne!(&seq, p);
                        debug_println!("  {} contains a prime {}", seq, p);

                        // Split into n families, x (L-y) (y (L-y))^i z for i in 0..n
                        let yless_core: Vec<_> = core.iter().copied().filter(|x| *x != d).collect();
                        let mut first_child = family.clone();
                        first_child.cores[i] = yless_core.clone();

                        let mut children = vec![first_child];

                        while children.len() < n {
                            let mut new = children.last().unwrap().clone();
                            new.digitseqs.insert(i + 1, d.into());
                            new.cores.insert(i + 1, yless_core.clone());
                            children.push(new);
                        }

                        // Simplify everything (don't do it while we're generating families),
                        // since that'd mess with indices).
                        for child in children.iter_mut() {
                            child.simplify();
                        }

                        debug_println!(
                            "  {} split into {}",
                            family,
                            children.iter().format(" and ")
                        );
                        return Some(children);
                    }
                }
            }
        }
        None
    }

    fn split_on_necessary_digit(&self, family: &Family) -> Option<Family> {
        // There's a case in base 11 (and probably others) where we have
        // just one core, where all the digits except one are even, and so
        // is the rest of the number.
        // This tells me that we are required to have at least one of that digit,
        // or else we'll forever be even.
        // This function detects that situation and splits the family accordingly.
        // TODO: generalize to multiple cores!
        // TODO: generalize to multiple digits?
        // TODO: does this belong in composite? not quite i think

        if family.cores.len() != 1 {
            return None;
        }

        if family.cores[0].len() <= 1 {
            return None;
        }

        let contracted = family.contract().value(self.base);

        for d in family.cores[0].iter().copied() {
            let g = gcd_reduce(
                // We want to try "no digits" and "all digits except d"
                std::iter::once(contracted.clone()).chain(
                    family.cores[0]
                        .iter()
                        .copied()
                        .filter(|d2| *d2 != d)
                        .map(|d2| family.substitute(0, d2).value(self.base)),
                ),
            );

            if g != BigUint::one() {
                // Got a match! Return xLyLz
                let mut new = family.clone();
                let d_less_core = family.cores[0]
                    .iter()
                    .copied()
                    .filter(|d2| *d2 != d)
                    .collect();

                new.digitseqs.insert(1, d.into());
                new.cores.insert(1, d_less_core);
                debug_println!("  {} must have a {}, transforming into {}", family, d, new);
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

fn is_substring_of_simple(needle: &DigitSeq, haystack: &SimpleFamily) -> bool {
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
