mod composite;
mod explore;
mod families;

use std::cell::RefCell;
use std::time::{Duration, Instant};

use itertools::Itertools;
use log::{debug, info, trace};
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use num_traits::One;

use self::composite::{find_even_odd_factor, find_perpetual_factor, shares_factor_with_base};
use self::explore::{Frontier, TreeTracer, Weight};
use crate::data_structures::{is_proper_substring, CandidateSequences};
use crate::digits::DigitSeq;
use crate::math::gcd_reduce;
use crate::search::explore::Explore;
use crate::SearchResults;

pub use self::families::{Family, SimpleFamily};

pub fn search_for_simple_families(
    base: u8,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    stop_when_simple: bool,
    tree_log: bool,
) -> SearchResults {
    if tree_log {
        search_for_simple_families_impl::<TreeTracer<SearchNode>>(
            base,
            max_weight,
            max_iter,
            stop_when_simple,
        )
    } else {
        search_for_simple_families_impl::<Frontier<SearchNode>>(
            base,
            max_weight,
            max_iter,
            stop_when_simple,
        )
    }
}

fn search_for_simple_families_impl<E: Explore>(
    base: u8,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    stop_when_simple: bool,
) -> SearchResults {
    let mut ctx = SearchContext::new(base);
    let mut explorer = E::start(SearchNode::Arbitrary(Family::any(base)));

    while let Some(weight) = explorer.min_weight() {
        if let Some(max) = max_weight {
            if weight > max {
                info!("Reached weight cutoff; stopping...");
                break;
            }
        }

        if let Some(max) = max_iter {
            if ctx.iter >= max {
                info!("Reached iteration cutoff; stopping...");
                break;
            }
        }

        if stop_when_simple && explorer.all_simple() {
            info!("All remaining families are simple; stopping...");
            break;
        }

        info!(
            "Iteration {} - Weight {} - {} branches",
            ctx.iter,
            weight,
            explorer.len()
        );

        explorer.explore_next(|node| ctx.explore_node(node));
        ctx.iter += 1;
    }

    ctx.primes.sort();
    explorer.print_tree_to_stdout();

    // Pull the unsolved branches and return them
    let mut ret = SearchResults {
        primes: ctx.primes,
        simple_families: vec![],
        other_families: vec![],
        stats: ctx.stats,
    };
    for family in explorer.iter().cloned() {
        match family {
            SearchNode::Arbitrary(family) => ret.other_families.push(family),
            SearchNode::Simple(simple_family) => ret.simple_families.push(simple_family),
        }
    }

    ret
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

    /// iteration counter; counts how many families we've looked at
    pub iter: usize,
    /// primes we've discovered so far, in two different formats
    /// depending on the order we discovered these, they may not be minimal!
    pub primes: CandidateSequences,

    /// For potentially getting insight into what's going on
    pub stats: RefCell<Stats>,
}

#[derive(Debug, Clone)]
pub enum SearchNode {
    Arbitrary(Family),
    Simple(SimpleFamily),
}

impl Weight for SearchNode {
    fn weight(&self) -> usize {
        match self {
            SearchNode::Arbitrary(x) => x.weight(),
            SearchNode::Simple(x) => x.before.0.len() + x.num_repeats + x.after.0.len(),
        }
    }
}

impl SearchContext {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            stats: RefCell::new(Stats::default()),
        }
    }

    fn explore_node(&mut self, family: SearchNode) -> (Vec<SearchNode>, String) {
        // Say our family is xL*z.
        // We want to explore all possible children with weight one more than this one.
        let children = match family {
            SearchNode::Arbitrary(family) => {
                debug!(" Exploring {}", family);
                self.explore_family(family)
            }
            SearchNode::Simple(family) => {
                debug!(" Exploring simple {}", family);
                self.explore_simple_family(family)
            }
        };
        self.stats.borrow_mut().num_branches_explored += 1;
        children
    }

    fn explore_family(&mut self, family: Family) -> (Vec<SearchNode>, String) {
        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = family.contract();
        if let Some(p) = self.test_for_contained_prime(&seq) {
            assert_ne!(&seq, p);
            let reason = format!("  Discarding {}, contains prime {}", family, p);
            debug!("{}", reason);
            return (vec![], reason);
        }

        trace!("  Testing for primality {}", seq);
        let value = seq.value(self.base);
        if self.test_for_prime(&value) {
            let reason = format!("  Saving {}, contracts to prime", family);
            self.primes.insert(seq);
            debug!("{}", reason);
            return (vec![], reason);
        }

        // Then, we try to reduce the cores.
        let mut family = self.reduce_cores(family);
        family.simplify();
        if family.cores.is_empty() {
            let reason = format!("  {} was reduced to trivial string", family);
            debug!("{}", reason);
            return (vec![], reason);
        }

        // Now, run some tests to see whether this family is guaranteed to
        // be composite.
        if self.test_for_perpetual_composite(&family) {
            let reason = format!("  Discarding {}, is always composite", family);
            debug!("{}", reason);
            return (vec![], reason);
        }

        // TODO: is this right?
        // Check if this family is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        let mut family = match SimpleFamily::try_from(family) {
            Ok(simple) => {
                return (
                    vec![SearchNode::Simple(simple)],
                    "reduced to simple".to_string(),
                );
            }
            // can't convert, put it back to normal
            Err(f) => f,
        };

        // Let's see if we can split it in an interesting way
        if let Some(children) = self.split_on_repeat(&family, 3) {
            return (
                children.into_iter().map(SearchNode::Arbitrary).collect(),
                "split on repeat".to_string(),
            );
        }

        if family.weight() >= 2 {
            if let Some(child) = self.split_on_necessary_digit(&family) {
                return (
                    vec![SearchNode::Arbitrary(child)],
                    "split on necessary digit".to_string(),
                );
            }
        }

        // If we couldn't eliminate the family, let's split it, left or right.
        // We can't split on a non-empty core, but after we simplify, we shouldn't
        // have to worry about that.
        let slot = family.weight() % family.cores.len();
        debug_assert!(!family.cores[slot].is_empty());
        let mut children = if family.weight() == 1 {
            debug!("  Splitting {} left", family);
            family.split_left(slot)
        } else {
            debug!("  Splitting {} right", family);
            family.split_right(slot)
        };

        // We also need to consider the case where the chosen core expands to
        // the empty string. However, in the case where there's one core, this
        // is pretty redundant with the work we're doing in reduce_core().
        // For example: if we reduce a[xyz]c, we test the primality of axc, ayc
        // and azc. So after we split, and get ax[xyz]c, there's no need to
        // test ax[]c again.
        if family.cores.len() > 1 {
            family.cores[slot].clear();
            family.simplify();
            children.push(family);
        }

        (
            children.into_iter().map(SearchNode::Arbitrary).collect(),
            "split".to_string(),
        )
    }

    fn explore_simple_family(&mut self, mut family: SimpleFamily) -> (Vec<SearchNode>, String) {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.
        // Even composite testing isn't very useful here, since we
        // can't get here without passing through a composite test.
        // All we can do is add another digit and see if it becomes
        // prime or not. Or contains another prime.

        // Test if it contains a prime
        let start = Instant::now();
        for prime in self.primes.iter() {
            let result = is_substring_of_simple(prime, &family);
            self.stats.borrow_mut().num_simple_substring_checks += 1;
            if let SubstringResult::Yes = result {
                let reason = format!("  Discarding {}, contains prime {}", family, prime);
                self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();
                debug!("{}", reason);
                return (vec![], reason);
            }
        }
        self.stats.borrow_mut().duration_simple_substring_checks += start.elapsed();

        // Test if it is a prime
        let value = family.value(self.base);

        if self.test_for_prime(&value) {
            let reason = format!("  Saving {}, is prime", family);
            debug!("{}", reason);
            let seq = family.sequence();
            self.primes.insert(seq);
            return (vec![], reason);
        }

        family.num_repeats += 1;
        (vec![SearchNode::Simple(family)], "increment".to_string())
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
                    debug!("  Discarding {}, contains prime {}", seq, p);
                    continue;
                }

                trace!("  Testing for primality {}", seq);
                let value = seq.value(self.base);
                if self.test_for_prime(&value) {
                    debug!("  Saving {}, is prime", seq);
                    self.primes.insert(seq);
                } else {
                    allowed_digits.push(digit);
                }
            }

            *core = allowed_digits;
        }
        // Now we've reduced the core, and have a new family.
        debug!("  Reducing {} to {}", old_family, family);
        family
    }

    fn test_for_contained_prime(&self, seq: &DigitSeq) -> Option<&DigitSeq> {
        let start = Instant::now();
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        let result = self.primes.iter().find(|subseq| {
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
            debug!("  {} has divisor {}", family, factor);
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect families like 4[6]9 being divisible by 7.
        for stride in 1..=2 {
            if let Some(factors) = find_perpetual_factor(self.base, family, stride) {
                debug!(
                    "  {} is divisible by {}",
                    family,
                    factors.iter().format(", ")
                );
                return true;
            }
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(self.base, family) {
            debug!(
                "  {} is divisible by either {} or {}",
                family, even_factor, odd_factor
            );
            return true;
        }

        false
    }

    fn split_on_repeat(&self, family: &Family, max_repeats: usize) -> Option<Vec<Family>> {
        debug!(" Trying to split {}", family);
        for (i, core) in family.cores.iter().enumerate() {
            for d in core.iter().copied() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = family.substitute_multiple(i, std::iter::repeat_n(d, n));
                    if let Some(p) = self.test_for_contained_prime(&seq) {
                        assert_ne!(&seq, p);
                        debug!("  {} contains a prime {}", seq, p);

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

                        debug!(
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
                debug!("  {} must have a {}, transforming into {}", family, d, new);
                return Some(new);
            }
        }
        None
    }
}

pub enum SubstringResult {
    Yes,
    Never,
    Eventually(usize),
}

pub fn is_substring_of_simple(needle: &DigitSeq, haystack: &SimpleFamily) -> SubstringResult {
    let mut needle_iter = needle.0.iter().copied().peekable();
    let mut repeats_required = 0;

    // Three stages: go through before, then center, then after.
    // Try to consume the whole needle.
    for d in haystack.before.0.iter().copied() {
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
            Some(d2) if haystack.center == *d2 => {
                repeats_required += 1;
                needle_iter.next();
            }
            // different digit, time to leave
            Some(_) => break,
            // done with the needle!
            None => break,
        }
    }

    for d in haystack.after.0.iter().copied() {
        match needle_iter.peek() {
            Some(d2) if d == *d2 => {
                needle_iter.next();
            }
            Some(_) => {}
            None => break,
        }
    }

    if needle_iter.peek().is_some() {
        SubstringResult::Never
    } else if repeats_required <= haystack.num_repeats {
        SubstringResult::Yes
    } else {
        SubstringResult::Eventually(repeats_required)
    }
}

impl std::fmt::Display for SearchNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SearchNode::Arbitrary(p) => p.fmt(f),
            SearchNode::Simple(p) => p.fmt(f),
        }
    }
}
