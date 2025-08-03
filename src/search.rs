mod composite;
mod explore;
mod split;

use std::ops::ControlFlow;
use std::time::{Duration, Instant};

use itertools::Itertools;
use log::{debug, info, trace};
use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use num_traits::One;

use self::composite::{find_even_odd_factor, find_perpetual_factor, shares_factor_with_base};
use self::explore::{Frontier, Weight};
use crate::data_structures::{
    is_proper_substring, AppendTree, AppendTreeNodeID, CandidateSequences,
};
use crate::digits::DigitSeq;
use crate::families::{Core, Family, Sequence, SimpleFamily};
use crate::math::gcd_reduce;
use crate::search::composite::composite_checks_for_simple;
use crate::SearchResults;

#[macro_export]
macro_rules! log_to_tree {
    ($tracer:expr, $lvl:expr, $($arg:tt)+) => {
        if log::log_enabled!($lvl) {
            $tracer.log(format!($($arg)+))
        }
    };
}

#[macro_export]
macro_rules! debug_to_tree {
    ($tracer:expr, $($arg:tt)+) => {
        crate::log_to_tree!($tracer, log::Level::Debug, $($arg)+)
    };
}

pub fn search_for_simple_families(
    base: u8,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    stop_when_simple: bool,
    tree_log: bool,
) -> SearchResults {
    let mut tree = SearchTree::new(base, tree_log);

    tree.explore_until(|tree| {
        let weight = match tree.frontier.min_weight() {
            Some(w) => w,
            // this means the frontier is empty!
            None => return ControlFlow::Break(()),
        };

        if let Some(max) = max_weight {
            if weight > max {
                info!("Reached weight cutoff; stopping...");
                return ControlFlow::Break(());
            }
        }

        if let Some(max) = max_iter {
            if tree.ctx.iter >= max {
                info!("Reached iteration cutoff; stopping...");
                return ControlFlow::Break(());
            }
        }

        if stop_when_simple && {
            tree.frontier
                .iter()
                .all(|node| matches!(node.family, NodeType::Simple(_)))
        } {
            info!("All remaining families are simple; stopping...");
            return ControlFlow::Break(());
        }

        info!(
            "Iteration {} - Weight {} - {} branches",
            tree.ctx.iter,
            weight,
            tree.frontier.len()
        );

        ControlFlow::Continue(())
    });

    tree.into_results()
}

pub struct SearchTree {
    /// Nodes that we haven't explored yet
    pub frontier: Frontier<SearchNode>,
    /// Everything else
    pub ctx: SearchContext,
}

impl SearchTree {
    pub fn new(base: u8, tree_log: bool) -> Self {
        let ctx = SearchContext::new(base, tree_log);
        let initial_node = SearchNode {
            family: NodeType::Arbitrary(Family::any(base)),
            id: ctx.tracer.root(),
        };
        let frontier = Frontier::start(initial_node);

        Self { frontier, ctx }
    }

    pub fn explore_until(
        &mut self,
        mut stop_condition: impl FnMut(&SearchTree) -> ControlFlow<()>,
    ) {
        loop {
            if stop_condition(self).is_break() {
                break;
            }

            // Explore and break if nothing happened
            if self.explore_once().is_break() {
                break;
            }
        }
    }

    fn explore_once(&mut self) -> ControlFlow<()> {
        self.frontier
            .explore_next(|node| self.ctx.explore_node(node))?;

        self.ctx.iter += 1;
        ControlFlow::Continue(())
    }

    pub fn into_results(mut self) -> SearchResults {
        self.ctx.primes.sort();

        // print the tree to stdout if we're tracing
        match self.ctx.tracer {
            Tracer::Real(t, _) => t.pretty_print_to_stdout(),
            Tracer::Dummy(_) => {}
        }

        // Pull the unsolved branches and return them
        let mut ret = SearchResults {
            primes: self.ctx.primes,
            simple_families: vec![],
            other_families: vec![],
            stats: self.ctx.stats,
        };
        for node in self.frontier.iter().cloned() {
            match node.family {
                NodeType::Arbitrary(family) => ret.other_families.push(family),
                // TODO: return the whole node, so the other stages can benefit here!
                NodeType::Simple(node) => ret.simple_families.push(node.family),
            }
        }

        ret
    }
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
    pub stats: Stats,
    /// For tracking our paths through the search space in a more
    /// understandable format.
    tracer: Tracer,
}

#[derive(Debug, Clone)]
struct SearchNode {
    family: NodeType,
    id: AppendTreeNodeID,
}

#[derive(Debug, Clone)]
enum NodeType {
    Arbitrary(Family),
    Simple(SimpleNode),
}

#[derive(Debug, Clone)]
struct SimpleNode {
    family: SimpleFamily,
    // Some simple families are too large to fit into
    // native integer types :/ Hopefully they get killed
    // by other means in the first stage, since otherwise
    // they'll fail in the second stage since we can't
    // make a sequence for them.
    sequence: Option<Sequence>,
    composite_tested: bool,
}

enum Never {}

enum Tracer {
    // The tree and the current index we're looking at
    Real(AppendTree<String>, AppendTreeNodeID),
    // we need to get an id sometimes, but we can't make children or items.
    // using an empty type forces this :)
    Dummy(AppendTree<Never>),
}

impl SearchContext {
    pub fn new(base: u8, tree_log: bool) -> Self {
        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            stats: Stats::default(),
            tracer: if tree_log {
                Tracer::new()
            } else {
                Tracer::dummy()
            },
        }
    }

    fn explore_node(&mut self, node: SearchNode) -> Vec<SearchNode> {
        let node_id = node.id;
        self.tracer.set_id(node_id);

        // Say our family is xL*z.
        // We want to explore all possible children with weight one more than this one.
        let children = match node.family {
            NodeType::Arbitrary(family) => {
                debug!(" Exploring {}", family);
                self.explore_family(family)
            }
            NodeType::Simple(node) => {
                debug!(" Exploring simple {}", node.family);
                self.explore_simple_family(node)
            }
        };
        self.stats.num_branches_explored += 1;
        let children = children
            .into_iter()
            .map(|family| {
                let child_id = self
                    .tracer
                    .make_child(node_id, family.to_string())
                    .expect("node id must be in tree");
                SearchNode {
                    family,
                    id: child_id,
                }
            })
            .collect();
        children
    }

    fn explore_family(&mut self, family: Family) -> Vec<NodeType> {
        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = family.contract();
        // TODO: this borrows self, preventing us from using self.tracer later.
        // can that be improved?
        if let Some(p) = self.test_for_contained_prime(&seq).cloned() {
            assert_ne!(seq, p);
            debug!("  Discarding {}, contains prime {}", family, p);
            debug_to_tree!(self.tracer, "Discarding, contains prime {}", p);
            return vec![];
        }

        trace!("  Testing for primality {}", seq);
        let value = seq.value(self.base);
        if self.test_for_prime(&value) {
            debug!("  Saving {}, contracts to prime", family);
            debug_to_tree!(self.tracer, "Saving, contracts to prime");
            self.primes.insert(seq);
            return vec![];
        }

        // Then, we try to reduce the cores.
        let mut family = self.reduce_cores(family);
        family.simplify();
        if family.cores.is_empty() {
            debug!("  {} was reduced to trivial string", family);
            debug_to_tree!(self.tracer, "Reduced to trivial string");
            return vec![];
        }

        // Now, run some tests to see whether this family is guaranteed to
        // be composite.
        if self.test_for_perpetual_composite(&family) {
            debug!("  Discarding {}, is always composite", family);
            debug_to_tree!(self.tracer, "Discarding, is always composite");
            return vec![];
        }

        // TODO: is this right?
        // Check if this family is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        let mut family = match SimpleFamily::try_from(family) {
            Ok(family) => {
                let sequence = Sequence::try_from_family(&family, self.base).ok();
                return vec![NodeType::Simple(SimpleNode {
                    family,
                    sequence,
                    composite_tested: false,
                })];
            }
            // can't convert, put it back to normal
            Err(f) => f,
        };

        // Let's see if we can split it in an interesting way
        if let Some(children) = self.split_on_repeat(&family, 3) {
            return children.into_iter().map(NodeType::Arbitrary).collect();
        }

        if family.weight() >= 2 {
            if let Some(child) = self.split_on_necessary_digit(&family) {
                return vec![NodeType::Arbitrary(child)];
            }
        }

        if family.weight() >= 2 {
            if let Some(children) = self.split_on_incompatible_digits(&family) {
                return children.into_iter().map(NodeType::Arbitrary).collect();
            }
        }

        // If we couldn't eliminate the family, let's split it, left or right.
        // We can't split on a non-empty core, but after we simplify, we shouldn't
        // have to worry about that.
        let weight = family.weight();
        let slot = (2 * weight) % family.cores.len();
        debug_assert!(!family.cores[slot].is_empty());
        let mut children = if family.weight() == 1 {
            debug!("  Splitting {} left on core {}", family, slot);
            debug_to_tree!(self.tracer, "Splitting left on core {}", slot);
            family.split_left(slot)
        } else {
            debug!("  Splitting {} right on core {}", family, slot);
            debug_to_tree!(self.tracer, "Splitting right on core {}", slot);
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

        children.into_iter().map(NodeType::Arbitrary).collect()
    }

    fn explore_simple_family(&mut self, mut node: SimpleNode) -> Vec<NodeType> {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.

        // We should do a composite test though; there's some specialized
        // composite tests that we can only do on simple families. We should
        // only do them once though.
        if !node.composite_tested {
            if composite_checks_for_simple(self.base, &node) {
                debug!("  Discarding {}, is always composite", node.family);
                debug_to_tree!(self.tracer, "Discarding, is always composite");
                return vec![];
            }
            node.composite_tested = true;
        }

        // Other that that, all we can do is add another digit and see if it
        // becomes prime or not. Or contains another prime.

        let start = Instant::now();
        for prime in self.primes.iter() {
            self.stats.num_simple_substring_checks += 1;
            if node
                .family
                .will_contain_at(prime)
                .is_some_and(|n| n <= node.family.min_repeats)
            {
                debug!("  Discarding {}, contains prime {}", node.family, prime);
                debug_to_tree!(self.tracer, "Discarding, contains prime {}", prime);
                self.stats.duration_simple_substring_checks += start.elapsed();
                return vec![];
            }
        }
        self.stats.duration_simple_substring_checks += start.elapsed();

        // Test if it is a prime
        let value = node.family.value(self.base);

        if self.test_for_prime(&value) {
            debug!("  Saving {}, is prime", node.family);
            debug_to_tree!(self.tracer, "Saving, is prime");
            let seq = node.family.sequence();
            self.primes.insert(seq);
            return vec![];
        }

        node.family.min_repeats += 1;
        vec![NodeType::Simple(node)]
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

                if let Some(p) = self.test_for_contained_prime(&seq).cloned() {
                    assert_ne!(seq, p);
                    debug!("  Discarding {}, contains prime {}", seq, p);
                    debug_to_tree!(self.tracer, "Discarding {}, contains prime {}", seq, p);
                    continue;
                }

                trace!("  Testing for primality {}", seq);
                let value = seq.value(self.base);
                if self.test_for_prime(&value) {
                    debug!("  Saving {}, is prime", seq);
                    debug_to_tree!(self.tracer, "Saving {}, is prime", seq);
                    self.primes.insert(seq);
                } else {
                    allowed_digits.push(digit);
                }
            }

            *core = Core::new(allowed_digits);
        }
        // Now we've reduced the core, and have a new family.
        debug!("  Reducing {} to {}", old_family, family);
        debug_to_tree!(self.tracer, "Reducing to {}", family);
        family
    }

    fn test_for_contained_prime(&mut self, seq: &DigitSeq) -> Option<&DigitSeq> {
        let start = Instant::now();
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        let result = self.primes.iter().find(|subseq| {
            self.stats.num_substring_checks += 1;
            is_proper_substring(subseq, seq)
        });

        self.stats.duration_substring_checks += start.elapsed();
        result
    }

    fn test_for_prime(&mut self, value: &BigUint) -> bool {
        let start = Instant::now();
        let result = is_prime(value, None).probably();
        self.stats.duration_primality_checks += start.elapsed();
        self.stats.num_primality_checks += 1;
        result
    }

    fn test_for_perpetual_composite(&mut self, family: &Family) -> bool {
        // This function is used to eliminate families that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity, unless specified otherwise.

        // p divides BASE (e.g., 2, 5)
        if let Some(factor) = shares_factor_with_base(self.base, family) {
            debug!("  {} has divisor {}", family, factor);
            debug_to_tree!(self.tracer, "Has divisor {}", factor);
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect families like 4[6]9 being divisible by 7.
        for stride in 1..=3 {
            if let Some(factors) = find_perpetual_factor(self.base, family, stride) {
                debug!(
                    "  {} is divisible by {}",
                    family,
                    factors.iter().format(", ")
                );
                debug_to_tree!(self.tracer, "Divisible by {}", factors.iter().format(","));
                return true;
            }
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(self.base, family) {
            debug!(
                "  {} is divisible by either {} or {}",
                family, even_factor, odd_factor
            );
            debug_to_tree!(
                self.tracer,
                "Divisible by either {} or {}",
                even_factor,
                odd_factor
            );
            return true;
        }

        false
    }
}

impl Tracer {
    fn new() -> Self {
        let tree = AppendTree::new();
        let root = tree.root();
        Self::Real(tree, root)
    }

    fn dummy() -> Self {
        Self::Dummy(AppendTree::new())
    }

    fn root(&self) -> AppendTreeNodeID {
        match self {
            Tracer::Real(t, _) => t.root(),
            Tracer::Dummy(t) => t.root(),
        }
    }

    fn make_child(
        &mut self,
        node_id: AppendTreeNodeID,
        tag: String,
    ) -> Result<AppendTreeNodeID, String> {
        match self {
            Tracer::Real(t, _) => t.make_child(node_id, tag),
            // we're never going to log anything, so just keep returning the root
            // to satisfy the type system
            Tracer::Dummy(t) => Ok(t.root()),
        }
    }

    fn set_id(&mut self, node_id: AppendTreeNodeID) {
        match self {
            Tracer::Real(_, id) => *id = node_id,
            Tracer::Dummy(_) => {}
        }
    }

    fn log(&mut self, item: String) {
        match self {
            Tracer::Real(t, id) => {
                t.append(*id, item).expect("logging to nonexistent id");
            }
            Tracer::Dummy(_) => {}
        }
    }
}

impl std::fmt::Display for NodeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            NodeType::Arbitrary(p) => p.fmt(f),
            NodeType::Simple(node) => node.family.fmt(f),
        }
    }
}

impl std::fmt::Display for SearchNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.family {
            NodeType::Arbitrary(p) => p.fmt(f),
            NodeType::Simple(node) => node.family.fmt(f),
        }
    }
}

impl Weight for SearchNode {
    fn weight(&self) -> usize {
        match &self.family {
            NodeType::Arbitrary(x) => x.weight(),
            NodeType::Simple(node) => {
                node.family.before.0.len() + node.family.min_repeats + node.family.after.0.len()
            }
        }
    }
}
