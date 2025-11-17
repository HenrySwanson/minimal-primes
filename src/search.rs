mod composite;
mod frontier;
mod gcd;
mod split;

use std::ops::ControlFlow;
use std::time::Instant;

use itertools::Itertools;
use log::{debug, trace};
use num_bigint::BigUint;
use num_prime::buffer::PrimeBufferExt;

use self::composite::{find_even_odd_factor, find_periodic_factor, shares_factor_with_base};
use self::frontier::{Frontier, Weight};
use crate::candidates::CandidateIndices;
use crate::context::SearchContext;
use crate::digits::DigitSeq;
use crate::families::{Core, Family, SimpleFamily};
use crate::logging::AppendTreeNodeID;
use crate::search::composite::{
    check_residues_mod_30, composite_checks_for_simple, find_common_factor, find_two_factors,
};
use crate::RemainingNodes;

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
        $crate::log_to_tree!($tracer, log::Level::Debug, $($arg)+)
    };
}

pub struct SearchTree {
    pub nodes: Frontier<SearchNode>,
}

impl SearchTree {
    pub fn new(ctx: &SearchContext) -> Self {
        let initial_node = SearchNode {
            family: NodeType::Arbitrary(Family::any(ctx.base)),
            possible_contained_primes: ctx.primes.indices_all(),
            id: ctx.tracer.root(),
        };
        let frontier = Frontier::start(initial_node);

        Self { nodes: frontier }
    }

    pub fn num_nodes_to_solve(&self) -> usize {
        self.nodes
            .iter()
            .filter(|node| match &node.family {
                NodeType::Arbitrary(_) => true,
                NodeType::Simple(simple_node) => !simple_node.composite_tested,
            })
            .count()
    }

    pub fn any_nodes_to_solve(&self) -> bool {
        self.nodes.iter().any(|node| match &node.family {
            NodeType::Arbitrary(_) => true,
            NodeType::Simple(simple_node) => !simple_node.composite_tested,
        })
    }

    pub fn explore_once(&mut self, ctx: &mut SearchContext) -> ControlFlow<()> {
        self.nodes.explore_next(|node| node.explore(ctx))?;

        ctx.iter += 1;
        ControlFlow::Continue(())
    }

    pub fn into_results(self) -> RemainingNodes {
        // Pull the unsolved branches and return them
        let mut ret = RemainingNodes {
            simple_families: vec![],
            other_families: vec![],
        };
        for node in self.nodes.into_iter() {
            match node.family {
                NodeType::Arbitrary(family) => ret.other_families.push(family),
                // TODO: return the whole node, so the other stages can benefit here!
                NodeType::Simple(node) => ret.simple_families.push(node.family),
            }
        }

        ret
    }
}

#[derive(Debug, Clone)]
pub struct SearchNode {
    family: NodeType,
    // TODO: put this in simplenode too
    possible_contained_primes: CandidateIndices,
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
    composite_tested: bool,
    /// When min_repeats equals this number, we can delete this
    /// family, because it contains this prime.
    dies_at: Option<(usize, DigitSeq)>,
}

impl SearchNode {
    fn explore(self, ctx: &mut SearchContext) -> Vec<SearchNode> {
        let node_id = self.id;
        ctx.tracer.set_id(node_id);

        // Say our family is xL*z.
        // We want to explore all possible children with weight one more than this one.
        let mut pcp = self.possible_contained_primes;
        let children = match self.family {
            NodeType::Arbitrary(family) => {
                debug!(" Exploring {family}");
                family.explore(&mut pcp, ctx)
            }
            NodeType::Simple(node) => {
                debug!(" Exploring simple {}", node.family);
                node.explore(&mut pcp, ctx)
            }
        };
        ctx.stats.num_branches_explored += 1;

        children
            .into_iter()
            .map(|family| {
                let child_id = ctx
                    .tracer
                    .make_child(node_id, family.to_string())
                    .expect("node id must be in tree");
                SearchNode {
                    family,
                    possible_contained_primes: pcp.clone(),
                    id: child_id,
                }
            })
            .collect()
    }
}

impl Family {
    fn explore(
        mut self,
        possible_contained_primes: &mut CandidateIndices,
        ctx: &mut SearchContext,
    ) -> Vec<NodeType> {
        // We have to be careful about not generating leading zeros. There's
        // a lot of different ways we could do that, but one possible way is
        // just to rig things so that on the first round, we split left (and
        // check for primes).
        // TODO: can we do this elsewhere? maybe some post-processing step in
        // the explore_node? this feels iffy

        // Now is a good time for us to narrow down the potential primes this
        // family could contain.
        let mut new_contained_primes = ctx.primes.indices_none();
        for (i, prime) in ctx.primes.get_many(possible_contained_primes) {
            let start = Instant::now();
            if self.could_contain(prime) {
                new_contained_primes.add(i);
            }
            ctx.stats.num_could_contains += 1;
            ctx.stats.duration_could_contains += start.elapsed();
        }
        *possible_contained_primes = new_contained_primes;

        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = self.contract();
        // TODO: this borrows self, preventing us from using self.tracer later.
        // can that be improved?
        if let Some(p) = ctx
            .test_for_contained_prime(&seq, possible_contained_primes)
            .cloned()
        {
            assert_ne!(seq, p);
            debug!("  Discarding {self}, contains prime {p}");
            debug_to_tree!(ctx.tracer, "Discarding, contains prime {}", p);
            ctx.stats.branch_stats.contains_prime += 1;
            return vec![];
        }

        trace!("  Testing for primality {seq}");
        let value = seq.value(ctx.base);
        if ctx.test_for_prime(&value) {
            debug!("  Saving {self}, contracts to prime");
            debug_to_tree!(ctx.tracer, "Saving, contracts to prime");
            ctx.primes.insert(seq);
            ctx.stats.branch_stats.is_new_prime += 1;
            return vec![];
        }

        // Then, we try to reduce the cores.
        self.reduce_cores(possible_contained_primes, ctx);
        self.simplify();
        if self.cores.is_empty() {
            debug!("  {self} was reduced to trivial string");
            debug_to_tree!(ctx.tracer, "Reduced to trivial string");
            ctx.stats.branch_stats.is_trivial_string += 1;
            return vec![];
        }

        // Now, run some tests to see whether this family is guaranteed to
        // be composite.
        if self.test_for_perpetual_composite(ctx) {
            debug!("  Discarding {self}, is always composite");
            debug_to_tree!(ctx.tracer, "Discarding, is always composite");
            ctx.stats.branch_stats.detected_composite += 1;
            return vec![];
        }

        // TODO: is this right?
        // Check if this family is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        if let Ok(family) = SimpleFamily::try_from(self.clone()) {
            ctx.stats.branch_stats.simplified += 1;
            return vec![NodeType::Simple(SimpleNode {
                family,
                composite_tested: false,
                dies_at: None,
            })];
        }

        // Let's see if we can split it in an interesting way
        // TODO: context-ify the splitting functions too!
        if self.weight() >= 2 {
            if let Some(children) = ctx.split_on_limited_digit(&self, 3, possible_contained_primes)
            {
                ctx.stats.branch_stats.split_on_limited_digit += 1;
                return children.into_iter().map(NodeType::Arbitrary).collect();
            }
        }

        if self.weight() >= 4 {
            if let Some(children) =
                ctx.split_on_incompatible_digits_different_cores(&self, possible_contained_primes)
            {
                ctx.stats.branch_stats.split_on_incompatible_different_cores += 1;
                return children.into_iter().map(NodeType::Arbitrary).collect();
            }

            if let Some(children) =
                ctx.split_on_incompatible_digits(&self, possible_contained_primes)
            {
                ctx.stats.branch_stats.split_on_incompatible_same_core += 1;
                return children.into_iter().map(NodeType::Arbitrary).collect();
            }
        }

        // this was introduced to kill long derivation chains of the form x[ab]*y that
        // we have trouble with otherwise. only invoke it when we are really stuck on
        // something.
        if self.weight() >= 10 {
            if let Some(children) =
                ctx.split_on_forbidden_sandwich(&self, possible_contained_primes)
            {
                ctx.stats.branch_stats.split_on_forbidden_sandwich += 1;
                return children.into_iter().map(NodeType::Arbitrary).collect();
            }
        }

        if self.weight() >= 5 {
            // This one doesn't actually simplify any cores, in fact, it'll make the
            // branch more complicated!. However, it might kickstart some branch elimination
            // by making us intersect another prime. So this check should always be last.
            if let Some(child) = ctx.split_on_necessary_digit(&self) {
                ctx.stats.branch_stats.split_on_necessary_digit += 1;
                return vec![NodeType::Arbitrary(child)];
            }
        }

        // If we couldn't eliminate the family, let's split it, left or right.
        // We can't split on a non-empty core, but after we simplify, we shouldn't
        // have to worry about that.

        // Which core do we split on? And should we split left or right?
        // We want to pick arbitrarily and get a good distribution over the cores.
        // Previously, we just used the weight of the family as a counter, and
        // looked at it mod (#cores) and mod 2 to make that decision. But that
        // resulted in occasional "resonance", where we'd keep making the same
        // decision over and over again. So instead, we still use the family weight,
        // but we jumble it up with some silly nonsense in order to get a more
        // unpredictable sequence.
        fn bit_mixer(h: usize) -> usize {
            // stolen from murmur3's finalizer
            let mut h = h as u64;
            h ^= h >> 33;
            h = h.wrapping_mul(0xff51afd7ed558ccd);
            h ^= h >> 33;
            h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
            h ^= h >> 33;
            h as usize
        }
        let magic = match self.weight() {
            0 => 0,
            1 => 1,
            w => bit_mixer(w),
        };

        let slot = (magic >> 1) % self.cores.len();
        debug_assert!(!self.cores[slot].is_empty());
        let mut children = if magic % 2 == 1 {
            debug!("  Splitting {self} left on core {slot}");
            debug_to_tree!(ctx.tracer, "Splitting left on core {}", slot);
            self.expand(slot)
        } else {
            debug!("  Splitting {self} right on core {slot}");
            debug_to_tree!(ctx.tracer, "Splitting right on core {}", slot);
            self.expand_right(slot)
        };

        // We also need to consider the case where the chosen core expands to
        // the empty string. However, in the case where there's one core, this
        // is pretty redundant with the work we're doing in reduce_core().
        // For example: if we reduce a[xyz]c, we test the primality of axc, ayc
        // and azc. So after we split, and get ax[xyz]c, there's no need to
        // test ax[]c again.
        if self.cores.len() > 1 {
            self.cores[slot].clear();
            self.simplify();
            children.push(self);
        }

        ctx.stats.branch_stats.explored_generically += 1;
        children.into_iter().map(NodeType::Arbitrary).collect()
    }
}

impl SimpleNode {
    fn explore(
        mut self,
        possible_contained_primes: &mut CandidateIndices,
        ctx: &mut SearchContext,
    ) -> Vec<NodeType> {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.

        // We should do a composite test though; there's some specialized
        // composite tests that we can only do on simple families. We should
        // only do them once though.
        if !self.composite_tested {
            if composite_checks_for_simple(ctx.base, &self) {
                debug!("  Discarding {}, is always composite", self.family);
                debug_to_tree!(ctx.tracer, "Discarding, is always composite");
                ctx.stats.branch_stats.detected_composite += 1;
                return vec![];
            }
            self.composite_tested = true;
        }

        // Other that that, all we can do is add another digit and see if it
        // becomes prime or not. Or contains another prime.

        // On a previous loop, we may have established when this family contains
        // a prime. Check it.
        if let Some((dies_at, prime)) = &self.dies_at {
            if self.family.min_repeats >= *dies_at {
                debug!("  Discarding {}, contains prime {}", self.family, prime);
                debug_to_tree!(ctx.tracer, "Discarding, contains prime {}", prime);
                ctx.stats.branch_stats.contains_prime += 1;
                return vec![];
            }
        }

        // Now check any new primes.
        for (_, prime) in ctx.primes.get_many(possible_contained_primes) {
            let start = Instant::now();
            if let Some(n) = self.family.will_contain_at(prime) {
                if n <= self.family.min_repeats {
                    debug!("  Discarding {}, contains prime {}", self.family, prime);
                    debug_to_tree!(ctx.tracer, "Discarding, contains prime {}", prime);
                    ctx.stats.branch_stats.contains_prime += 1;
                    return vec![];
                }

                // otherwise, we should incorporate this into dies_at
                match self.dies_at {
                    // its better to have smaller n; skip this if so
                    Some((old_n, _)) if old_n <= n => {}
                    // otherwise, update
                    Some(_) | None => self.dies_at = Some((n, prime.clone())),
                }
            }
            ctx.stats.num_simple_substring_checks += 1;
            ctx.stats.duration_simple_substring_checks += start.elapsed();
        }
        *possible_contained_primes = ctx.primes.indices_none(); // resets our collection

        // Test if it is a prime
        let value = self.family.value(ctx.base);

        if ctx.test_for_prime(&value) {
            debug!("  Saving {}, is prime", self.family);
            debug_to_tree!(ctx.tracer, "Saving, is prime");
            ctx.stats.branch_stats.is_new_prime += 1;
            let seq = self.family.contract();
            ctx.primes.insert(seq);
            return vec![];
        }

        self.family.min_repeats += 1;
        ctx.stats.branch_stats.explored_generically += 1;
        vec![NodeType::Simple(self)]
    }
}

impl Family {
    fn reduce_cores(
        &mut self,
        possible_contained_primes: &CandidateIndices,
        ctx: &mut SearchContext,
    ) {
        let old_family = self.clone();
        for (i, core) in self.cores.iter_mut().enumerate() {
            // Substitute elements from the core into the string to see if any
            // of them contain or are a prime.
            // NOTE: this is where we generate minimal primes of (weight + 1), so
            // next loop, those should all be available.
            let mut allowed_digits = vec![];
            for digit in core.iter() {
                let seq = old_family.substitute(i, digit);

                if let Some(p) = ctx
                    .test_for_contained_prime(&seq, possible_contained_primes)
                    .cloned()
                {
                    assert_ne!(seq, p);
                    debug!("  Discarding {seq}, contains prime {p}");
                    debug_to_tree!(ctx.tracer, "Discarding {}, contains prime {}", seq, p);
                    continue;
                }

                trace!("  Testing for primality {seq}");
                let value = seq.value(ctx.base);
                if ctx.test_for_prime(&value) {
                    debug!("  Saving {seq}, is prime");
                    debug_to_tree!(ctx.tracer, "Saving {}, is prime", seq);
                    ctx.primes.insert(seq);
                } else {
                    allowed_digits.push(digit);
                }
            }

            *core = Core::new(allowed_digits);
        }
        // Now we've reduced the core, and have a new family.
        debug!("  Reducing {old_family} to {self}");
        debug_to_tree!(ctx.tracer, "Reducing to {}", self);
    }
}

impl SearchContext {
    /// Checks whether any of the candidate primes are (properly) contained in the
    /// given sequence. If so, return a reference to that prime.
    fn test_for_contained_prime(
        &mut self,
        seq: &DigitSeq,
        possible_contained_primes: &CandidateIndices,
    ) -> Option<&DigitSeq> {
        let start = Instant::now();
        // We don't need to search for *all* possible primes, just the minimal
        // ones. And if we've been doing our job right, we should have a complete
        // list of them (up to a length limit).
        let result = self
            .primes
            .get_many(possible_contained_primes)
            .find(|(_, subseq)| {
                self.stats.num_substring_checks += 1;
                seq.properly_contains(subseq)
            });

        self.stats.duration_substring_checks += start.elapsed();
        let (_, seq) = result?;
        Some(seq)
    }

    fn test_for_prime(&mut self, value: &BigUint) -> bool {
        let start = Instant::now();
        let result = self.prime_buffer.is_prime(value, None).probably();
        self.stats.duration_primality_checks += start.elapsed();
        self.stats.num_primality_checks += 1;
        result
    }
}

impl Family {
    fn test_for_perpetual_composite(&self, ctx: &mut SearchContext) -> bool {
        // This function is used to eliminate families that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity, unless specified otherwise.

        // p divides BASE (e.g., 2, 5)
        if let Some(factor) = shares_factor_with_base(ctx.base, self) {
            debug!("  {self} has divisor {factor}");
            debug_to_tree!(ctx.tracer, "Has divisor {}", factor);
            return true;
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect families like 4[6]9 being divisible by 7.
        if let Some(divisor) = find_common_factor(ctx.base, self) {
            debug!("  {self} is divisible by {divisor}");
            debug_to_tree!(ctx.tracer, "Divisible by {}", divisor);
            return true;
        }
        // start at stride 2; `find_guaranteed_factor` effectively handles stride 1
        for stride in 2..=4 {
            if let Some(factors) = find_periodic_factor(ctx.base, self, stride) {
                debug!("  {} is divisible by {}", self, factors.iter().format(", "));
                debug_to_tree!(ctx.tracer, "Divisible by {}", factors.iter().format(","));
                return true;
            }
        }
        if let Some((even_factor, odd_factor)) = find_two_factors(ctx.base, self) {
            debug!("  {self} is divisible by either {even_factor} or {odd_factor} (#1)");
            debug_to_tree!(
                ctx.tracer,
                "Divisible by either {} or {} (#1)",
                even_factor,
                odd_factor
            );
            return true;
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(ctx.base, self) {
            debug!("  {self} is divisible by either {even_factor} or {odd_factor} (#2)");
            debug_to_tree!(
                ctx.tracer,
                "Divisible by either {} or {} (#2)",
                even_factor,
                odd_factor
            );
            return true;
        }

        if check_residues_mod_30(ctx.base, self) {
            debug!("  {self} always shares a factor with 30");
            debug_to_tree!(ctx.tracer, "Always shares a factor with 30");
            return true;
        }

        false
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
