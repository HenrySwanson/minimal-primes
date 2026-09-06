mod composite;
mod context;
mod frontier;
mod gcd;
mod split;
mod trace;

use std::ops::ControlFlow;

use itertools::Itertools;
use log::{debug, trace};
use num_bigint::BigUint;
use num_prime::buffer::PrimeBufferExt;

use self::composite::{find_even_odd_factor, find_periodic_factor, shares_factor_with_base};
pub use self::context::{
    CompositeReason, ExploreEvent, SearchContext, SplitDirection, print_stats,
};
use self::frontier::{Frontier, Weight};
pub use self::trace::TraceRecord;
use crate::RemainingNodes;
use crate::candidates::CandidateIndices;
use crate::digits::DigitSeq;
use crate::families::{Core, Family, SimpleFamily};
use crate::search::composite::{
    check_residues_mod_30, composite_checks_for_simple, find_common_factor, find_two_factors,
};

pub struct SearchTree {
    pub nodes: Frontier<SearchNode>,
}

impl SearchTree {
    pub fn new(ctx: &mut SearchContext) -> Self {
        let initial_node = SearchNode {
            id: ctx.alloc_node_id(),
            parent_id: None,
            node_type: NodeType::Arbitrary(FamilyNode {
                family: Family::any(ctx.base),
                possible_contained_primes: ctx.primes.indices_all(),
            }),
        };
        let frontier = Frontier::start(initial_node);

        Self { nodes: frontier }
    }

    pub fn num_nodes_to_solve(&self) -> usize {
        self.nodes
            .iter()
            .filter(|node| match &node.node_type {
                NodeType::Arbitrary(_) => true,
                NodeType::Simple(simple_node) => !simple_node.composite_tested,
            })
            .count()
    }

    pub fn any_nodes_to_solve(&self) -> bool {
        self.nodes.iter().any(|node| match &node.node_type {
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
            match node.node_type {
                NodeType::Arbitrary(node) => ret.other_families.push(node.family),
                // TODO: return the whole node, so the other stages can benefit here!
                NodeType::Simple(node) => ret.simple_families.push(node.family),
            }
        }

        ret
    }
}

#[derive(Debug, Clone)]
pub struct SearchNode {
    id: u64,
    parent_id: Option<u64>,
    node_type: NodeType,
}

#[derive(Debug, Clone)]
enum NodeType {
    Arbitrary(FamilyNode),
    Simple(SimpleNode),
}

#[derive(Debug, Clone)]
struct FamilyNode {
    family: Family,
    possible_contained_primes: CandidateIndices,
}

#[derive(Debug, Clone)]
struct SimpleNode {
    family: SimpleFamily,
    composite_tested: bool,
    /// When min_repeats equals this number, we can delete this
    /// family, because it contains this prime.
    dies_at: DiesAt,
    /// What the index of the first prime we've never seen is
    start_unknown_primes: usize,
}

impl SearchNode {
    fn explore(self, ctx: &mut SearchContext) -> Vec<SearchNode> {
        let node_id = self.id;
        let parent_id = self.parent_id;

        // We're about to relinquish our ownership of "family", so if tracing
        // is enabled, we should call our display method now before we lose
        // the opportunity.
        let family_display = ctx.trace.is_some().then(|| self.node_type.to_string());

        // Say our family is xL*z.
        // We want to explore all possible children with weight one more than this one.
        let (children, event) = match self.node_type {
            NodeType::Arbitrary(node) => {
                debug!(" Exploring {}", node.family);
                node.explore(ctx)
            }
            NodeType::Simple(node) => {
                debug!(" Exploring simple {}", node.family);
                node.explore(ctx)
            }
        };

        // Log what kind of split or discard occured.
        ctx.stats.num_branches_explored += 1;
        ctx.stats.branch_stats.record(&event);

        if let (Some(trace), Some(family)) = (&mut ctx.trace, family_display) {
            trace.record(node_id, parent_id, family, event);
        }

        // Wrap the child nodes into proper SearchNodes and return them
        children
            .into_iter()
            .map(|node_type| SearchNode {
                id: ctx.alloc_node_id(),
                parent_id: Some(node_id),
                node_type,
            })
            .collect()
    }
}

impl FamilyNode {
    fn explore(mut self, ctx: &mut SearchContext) -> (Vec<NodeType>, ExploreEvent) {
        // We have to be careful about not generating leading zeros. There's
        // a lot of different ways we could do that, but one possible way is
        // just to rig things so that on the first round, we split left (and
        // check for primes).
        // TODO: can we do this elsewhere? maybe some post-processing step in
        // the explore_node? this feels iffy

        // Now is a good time for us to narrow down the potential primes this
        // family could contain.
        let mut new_contained_primes = ctx.primes.indices_none();
        for (i, prime) in ctx.primes.get_many(&self.possible_contained_primes) {
            if self.family.could_contain(prime) {
                new_contained_primes.add(i);
            }
            ctx.stats.num_could_contains += 1;
        }
        self.possible_contained_primes = new_contained_primes;

        // Test this for primality
        // TODO: normally we've tested this already, in reduce_cores,
        // but split_on_repeat can produce strings we've never tested :/
        // What's a better way to avoid this redundancy?
        let seq = self.family.contract();
        if let Some(p) = ctx.test_for_contained_prime(&seq, &self.possible_contained_primes) {
            assert_ne!(seq, *p);
            debug!("  Discarding {}, contains prime {}", self.family, p);
            return (vec![], ExploreEvent::ContainsPrime(p.clone()));
        }

        trace!("  Testing for primality {seq}");
        let value = seq.value(ctx.base);
        if ctx.test_for_prime(&value) {
            debug!("  Saving {}, contracts to prime", self.family);
            ctx.primes.insert(seq);
            return (vec![], ExploreEvent::IsNewPrime);
        }

        // Then, we try to reduce the cores, notifying the tracer if that changed
        // anything.
        let mut anything_changed = false;
        anything_changed |= self.reduce_cores(ctx);
        anything_changed |= self.family.simplify();
        if anything_changed && let Some(trace) = &mut ctx.trace {
            trace.set_reduction(self.family.to_string());
        }
        if self.family.cores.is_empty() {
            debug!("  {} was reduced to trivial string", self.family);
            return (vec![], ExploreEvent::NoCoresRemaining);
        }

        // Now, run some tests to see whether this family is guaranteed to
        // be composite.
        if let Some(reason) = self.family.test_for_perpetual_composite(ctx) {
            debug!("  Discarding {}, is always composite", self.family);
            return (vec![], ExploreEvent::DetectedComposite(reason));
        }

        // TODO: is this right?
        // Check if this family is simple or not. If it is, we should
        // re-enqueue it as such. (Note: this is after composite check!)
        if let Ok(family) = SimpleFamily::try_from(self.family.clone()) {
            // Take all the primes we know could be contained in this family,
            // and check exactly when this family meets them.
            let mut dies_at = DiesAt::Unknown;
            for (_, p) in ctx.primes.get_many(&self.possible_contained_primes) {
                if let Some(n) = family.will_contain_at(p) {
                    dies_at.update(n, p);
                }
            }

            return (
                vec![NodeType::Simple(SimpleNode {
                    family,
                    composite_tested: false,
                    dies_at,
                    start_unknown_primes: ctx.primes.upper_bound(),
                })],
                ExploreEvent::Simplified,
            );
        }

        // Let's see if we can split it in an interesting way
        // TODO: context-ify the splitting functions too!
        if self.family.weight() >= 2
            && let Some((children, event)) =
                ctx.split_on_limited_digit(&self.family, 3, &self.possible_contained_primes)
        {
            let children = children
                .into_iter()
                .map(|family| {
                    NodeType::Arbitrary(FamilyNode {
                        family,
                        possible_contained_primes: self.possible_contained_primes.clone(),
                    })
                })
                .collect();
            return (children, event);
        }

        if self.family.weight() >= 4 {
            if let Some((children, event)) = ctx.split_on_incompatible_digits_different_cores(
                &self.family,
                &self.possible_contained_primes,
            ) {
                let children = children
                    .into_iter()
                    .map(|family| {
                        NodeType::Arbitrary(FamilyNode {
                            family,
                            possible_contained_primes: self.possible_contained_primes.clone(),
                        })
                    })
                    .collect();
                return (children, event);
            }

            if let Some((children, event)) =
                ctx.split_on_incompatible_digits(&self.family, &self.possible_contained_primes)
            {
                let children = children
                    .into_iter()
                    .map(|family| {
                        NodeType::Arbitrary(FamilyNode {
                            family,
                            possible_contained_primes: self.possible_contained_primes.clone(),
                        })
                    })
                    .collect();
                return (children, event);
            }
        }

        // this was introduced to kill long derivation chains of the form x[ab]*y that
        // we have trouble with otherwise. only invoke it when we are really stuck on
        // something.
        if self.family.weight() >= 10
            && let Some((children, event)) =
                ctx.split_on_forbidden_sandwich(&self.family, &self.possible_contained_primes)
        {
            let children = children
                .into_iter()
                .map(|family| {
                    NodeType::Arbitrary(FamilyNode {
                        family,
                        possible_contained_primes: self.possible_contained_primes.clone(),
                    })
                })
                .collect();
            return (children, event);
        }

        if self.family.weight() >= 5 {
            // This one doesn't actually simplify any cores, in fact, it'll make the
            // branch more complicated!. However, it might kickstart some branch elimination
            // by making us intersect another prime. So this check should always be last.
            if let Some((child, event)) = ctx.split_on_necessary_digit(&self.family) {
                let children = vec![NodeType::Arbitrary(FamilyNode {
                    family: child,
                    possible_contained_primes: self.possible_contained_primes,
                })];
                return (children, event);
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
        let magic = match self.family.weight() {
            0 => 0,
            1 => 1,
            w => bit_mixer(w),
        };

        let slot = (magic >> 1) % self.family.cores.len();
        debug_assert!(!self.family.cores[slot].is_empty());
        let direction = if magic % 2 == 1 {
            SplitDirection::Left
        } else {
            SplitDirection::Right
        };
        let mut children = match direction {
            SplitDirection::Left => {
                debug!("  Splitting {} left on core {slot}", self.family);
                self.family.expand_left(slot)
            }
            SplitDirection::Right => {
                debug!("  Splitting {} right on core {slot}", self.family);
                self.family.expand_right(slot)
            }
        };

        // We also need to consider the case where the chosen core expands to
        // the empty string. However, in the case where there's one core, this
        // is pretty redundant with the work we're doing in reduce_core().
        // For example: if we reduce a[xyz]c, we test the primality of axc, ayc
        // and azc. So after we split, and get ax[xyz]c, there's no need to
        // test ax[]c again.
        if self.family.cores.len() > 1 {
            self.family.cores[slot].clear();
            self.family.simplify();
            children.push(self.family);
        }

        let children = children
            .into_iter()
            .map(|family| {
                NodeType::Arbitrary(FamilyNode {
                    family,
                    possible_contained_primes: self.possible_contained_primes.clone(),
                })
            })
            .collect();
        (
            children,
            ExploreEvent::SplitGenerically {
                core_idx: slot,
                direction,
            },
        )
    }
}

impl SimpleNode {
    fn explore(mut self, ctx: &mut SearchContext) -> (Vec<NodeType>, ExploreEvent) {
        // There's a lot less we can do here! We can't split anything,
        // we can't reduce cores, etc, etc.

        // We should do a composite test though; there's some specialized
        // composite tests that we can only do on simple families. We should
        // only do them once though.
        if !self.composite_tested {
            if composite_checks_for_simple(ctx.base, &self) {
                debug!("  Discarding {}, is always composite", self.family);
                return (
                    vec![],
                    ExploreEvent::DetectedComposite(CompositeReason::FactorsAlgebraically),
                );
            }
            self.composite_tested = true;
        }

        // Other that that, all we can do is add another digit and see if it
        // becomes prime or not. Or contains another prime.

        // On a previous loop, we may have established when this family contains
        // a prime. Check it.
        if let DiesAt::KilledBy(dies_at, prime) = &self.dies_at
            && self.family.min_repeats >= *dies_at
        {
            debug!("  Discarding {}, contains prime {}", self.family, prime);
            return (vec![], ExploreEvent::ContainsPrime(prime.clone()));
        }

        // Now check any new primes.
        for (_, prime) in ctx.primes.get_tail(self.start_unknown_primes) {
            if let Some(n) = self.family.will_contain_at(prime) {
                if n <= self.family.min_repeats {
                    debug!("  Discarding {}, contains prime {}", self.family, prime);
                    return (vec![], ExploreEvent::ContainsPrime(prime.clone()));
                }

                // otherwise, we should incorporate this into dies_at
                self.dies_at.update(n, prime);
            }
            ctx.stats.num_simple_substring_checks += 1;
        }
        self.start_unknown_primes = ctx.primes.upper_bound(); // resets our collection

        // Test if it is a prime
        let value = self.family.value(ctx.base);

        if ctx.test_for_prime(&value) {
            debug!("  Saving {}, is prime", self.family);
            let seq = self.family.contract();
            ctx.primes.insert(seq);
            return (vec![], ExploreEvent::IsNewPrime);
        }

        self.family.min_repeats += 1;
        (
            vec![NodeType::Simple(self)],
            ExploreEvent::IncrementedRepeat,
        )
    }
}

impl FamilyNode {
    /// Removes digits from the family's cores that would cause the family to
    /// contain a known prime.
    ///
    /// Returns `true` if any cores were modified, otherwise `false`.
    fn reduce_cores(&mut self, ctx: &mut SearchContext) -> bool {
        let mut anything_changed = false;
        let old_family = self.family.clone();

        for (i, core) in self.family.cores.iter_mut().enumerate() {
            // Substitute elements from the core into the string to see if any
            // of them contain or are a prime.
            // NOTE: this is where we generate minimal primes of (weight + 1), so
            // next loop, those should all be available.
            let mut allowed_digits = vec![];
            for digit in core.iter() {
                let seq = old_family.substitute(i, digit);

                if let Some(p) = ctx.test_for_contained_prime(&seq, &self.possible_contained_primes)
                {
                    assert_ne!(seq, *p);
                    debug!("  Discarding {seq}, contains prime {p}");
                    anything_changed = true;
                    continue;
                }

                trace!("  Testing for primality {seq}");
                let value = seq.value(ctx.base);
                if ctx.test_for_prime(&value) {
                    debug!("  Saving {seq}, is prime");
                    ctx.primes.insert(seq);
                    anything_changed = true;
                } else {
                    allowed_digits.push(digit);
                }
            }

            *core = Core::new(allowed_digits);
        }
        // Now we've reduced the core, and have a new family.
        debug!("  Reducing {} to {}", old_family, self.family);
        anything_changed
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

        let (_, seq) = result?;
        Some(seq)
    }

    fn test_for_prime(&mut self, value: &BigUint) -> bool {
        let result = self.prime_buffer.is_prime(value, None).probably();
        self.stats.num_primality_checks += 1;
        result
    }
}

impl Family {
    fn test_for_perpetual_composite(&self, ctx: &mut SearchContext) -> Option<CompositeReason> {
        // This function is used to eliminate families that will always result
        // in composite numbers, letting us cut off infinite branches of the
        // search space.
        // There are a few possible ways this can happen. We'll use base 10
        // in the comments for familiarity, unless specified otherwise.

        // p divides BASE (e.g., 2, 5)
        if let Some(factor) = shares_factor_with_base(ctx.base, self) {
            debug!("  {self} has divisor {factor}");
            return Some(CompositeReason::SharesFactorWithBase(factor));
        }
        // p does not divide BASE (e.g. 7)
        // -------------------------------
        // This is how we detect families like 4[6]9 being divisible by 7.
        if let Some(factor) = find_common_factor(ctx.base, self) {
            debug!("  {self} is divisible by {factor}");
            return Some(CompositeReason::CommonFactor(factor));
        }
        // start at stride 2; `find_guaranteed_factor` effectively handles stride 1
        for stride in 2..=4 {
            if let Some(factors) = find_periodic_factor(ctx.base, self, stride) {
                debug!("  {} is divisible by {}", self, factors.iter().format(", "));
                return Some(CompositeReason::PeriodicFactors(factors));
            }
        }
        if let Some((core_idx, even_factor, odd_factor)) = find_two_factors(ctx.base, self) {
            debug!("  {self} is divisible by either {even_factor} or {odd_factor} (#1)");
            return Some(CompositeReason::LocalAlternatingFactors {
                core_idx,
                even_factor,
                odd_factor,
            });
        }

        if let Some((even_factor, odd_factor)) = find_even_odd_factor(ctx.base, self) {
            debug!("  {self} is divisible by either {even_factor} or {odd_factor} (#2)");
            return Some(CompositeReason::GlobalAlternatingFactors {
                even_factor,
                odd_factor,
            });
        }

        if check_residues_mod_30(ctx.base, self) {
            debug!("  {self} always shares a factor with 30");
            return Some(CompositeReason::NeverCoprimeTo30);
        }

        None
    }
}

impl std::fmt::Display for NodeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            NodeType::Arbitrary(node) => node.family.fmt(f),
            NodeType::Simple(node) => node.family.fmt(f),
        }
    }
}

impl std::fmt::Display for SearchNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.node_type.fmt(f)
    }
}

impl Weight for SearchNode {
    fn weight(&self) -> usize {
        match &self.node_type {
            NodeType::Arbitrary(node) => node.family.weight(),
            NodeType::Simple(node) => {
                node.family.bare.before.0.len()
                    + node.family.min_repeats
                    + node.family.bare.after.0.len()
            }
        }
    }
}

#[derive(Debug, Clone)]
pub enum DiesAt {
    Unknown,
    KilledBy(usize, DigitSeq),
}

impl DiesAt {
    pub fn update(&mut self, n: usize, p: &DigitSeq) {
        match self {
            Self::KilledBy(old_n, _) if *old_n <= n => {}
            Self::KilledBy(..) | Self::Unknown => *self = Self::KilledBy(n, p.clone()),
        }
    }
}
