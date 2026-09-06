use num_bigint::BigUint;
use num_prime::buffer::NaiveBuffer;
use serde::{Deserialize, Serialize};

use crate::candidates::CandidateSequences;
use crate::digits::{Digit, DigitSeq};
use crate::search::trace::TraceWriter;

pub struct SearchContext {
    pub base: u8,

    /// iteration counter; counts how many families we've looked at
    pub iter: usize,
    /// primes we've discovered so far, in two different formats
    /// depending on the order we discovered these, they may not be minimal!
    pub primes: CandidateSequences,

    /// For primality testing. Re-using this avoids unnecessary re-computation
    /// of primes.
    pub prime_buffer: NaiveBuffer,
    /// For potentially getting insight into what's going on
    pub stats: SearchStats,

    /// Used for assigning unique IDs to search nodes. This is necessary for
    /// tree-tracing to work.
    next_node_id: u64,
    /// Where to log the tree-tracing events to, if at all.
    pub(super) trace: Option<TraceWriter>,
}

#[derive(Debug, Default)]
pub struct SearchStats {
    /// How many primality checks we ran.
    pub num_primality_checks: usize,
    /// How many primality checks were settled by the small-factor screen,
    /// without running a Miller-Rabin round.
    pub num_screened_out: usize,
    /// How many times we checked if a prime was contained in a digit sequence.
    pub num_substring_checks: usize,
    /// How many times we checked if a prime was contained in a `SimpleFamily`.
    pub num_simple_substring_checks: usize,
    /// How many times we checked if a prime could be contained in `Family`.
    pub num_could_contains: usize,
    /// How many branches of the search tree we explored.
    pub num_branches_explored: usize,
    /// What happened to those branches.
    pub branch_stats: BranchStats,
}

#[derive(Debug, Default)]
pub struct BranchStats {
    pub leading_zeros: usize,
    pub contains_prime: usize,
    pub is_new_prime: usize,
    pub is_trivial_string: usize,
    pub detected_composite: usize,
    pub simplified: usize,
    pub split_on_limited_digit: usize,
    pub split_on_incompatible_same_core: usize,
    pub split_on_incompatible_different_cores: usize,
    pub split_on_forbidden_sandwich: usize,
    pub split_on_necessary_digit: usize,
    pub explored_generically: usize,
}

impl BranchStats {
    pub fn record(&mut self, event: &ExploreEvent) {
        match event {
            ExploreEvent::ContainsPrime(_) => self.contains_prime += 1,
            ExploreEvent::IsNewPrime => self.is_new_prime += 1,
            ExploreEvent::NoCoresRemaining => self.is_trivial_string += 1,
            ExploreEvent::DetectedComposite(_) => self.detected_composite += 1,
            ExploreEvent::Simplified => self.simplified += 1,
            ExploreEvent::SplitOnLimitedDigit { .. } => self.split_on_limited_digit += 1,
            ExploreEvent::SplitOnIncompatibleDifferentCores { .. } => {
                self.split_on_incompatible_different_cores += 1
            }
            ExploreEvent::SplitOnIncompatibleSameCore { .. } => {
                self.split_on_incompatible_same_core += 1
            }
            ExploreEvent::SplitOnForbiddenSandwich { .. } => self.split_on_forbidden_sandwich += 1,
            ExploreEvent::SplitOnNecessaryDigit { .. } => self.split_on_necessary_digit += 1,
            ExploreEvent::SplitGenerically { .. } | ExploreEvent::IncrementedRepeat => {
                self.explored_generically += 1
            }
        }
    }
}

/// Which way [ExploreEvent::SplitArbitrarily] expanded the chosen core.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum SplitDirection {
    Left,
    Right,
}

/// Describes what happened to a family (simple or otherwise) when we explored
/// it (i.e., expanded it into a (possibly empty) list of children).
///
/// This is what powers the [BranchStats] and also the tree-based trace
/// logging in [crate::search::trace] that lets us reconstruct the search
/// tree after the fact.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExploreEvent {
    /// This family's smallest member contains an already-known prime, so every
    /// member of the family does too. Eliminates the branch.
    ContainsPrime(DigitSeq),
    /// The family's smallest member is itself a new minimal prime. Eliminates
    /// the branch.
    IsNewPrime,
    /// All cores reduced away, leaving a fixed string with no possible children.
    /// Eliminates the branch.
    NoCoresRemaining,
    /// One of the compositeness lemmas proved every member is composite. Eliminates
    /// the branch.
    DetectedComposite(CompositeReason),
    /// The family has only one core with only one digit, so it was converted into
    /// a [crate::families::SimpleFamily].
    Simplified,
    /// Lemma 21: `digit`, repeated `n` times in `core_idx`, is forbidden.
    SplitOnLimitedDigit {
        core_idx: usize,
        digit: Digit,
        n: usize,
    },
    /// Lemma 27: `a` (from `core_i`) and `b` (from `core_j`) can't co-occur.
    SplitOnIncompatibleDifferentCores {
        core_i: usize,
        core_j: usize,
        a: Digit,
        b: Digit,
    },
    /// Lemmas 23/25: `ab` can't appear in `core_idx`. If `reverse_also_forbidden`
    /// is set, then `ba` also is forbidden (in other words, `a` and `b` can't
    /// co-occur at all in that core).
    SplitOnIncompatibleSameCore {
        core_idx: usize,
        first: Digit,
        second: Digit,
        reverse_also_forbidden: bool,
    },
    /// Lemma 31: the "sandwich" pattern `aba` is forbidden in `core_idx`.
    SplitOnForbiddenSandwich { core_idx: usize, a: Digit, b: Digit },
    /// Lemma 29: `digit` is required to appear at least once in `core_idx`.
    SplitOnNecessaryDigit { core_idx: usize, digit: Digit },
    /// No lemma applied; split arbitrarily left/right on `core_idx`.
    SplitGenerically {
        core_idx: usize,
        direction: SplitDirection,
    },
    /// A [crate::families::SimpleFamily] didn't die and no lemma applied;
    /// just incremented `min_repeats` and kept going.
    IncrementedRepeat,
}

/// Describes why a family was proven to be composite.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CompositeReason {
    /// The family shares a factor with the base.
    SharesFactorWithBase(u8),
    /// All members of the family are divisible by some shared factor.
    CommonFactor(BigUint),
    /// The one-core family has a periodic factor sequence; the nth member
    /// of the family is divisible by the nth element in the sequence.
    PeriodicFactors(Vec<BigUint>),
    /// This one's a tricky one to explain, but given some core `m`, members
    /// of the family where `m` is expanded an even number of times all have
    /// one factor, and members where it's expanded an odd number of times
    /// have another.
    LocalAlternatingFactors {
        core_idx: usize,
        even_factor: BigUint,
        odd_factor: BigUint,
    },
    /// Members of the family are divisible by one of the given factors,
    /// depending on whether the total number of digits contributed by the
    /// cores is odd or even.
    GlobalAlternatingFactors {
        even_factor: BigUint,
        odd_factor: BigUint,
    },
    /// Members of the family are never coprime to 30.
    NeverCoprimeTo30,
    /// Each member of the family factors as something like a sum of cubes,
    /// difference of squares, etc.
    FactorsAlgebraically,
}

impl SearchContext {
    /// Creates some fresh context for the given `base`.
    ///
    /// If `trace` is set, creates a fresh trace file under `results/` and
    /// logs every explored node's event to it.
    pub fn new(base: u8, trace: bool) -> Self {
        let trace = trace.then(|| TraceWriter::create(base).expect("failed to create trace file"));

        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            prime_buffer: NaiveBuffer::new(),
            stats: SearchStats::default(),
            next_node_id: 0,
            trace,
        }
    }

    /// Allocates a fresh id for a newly created search node.
    pub fn alloc_node_id(&mut self) -> u64 {
        let id = self.next_node_id;
        self.next_node_id += 1;
        id
    }
}

pub fn print_stats(stats: &SearchStats) {
    println!("{} branches explored", stats.num_branches_explored);
    println!(
        "{} primality tests ({} settled by the small-factor screen)",
        stats.num_primality_checks, stats.num_screened_out,
    );
    println!("{} calls Family::could_contain", stats.num_could_contains,);
    println!("{} substring tests", stats.num_substring_checks,);
    println!(
        "{} simple substring tests",
        stats.num_simple_substring_checks,
    );
    let branch_stats = &stats.branch_stats;
    println!(
        "{} branches eliminated with leading zeros",
        branch_stats.leading_zeros
    );
    println!(
        "{} branches eliminated for containing a prime",
        branch_stats.contains_prime
    );
    println!(
        "{} branches eliminated by discovering a new prime",
        branch_stats.is_new_prime
    );
    println!(
        "{} branches eliminated by reducing to trivial string",
        branch_stats.is_trivial_string
    );
    println!(
        "{} branches eliminated for compositeness",
        branch_stats.detected_composite
    );
    println!(
        "{} branches simplified into simple families",
        branch_stats.simplified
    );
    println!(
        "{} branches split on a limited digit",
        branch_stats.split_on_limited_digit
    );
    println!(
        "{} branches split on incompatible digits (same core)",
        branch_stats.split_on_incompatible_same_core
    );
    println!(
        "{} branches split on incompatible digits (different cores)",
        branch_stats.split_on_incompatible_different_cores
    );
    println!(
        "{} branches split on a forbidden sandwich",
        branch_stats.split_on_forbidden_sandwich
    );
    println!(
        "{} branches split on a necessary digit",
        branch_stats.split_on_necessary_digit
    );
    println!(
        "{} branches explored generically",
        branch_stats.explored_generically
    );
}
