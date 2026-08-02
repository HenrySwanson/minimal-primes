use num_prime::buffer::NaiveBuffer;

use crate::candidates::CandidateSequences;
use crate::digits::{Digit, DigitSeq};

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
}

#[derive(Debug, Default)]
pub struct SearchStats {
    pub num_primality_checks: usize,
    pub num_substring_checks: usize,
    pub num_simple_substring_checks: usize,
    pub num_could_contains: usize,
    pub num_branches_explored: usize,
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
            ExploreEvent::DetectedComposite => self.detected_composite += 1,
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
#[derive(Debug, Clone, Copy)]
pub enum SplitDirection {
    Left,
    Right,
}

/// Describes what happened to a family (simple or otherwise) when we explored
/// it (i.e., expanded it into a (possibly empty) list of children).
///
/// This is what powers the [BranchStats] and also eventually the tree-based
/// logging that'll let me trace through the search tree after-the-fact.
#[derive(Debug, Clone)]
#[expect(dead_code)] // TODO: remove this when we use the bodies in the tree tracer thing
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
    // TODO: record the reason for compositeness here too?
    DetectedComposite,
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

impl SearchContext {
    pub fn new(base: u8) -> Self {
        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            prime_buffer: NaiveBuffer::new(),
            stats: SearchStats::default(),
        }
    }
}

pub fn print_stats(stats: &SearchStats) {
    println!("{} branches explored", stats.num_branches_explored);
    println!("{} primality tests", stats.num_primality_checks,);
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
