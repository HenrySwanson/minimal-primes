use std::time::Duration;

use num_prime::buffer::NaiveBuffer;

use crate::candidates::CandidateSequences;
use crate::logging::Tracer;

pub struct SearchContext {
    pub base: u8,

    /// iteration counter; counts how many families we've looked at
    pub iter: usize,
    /// primes we've discovered so far, in two different formats
    /// depending on the order we discovered these, they may not be minimal!
    pub primes: CandidateSequences,

    /// For primality testing. Re-using this avoids unnecessary re-computation
    /// of primes.
    /// TODO: re-use in sieving too!
    pub prime_buffer: NaiveBuffer,
    /// For potentially getting insight into what's going on
    pub stats: Stats,
    /// For tracking our paths through the search space in a more
    /// understandable format.
    pub tracer: Tracer,
}

#[derive(Debug, Default)]
pub struct Stats {
    pub num_primality_checks: usize,
    pub duration_primality_checks: Duration,
    pub num_substring_checks: usize,
    pub duration_substring_checks: Duration,
    pub num_simple_substring_checks: usize,
    pub duration_simple_substring_checks: Duration,
    pub num_could_contains: usize,
    pub duration_could_contains: Duration,
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

impl SearchContext {
    pub fn new(base: u8, tree_log: bool) -> Self {
        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            prime_buffer: NaiveBuffer::new(),
            stats: Stats::default(),
            tracer: if tree_log {
                Tracer::new()
            } else {
                Tracer::dummy()
            },
        }
    }
}
