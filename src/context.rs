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
    pub stats: SearchStats,
    /// For tracking our paths through the search space in a more
    /// understandable format.
    pub tracer: Tracer,
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

impl SearchContext {
    pub fn new(base: u8, tree_log: bool) -> Self {
        Self {
            base,
            iter: 0,
            primes: CandidateSequences::new(),
            prime_buffer: NaiveBuffer::new(),
            stats: SearchStats::default(),
            tracer: if tree_log {
                Tracer::new()
            } else {
                Tracer::dummy()
            },
        }
    }
}

pub fn print_stats(stats: &SearchStats) {
    println!("{} branches explored", stats.num_branches_explored);
    println!(
        "{} primality tests",
        stats.num_primality_checks,
    );
    println!(
        "{} calls Family::could_contain",
        stats.num_could_contains,
    );
    println!(
        "{} substring tests",
        stats.num_substring_checks,
    );
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
