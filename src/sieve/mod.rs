use std::ops::Range;
use std::time::Instant;

use log::debug;
use num_bigint::BigUint;
use num_prime::buffer::{NaiveBuffer, PrimeBufferExt};

use self::bsgs::sieve;
pub use self::context::SieveContext;
use self::sequence_slice::SequenceSlice;
pub use self::stats::{suggest_next_p_max, SieveStats};
use crate::digits::{Digit, DigitSeq};
use crate::families::SimpleFamily;
use crate::sequence::Sequence;

mod bsgs;
mod context;
mod sequence_slice;
mod stats;

/// Entry point for eliminating simple families through sieving.
pub fn do_one_round(
    ctx: &mut SieveContext,
    remaining_branches: &mut Vec<(SimpleFamily, Sequence)>,
    n_range: &Range<usize>,
    p_max: u64,
) -> SieveStats {
    let base = ctx.base;
    let mut stats = SieveStats::default();

    let mut sequences_to_sieve = vec![];
    let mut slices_to_sieve = vec![];

    // Real quick, check if this can be eliminated via a minimal prime
    // TODO: shouldn't this come _after_ sieving?
    for (simple, seq) in std::mem::take(remaining_branches) {
        if let Some(p) = ctx
            .primes
            .iter()
            .find(|p| simple.will_contain_at(p).is_some_and(|n| n < n_range.start))
        {
            println!("{simple} can be eliminated, since it contains {p}");
            continue;
        }

        sequences_to_sieve.push(simple);
        slices_to_sieve.push(SequenceSlice::new(seq, n_range.clone()))
    }

    stats.num_slices = slices_to_sieve.len();
    let candidates_before: usize = slices_to_sieve.iter().map(|s| s.num_remaining()).sum();

    // Now sieve all these slices at once
    println!(
        "Sieving {} families for n from {} to {}",
        slices_to_sieve.len(),
        n_range.start,
        n_range.end,
    );
    sieve(
        base,
        &mut slices_to_sieve,
        p_max,
        &mut ctx.prime_buffer,
        &mut stats,
    );

    let candidates_after: usize = slices_to_sieve.iter().map(|s| s.num_remaining()).sum();
    stats.num_eliminated_by_sieve = candidates_before - candidates_after;
    stats.num_candidates_after_sieve = candidates_after;

    for (simple, slice) in std::iter::zip(sequences_to_sieve, slices_to_sieve) {
        // Iterate through the unmarked n and manually check primality
        println!(
            "Investigating the {}/{} terms remaining in {}",
            slice.num_remaining(),
            n_range.len(),
            simple
        );

        match last_resort(base, &slice, &mut ctx.prime_buffer, &mut stats) {
            Some((i, p)) => {
                let digitseq =
                    DigitSeq(p.to_radix_be(base.into()).into_iter().map(Digit).collect());
                println!("Found prime at exponent {i}: {digitseq}");
                ctx.primes.insert(digitseq);
                stats.num_primes_found += 1;
            }
            None => {
                println!("Unable to find prime in the given range: {simple}");
                remaining_branches.push((simple, slice.seq))
            }
        }
    }

    stats.print(n_range, p_max);
    stats
}

/// Small convenience function for sieving a single sequence.
pub fn find_first_prime(
    base: u8,
    k: u64,
    c: i64,
    d: u64,
    n_lo: usize,
    n_hi: usize,
    p_max: u64,
) -> Option<(usize, BigUint)> {
    let seq = Sequence::new(k, c, d);
    let slice = SequenceSlice::new(seq, n_lo..n_hi);

    let mut slices = [slice];
    let mut prime_buffer = NaiveBuffer::new();
    let mut stats = SieveStats::default();
    sieve(base, &mut slices, p_max, &mut prime_buffer, &mut stats);
    last_resort(base, &slices[0], &mut prime_buffer, &mut stats)
}

fn last_resort(
    base: u8,
    slice: &SequenceSlice,
    prime_buffer: &mut NaiveBuffer,
    stats: &mut SieveStats,
) -> Option<(usize, BigUint)> {
    for exponent in slice.iter_remaining() {
        let value = slice.seq.compute_term(exponent as u32, base.into());
        debug!("  Check {} at n={}", slice.seq, exponent);

        let start = Instant::now();
        let is_prime = prime_buffer.is_prime(&value, None).probably();
        stats.record_primality_test(start.elapsed());

        if is_prime {
            return Some((exponent, value));
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prime_finding() {
        // Let's test some sequences where we already know the answer :)
        // We're looking for n = (# digits - # of digits in the first part).

        // Base 17: A0*1 is first prime at 1357 digits.
        // Sequence is 10*17^n+1, n>=1
        let x = find_first_prime(17, 10, 1, 1, 1, 2000, 100_000);
        assert_eq!(x.unwrap().0, 1357 - 1);

        // Base 23: E0*KLE is first prime at 1658 digits.
        // Sequence is 14*23^n+11077, n>=3
        let x = find_first_prime(23, 14, 11077, 1, 3, 2000, 100_000);
        assert_eq!(x.unwrap().0, 1658 - 1);

        // Base 11: 44*1 is first prime at 45 digits.
        // Sequence is (44*b^n - 34)/d, n>=1
        let x = find_first_prime(11, 44, -34, 10, 1, 100, 1000);
        assert_eq!(x.unwrap().0, 45 - 1);

        // Base 13: 80*111 is first prime at at 32021 digits.
        // Sequence is 8*13^n+183, n>=3
        // This takes too long.
        // let x = find_first_prime(13, 8, 183, 3, 40000).unwrap();
        // assert_eq!(x.0, 32021 - 1);
    }
}
