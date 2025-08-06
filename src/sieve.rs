use std::collections::HashMap;

use bitvec::prelude::BitVec;
use log::{debug, info};
use num_bigint::BigUint;
use num_modular::{ModularCoreOps, ModularPow, ModularUnaryOps};
use num_prime::buffer::NaiveBuffer;
use num_prime::nt_funcs::is_prime;

use crate::families::Sequence;

#[derive(Debug)]
pub struct SequenceSlice {
    pub seq: Sequence,
    n_lo: usize,
    n_bitvec: BitVec,
}

impl SequenceSlice {
    pub fn new(seq: Sequence, n_lo: usize, n_hi: usize) -> Self {
        assert!(n_lo <= n_hi);
        let n_range = n_hi - n_lo;

        Self {
            seq,
            n_lo,
            n_bitvec: BitVec::repeat(true, n_range),
        }
    }

    pub fn n_lo(&self) -> usize {
        self.n_lo
    }

    pub fn n_hi(&self) -> usize {
        self.n_lo + self.n_bitvec.len()
    }

    pub fn check_n(&self, n: usize) -> bool {
        self.n_bitvec[n - self.n_lo]
    }

    pub fn eliminate_n(&mut self, n: usize) {
        self.n_bitvec.set(n - self.n_lo, false);
    }

    pub fn eliminate_multiple(&mut self, p: u64, base: u64, start: usize, spacing: usize) {
        let mut idx = start - self.n_lo;

        // It's possible the term we're about to eliminate is actually p itself.
        // Let's avoid that, if so.
        // TODO: actually, this should report a prime immediately! not that that'll
        // happen in the interesting cases, but still!
        if self.seq.check_term_equal(base, p, start) {
            idx += spacing;
        }

        // Insane edge case: it could also be zero! In that case, bump it up twice.
        if self.seq.check_term_equal(base, 0, start) {
            idx += 2 * spacing;
        }

        while let Some(mut slot) = self.n_bitvec.get_mut(idx) {
            slot.set(false);
            idx += spacing;
        }
    }
}

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
    let slice = SequenceSlice::new(seq, n_lo, n_hi);

    let mut slices = [slice];
    sieve(base, &mut slices, p_max);
    last_resort(base, &slices[0])
}

pub fn sieve(
    base: u8,
    slices: &mut [SequenceSlice],
    // TODO: how many? can i decide from "outside"?
    p_max: u64,
) {
    if slices.is_empty() {
        return;
    }

    // Decide how many steps for baby-step giant-step
    // TODO: i think i can make it all the same range?
    let n_range = slices
        .iter()
        .map(|slice| slice.n_bitvec.len())
        .max()
        .unwrap();
    let num_baby_steps = (n_range as f64).sqrt() as usize;
    let num_giant_steps = n_range.div_ceil(num_baby_steps);

    // Now go and eliminate a bunch of terms
    let mut prime_buffer = NaiveBuffer::new();
    for p in prime_buffer.primes(p_max) {
        baby_step_giant_step(base.into(), *p, num_baby_steps, num_giant_steps, slices);
    }
}

pub fn last_resort(base: u8, slice: &SequenceSlice) -> Option<(usize, BigUint)> {
    // Lastly, iterate through the remaining numbers and see if they're prime
    info!(
        "{} has {}/{} terms remaining",
        slice.seq,
        slice.n_bitvec.count_ones(),
        slice.n_bitvec.len()
    );
    for i in slice.n_bitvec.iter_ones() {
        let exponent = slice.n_lo + i;

        let value = slice.seq.compute_term(exponent as u32, base.into());
        debug!("  Check {} at n={}", slice.seq, exponent);

        if is_prime(&value, None).probably() {
            return Some((exponent, value));
        }
    }

    None
}

fn baby_step_giant_step(
    base: u64,
    p: u64,
    num_baby_steps: usize,
    num_giant_steps: usize,
    slices: &mut [SequenceSlice],
) {
    // Now that we're dealing with multiple simultaneous sequences, we may need
    // to skip over some of them. We do so with this vector.
    let mut skip = vec![false; slices.len()];

    // This works by solving b^n = (-c/k) mod p for n, using baby-step-giant-step.

    // What about d?
    // If p doesn't divide d, then we don't have to worry; p divides (kb^n+c)/d exactly
    // when it divides kb^n+c.
    // If p does divide d, then we have to be more careful, and count the number of ps.
    // For now though, we just skip that prime for that sequence! (TODO)

    // Compute some inverses!
    let binv = match base.invm(&p) {
        Some(x) => x,
        None => {
            // If p divides b, then the term will be equivalent to c mod p.
            // If c is zero, this is always divisible by p, otherwise, it never
            // is.
            // The former situation should never arise in practice.
            debug!("Completely skipping prime {p}, it divides the base b={base}");
            for slice in slices {
                assert_ne!(
                    slice.seq.c.unsigned_abs() % p,
                    0,
                    "Sequence {:?} is always divisible by {}",
                    slice.seq,
                    p
                );
            }
            return;
        }
    };

    // We want to compute -c/k for each of our sequences.
    let ck: Vec<_> = slices
        .iter()
        .enumerate()
        .map(|(i, slice)| {
            // Here is a convenient place to check d
            if
            /* p < seq.d && */
            slice.seq.d.is_multiple_of(p) {
                // TODO: log something
                skip[i] = true;
                return 0;
            }

            // If k is zero mod p, we have the same situation as with b, but
            // we only need to check one c.
            let kinv = match slice.seq.k.invm(&p) {
                Some(x) => x,
                None => {
                    assert_ne!(
                        slice.seq.c.unsigned_abs() % p,
                        0,
                        "Sequence {:?} is always divisible by {}",
                        slice.seq,
                        p
                    );
                    skip[i] = true;
                    return 0;
                }
            };

            // We first have to get c as a u64 before we can do mod-p math with it.
            let neg_c_mod_p = if slice.seq.c >= 0 {
                slice.seq.c.unsigned_abs().negm(&p)
            } else {
                slice.seq.c.unsigned_abs()
            };

            // -c/k mod p
            neg_c_mod_p.mulm(kinv, &p)
        })
        .collect();

    // Take some baby steps
    let (baby_table, order) = baby_steps(base, p, num_baby_steps, slices[0].n_lo);
    debug_assert!(
        baby_table.keys().all(|x| *x != 0),
        "should never see 0s in baby_table when b has an inverse"
    );

    // If we know the order, our baby table contains all powers of b.
    // If ck is in there, we can do some elimination.
    if let Some(order) = order {
        for (idx, slice) in slices.iter_mut().enumerate() {
            if skip[idx] {
                continue;
            }

            if let Some(i) = baby_table.get(&ck[idx]) {
                // eliminate L + i, and all multiples of order afterward
                slice.eliminate_multiple(p, base, slice.n_lo + i, order);
            }
        }

        return;
    }

    // Otherwise, we'll do some giant steps
    let m: u64 = num_baby_steps.try_into().unwrap();
    let bm = binv.powm(m, &p);
    let mut ckb = ck;

    // Along the way, we'll try to find the order.
    let mut order = None;
    let mut idx_of_first_solution: Option<usize> = None;

    // Also set an array for tracking the actual solutions we find.
    let mut solutions = vec![None; slices.len()];
    let mut n_solutions = 0;

    for i in 0..num_giant_steps {
        // One extra for the repeat solution
        if n_solutions == slices.len() + 1 {
            break;
        }

        // Check whether this step gave a solution for any of the slices
        for (idx, slice) in slices.iter().enumerate() {
            // Ignore this slice if we've already solved it (unless we're still looking
            // for the order). Or if we're just skipping it outright.
            if skip[idx] || (solutions[idx].is_some() && idx_of_first_solution != Some(idx)) {
                continue;
            }

            // See if we got a hit on ckb
            if let Some(j) = baby_table.get(&ckb[idx]) {
                // Found a solution! (-c/k)b^(im) = b^(L+j), so we eliminate L+im+j
                let exp = slice.n_lo + i * num_baby_steps + j;

                n_solutions += 1;

                // Is this a repeat solution? We have to be slightly different if so.
                if idx_of_first_solution == Some(idx) {
                    let old_soln = solutions[idx].expect("first solution");
                    assert!(exp > old_soln);
                    order = Some(exp - old_soln);
                    idx_of_first_solution = None; // don't need it anymore
                } else {
                    // Save it normally
                    solutions[idx] = Some(exp);
                    if n_solutions == 1 {
                        idx_of_first_solution = Some(idx);
                    }
                }
            }

            // Either way, bump ckb
            ckb[idx] = ckb[idx].mulm(bm, &p);
        }
    }

    // Okay, after all that, what do we have?
    // We have some solutions to some of the sequences, and hopefully we have an order.
    // If we don't, just use p-1.
    let order = order.unwrap_or((p - 1) as usize);
    for (idx, slice) in slices.iter_mut().enumerate() {
        if skip[idx] {
            continue;
        }

        if let Some(exp) = solutions[idx] {
            slice.eliminate_multiple(p, base, exp, order);
        }
    }
}

fn baby_steps(
    base: u64,
    p: u64,
    num_baby_steps: usize,
    start_exp: usize,
) -> (HashMap<u64, usize>, Option<usize>) {
    let start_exp: u64 = start_exp.try_into().unwrap();

    let mut map = HashMap::with_capacity(num_baby_steps);
    let initial_value = base.powm(start_exp, &p);
    let mut value = initial_value;
    for i in 0..num_baby_steps {
        map.insert(value, i);
        value = value.mulm(base, &p);

        // We've looped all the way around! No need to insert any more entries,
        // we can return with knowledge of the order.
        if value == initial_value {
            return (map, Some(i + 1));
        }
    }

    (map, None)
}

#[cfg(test)]
mod tests {
    use num_prime::buffer::PrimeBufferExt;

    use super::*;

    #[test]
    fn test_bsgs() {
        // Start at zero for ease of understanding
        // 5*2^n+1
        let base = 2;
        let n_range = 100;
        let max_p = 100;

        let seq = Sequence::new(5, 1, 1);
        let slice = SequenceSlice::new(seq, 0, n_range);
        let mut slices = [slice];
        let mut prime_buffer = NaiveBuffer::new();
        for p in prime_buffer.primes(max_p) {
            baby_step_giant_step(base, *p, 10, 10, &mut slices);
        }

        let slice = &slices[0];
        for i in 0..n_range {
            // Check every remaining element in the sequence
            let exp = slice.n_lo + i;
            let elt = slice.seq.compute_term(exp as u32, base);
            let is_remaining = slice.check_n(i);
            let is_prime = prime_buffer.is_prime(&elt, None).probably();

            // If it's prime, we must not eliminate it.
            if is_prime {
                assert!(
                    is_remaining,
                    "{elt} = {seq} is prime at n={exp}, but was removed from the list"
                );
            // If it's composite, we might have eliminated it. Specifically,
            // if it has small factors, we should have been able to eliminate it.
            } else if is_remaining {
                // Conversely, if we didn't eliminate it, it should not have small factors.
                let factors = prime_buffer.factorize(elt.clone());
                let min_factor = factors.keys().min().expect("at least one factor");
                assert!(
                    min_factor >= &max_p.into(),
                    "{} = {} at n={} has unexpected small factor {}",
                    elt,
                    seq,
                    slice.n_lo + i,
                    min_factor
                );
            }
        }

        // Also check we eliminated a substantial number of them
        let num_remaining = slice.n_bitvec.count_ones();
        assert!(
            num_remaining < 20,
            "Expected to eliminate more options, there are {num_remaining} remaining"
        );
    }

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
