use std::ops::Range;

use bitvec::prelude::BitVec;
use log::debug;
use num_bigint::BigUint;
use num_modular::{ModularCoreOps, ModularPow, ModularUnaryOps};
use num_prime::buffer::{NaiveBuffer, PrimeBufferExt};

use crate::context::SearchContext;
use crate::digits::{Digit, DigitSeq};
use crate::families::SimpleFamily;
use crate::sequence::Sequence;

/// Entry point for eliminating simple families through sieving.
pub fn do_one_round(
    ctx: &mut SearchContext,
    remaining_branches: &mut Vec<(SimpleFamily, Sequence)>,
    n_range: &Range<usize>,
    p_max: u64,
) {
    let base = ctx.base;

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

    // Now sieve all these slices at once
    println!(
        "Sieving {} families for n from {} to {}",
        slices_to_sieve.len(),
        n_range.start,
        n_range.end,
    );
    sieve(base, &mut slices_to_sieve, p_max, &mut ctx.prime_buffer);

    for (simple, slice) in std::iter::zip(sequences_to_sieve, slices_to_sieve) {
        // Iterate through the unmarked n and manually check primality
        println!(
            "Investigating the {}/{} terms remaining in {}",
            slice.num_remaining(),
            n_range.len(),
            simple
        );

        match last_resort(base, &slice, &mut ctx.prime_buffer) {
            Some((i, p)) => {
                let digitseq =
                    DigitSeq(p.to_radix_be(base.into()).into_iter().map(Digit).collect());
                println!("Found prime at exponent {i}: {digitseq}");
                ctx.primes.insert(digitseq);
            }
            None => {
                println!("Unable to find prime in the given range: {simple}");
                remaining_branches.push((simple, slice.seq))
            }
        }
    }
}

#[derive(Debug)]
pub struct SequenceSlice {
    pub seq: Sequence,
    n_lo: usize,
    n_bitvec: BitVec,
}

impl SequenceSlice {
    pub fn new(seq: Sequence, range: Range<usize>) -> Self {
        Self {
            seq,
            n_lo: range.start,
            n_bitvec: BitVec::repeat(true, range.len()),
        }
    }

    #[cfg(test)]
    pub fn check_n(&self, n: usize) -> bool {
        self.n_bitvec[n - self.n_lo]
    }

    pub fn num_remaining(&self) -> usize {
        self.n_bitvec.count_ones()
    }

    pub fn iter_remaining(&self) -> impl Iterator<Item = usize> + use<'_> {
        self.n_bitvec.iter_ones().map(|i| self.n_lo + i)
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
        else if self.seq.check_term_equal(base, 0, start) {
            idx += 2 * spacing;
        }
        // Can it be negative? No, that would not be meaningful for the kinds
        // of sequences we're considering.
        if self.seq.c < 0 {
            // 0th term is (k*1+c) / d, which can only go negative if c is large and negative
            debug_assert!(self.seq.k > self.seq.c.unsigned_abs())
        }

        while let Some(mut slot) = self.n_bitvec.get_mut(idx) {
            slot.set(false);
            idx += spacing;
        }
    }
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
    sieve(base, &mut slices, p_max, &mut prime_buffer);
    last_resort(base, &slices[0], &mut prime_buffer)
}

fn sieve(
    base: u8,
    slices: &mut [SequenceSlice],
    // TODO: how many? can i decide from "outside"?
    p_max: u64,
    prime_buffer: &mut NaiveBuffer,
) {
    // The modular arithmetic in baby_step_giant_step reduces every value mod
    // p before operating on it, which lets it stay in u32 (doubling to u64
    // internally) instead of num_modular's default u64->u128 widening (which
    // triggers a slow software 128-bit division on every step). That only
    // works if p itself fits in a u32, which holds for any realistic p_max.
    assert!(
        p_max <= u32::MAX as u64,
        "sieve() requires p_max to fit in a u32 (got {p_max})"
    );

    // Decide how many steps for baby-step giant-step
    // TODO: will the input slices have different sizes?
    let Some(n_range) = slices.iter().map(|slice| slice.n_bitvec.len()).max() else {
        // no slices; return immediately
        return;
    };

    // TODO: with multiple sequences, explore different ways to divvy up
    // baby and giant steps
    let num_baby_steps = n_range.isqrt();
    let num_giant_steps = n_range.div_ceil(num_baby_steps);

    // The baby-step table is effectively a hashmap, but actually using one is
    // not the fastest choice. These are only sqrt(n_range) in size, so they're
    // pretty small, and since we're populating it for each prime up to p_max,
    // we really want to re-use our storage. Let's just use a sorted vector.
    let mut baby_table: Vec<(u32, usize)> = Vec::with_capacity(num_baby_steps);

    // Now go and eliminate a bunch of terms
    for p in prime_buffer.primes(p_max) {
        baby_step_giant_step(
            base.into(),
            *p as u32,
            num_baby_steps,
            num_giant_steps,
            slices,
            &mut baby_table,
        );
    }
}

/// Looks up `key` in a baby-step table sorted by [baby_steps].
fn lookup_baby_step(table: &[(u32, usize)], key: u32) -> Option<usize> {
    let index = table.binary_search_by_key(&key, |&(k, _)| k).ok()?;

    Some(table[index].1)
}

/// Inverts every value in `values` mod `p`, using just one modular inversion,
/// using Montgomery's batch inversion trick.
///
/// If a value is 0 it is skipped entirely and will not be modified.
fn batch_invm(values: Vec<u32>, p: u32) -> Vec<u32> {
    // prefix[i] is the product of values[0..i] mod p, so it starts with 1 and
    // ends with the full product (stuttering on any zeros along the way)
    let mut prefix = Vec::with_capacity(values.len() + 1);
    prefix.push(1);
    for v in values.iter().copied() {
        let last = *prefix.last().unwrap();
        if v == 0 {
            prefix.push(last);
        } else {
            prefix.push(last.mulm(v, &p));
        }
    }

    // Now do our one and only modular inversion. This gives us the inverse of
    // "everything multiplied together", i.e., (v_0 ... v_{n-1})^-1
    let Some(mut big_inverse) = prefix
        .last()
        .copied()
        .expect("prefix should have length > 0")
        .invm(&p)
    else {
        panic!("Somehow got something non-invertible; was {p} not prime? Values were: {values:?}")
    };

    // Now backtrack, filling in our `values` vector in reverse. (We can re-use that
    // buffer!)
    let mut result = values;
    for i in (0..result.len()).rev() {
        // loop invariant: big_inverse is the inverse of (v_0 ... v_i)
        
        let v_i = result[i];
        if v_i == 0 {
            continue; // don't do anything for zeros
        }
        
        // (v_0 ... v_i)^-1 * (v_0 ... v_{i-1}) = v_i^-1
        result[i] = big_inverse.mulm(prefix[i], &p);
        // maintain the invariant
        big_inverse = big_inverse.mulm(v_i, &p);
    }

    result
}

fn last_resort(
    base: u8,
    slice: &SequenceSlice,
    prime_buffer: &mut NaiveBuffer,
) -> Option<(usize, BigUint)> {
    for exponent in slice.iter_remaining() {
        let value = slice.seq.compute_term(exponent as u32, base.into());
        debug!("  Check {} at n={}", slice.seq, exponent);

        if prime_buffer.is_prime(&value, None).probably() {
            return Some((exponent, value));
        }
    }

    None
}

fn baby_step_giant_step(
    base: u64,
    p: u32,
    num_baby_steps: usize,
    num_giant_steps: usize,
    slices: &mut [SequenceSlice],
    baby_table: &mut Vec<(u32, usize)>,
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

    // All the modular arithmetic below operates on values already reduced
    // mod p, so it fits in u32 and doubles to u64 internally instead of
    // num_modular's default u64->u128 widening for u64 operands (which
    // triggers a slow software 128-bit division on every step).
    let p64 = u64::from(p);
    let base_mod_p = (base % p64) as u32;

    // Compute some inverses!
    let binv = match base_mod_p.invm(&p) {
        Some(x) => x,
        None => {
            // If p divides b, then the term will be equivalent to c mod p.
            // If c is zero, this is always divisible by p, otherwise, it never
            // is.
            // The former situation should never arise in practice.
            debug!("Completely skipping prime {p}, it divides the base b={base}");
            for slice in slices {
                assert_ne!(
                    slice.seq.c.unsigned_abs() % p64,
                    0,
                    "Sequence {:?} is always divisible by {}",
                    slice.seq,
                    p
                );
            }
            return;
        }
    };

    // Next, we need to compute -c/k for each of our sequences. However,
    // modular inversion is expensive, and we're doing a lot of it! We can
    // speed things up quite a bit with Montgomery's "batch-inversion" trick.

    // We want to compute -c/k for each of our sequences. Rather than doing a
    // full extended-Euclidean inversion of k for every slice (the dominant
    // cost of this whole function, since it happens once per (prime, slice)
    // pair), reduce k and c mod p for every eligible slice first, then invert
    // all the k's in one shot with Montgomery's batch-inversion trick: one
    // true modular inverse (of the product) plus O(#slices) multiplications,
    // instead of one inverse per slice.
    let mut k_mods = vec![0u32; slices.len()];
    let mut neg_c_mods = vec![0u32; slices.len()];
    for (i, slice) in slices.iter().enumerate() {
        // Here is a convenient place to check d
        if slice.seq.d.is_multiple_of(p64) {
            // TODO: log something
            skip[i] = true;
            continue;
        }

        k_mods[i] = (slice.seq.k % p64) as u32;
        neg_c_mods[i] = (slice.seq.c.unsigned_abs() % p64) as u32;
        if slice.seq.c > 0 {
            // note that c_mods is not fully reduced into [0, p)
            neg_c_mods[i] = neg_c_mods[i].negm(&p);
        }

        // If p divides k, then the term will always be equivalent to c mod p,
        // so we need to check if c is also divisible by p. If so, we should skip
        // this prime.
        if k_mods[i] == 0 {
            assert_ne!(
                neg_c_mods[i], 0,
                "Sequence {:?} is always divisible by {}",
                slice.seq, p
            );
            skip[i] = true;
        }
    }

    let k_invs = batch_invm(k_mods, p);
    let ck: Vec<u32> = k_invs
        .iter()
        .zip(neg_c_mods)
        .map(|(k_inv, neg_c)| neg_c.mulm(k_inv, &p))
        .collect();

    // Take some baby steps
    let order = baby_steps(base_mod_p, p, num_baby_steps, slices[0].n_lo, baby_table);
    debug_assert!(
        baby_table.iter().all(|&(x, _)| x != 0),
        "should never see 0s in baby_table when b has an inverse"
    );

    // If we know the order, our baby table contains all powers of b.
    // If ck is in there, we can do some elimination.
    if let Some(order) = order {
        for (idx, slice) in slices.iter_mut().enumerate() {
            if skip[idx] {
                continue;
            }

            if let Some(i) = lookup_baby_step(baby_table, ck[idx]) {
                // eliminate L + i, and all multiples of order afterward
                slice.eliminate_multiple(p64, base, slice.n_lo + i, order);
            }
        }

        return;
    }

    // Otherwise, we'll do some giant steps
    let m: u32 = num_baby_steps.try_into().unwrap();
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
            if let Some(j) = lookup_baby_step(baby_table, ckb[idx]) {
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
            slice.eliminate_multiple(p64, base, exp, order);
        }
    }
}

/// Populates `table` with the baby steps.
///
/// Specifically, table[i] = base^(start_exp + i) mod p.
///
/// If `num_baby_steps` is less than or equal to the order of base mod p,
/// return Some(order), and `table` will have length of that order. Otherwise,
/// the table will be filled up to `num_baby_steps` and we will return `None`.
fn baby_steps(
    base_mod_p: u32,
    p: u32,
    num_baby_steps: usize,
    start_exp: usize,
    table: &mut Vec<(u32, usize)>,
) -> Option<usize> {
    let start_exp: u32 = start_exp
        .try_into()
        .expect("sieve exponent range should fit in a u32");

    table.clear();
    let initial_value = base_mod_p.powm(start_exp, &p);
    let mut value = initial_value;
    for i in 0..num_baby_steps {
        table.push((value, i));
        value = value.mulm(base_mod_p, &p);

        // We've looped all the way around! No need to insert any more entries,
        // we can return with knowledge of the order.
        if value == initial_value {
            table.sort_unstable_by_key(|&(k, _)| k);
            return Some(i + 1);
        }
    }

    table.sort_unstable_by_key(|&(k, _)| k);
    None
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
        let slice = SequenceSlice::new(seq, 0..n_range);
        let mut slices = [slice];
        let mut prime_buffer = NaiveBuffer::new();
        let mut baby_table = Vec::new();
        for p in prime_buffer.primes(max_p) {
            baby_step_giant_step(base, *p as u32, 10, 10, &mut slices, &mut baby_table);
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

    #[test]
    fn test_batch_invm() {
        // Let's pick some numbers to invert, including a zero to make sure
        // it's skipped rather than breaking the rest of the batch.
        let p = 101;
        let values = vec![1, 2, 3, 5, 7, 100, 50, 0, 99];
        let inverses = batch_invm(values.clone(), p);

        assert_eq!(inverses.len(), values.len());
        for (&v, &inv) in values.iter().zip(inverses.iter()) {
            if v == 0 {
                assert_eq!(inv, 0, "zero should be left unmodified");
            } else {
                assert_eq!(
                    v.mulm(inv, &p),
                    1,
                    "{v} * {inv} should be 1 mod {p}, got {}",
                    v.mulm(inv, &p)
                );
            }

            assert!(inv < p, "{v} should be between 0 and {p}");
        }
    }

    #[test]
    fn test_batch_invm_empty() {
        let inverses = batch_invm(Vec::new(), 101);
        assert!(inverses.is_empty());
    }

    #[test]
    fn test_batch_invm_single() {
        let p = 13;
        let inverses = batch_invm(vec![5], p);
        // 5 * 8 = 40 = 1 mod 13
        assert_eq!(inverses, vec![8]);
    }

    #[test]
    fn test_batch_invm_all_zero() {
        let values = vec![0, 0, 0];
        let inverses = batch_invm(values.clone(), 101);
        assert_eq!(inverses, values);
    }
}
