use std::ops::Range;
use std::time::{Duration, Instant};

use bitvec::prelude::BitVec;
use log::debug;
use num_bigint::BigUint;
use num_modular::{ModularCoreOps, ModularPow, ModularUnaryOps};
use num_prime::buffer::{NaiveBuffer, PrimeBufferExt};

use crate::context::SearchContext;
use crate::digits::{Digit, DigitSeq};
use crate::families::SimpleFamily;
use crate::sequence::Sequence;

/// Stats for a single round of sieving (one n_range). A fresh one is created
/// (and thus implicitly reset) at the start of every [do_one_round] call, and
/// printed at the end of it, so these numbers are always specific to one
/// round rather than cumulative across the whole run.
#[derive(Debug, Default)]
pub struct SieveStats {
    /// number of sequences being sieved this round
    pub num_slices: usize,

    /// number of primes iterated over while running BSGS this round
    pub num_bsgs_primes: usize,
    /// total time spent inside baby_step_giant_step this round
    pub duration_bsgs: Duration,

    /// number of (slice, n) candidates the BSGS sieve eliminated this round
    pub num_eliminated_by_sieve: usize,
    /// number of (slice, n) candidates that survived sieving and had to be
    /// checked individually
    pub num_candidates_after_sieve: usize,

    /// number of individual candidates actually run through a primality test
    /// this round (<= num_candidates_after_sieve, since we stop checking a
    /// slice as soon as we find a prime in it)
    pub num_primality_tests: usize,
    /// total time spent on those primality tests
    pub duration_primality_tests: Duration,

    /// number of new primes found this round
    pub num_primes_found: usize,
}

impl SieveStats {
    fn record_bsgs(&mut self, count: usize, elapsed: Duration) {
        self.num_bsgs_primes += count;
        self.duration_bsgs += elapsed;
    }

    fn record_primality_test(&mut self, elapsed: Duration) {
        self.num_primality_tests += 1;
        self.duration_primality_tests += elapsed;
    }

    pub fn print(&self, n_range: &Range<usize>, p_max: u64) {
        println!(
            "---- SIEVE STATS (n = {}..{}, p_max = {}) ----",
            n_range.start, n_range.end, p_max
        );
        println!(
            "{} sequences sieved over {} values of n",
            self.num_slices,
            n_range.len()
        );
        println!(
            "{} primes used in BSGS ({}ms total, {}/prime, {}/prime/seq)",
            self.num_bsgs_primes,
            self.duration_bsgs.as_millis(),
            format_avg(self.duration_bsgs, self.num_bsgs_primes),
            format_avg(self.duration_bsgs, self.num_bsgs_primes * self.num_slices),
        );
        println!(
            "{} candidates eliminated by sieving, {} left over for primality testing",
            self.num_eliminated_by_sieve, self.num_candidates_after_sieve,
        );
        println!(
            "{} primality tests performed ({}ms total, {}/test)",
            self.num_primality_tests,
            self.duration_primality_tests.as_millis(),
            format_avg(self.duration_primality_tests, self.num_primality_tests),
        );
        println!("{} new primes found", self.num_primes_found);
    }
}

/// Formats `d / n` as a human-readable per-unit duration, or "n/a" if n is 0.
fn format_avg(d: Duration, n: usize) -> String {
    if n == 0 {
        return "n/a".to_string();
    }
    let micros_per = d.as_micros() as f64 / n as f64;
    if micros_per >= 1000.0 {
        format!("{:.3}ms", micros_per / 1000.0)
    } else {
        format!("{:.3}us", micros_per)
    }
}

/// e^-gamma, from Mertens' third theorem: the fraction of integers surviving
/// sieving by every prime up to P is asymptotically e^-gamma / ln(P).
/// e^-γ, where γ is the Euler-Mascheroni constant. Computed with WolframAlpha.
const E_NEG_GAMMA: f64 = 0.561_459_483_566_885_1;

/// Suggests the right p_max for the next round of sieving, based on the
/// timing statistics from the previous round.
///
/// The idea is that the more we sieve, the more numbers we can eliminate,
/// which means the fewer (expensive) primality tests we have to run, but
/// the sieving itself takes time. The tricky part is that it's not as easy
/// as "one more prime in the sieve => M fewer primality tests".
///
/// Mertens' third theorem tells us the probability of a number surviving a
/// sieve of primes up to P -- it is approximately e^-γ / ln(P), where γ is
/// the Euler-Mascheroni constant.
///
/// Meanwhile, the cost of doing that sieve is the cost of one round of BSGS,
/// times the number of primes less than P (approximately P/ln(P)).
///
/// So, if N is the number of candidates, cost_bsgs is the cost-per-prime, and
/// cost_test is the cost of a primality test, the total cost is:
///
///   cost_bsgs * P / ln(P) + cost_test * N * e^-gamma / ln(P)
///
/// Taking the derivative with respect to P, and setting it equal to zero, we get:
///
///   cost_bsgs * (ln(P) - 1) / ln(P)^2 = cost_test * N * e^-gamma / P / ln(P)^2
///
///   P * (ln(P) - 1) = cost_test * N * e^-gamma / cost_bsgs
///
/// If this calculation can't be done for some reason (empty stats?), returns
/// `None` instead.
pub fn suggest_next_p_max(
    prev: &SieveStats,
    prev_range_len: usize,
    next_range_len: usize,
    next_num_slices: usize,
) -> Option<u64> {
    if prev.num_bsgs_primes == 0 || prev.num_primality_tests == 0 || prev.num_slices == 0 {
        return None;
    }

    // The cost of BSGS scales up with sqrt(range * num_slices) (see the BSGS
    // code above), so take the previous cost and scale it accordingly.
    let scale = ((next_range_len * next_num_slices) as f64
        / (prev_range_len * prev.num_slices) as f64)
        .sqrt();
    let cost_bsgs = scale * (prev.duration_bsgs.as_secs_f64() / prev.num_bsgs_primes as f64);

    // The cost of primality testing grows as the number of digits gets larger.
    // I don't have a hard-and-fast rule here, but it seems to grow about 3x
    // each time the number of digits doubles.
    let scale = 3.0;
    let cost_test =
        scale * prev.duration_primality_tests.as_secs_f64() / prev.num_primality_tests as f64;

    // Remember, we want to solve P * (ln(P) - 1) = cost_test * N * e^-gamma / cost_bsgs.
    let total_candidates = (next_range_len * next_num_slices) as f64;
    let rhs = total_candidates * E_NEG_GAMMA * cost_test / cost_bsgs;

    // We approximate a solution by doing some fixed-point iteration on it.
    let mut p = rhs;
    for _ in 0..10 {
        // If p < e, then the denominator below is negative. I don't see
        // how this could happen, but bail out and return something small.
        if p < std::f64::consts::E {
            return Some(2);
        }
        p = rhs / (p.ln() - 1.0);
    }

    Some(p.round() as u64)
}

/// Entry point for eliminating simple families through sieving.
pub fn do_one_round(
    ctx: &mut SearchContext,
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
    let mut stats = SieveStats::default();
    sieve(base, &mut slices, p_max, &mut prime_buffer, &mut stats);
    last_resort(base, &slices[0], &mut prime_buffer, &mut stats)
}

fn sieve(
    base: u8,
    slices: &mut [SequenceSlice],
    // TODO: how many? can i decide from "outside"?
    p_max: u64,
    prime_buffer: &mut NaiveBuffer,
    stats: &mut SieveStats,
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

    // Traditionally, baby-step-giant-step uses sqrt(N) baby steps and sqrt(N) giant
    // steps, since this is what minimizes the sum.
    // But we're sieving multiple sequences at once, which changes the tradeoffs
    // quite a bit!
    //
    // Let's say we have m baby-steps, M giant-steps, and S slices. We build the
    // baby-step table only once, so it costs O(m). However, we iterate over the
    // giant-steps for each slice, so it costs O(M * S). (There's some log m stuff
    // in there because we have to sort the baby-table but it's probably fine to
    // ignore...)
    //
    // Since m*M has to be at least N, this gives us a total cost of m + (N/m)S,
    // which is minimized at m = sqrt(N*S).
    //
    // Basically, having more slices means that we can push more cost into the
    // (shared) baby table so we don't have to pay it for the (not-shared) giant
    // steps.
    let num_baby_steps = ((n_range * slices.len()) as f64).sqrt();
    let num_baby_steps = (num_baby_steps.round() as usize).clamp(1, n_range);
    let num_giant_steps = n_range.div_ceil(num_baby_steps);

    // The baby-step table is effectively a hashmap, but actually using one is
    // not the fastest choice. These are only ~sqrt(N*S) in size, so
    // they're pretty small, and since we're populating it for each prime up
    // to p_max, we really want to re-use our storage. Let's just use a sorted
    // vector.
    let mut baby_table: Vec<(u32, usize)> = Vec::with_capacity(num_baby_steps);

    // Now go and eliminate a bunch of terms
    let start = Instant::now();
    let mut num_primes = 0;
    for p in prime_buffer.primes(p_max) {
        baby_step_giant_step(
            base.into(),
            *p as u32,
            num_baby_steps,
            num_giant_steps,
            slices,
            &mut baby_table,
        );
        num_primes += 1;
    }
    stats.record_bsgs(num_primes, start.elapsed());
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
    fn test_suggest_next_p_max_no_history() {
        let stats = SieveStats::default();
        assert_eq!(suggest_next_p_max(&stats, 16, 32, 100), None);
    }

    #[test]
    fn test_suggest_next_p_max_solves_the_balance_condition() {
        // Made-up but self-consistent stats: 1000 primes cost 1ms total in
        // BSGS (1us/prime), and 500 primality tests cost 100ms total
        // (200us/test).
        let stats = SieveStats {
            num_slices: 10,
            num_bsgs_primes: 1000,
            duration_bsgs: Duration::from_millis(1),
            num_primality_tests: 500,
            duration_primality_tests: Duration::from_millis(100),
            ..Default::default()
        };

        let p = suggest_next_p_max(&stats, 16, 16, 10).expect("should have a suggestion");

        // Check it actually (approximately) solves p*(ln(p)-1) = rhs for the
        // rhs this scenario implies, rather than just trusting the formula
        // was transcribed correctly. Mirrors suggest_next_p_max's own math,
        // including the flat 3x fudge factor for primality-test cost growth.
        let t_bsgs = stats.duration_bsgs.as_secs_f64() / stats.num_bsgs_primes as f64;
        let t_test =
            3.0 * stats.duration_primality_tests.as_secs_f64() / stats.num_primality_tests as f64;
        let n = 16.0 * 10.0;
        let rhs = n * E_NEG_GAMMA * t_test / t_bsgs;

        let p = p as f64;
        let residual = (p * (p.ln() - 1.0) - rhs).abs() / rhs;
        assert!(
            residual < 0.01,
            "p={p} doesn't satisfy p*(ln(p)-1)={rhs} closely enough (residual {residual})"
        );
    }

    #[test]
    fn test_suggest_next_p_max_grows_with_more_candidates() {
        // Same per-unit costs, but the next round covers more ground (more
        // candidates to potentially save primality tests on) -- should want
        // to sieve deeper.
        let stats = SieveStats {
            num_slices: 10,
            num_bsgs_primes: 1000,
            duration_bsgs: Duration::from_millis(1),
            num_primality_tests: 500,
            duration_primality_tests: Duration::from_millis(100),
            ..Default::default()
        };

        let p_small = suggest_next_p_max(&stats, 16, 16, 10).unwrap();
        let p_large = suggest_next_p_max(&stats, 16, 16, 1000).unwrap();
        assert!(
            p_large > p_small,
            "expected deeper sieving for a round with more candidates: {p_large} <= {p_small}"
        );
    }

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
