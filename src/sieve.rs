use std::collections::HashMap;

use bitvec::prelude::BitVec;
use num_bigint::BigUint;
use num_integer::Integer;
use num_modular::{ModularCoreOps, ModularPow, ModularUnaryOps};
use num_prime::buffer::{NaiveBuffer, PrimeBufferExt};

use crate::math::gcd_reduce;

#[derive(Debug)]
pub struct Sequence {
    pub k: u64,
    pub c: i64,
    pub d: u64,
    n_lo: usize,
    n_bitvec: BitVec,
}

impl Sequence {
    pub fn new(k: u64, c: i64, d: u64, n_lo: usize, n_hi: usize) -> Self {
        assert!(n_lo <= n_hi);
        let n_range = n_hi - n_lo;

        assert_ne!(k, 0);
        assert_ne!(c, 0);
        assert_ne!(d, 0);
        // k and c have to be opposites mod d
        assert_eq!(k.checked_add_signed(c).unwrap() % d, 0);

        // Do some quick reduction to put it in lowest terms
        let gcd = gcd_reduce([k, c.unsigned_abs(), d]);

        Self {
            k: k / gcd,
            // casting is okay because 0 < gcd <= |c|
            c: c / (gcd as i64),
            d: d / gcd,
            n_lo,
            n_bitvec: BitVec::repeat(true, n_range),
        }
    }

    pub fn compute_term(&self, n: u32, base: u64) -> BigUint {
        let bn = BigUint::from(base).pow(n);
        let kbnc = if self.c > 0 {
            self.k * bn + self.c.unsigned_abs()
        } else {
            self.k * bn - self.c.unsigned_abs()
        };
        let (q, r) = kbnc.div_rem(&self.d.into());
        debug_assert_eq!(r, BigUint::ZERO);
        q
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
        if self.check_term_equal(base, p, start) {
            idx += spacing;
        }

        while let Some(mut slot) = self.n_bitvec.get_mut(idx) {
            slot.set(false);
            idx += spacing;
        }
    }

    fn check_term_equal(&self, base: u64, p: u64, n: usize) -> bool {
        let mut x = match p.checked_add_signed(-self.c) {
            Some(x) => x,
            None => return false,
        };

        if x % self.k != 0 {
            return false;
        }
        x /= self.k;

        for _ in 0..n {
            if x % base != 0 {
                return false;
            }
            x /= base;
        }

        x == 1
    }
}

pub fn find_first_prime(
    base: u8,
    k: u64,
    c: i64,
    d: u64,
    n_lo: usize,
    n_hi: usize,
    // TODO: how many? can i decide from "outside"?
    p_max: u64,
) -> Option<(usize, BigUint)> {
    let n_range = n_hi - n_lo;
    let mut seq = Sequence::new(k, c, d, n_lo, n_hi);

    // Decide how many steps for baby-step giant-step
    let num_baby_steps = (n_range as f64).sqrt() as usize;
    let num_giant_steps = n_range.div_ceil(num_baby_steps);

    // Now go and eliminate a bunch of terms
    let mut prime_buffer = NaiveBuffer::new();
    for p in prime_buffer.primes(p_max) {
        baby_step_giant_step(base.into(), *p, num_baby_steps, num_giant_steps, &mut seq);
    }

    // Lastly, iterate through the remaining numbers and see if they're prime
    println!(
        "{}/{} remaining",
        seq.n_bitvec.count_ones(),
        seq.n_bitvec.len()
    );
    for i in seq.n_bitvec.iter_ones() {
        let exponent = seq.n_lo + i;
        println!("Start computing #{}", exponent);

        // TODO: re-use the previous computation?
        let value = seq.compute_term(exponent as u32, base.into());
        println!("Start checking #{}", exponent);

        if prime_buffer.is_prime(&value, None).probably() {
            return Some((exponent, value));
        }

        println!("Done checking #{}", exponent);
    }

    None
}

fn baby_step_giant_step(
    base: u64,
    p: u64,
    num_baby_steps: usize,
    num_giant_steps: usize,
    seq: &mut Sequence,
) {
    // This works by solving b^n = (-c/k) mod p for n, using baby-step-giant-step.

    // What about d?
    // If p doesn't divide d, then we don't have to worry; p divides (kb^n+c)/d exactly
    // when it divides kb^n+c.
    // If p does divide d, then we have to be more careful, and count the number of ps.
    // For now though, skip it! (TODO)
    if
    /* p < seq.d && */
    seq.d.is_multiple_of(&p) {
        println!(
            "Skipping prime {}, it divides the denominator d={}",
            p, seq.d
        );
        return;
    }

    // Compute some inverses. If k or b is zero mod p, we need to check whether c is
    // zero mod p.
    let (kinv, binv) = match (seq.k.invm(&p), base.invm(&p)) {
        (Some(x), Some(y)) => (x, y),
        (_, _) => {
            // This'll always be equivalent to c mod p. If we're using this right,
            // we'll have non-zero c%p, so we can't eliminate anything.
            assert_ne!(
                seq.c.unsigned_abs() % p,
                0,
                "Sequence {:?} is always divisible by {}",
                seq,
                p
            );
            return;
        }
    };

    // We first have to get c as a u64 before we can do mod-p math with it.
    let neg_c_mod_p = if seq.c >= 0 {
        seq.c.unsigned_abs().negm(&p)
    } else {
        seq.c.unsigned_abs()
    };

    let ck = neg_c_mod_p.mulm(kinv, &p);
    let (baby_table, order) = baby_steps(base, p, num_baby_steps, seq.n_lo);

    // If we know the order, our baby table contains all powers of b.
    // If ck is in there, we can do some elimination.
    if let Some(order) = order {
        if let Some(i) = baby_table.get(&ck) {
            // eliminate L + i, and all multiples of order afterward
            seq.eliminate_multiple(p, base, seq.n_lo + i, order);
        }
        return;
    }

    // Otherwise, we'll do some giant steps
    let m: u64 = num_baby_steps.try_into().unwrap();
    let bm = binv.powm(m, &p);
    let mut ckb = ck;

    // We know that the order divides p-1, so worst-case scenario,
    // we can use that as the order.
    let mut order: usize = (p - 1).try_into().unwrap();
    let mut first_solution: Option<usize> = None;

    for i in 0..num_giant_steps {
        if let Some(j) = baby_table.get(&ckb) {
            // Found a solution! (-c/k)b^(im) = b^(L+j), so we eliminate L+im+j
            let exp = seq.n_lo + i * num_baby_steps + j;

            // Is this the first or second solution we've found?
            match first_solution {
                Some(n) => {
                    // Great! Two solutions tell us the order of b mod p.
                    assert!(exp > n);
                    order = exp - n;
                    break;
                }
                None => first_solution = Some(exp),
            }
        }

        // Bump ckb
        ckb = ckb.mulm(bm, &p);
    }

    if let Some(exp) = first_solution {
        seq.eliminate_multiple(p, base, exp, order);
    }
}

fn baby_steps(
    base: u64,
    p: u64,
    num_baby_steps: usize,
    start_exp: usize,
) -> (HashMap<u64, usize>, Option<usize>) {
    let start_exp: u64 = start_exp.try_into().unwrap();

    let mut map = HashMap::new();
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
    use super::*;

    #[test]
    fn test_bsgs() {
        // Start at zero for ease of understanding
        // 5*2^n+1
        let base = 2;
        let n_range = 100;
        let max_p = 100;

        let mut seq = Sequence::new(5, 1, 1, 0, n_range);
        let mut prime_buffer = NaiveBuffer::new();
        for p in prime_buffer.primes(max_p) {
            baby_step_giant_step(base, *p, 10, 10, &mut seq);
        }

        for i in 0..n_range {
            // Check every remaining element in the sequence
            let exp = seq.n_lo + i;
            let elt = seq.compute_term(exp as u32, base);
            let is_remaining = seq.check_n(i);
            let is_prime = prime_buffer.is_prime(&elt, None).probably();

            // If it's prime, we must not eliminate it.
            if is_prime {
                assert!(
                    is_remaining,
                    "{} = {}*{}^{}+{} is prime, but was removed from the list",
                    elt, seq.k, base, exp, seq.c
                );
            // If it's composite, we might have eliminated it. Specifically,
            // if it has small factors, we should have been able to eliminate it.
            } else if is_remaining {
                // Conversely, if we didn't eliminate it, it should not have small factors.
                let factors = prime_buffer.factorize(elt.clone());
                let min_factor = factors.keys().min().expect("at least one factor");
                assert!(
                    min_factor >= &max_p.into(),
                    "{} = {}*{}^{}+{} has unexpected small factor {}",
                    elt,
                    seq.k,
                    base,
                    seq.n_lo + i,
                    seq.c,
                    min_factor
                );
            }
        }

        // Also check we eliminated a substantial number of them
        let num_remaining = seq.n_bitvec.count_ones();
        assert!(
            num_remaining < 20,
            "Expected to eliminate more options, there are {} remaining",
            num_remaining
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
