//! This module is used to check whether patterns are always composite,
//! letting us discard the whole branch.

use itertools::Itertools;
use num_bigint::BigUint;

use crate::data::Digit;
use crate::data::DigitSeq;
use crate::data::Pattern;

/// Checks whether this pattern shares a factor with the base.
/// Basically just checks the last digit.
pub fn shares_factor_with_base(base: u8, pattern: &Pattern) -> Option<BigUint> {
    // Get the last digit of the pattern
    let last_seq = pattern.digitseqs.last().expect("digitseqs nonempty");
    let d = last_seq.0.last()?;

    let gcd = gcd(d.0.into(), base.into());
    debug_assert!(gcd != BigUint::ZERO);
    if gcd != big_one() {
        Some(gcd)
    } else {
        None
    }
}

/// Given a pattern of the shape xLz, checks whether there are any "periodic"
/// factors, up to the given limit.
///
/// For a period n, a periodic factor sequence is a list of numbers f_1, ..., f_n
/// such that f_i divides xL^(i+kn)z for all k.
///
/// For example, a period-1 factor is a single number that divides xz, xLz, xLLz, ...,
/// and a period-2 factor is two numbers N and M such that N divides xz, xLLz, xL^4z, ...,
/// and M divides xLz, xLLLz, xL^5z, ....
///
/// TODO: period one but allow multiple cores
pub fn find_perpetual_factor(base: u8, pattern: &Pattern, stride: usize) -> Option<Vec<BigUint>> {
    let one = BigUint::from(1_u32);

    // TODO: generalize
    if pattern.cores.len() != 1 {
        return None;
    }

    let mut gcds = vec![BigUint::ZERO; stride];
    for (i, gcd_i) in gcds.iter_mut().enumerate() {
        // The smaller of the two sets: xL^iz
        for center in pattern.cores[0]
            .iter()
            .copied()
            .combinations_with_replacement(i)
        {
            let value = DigitSeq::concat_value(
                [&pattern.digitseqs[0], &center.into(), &pattern.digitseqs[1]],
                base,
            );
            // Update the GCD. If we ever see a 1, it's always going to
            // be that way, so bail out instantly.
            *gcd_i = gcd(gcd_i.clone(), value);
            if *gcd_i == one {
                return None;
            }
        }
        // The larger of the two sets: xL^(i+stride)z
        for center in pattern.cores[0]
            .iter()
            .copied()
            .combinations_with_replacement(i + stride)
        {
            let value = DigitSeq::concat_value(
                [&pattern.digitseqs[0], &center.into(), &pattern.digitseqs[1]],
                base,
            );
            *gcd_i = gcd(gcd_i.clone(), value);
            if *gcd_i == one {
                return None;
            }
        }
    }
    for g in &gcds {
        debug_assert_ne!(*g, one);
    }
    Some(gcds)
}

pub fn find_even_odd_factor(base: u8, pattern: &Pattern) -> Option<(BigUint, BigUint)> {
    // Similar to find_perpetual_factor, but checks strings of even length and odd length separately.

    // TODO: can this be generalized to more than two cores?
    if pattern.cores.len() != 2 {
        return None;
    }

    // TODO: generalize this!
    fn foo(
        initial_gcd: BigUint,
        base: u8,
        x: &DigitSeq,
        a: &[Digit],
        y: &DigitSeq,
        b: &[Digit],
        z: &DigitSeq,
        a_repeat: usize,
        b_repeat: usize,
    ) -> BigUint {
        let mut gcd_accum = initial_gcd;

        for a_choices in a.iter().copied().combinations_with_replacement(a_repeat) {
            for b_choices in b.iter().copied().combinations_with_replacement(b_repeat) {
                let mut seq = x.clone();
                seq += DigitSeq(a_choices.clone());
                seq += y;
                seq += DigitSeq(b_choices);
                seq += z;

                let n = seq.value(base);
                gcd_accum = gcd(gcd_accum, n);
                if gcd_accum == BigUint::from(1_u32) {
                    return gcd_accum;
                }
            }
        }
        gcd_accum
    }

    // Take the pattern xA*yB*z
    // - even number of A+B: x(AA)*y(BB)*z, xA(AA)*yB(BB)*z
    // -  odd number of A+B: xA(AA)*y(BB)*z, z(AA)*yB(BB)*z
    let mut even_gcd = BigUint::ZERO;
    let mut odd_gcd = BigUint::ZERO;

    let bar = |gcd, repeat_a, repeat_b| {
        foo(
            gcd,
            base,
            &pattern.digitseqs[0],
            &pattern.cores[0],
            &pattern.digitseqs[1],
            &pattern.cores[1],
            &pattern.digitseqs[2],
            repeat_a,
            repeat_b,
        )
    };

    even_gcd = bar(even_gcd, 0, 0);
    even_gcd = bar(even_gcd, 2, 0);
    even_gcd = bar(even_gcd, 0, 2);
    even_gcd = bar(even_gcd, 1, 1);
    even_gcd = bar(even_gcd, 1, 3);
    even_gcd = bar(even_gcd, 3, 1);

    odd_gcd = bar(odd_gcd, 1, 0);
    odd_gcd = bar(odd_gcd, 3, 0);
    odd_gcd = bar(odd_gcd, 1, 2);
    odd_gcd = bar(odd_gcd, 0, 1);
    odd_gcd = bar(odd_gcd, 2, 1);
    odd_gcd = bar(odd_gcd, 0, 3);

    if even_gcd > big_one() && odd_gcd > big_one() {
        return Some((even_gcd, odd_gcd));
    }

    None
}

pub fn big_one() -> BigUint {
    BigUint::from(1_u8)
}

pub fn gcd(mut a: BigUint, mut b: BigUint) -> BigUint {
    // Since we're working with a binary representation, it's more efficient
    // to use this algorithm than the standard Euclidean one.

    // First, pull out all the factors of 2.
    // gcd(2^i a, 2^j b) = 2^k gcd(a, b) when a, b odd and k = min(i, j);
    // But since x.trailing_zeros() is None when x is 0, take this opportunity
    // to check for the edge case of a = 0 or b = 0.
    let i = match a.trailing_zeros() {
        Some(i) => i,
        None => return b,
    };
    let j = match b.trailing_zeros() {
        Some(j) => j,
        None => return a,
    };

    // How many of those 2s do we want to keep?
    let k = i.min(j);
    a >>= i;
    b >>= j;

    loop {
        // Now they're both odd
        debug_assert!(a.bit(0));
        debug_assert!(b.bit(0));

        // Swap so a is larger
        if a < b {
            std::mem::swap(&mut a, &mut b);
        }

        // Subtract, just like with Euclid
        // gcd(u, v) = gcd(u, v-u)
        a -= &b;

        // Now a is even; remove its 2s (again checking for the case of a = 0).
        // gcd(2^i a, b) = gcd(a, b) when b is odd
        match a.trailing_zeros() {
            Some(i) => {
                a >>= i;
            }
            None => {
                // gcd(0, b) = b, then add in the 2^k we removed earlier
                return b << k;
            }
        }
    }
}
