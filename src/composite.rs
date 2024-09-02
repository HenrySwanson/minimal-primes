//! This module is used to check whether families are always composite,
//! letting us discard the whole branch.

use itertools::Itertools;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::identities::One;

use crate::data::Digit;
use crate::data::DigitSeq;
use crate::data::Family;
use crate::math::gcd_reduce;

/// Checks whether this family shares a factor with the base.
/// Basically just checks the last digit.
pub fn shares_factor_with_base(base: u8, family: &Family) -> Option<u8> {
    // Get the last digit of the family
    let last_seq = family.digitseqs.last().expect("digitseqs nonempty");
    let d = last_seq.0.last()?;

    let gcd = d.0.gcd(&base);
    debug_assert!(gcd != 0);
    if gcd != 1 {
        Some(gcd)
    } else {
        None
    }
}

/// Given a family of the shape xLz, checks whether there are any "periodic"
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
pub fn find_perpetual_factor(base: u8, family: &Family, stride: usize) -> Option<Vec<BigUint>> {
    let one = BigUint::from(1_u32);

    // TODO: generalize
    if family.cores.len() != 1 {
        return None;
    }
    let core = &family.cores[0];

    let mut gcds = vec![];
    for i in 0..stride {
        let g = gcd_reduce(
            // The smaller of the two sets: xL^iz
            core.iter()
                .copied()
                .combinations_with_replacement(i)
                .chain(
                    // The larger of the two sets: xL^(i+stride)z
                    core.iter()
                        .copied()
                        .combinations_with_replacement(i + stride),
                )
                .map(|center| {
                    DigitSeq::concat_value(
                        [&family.digitseqs[0], &center.into(), &family.digitseqs[1]],
                        base,
                    )
                }),
        );

        // If any of the GCDs are 1, bail out immediately.
        if g == BigUint::one() {
            return None;
        }

        gcds.push(g);
    }
    for g in &gcds {
        debug_assert_ne!(*g, one);
    }
    Some(gcds)
}

pub fn find_even_odd_factor(base: u8, family: &Family) -> Option<(BigUint, BigUint)> {
    // Similar to find_perpetual_factor, but checks strings of even length and odd length separately.

    // TODO: can this be generalized to more than two cores?
    if family.cores.len() != 2 {
        return None;
    }

    // Complicated, but actually does help quite a bit.
    fn family_iter_helper<'f>(
        base: u8,
        x: &'f DigitSeq,
        a: &'f [Digit],
        y: &'f DigitSeq,
        b: &'f [Digit],
        z: &'f DigitSeq,
        a_repeat: usize,
        b_repeat: usize,
    ) -> impl Iterator<Item = BigUint> + 'f {
        a.iter()
            .copied()
            .combinations_with_replacement(a_repeat)
            .cartesian_product(b.iter().copied().combinations_with_replacement(b_repeat))
            .map(move |(a_choices, b_choices)| {
                let mut seq = x.clone();
                seq += DigitSeq(a_choices.clone());
                seq += y;
                seq += DigitSeq(b_choices);
                seq += z;

                seq.value(base)
            })
    }

    let bar = |&(repeat_a, repeat_b)| {
        family_iter_helper(
            base,
            &family.digitseqs[0],
            &family.cores[0],
            &family.digitseqs[1],
            &family.cores[1],
            &family.digitseqs[2],
            repeat_a,
            repeat_b,
        )
    };

    // Take the family xA*yB*z
    // - even number of A+B: x(AA)*y(BB)*z, xA(AA)*yB(BB)*z
    // -  odd number of A+B: xA(AA)*y(BB)*z, z(AA)*yB(BB)*z
    let even_repeats: [(usize, usize); 6] = [(0, 0), (2, 0), (0, 2), (1, 1), (1, 3), (3, 1)];
    let even_gcd = gcd_reduce(even_repeats.iter().flat_map(bar));

    if even_gcd == BigUint::one() {
        return None;
    }

    let odd_repeats: [(usize, usize); 6] = [(1, 0), (3, 0), (1, 2), (0, 1), (2, 1), (0, 3)];
    let odd_gcd = gcd_reduce(odd_repeats.iter().flat_map(bar));

    if odd_gcd == BigUint::one() {
        return None;
    }

    debug_assert_ne!(even_gcd, BigUint::ZERO);
    debug_assert_ne!(odd_gcd, BigUint::ZERO);

    Some((even_gcd, odd_gcd))
}
