//! This module is used to check whether families are always composite,
//! letting us discard the whole branch.

use itertools::Itertools;
use num_bigint::BigUint;
use num_prime::ExactRoots;

use crate::digits::DigitSeq;
use crate::families::Core;
use crate::families::Family;
use crate::families::Sequence;
use crate::search::gcd::nontrivial_gcd;
use crate::search::SimpleNode;

/// Checks whether this family shares a factor with the base.
/// Basically just checks the last digit.
pub fn shares_factor_with_base(base: u8, family: &Family) -> Option<u8> {
    // Get the last digit of the family
    let last_seq = family.digitseqs.last().expect("digitseqs nonempty");
    let d = last_seq.0.last()?;

    nontrivial_gcd(&d.0, &base)
}

/// Given a family of the shape xLyMz..., checks whether it has a single factor
/// dividing every member of the family.
pub fn find_guaranteed_factor(base: u8, family: &Family) -> Option<BigUint> {
    let mut gcd = family.contract().value(base);

    for (i, core) in family.cores.iter().enumerate() {
        for d in core.iter() {
            gcd = nontrivial_gcd(&gcd, &family.substitute(i, d).value(base))?;
        }
    }

    Some(gcd)
}

/// Given a family of the shape xLz, checks whether there are any "periodic"
/// factors, up to the given limit. Only works on families with one core.
///
/// For a period n, a periodic factor sequence is a list of numbers f_1, ..., f_n
/// such that f_i divides xL^(i+kn)z for all k.
///
/// For example, a period-1 factor is a single number that divides xz, xLz, xLLz, ...,
/// and a period-2 factor is two numbers N and M such that N divides xz, xLLz, xL^4z, ...,
/// and M divides xLz, xLLLz, xL^5z, ....
pub fn find_periodic_factor(base: u8, family: &Family, stride: usize) -> Option<Vec<BigUint>> {
    let one = BigUint::from(1_u32);

    // TODO: generalize
    if family.cores.len() != 1 {
        return None;
    }
    let core = &family.cores[0];

    let mut gcds = vec![];
    for k in 0..stride {
        let mut g = BigUint::ZERO;

        // TODO: shouldn't this be cartesian product??
        for center in core.iter().combinations_with_replacement(k) {
            let value = DigitSeq::concat_value(
                [&family.digitseqs[0], &center.into(), &family.digitseqs[1]],
                base,
            );
            g = nontrivial_gcd(&g, &value)?;
        }

        for center in core.iter().combinations_with_replacement(k + stride) {
            let value = DigitSeq::concat_value(
                [&family.digitseqs[0], &center.into(), &family.digitseqs[1]],
                base,
            );
            g = nontrivial_gcd(&g, &value)?;
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
        a: &'f Core,
        y: &'f DigitSeq,
        b: &'f Core,
        z: &'f DigitSeq,
        a_repeat: usize,
        b_repeat: usize,
    ) -> impl Iterator<Item = BigUint> + 'f {
        a.iter()
            .combinations_with_replacement(a_repeat)
            .cartesian_product(b.iter().combinations_with_replacement(b_repeat))
            .map(move |(a_choices, b_choices)| {
                let mut seq = x.clone();
                seq += DigitSeq(a_choices.clone());
                seq += y;
                seq += DigitSeq(b_choices);
                seq += z;

                seq.value(base)
            })
    }

    let bar = |repeat_a, repeat_b| {
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
    let mut even_gcd = BigUint::ZERO;
    for (repeat_a, repeat_b) in even_repeats {
        for value in bar(repeat_a, repeat_b) {
            even_gcd = nontrivial_gcd(&even_gcd, &value)?;
        }
    }

    let odd_repeats: [(usize, usize); 6] = [(1, 0), (3, 0), (1, 2), (0, 1), (2, 1), (0, 3)];
    let mut odd_gcd = BigUint::ZERO;
    for (repeat_a, repeat_b) in odd_repeats {
        for value in bar(repeat_a, repeat_b) {
            odd_gcd = nontrivial_gcd(&odd_gcd, &value)?;
        }
    }

    debug_assert_ne!(even_gcd, BigUint::ZERO);
    debug_assert_ne!(odd_gcd, BigUint::ZERO);

    Some((even_gcd, odd_gcd))
}

pub fn composite_checks_for_simple(base: u8, node: &SimpleNode) -> bool {
    // Get the sequence for this. As a reminder, it looks like:
    // (k B^n + c) / d, where d = B-1
    let sequence = match node.sequence {
        Some(s) => s,
        None => return false, // nothing we can check :(
    };

    check_sum_diff_of_cubes(base, &sequence) || check_diff_of_squares(base, &sequence)
}

macro_rules! bail_if_none {
    ($e:expr) => {
        match $e {
            Some(x) => x,
            None => return false,
        }
    };
}

fn check_sum_diff_of_cubes(base: u8, sequence: &Sequence) -> bool {
    // We need B to be a cube, and also k and c.
    // What about d? We can mostly ignore it, but we do have to be careful
    // that neither of the factors is smaller than d, otherwise we could
    // unwittingly be factoring into p * 1.
    let _cbrt_base = bail_if_none!(base.cbrt_exact());
    let cbrt_k: i64 = bail_if_none!(sequence.k.cbrt_exact())
        .try_into()
        .expect("cube root should be small enough to fit in i64");
    let cbrt_c = bail_if_none!(sequence.c.cbrt_exact());

    // It's a possibility! Now check the factor size.
    // a^3 + b^3 = (a + b)(a^2 - ab + b^2)
    let d = sequence.d.try_into().expect("d should be small");
    let x = cbrt_k + cbrt_c;
    let y = cbrt_k * cbrt_k - cbrt_k * cbrt_c + cbrt_c * cbrt_c;
    x.abs() > d && y.abs() > d

    // TODO: if that check fails, we should retry with higher n! we might
    // still be composite, and we should check another time
    // TODO: also log (to tree) what the specific reason is!
}

fn check_diff_of_squares(base: u8, sequence: &Sequence) -> bool {
    // we need the base, k and -c to be squares
    let _sqrt_base = bail_if_none!(base.sqrt_exact());
    let sqrt_k: i64 = bail_if_none!(sequence.k.sqrt_exact())
        .try_into()
        .expect("square root should be small enough to fit in i64");
    let sqrt_neg_c = bail_if_none!((-sequence.c).sqrt_exact());

    // Check the factor size
    // TODO: for 1* in base 9, this won't work! you need to actually implement
    // the fancy "definitely composite after n but maybe prime before that"
    // tracking.
    true
}
