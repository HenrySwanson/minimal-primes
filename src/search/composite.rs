//! This module is used to check whether families are always composite,
//! letting us discard the whole branch.

use itertools::Itertools;
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_prime::ExactRoots;

use crate::digits::DigitSeq;
use crate::families::{BigSequence, Core, Family};
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

/// This is the most complex check we have, I think.
///
/// Given a family x_1 L_1 ... x_n L_n, and some index m, partitions the elements
/// in that family into two parts:
/// - those with an even number of digits contributed from L_m
/// - those with an odd number of digits
///
/// If those families each have common factors, returns them.
pub fn find_two_factors(base: u8, family: &Family) -> Option<(BigUint, BigUint)> {
    for i in 0..family.cores.len() {
        if let Some(factors) = find_two_factors_helper(base, family, i) {
            return Some(factors);
        }
    }
    None
}

fn find_two_factors_helper(base: u8, family: &Family, i: usize) -> Option<(BigUint, BigUint)> {
    let core_i = &family.cores[i];

    // gcd(∅, j, ii) for all j≠i
    let mut even_gcd = family.contract().value(base);
    for (j, core_j) in family.cores.iter().enumerate() {
        if i == j {
            continue;
        }

        for d in core_j.iter() {
            let val = family.substitute(j, d).value(base);
            even_gcd = nontrivial_gcd(&even_gcd, &val)?;
        }
    }
    for d1 in core_i.iter() {
        for d2 in core_i.iter() {
            let val = family.substitute_multiple(i, [d1, d2]).value(base);
            even_gcd = nontrivial_gcd(&even_gcd, &val)?;
        }
    }

    // gcd(i, ij, iii) for all j≠i
    let mut odd_gcd = BigUint::ZERO;
    for d1 in core_i.iter() {
        let val = family.substitute(i, d1).value(base);
        odd_gcd = nontrivial_gcd(&odd_gcd, &val)?;

        for (j, core_j) in family.cores.iter().enumerate() {
            if i == j {
                continue;
            }

            for dj in core_j.iter() {
                let val = family.substitute_two(i, d1, j, dj).value(base);
                odd_gcd = nontrivial_gcd(&odd_gcd, &val)?;
            }
        }

        for d2 in core_i.iter() {
            for d3 in core_i.iter() {
                let val = family.substitute_multiple(i, [d1, d2, d3]).value(base);
                odd_gcd = nontrivial_gcd(&odd_gcd, &val)?;
            }
        }
    }

    Some((even_gcd, odd_gcd))
}

pub fn find_even_odd_factor(base: u8, family: &Family) -> Option<(BigUint, BigUint)> {
    // Similar to find_perpetual_factor, but checks strings of even length and odd length separately.

    // TODO: can this be generalized to more than two cores?
    if family.cores.len() != 2 {
        return None;
    }

    // Complicated, but actually does help quite a bit.
    let family_iter_helper = |repeat_a, repeat_b| {
        let x = &family.digitseqs[0];
        let a = &family.cores[0];
        let y = &family.digitseqs[1];
        let b = &family.cores[1];
        let z = &family.digitseqs[2];
        a.iter()
            .combinations_with_replacement(repeat_a)
            .cartesian_product(b.iter().combinations_with_replacement(repeat_b))
            .map(move |(a_choices, b_choices)| {
                let mut seq = x.clone();
                seq += DigitSeq(a_choices.clone());
                seq += y;
                seq += DigitSeq(b_choices);
                seq += z;

                seq.value(base)
            })
    };

    // Take the family xA*yB*z
    // - even number of A+B: x(AA)*y(BB)*z, xA(AA)*yB(BB)*z
    // -  odd number of A+B: xA(AA)*y(BB)*z, z(AA)*yB(BB)*z
    let even_repeats: [(usize, usize); 6] = [(0, 0), (2, 0), (0, 2), (1, 1), (1, 3), (3, 1)];
    let mut even_gcd = BigUint::ZERO;
    for (repeat_a, repeat_b) in even_repeats {
        for value in family_iter_helper(repeat_a, repeat_b) {
            even_gcd = nontrivial_gcd(&even_gcd, &value)?;
        }
    }

    let odd_repeats: [(usize, usize); 6] = [(1, 0), (3, 0), (1, 2), (0, 1), (2, 1), (0, 3)];
    let mut odd_gcd = BigUint::ZERO;
    for (repeat_a, repeat_b) in odd_repeats {
        for value in family_iter_helper(repeat_a, repeat_b) {
            odd_gcd = nontrivial_gcd(&odd_gcd, &value)?;
        }
    }

    debug_assert_ne!(even_gcd, BigUint::ZERO);
    debug_assert_ne!(odd_gcd, BigUint::ZERO);

    Some((even_gcd, odd_gcd))
}

pub fn check_residues_mod_30(base: u8, family: &Family) -> bool {
    // This is a weird one. It's Lemma 34 in Curtis Bright's paper.
    // Basically, for a given N, we can compute all the residues of a given family mod N.
    // If these are all > 1, then we know the family has a non-trivial factor.
    // Denote [L] as the set of residues mod N.
    // Like Bright, we use N = 30.

    // Small caveat: if L contains 2, 3, or 5, then this test will fail. The small
    // prime will not be detected, and we will claim this family is composite.
    // But we can just check the length of the family (if the base is not tiny)
    // to bail out of this.
    if base < 5 || family.weight() <= 1 {
        return false;
    }

    let residues = get_residues_mod_30(base, family);
    let has_bad_residue = residues
        .iter()
        .enumerate()
        .any(|(i, is_residue)| *is_residue && i.gcd(&30) == 1);

    !has_bad_residue
}

fn get_residues_mod_30(base: u8, family: &Family) -> [bool; 30] {
    let mut residues = [false; 30];

    fn process_seq(base: u8, seq: &DigitSeq, residues: &mut [bool; 30]) {
        // [Lx] = [L] * base^|x| + (x mod N)
        let mut base_pow: usize = 1;
        for _ in &seq.0 {
            base_pow = (base_pow * base as usize) % 30;
        }

        let x_mod_30: usize = (seq.value(base) % 30_usize)
            .try_into()
            .expect("30 fits in a usize");
        let mut new_residues = [false; 30];
        for (i, is_residue) in residues.iter().enumerate() {
            if *is_residue {
                let j = (i * base_pow + x_mod_30) % 30;
                new_residues[j] = true;
            }
        }

        *residues = new_residues;
    }

    fn process_core(base: u8, core: &Core, residues: &mut [bool; 30]) {
        // [LX*] = [L] U [L x_i] U [L x_i x_j] U ...
        // Eventually this union chain will stabilize, so let's just
        // compute it until it does.
        loop {
            let mut has_changed = false;

            for d in core.iter() {
                for i in 0..30 {
                    if residues[i] {
                        let j = (i * base as usize + d.0 as usize) % 30;
                        if !residues[j] {
                            residues[j] = true;
                            has_changed = true;
                        }
                    }
                }
            }

            if !has_changed {
                return;
            }
        }
    }

    residues[0] = true; // starting state
    for (digitseq, core) in family.digitseqs.iter().zip(&family.cores) {
        process_seq(base, digitseq, &mut residues);
        process_core(base, core, &mut residues);
    }
    process_seq(base, family.digitseqs.last().unwrap(), &mut residues);
    residues
}

pub fn composite_checks_for_simple(base: u8, node: &SimpleNode) -> bool {
    // Get the sequence for this. As a reminder, it looks like:
    // (k B^n + c) / d, where d = B-1
    let sequence = BigSequence::from_family(&node.family, base);

    check_sum_diff_of_cubes(base, &sequence)
        || check_diff_of_squares(base, &sequence)
        || check_diff_of_squares_or_divisor(base, &sequence)
}

macro_rules! bail_if_none {
    ($e:expr) => {
        match $e {
            Some(x) => x,
            None => return false,
        }
    };
}

fn check_sum_diff_of_cubes(base: u8, sequence: &BigSequence) -> bool {
    // We need B to be a cube, and also k and c.
    // What about d? We can mostly ignore it, but we do have to be careful
    // that neither of the factors is smaller than d, otherwise we could
    // unwittingly be factoring into p * 1.
    let _cbrt_base = bail_if_none!(base.cbrt_exact());
    let cbrt_k: BigInt = bail_if_none!(sequence.k.cbrt_exact()).into();
    let cbrt_c = bail_if_none!(sequence.c.cbrt_exact());

    // It's a possibility! Now check the factor size.
    // a^3 + b^3 = (a + b)(a^2 - ab + b^2)
    let d: BigUint = sequence.d.into();
    let x = &cbrt_k + &cbrt_c;
    let y = &cbrt_k * &cbrt_k - &cbrt_k * &cbrt_c + &cbrt_c * &cbrt_c;
    x.magnitude() > &d && y.magnitude() > &d

    // TODO: if that check fails, we should retry with higher n! we might
    // still be composite, and we should check another time
    // TODO: also log (to tree) what the specific reason is!
}

fn check_diff_of_squares(base: u8, sequence: &BigSequence) -> bool {
    // we need the base, k and -c to be squares
    let _sqrt_base = bail_if_none!(base.sqrt_exact());
    let _sqrt_k: BigInt = bail_if_none!(sequence.k.sqrt_exact()).into();
    let _sqrt_neg_c = bail_if_none!((-&sequence.c).sqrt_exact());

    // Check the factor size
    // TODO: for 1* in base 9, this won't work! you need to actually implement
    // the fancy "definitely composite after n but maybe prime before that"
    // tracking.
    true
}

// TODO: i feel like this could be simplified
fn check_diff_of_squares_or_divisor(base: u8, sequence: &BigSequence) -> bool {
    // maybe kB^n+c factors in alternating ways
    // - difference of squares
    // - common factor
    //
    // this can happen when k is a square or when kB is a square. -c has to
    // be a square either way though.
    bail_if_none!((-&sequence.c).sqrt_exact());

    let base: u64 = base.into();
    let k: BigInt = sequence.k.clone().into();
    let kb = &k * base;
    if sequence.k.is_square() {
        // Factors as a difference of squares for even n. Now we just
        // need to see if the odd n terms have a shared factor.
        // It suffices to check n=1 and n=3.
        let a = (&kb + &sequence.c) / sequence.d;
        let b = (&kb * base * base + &sequence.c) / sequence.d;
        if nontrivial_gcd(&a, &b).is_some() {
            println!("{sequence} is composite ({a} and {b})");
            return true;
        }
    }

    if kb.is_square() {
        // Factors as a difference of squares for odd n. Check
        // the even terms for common factors.
        let a = (k + &sequence.c) / sequence.d;
        let b = (kb * base + &sequence.c) / sequence.d;
        if nontrivial_gcd(&a, &b).is_some() {
            println!("{sequence} is composite ({a} and {b})");
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digits::Digit;

    fn parse_family(input: &str, base: u8) -> Family {
        let base: u32 = base.into();
        let mut digitseqs = vec![];
        let mut cores = vec![];

        fn string_to_vec_digit(s: &str, base: u32) -> Vec<Digit> {
            s.split_inclusive(|_| true)
                .map(|ch| Digit(u8::from_str_radix(ch, base).unwrap()))
                .collect()
        }

        for s in input.split('*') {
            // s should be `x[y]` (the * is removed), or `z` for the last one
            if let Some(s) = s.strip_suffix(']') {
                let (x, y) = s.split_once('[').unwrap();
                digitseqs.push(DigitSeq(string_to_vec_digit(x, base)));
                cores.push(Core::new(string_to_vec_digit(y, base)));
            } else {
                digitseqs.push(DigitSeq(string_to_vec_digit(s, base)));
            }
        }

        assert_eq!(digitseqs.len(), cores.len() + 1);
        Family { digitseqs, cores }
    }

    #[test]
    fn example_5() {
        // Base 10: 46*9 should be divisible by 7 always
        let base = 10;
        let family = parse_family("4[6]*9", base);
        assert_eq!(
            Some(BigUint::from(7_u32)),
            find_guaranteed_factor(base, &family)
        );
    }

    #[test]
    fn example_7() {
        // Base 9: 61* should be divisible by 2 and 5 alternately
        let base = 9;
        let family = parse_family("6[1]*", base);
        assert_eq!(
            Some(vec![BigUint::from(2_u32), BigUint::from(5_u32)]),
            find_periodic_factor(base, &family, 2)
        );
    }

    #[test]
    fn example_8() {
        // Base 16: 8A0A*1 should be divisible by 7, 13, 3
        let base = 16;
        let family = parse_family("8A0[A]*1", base);
        assert_eq!(
            Some(vec![
                BigUint::from(7_u32),
                BigUint::from(13_u32),
                BigUint::from(3_u32)
            ]),
            find_periodic_factor(base, &family, 3)
        );
    }

    #[test]
    fn example_10() {
        // Base 10: 90*80*1 is always divisible by 9
        let base = 10;
        let family = parse_family("9[0]*8[0]*1", base);
        assert_eq!(
            Some(BigUint::from(9_u32)),
            find_guaranteed_factor(base, &family)
        );
    }

    #[test]
    fn example_11() {
        // Base 11: 44[0]*A[1]*1 should be composite by looking at slot 2
        let base = 11;
        let family = parse_family("44[0]*A[1]*1", base);
        assert_eq!(
            Some((BigUint::from(3_u32), BigUint::from(2_u32))),
            find_two_factors_helper(base, &family, 1)
        );
    }

    #[test]
    fn example_12() {
        // Base 9: 1*61*
        let base = 9;
        let family = parse_family("[1]*6[1]*", base);
        assert_eq!(
            Some((BigUint::from(2_u32), BigUint::from(5_u32))),
            find_even_odd_factor(base, &family)
        );
    }

    fn get_residues_as_string(base: u8, family: &Family) -> String {
        get_residues_mod_30(base, family)
            .iter()
            .enumerate()
            .filter_map(|(i, is_residue)| is_residue.then_some(format!("{i}")))
            .join(",")
    }

    #[test]
    fn example_35() {
        // Base 29: L1*61*LK*K
        let base = 29;
        let family = parse_family("L[1]*6[1]*L[K]*K", base);
        let residues_as_string = get_residues_as_string(base, &family);
        assert_eq!(residues_as_string, "4,5,6,14,15,16,25");
        assert!(check_residues_mod_30(base, &family));
    }

    #[test]
    fn example_custom() {
        // Found this one when I was experimenting on base 29
        let base = 29;
        let family = parse_family("[P]*MMMMM[R]*", base);
        let residues_as_string = get_residues_as_string(base, &family);
        assert_eq!(residues_as_string, "0,5,22,27");
        assert!(check_residues_mod_30(base, &family));
    }
}
