use itertools::Itertools;
use log::debug;

use crate::debug_to_tree;
use crate::digits::{Digit, DigitSeq};
use crate::families::Family;
use crate::search::gcd::nontrivial_gcd;
use crate::search::SearchContext;

// TODO: this probably shouldn't be searchcontext, but this works well now
impl SearchContext {
    /// Given a family `xLz`, checks if there's any y in L such that `x y^n z` is forbidden.
    /// If so, we split the family into `x (L-y) (y (L-y))^i z`.
    ///
    /// This is extended to multiple cores in a straightforward way.
    ///
    /// We check n from 1 to `max_repeats`.
    pub fn split_on_limited_digit(
        &mut self,
        family: &Family,
        max_repeats: usize,
    ) -> Option<Vec<Family>> {
        for (i, core) in family.cores.iter().enumerate() {
            for d in core.iter() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = family.substitute_multiple(i, std::iter::repeat_n(d, n));
                    if let Some(p) = self.test_for_contained_prime(&seq) {
                        assert_ne!(&seq, p);
                        debug!("  {seq} contains a prime {p}");

                        // Split into n families, x (L-y) (y (L-y))^i z for i in 0..n
                        let yless_core = core.clone().without(d);
                        let mut first_child = family.clone();
                        first_child.cores[i] = yless_core.clone();

                        let mut children = vec![first_child];

                        while children.len() < n {
                            let mut new = children.last().unwrap().clone();
                            new.digitseqs.insert(i + 1, d.into());
                            new.cores.insert(i + 1, yless_core.clone());
                            children.push(new);
                        }

                        // Simplify everything (don't do it while we're generating families),
                        // since that'd mess with indices).
                        for child in children.iter_mut() {
                            child.simplify();
                        }

                        debug!(
                            "  {} split into {}",
                            family,
                            children.iter().format(" and ")
                        );
                        debug_to_tree!(
                            self.tracer,
                            "Splitting into {}",
                            children.iter().format(" and ")
                        );
                        return Some(children);
                    }
                }
            }
        }
        None
    }

    /// Given a family `xLz`, if there's some y in L for which `x (L-y) z` is always composite,
    /// then we can split the family as `x L y (L-y) z`
    pub fn split_on_necessary_digit(&mut self, family: &Family) -> Option<Family> {
        // There's a case in base 11 (and probably others) where we have
        // just one core, where all the digits except one are even, and so
        // is the rest of the number.
        // This tells me that we are required to have at least one of that digit,
        // or else we'll forever be even.
        // This function detects that situation and splits the family accordingly.
        // TODO: does this belong in composite? not quite i think

        let contracted = family.contract().value(self.base);

        'cores: for (i, core) in family.cores.iter().enumerate() {
            // Try substituting everything from all other cores (saves some repeat
            // work in the try_digit closure.)
            let mut gcd_other_cores = contracted.clone();
            for (j, other_core) in family.cores.iter().enumerate() {
                if i == j {
                    continue;
                }

                for d in other_core.iter() {
                    gcd_other_cores = match nontrivial_gcd(
                        &gcd_other_cores,
                        &family.substitute(j, d).value(self.base),
                    ) {
                        Some(g) => g,
                        None => {
                            // there is no hope of finding a necessary digit in core i,
                            // move to the next one
                            continue 'cores;
                        }
                    };
                }
            }

            // Now try each digit in `core` individually and see if any of them work
            'digits: for d in core.iter() {
                let mut g = gcd_other_cores.clone();
                // We want to check substitution by all digits except d
                for d2 in core.iter() {
                    if d2 == d {
                        continue;
                    }

                    g = match nontrivial_gcd(&g, &family.substitute(i, d2).value(self.base)) {
                        Some(g) => g,
                        None => {
                            // nope, give up on this digit (but not this core!)
                            continue 'digits;
                        }
                    }
                }

                // If we got here, then g is a nontrivial common divisor of "substituting anything
                // except (i, d)", and so we must have at least one substitution of (i, d).
                // Loosely speaking, we now return `x L y (L-y) z`
                let mut new = family.clone();
                let d_less_core = core.clone().without(d);

                new.digitseqs.insert(i + 1, d.into());
                new.cores.insert(i + 1, d_less_core);
                debug!("  {family} must have a {d}, transforming into {new}");
                debug_to_tree!(self.tracer, "Must have a {}, transforming to {}", d, new);
                return Some(new);
            }
        }

        None
    }

    /// Given a family `xLz`, with a, b in L, if `xabz` or `xbaz` is forbidden. If so,
    /// we can reduce the family a bit.
    pub fn split_on_incompatible_digits(&mut self, family: &Family) -> Option<Vec<Family>> {
        for (i, core) in family.cores.iter().enumerate() {
            for (a, b) in core.iter().tuple_combinations() {
                if a == b {
                    continue;
                }

                // Check whether we can substitute a and b in either order.
                let seq_ab = family.substitute_multiple(i, [a, b]);
                let seq_ba = family.substitute_multiple(i, [b, a]);

                // TODO: no need to clone this!
                match (
                    self.test_for_contained_prime(&seq_ab).cloned(),
                    self.test_for_contained_prime(&seq_ba).cloned(),
                ) {
                    (Some(p), Some(q)) => {
                        // We can't have both a and b in this core.
                        // We could split into X[aY]Z and X[bY]Z, but that
                        // has potential duplicates. So let's split into three:
                        // X[Y]Z, X[Y]a[aY]Z, X[Y]b[bY]Z
                        assert_ne!(seq_ab, p);
                        assert_ne!(seq_ba, q);
                        debug!("  {seq_ab} contains a prime {p} and {seq_ba} contains a prime {q}");
                        debug_to_tree!(
                            self.tracer,
                            "digits {a} and {b} are incompatible in core {i}"
                        );

                        // Make three children
                        let mut with_neither = family.clone();
                        let mut with_a = family.clone();
                        let mut with_b = family.clone();
                        with_neither.cores[i].remove(a);
                        with_neither.cores[i].remove(b);
                        let neither_core = &with_neither.cores[i];
                        with_a.cores.insert(i, neither_core.clone());
                        with_a.digitseqs.insert(i + 1, DigitSeq(vec![a]));
                        with_a.cores[i + 1].remove(b);
                        with_b.cores.insert(i, neither_core.clone());
                        with_b.digitseqs.insert(i + 1, DigitSeq(vec![b]));
                        with_b.cores[i + 1].remove(a);

                        return Some(vec![with_neither, with_a, with_b]);
                    }
                    (Some(p), None) => {
                        // a can't occur before b
                        assert_ne!(seq_ab, p);
                        debug!("  {seq_ab} contains a prime {p}");
                        debug_to_tree!(
                            self.tracer,
                            "digits {a} and {b} are semi-incompatible in core {i}"
                        );

                        return Some(do_split_for_semi_incompatible(family, i, a, b));
                    }
                    (None, Some(q)) => {
                        // b can't occur before a; converse of the previous branch
                        assert_ne!(seq_ba, q);
                        debug!("  {seq_ba} contains a prime {q}");
                        debug_to_tree!(
                            self.tracer,
                            "digits {b} and {a} are semi-incompatible in core {i}"
                        );

                        // note that b and a are switched!
                        return Some(do_split_for_semi_incompatible(family, i, b, a));
                    }
                    (None, None) => {
                        // nope, nothing we can do
                    }
                }
            }
        }

        None
    }

    /// Given a family `xLyMz`, with a in L, and b in M, if `xaybz` is forbidden,
    /// then we could split the family into `x(L-a)yMz` and `xLy(M-b)z`.
    /// Because this would cause duplication issues though (consider strings with
    /// neither a nor b), we need to split it differently:
    /// - with no a: x(L-a)yMz
    /// - with an a: x(L-a)aLy(M-b)z
    pub fn split_on_incompatible_digits_different_cores(
        &mut self,
        family: &Family,
    ) -> Option<Vec<Family>> {
        // iterate over unordered pairs of cores
        for (j, core_j) in family.cores.iter().enumerate() {
            for (i, core_i) in family.cores.iter().enumerate() {
                if i >= j {
                    break;
                }
                // i < j

                for a in core_i.iter() {
                    for b in core_j.iter() {
                        // Check whether we can substitute in a and b
                        let seq = family.substitute_two(i, a, j, b);

                        if let Some(p) = self.test_for_contained_prime(&seq).cloned() {
                            assert_ne!(seq, p);

                            debug!("  {seq} contains a prime {p}");
                            debug_to_tree!(
                                self.tracer,
                                "digits {a} and {b} are incompatible in slots {i} and {j}"
                            );

                            // We can split the family into two:
                            // - with no a: x(L-a)yMz
                            // - with an a: x(L-a)aLy(M-b)z
                            let mut without_a = family.clone();
                            without_a.cores[i].remove(a);
                            let mut with_a = without_a.clone();
                            with_a.cores[j].remove(b);
                            with_a.digitseqs.insert(i + 1, DigitSeq(vec![a]));
                            with_a.cores.insert(i + 1, family.cores[i].clone());

                            return Some(vec![without_a, with_a]);
                        }
                    }
                }
            }
        }

        None
    }
}

/// Given a family X[abY]Z for which XabZ is forbidden,
/// splits it into:
/// - families with no a: X[bY]Z
/// - families with an a: X[bY]a[aY]Z
///
/// We could reduce to X[bY][aY]Z, but this is leads to huge families.
/// I think it's because XyyyZ can be parsed into that pattern multiple
/// ways.
fn do_split_for_semi_incompatible(family: &Family, i: usize, a: Digit, b: Digit) -> Vec<Family> {
    let mut without_a = family.clone();
    without_a.cores[i].remove(a);
    let mut with_a = without_a.clone();
    with_a.digitseqs.insert(i + 1, DigitSeq(vec![a]));
    with_a
        .cores
        .insert(i + 1, family.cores[i].clone().without(b));
    vec![without_a, with_a]
}
