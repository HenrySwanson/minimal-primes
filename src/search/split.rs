use itertools::Itertools;
use log::debug;

use crate::candidates::CandidateIndices;
use crate::digits::Digit;
use crate::families::Family;
use crate::search::context::ExploreEvent;
use crate::search::gcd::nontrivial_gcd;
use crate::search::SearchContext;

// TODO: this probably shouldn't be searchcontext, but this works well now
impl SearchContext {
    /// Given a family `xLz`, checks if there's any y in L such that `x y^n z`
    /// is forbidden. If so, we split the family into `x (L-y) (y (L-y))^i z`.
    ///
    /// Generalizes to multi-core families.
    ///
    /// This is Lemma 21 from Bright.
    ///
    /// We check n from 1 to `max_repeats` inclusive.
    pub fn split_on_limited_digit(
        &mut self,
        family: &Family,
        max_repeats: usize,
        possible_contained_primes: &CandidateIndices,
    ) -> Option<(Vec<Family>, ExploreEvent)> {
        for (i, core) in family.cores.iter().enumerate() {
            for d in core.iter() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = family.substitute_multiple(i, std::iter::repeat_n(d, n));
                    if let Some(p) = self.test_for_contained_prime(&seq, possible_contained_primes)
                    {
                        assert_ne!(&seq, p);
                        debug!("  {seq} contains a prime {p}");

                        // Split into n families, x (L-y) (y (L-y))^i z for i in 0..n
                        let yless_core = core.clone().without(d);
                        // xLz -> x(L-y)z
                        let mut first_child = family.clone();
                        first_child.cores[i] = yless_core.clone();

                        let mut children = vec![first_child];

                        while children.len() < n {
                            let mut new = children.last().unwrap().clone();
                            // x(L-y)z -> x(L-y)y(L-y)z
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
                        let event = ExploreEvent::SplitOnLimitedDigit {
                            core_idx: i,
                            digit: d,
                            n,
                        };
                        return Some((children, event));
                    }
                }
            }
        }
        None
    }

    /// Given a family `xLz`, if there's some y in L for which `x (L-y) z` is always composite,
    /// then we can split the family as `x L y (L-y) z`.
    ///
    /// This doesn't reduce the complexity of the cores, so its use should be limited. It
    /// does seem to help in small doses though.
    ///
    /// This is Lemma 29 in Bright.
    pub fn split_on_necessary_digit(&mut self, family: &Family) -> Option<(Family, ExploreEvent)> {
        // There's a case in base 11 (and probably others) where we have
        // just one core, where all the digits except one are even, and so
        // is the rest of the number.
        // This tells me that we are required to have at least one of that digit,
        // or else we'll forever be even.
        // This function detects that situation and splits the family accordingly.

        // TODO: should this call all the composite checks?

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
                let mut new = family.clone();
                let d_less_core = core.clone().without(d);

                // xLz -> xLy(L-y)z
                new.digitseqs.insert(i + 1, d.into());
                new.cores.insert(i + 1, d_less_core);
                debug!("  {family} must have a {d}, transforming into {new}");
                let event = ExploreEvent::SplitOnNecessaryDigit {
                    core_idx: i,
                    digit: d,
                };
                return Some((new, event));
            }
        }

        None
    }

    /// Given a family `xLz`, if there's some a, b in L such that `xabz` or `xbaz`
    /// (or both) is forbidden, we can split the family.
    ///
    /// If `xabz` is forbidden, we could reduce it to `x(L-a)(L-b)z`, but this leads
    /// to huge families. I think this happens because `xcz`, where c is some other
    /// digit in L, can be generated in multiple ways.
    ///
    /// Instead, we split into:
    /// - families with no a: `x(L-a)z`
    /// - families with an a: `x(L-a)a(L-b)z`
    ///
    /// If both `xabz` and `xbaz` are forbidden, we split it into:
    /// - families with neither: `x(L-a-b)z`
    /// - families with an a:    `x(L-a-b)a(L-b)z`
    /// - families with a b:     `x(L-a-b)b(L-a)z`
    ///
    /// This is similar to Lemmas 23 and 25 in Bright.
    pub fn split_on_incompatible_digits(
        &mut self,
        family: &Family,
        possible_contained_primes: &CandidateIndices,
    ) -> Option<(Vec<Family>, ExploreEvent)> {
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
                    self.test_for_contained_prime(&seq_ab, possible_contained_primes)
                        .cloned(),
                    self.test_for_contained_prime(&seq_ba, possible_contained_primes)
                        .cloned(),
                ) {
                    (Some(p), Some(q)) => {
                        // a and b can't co-occur in either order
                        assert_ne!(seq_ab, p);
                        assert_ne!(seq_ba, q);
                        debug!("  {seq_ab} contains a prime {p} and {seq_ba} contains a prime {q}");

                        // Make the family with neither a nor b
                        let mut with_neither = family.clone();
                        // xLz -> x(L-a-b)z
                        with_neither.cores[i].remove(a);
                        with_neither.cores[i].remove(b);
                        let neither_core = &with_neither.cores[i];

                        // Make the families with only a or b
                        let mut with_a = family.clone();
                        let mut with_b = family.clone();

                        // xLz -> x(L-a-b)aLz -> x(L-a-b)a(L-b)z
                        with_a.cores.insert(i, neither_core.clone());
                        with_a.digitseqs.insert(i + 1, a.into());
                        with_a.cores[i + 1].remove(b);

                        // converse
                        with_b.cores.insert(i, neither_core.clone());
                        with_b.digitseqs.insert(i + 1, b.into());
                        with_b.cores[i + 1].remove(a);

                        let event = ExploreEvent::SplitOnIncompatibleSameCore {
                            core_idx: i,
                            first: a,
                            second: b,
                            reverse_also_forbidden: true,
                        };
                        return Some((vec![with_neither, with_a, with_b], event));
                    }
                    (Some(p), None) => {
                        // a can't occur before b
                        assert_ne!(seq_ab, p);
                        debug!("  {seq_ab} contains a prime {p}");

                        return Some(do_split_for_semi_incompatible(family, i, a, b));
                    }
                    (None, Some(q)) => {
                        // b can't occur before a; converse of the previous branch
                        assert_ne!(seq_ba, q);
                        debug!("  {seq_ba} contains a prime {q}");

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
    ///
    /// However, this would cause us to consider strings in the family
    /// `x(L-a)y(M-b)z` twice, so instead, we split it differently:
    /// - with no a: `x(L-a)yMz`
    /// - with an a: `x(L-a)aLy(M-b)z`
    ///
    /// This is similar to Lemma 27 in Bright.
    pub fn split_on_incompatible_digits_different_cores(
        &mut self,
        family: &Family,
        possible_contained_primes: &CandidateIndices,
    ) -> Option<(Vec<Family>, ExploreEvent)> {
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

                        if let Some(p) = self
                            .test_for_contained_prime(&seq, possible_contained_primes)
                            .cloned()
                        {
                            assert_ne!(seq, p);

                            debug!("  {seq} contains a prime {p}");

                            // We can split the family into two:
                            // - with no a: x(L-a)yMz
                            // - with an a: x(L-a)aLy(M-b)z

                            // xLyMz -> x(L-a)yMz
                            let mut without_a = family.clone();
                            without_a.cores[i].remove(a);

                            // x(L-a)yMz -> x(L-a)y(M-b)z
                            let mut with_a = without_a.clone();
                            with_a.cores[j].remove(b);
                            // x(L-a)y(M-b)z -> x(L-a)aLy(M-b)z
                            with_a.digitseqs.insert(i + 1, a.into());
                            with_a.cores.insert(i + 1, family.cores[i].clone());

                            let event = ExploreEvent::SplitOnIncompatibleDifferentCores {
                                core_i: i,
                                core_j: j,
                                a,
                                b,
                            };
                            return Some((vec![without_a, with_a], event));
                        }
                    }
                }
            }
        }

        None
    }

    /// Given a family xLz, with a and b in L, if xabaz is forbidden, splits the family into:
    /// - no as: x(L-a)z
    /// - one a: x(L-a)a(L-a)z
    /// - 2+ as, but no bs between them: x(L-a)a(L-b)a(L-a)z
    ///   - this is unambiguous: the as must be the first and last ones
    ///
    /// This is similar to Lemma 31 in Bright.
    pub fn split_on_forbidden_sandwich(
        &mut self,
        family: &Family,
        possible_contained_primes: &CandidateIndices,
    ) -> Option<(Vec<Family>, ExploreEvent)> {
        // iterate over cores
        for (i, core) in family.cores.iter().enumerate() {
            for (a, b) in core.iter().cartesian_product(core.iter()) {
                if a == b {
                    continue;
                }

                // Check whether aba is forbidden
                let seq = family.substitute_multiple(i, [a, b, a]);
                if let Some(p) = self
                    .test_for_contained_prime(&seq, possible_contained_primes)
                    .cloned()
                {
                    assert_ne!(seq, p);

                    debug!("  {seq} contains a prime {p}");

                    // Split the family
                    // xLz -> x(L-a)z
                    let mut no_as = family.clone();
                    no_as.cores[i].remove(a);
                    let aless_core = &no_as.cores[i];
                    // x(L-a)z -> x(L-a)a(L-a)z
                    let mut one_a = no_as.clone();
                    one_a.digitseqs.insert(i + 1, a.into());
                    one_a.cores.insert(i + 1, aless_core.clone());
                    // x(L-a)a(L-a)z -> x(L-a)a(L-b)a(L-a)z
                    let mut more_as = one_a.clone();
                    more_as.digitseqs.insert(i + 1, a.into());
                    more_as
                        .cores
                        .insert(i + 1, family.cores[i].clone().without(b));

                    let event = ExploreEvent::SplitOnForbiddenSandwich { core_idx: i, a, b };
                    return Some((vec![no_as, one_a, more_as], event));
                }
            }
        }

        None
    }
}

/// Given a family `xLz` for which `xabz` is forbidden,
/// splits it into:
/// - families with no a: `x(L-a)z`
/// - families with an a: `x(L-a)a(L-b)z`
fn do_split_for_semi_incompatible(
    family: &Family,
    i: usize,
    a: Digit,
    b: Digit,
) -> (Vec<Family>, ExploreEvent) {
    // xLz -> x(L-a)z
    let mut without_a = family.clone();
    without_a.cores[i].remove(a);
    // x(L-a)a(L-b)z
    let mut with_a = without_a.clone();
    with_a.digitseqs.insert(i + 1, a.into());
    with_a
        .cores
        .insert(i + 1, family.cores[i].clone().without(b));

    let event = ExploreEvent::SplitOnIncompatibleSameCore {
        core_idx: i,
        first: a,
        second: b,
        reverse_also_forbidden: false,
    };
    (vec![without_a, with_a], event)
}
