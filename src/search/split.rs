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
    pub fn split_on_repeat(&mut self, family: &Family, max_repeats: usize) -> Option<Vec<Family>> {
        for (i, core) in family.cores.iter().enumerate() {
            for d in core.iter() {
                for n in 2..=max_repeats {
                    // Check whether x y^n z contains a prime subword
                    let seq = family.substitute_multiple(i, std::iter::repeat_n(d, n));
                    if let Some(p) = self.test_for_contained_prime(&seq) {
                        assert_ne!(&seq, p);
                        debug!("  {} contains a prime {}", seq, p);

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
        // TODO: generalize to multiple cores!
        // TODO: generalize to multiple digits?
        // TODO: does this belong in composite? not quite i think

        if family.cores.len() != 1 {
            return None;
        }
        let only_core = &family.cores[0];

        if only_core.len() <= 1 {
            return None;
        }

        let contracted = family.contract().value(self.base);

        let try_digit = |d: Digit| -> Option<()> {
            let mut g = contracted.clone();
            // We want to try "no digits" and "all digits except d"
            for d2 in only_core.iter() {
                if d2 == d {
                    continue;
                }

                g = nontrivial_gcd(&g, &family.substitute(0, d2).value(self.base))?;
            }

            Some(())
        };

        for d in family.cores[0].iter() {
            if let Some(()) = try_digit(d) {
                // Got a match! Return xLyLz
                let mut new = family.clone();
                let d_less_core = family.cores[0].clone().without(d);

                new.digitseqs.insert(1, d.into());
                new.cores.insert(1, d_less_core);
                debug!("  {} must have a {}, transforming into {}", family, d, new);
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
                        // Split into X[aY]Z and X[bY]Z
                        assert_ne!(seq_ab, p);
                        assert_ne!(seq_ba, q);
                        debug!(
                            "  {} contains a prime {} and {} contains a prime {}",
                            seq_ab, p, seq_ba, q
                        );
                        debug_to_tree!(
                            self.tracer,
                            "digits {a} and {b} are incompatible in core {i}"
                        );

                        // Best case scenario! Just make two copies, one without a and one
                        // without b.
                        let mut without_a = family.clone();
                        let mut without_b = family.clone();
                        without_a.cores[i].remove(a);
                        without_b.cores[i].remove(b);

                        return Some(vec![without_a, without_b]);
                    }
                    (Some(p), None) => {
                        // a can't occur before b
                        // Reduce to X[bY][aY]Z

                        // TODO: this makes things so much worse! we get families
                        // with oodles of cores, like
                        // DEBUG -  Exploring E[06EG]*[06CF]*[06CG]*[06CEG]*[0CF]*F[0]*[0]*[06A]*[06F]*[06]*[06]*6
                        // Fortunately, the branch above is fine, and actually downright helpful, but
                        // this one makes things pretty unpleasant. Gotta investigate that later.
                        continue;

                        assert_ne!(seq_ab, p);
                        debug!("  {} contains a prime {}", seq_ab, p);
                        debug_to_tree!(
                            self.tracer,
                            "digits {a} and {b} are semi-incompatible in core {i}"
                        );

                        let mut new = family.clone();
                        new.cores[i].remove(b);
                        new.digitseqs.insert(i + 1, DigitSeq::new());
                        new.cores.insert(i + 1, family.cores[i].without(a));

                        return Some(vec![new]);
                    }
                    (None, Some(q)) => {
                        // b can't occur before a
                        // Reduce to X[aY][bY]Z

                        continue; // see above

                        assert_ne!(seq_ba, q);
                        debug!("  {} contains a prime {}", seq_ba, q);
                        debug_to_tree!(
                            self.tracer,
                            "digits {b} and {a} are semi-incompatible in core {i}"
                        );

                        let mut new = family.clone();
                        new.cores[i].remove(a);
                        new.digitseqs.insert(i + 1, DigitSeq::new());
                        new.cores.insert(i + 1, family.cores[i].without(b));

                        return Some(vec![new]);
                    }
                    (None, None) => {
                        // nope, nothing we can do
                    }
                }
            }
        }

        None
    }
}
