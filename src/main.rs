use crate::data_structures::CandidateSequences;
use crate::digits::{Digit, DigitSeq};
use crate::families::{Family, Sequence, SimpleFamily};
use crate::search::{search_for_simple_families, Stats};
use crate::sieve::SequenceSlice;

use clap::Parser;
use itertools::Itertools;
use log::LevelFilter;
use num_bigint::{BigInt, BigUint};
use num_prime::nt_funcs::is_prime;

mod data_structures;
mod digits;
mod families;
mod logging;
mod math;
mod search;
mod sieve;

#[derive(Parser)]
struct Args {
    /// What to do
    #[command(subcommand)]
    command: Command,

    /// Log level
    #[arg(long, global = true, default_value_t = LevelFilter::Warn)]
    log_level: LevelFilter,
}

#[derive(clap::Subcommand)]
enum Command {
    /// Searches for minimal primes in the given base.
    Search(SearchArgs),
    /// Sieves through a sequence of the form k b^n + c.
    Sieve(SieveArgs),
}

#[derive(clap::Args)]
struct SearchArgs {
    /// Base, e.g., decimal, binary, etc.
    base: u8,

    /// Stop exploring when families get above this weight.
    #[arg(long)]
    max_weight: Option<usize>,

    /// Stop exploring after a specific number of iterations.
    #[arg(long)]
    max_iter: Option<usize>,

    /// Switch over to sieving after all remaining families are simple.
    #[arg(long)]
    with_sieve: bool,

    /// upper bound for n
    #[arg(long, default_value_t = 5_000)]
    n_hi: usize,

    /// max p to sieve with
    #[arg(long, default_value_t = 1_000_000)]
    p_max: u64,

    /// whether to log the whole search tree
    #[arg(long)]
    tree_log: bool,
}

#[derive(clap::Args)]
struct SieveArgs {
    /// Base, e.g., decimal, binary, etc.
    base: u8,

    /// k
    k: u64,

    /// c
    c: i64,

    /// lower bound for n
    n_lo: usize,

    /// upper bound for n
    n_hi: usize,

    /// max p to sieve with
    #[arg(default_value_t = 10_000_000)]
    p_max: u64,
}

fn main() {
    let args = Args::parse();

    // Set up logging
    log::set_boxed_logger(Box::new(logging::SimpleLogger)).expect("logger should be uninitialized");
    log::set_max_level(args.log_level);

    match args.command {
        Command::Search(cmd) => {
            do_search(&cmd);
        }
        Command::Sieve(cmd) => {
            do_sieve(&cmd);
        }
    }
}

pub struct SearchResults {
    pub primes: CandidateSequences,
    pub simple_families: Vec<SimpleFamily>,
    pub other_families: Vec<Family>,
    pub stats: Stats,
}

fn do_search(cmd: &SearchArgs) -> SearchResults {
    let mut results = first_stage(cmd);

    if !cmd.with_sieve {
        return results;
    }

    if !results.other_families.is_empty() {
        println!("Not all remaining branches are simple! Must bail out now.");
        return results;
    }

    let unsolved_families =
        intermediate_stage(cmd.base, results.simple_families, &mut results.primes);

    let (primes, unsolved) = second_stage(cmd, unsolved_families, results.primes);

    println!(
        "Final set of primes ({}): {}",
        primes.len(),
        primes.iter().format(", ")
    );
    println!("{} branches unsolved", unsolved.len());
    for x in &unsolved {
        println!("{}", x);
    }

    SearchResults {
        primes,
        simple_families: unsolved,
        other_families: vec![],
        stats: results.stats,
    }
}

fn first_stage(cmd: &SearchArgs) -> SearchResults {
    let results = search_for_simple_families(
        cmd.base,
        cmd.max_weight,
        cmd.max_iter,
        cmd.with_sieve,
        cmd.tree_log,
    );

    println!("---- BRANCHES REMAINING ----");
    for f in results.simple_families.iter() {
        println!("{}", f);
    }
    for f in results.other_families.iter() {
        println!("{}", f);
    }
    println!("---- MINIMAL PRIMES ----");
    println!("{}", results.primes.iter().format(", "));
    println!("------------");
    println!(
        "{} primes found, {} branches unresolved",
        results.primes.len(),
        results.simple_families.len() + results.other_families.len()
    );
    println!("---- STATS ----");
    println!("{} branches explored", results.stats.num_branches_explored);
    println!(
        "{} primality tests ({}ms)",
        results.stats.num_primality_checks,
        results.stats.duration_primality_checks.as_millis()
    );
    println!(
        "{} substring tests ({}ms)",
        results.stats.num_substring_checks,
        results.stats.duration_substring_checks.as_millis()
    );
    println!(
        "{} simple substring tests ({}ms)",
        results.stats.num_simple_substring_checks,
        results.stats.duration_simple_substring_checks.as_millis()
    );

    results
}

fn intermediate_stage(
    base: u8,
    unsolved_families: Vec<SimpleFamily>,
    primes: &mut CandidateSequences,
) -> Vec<SimpleFamily> {
    // It's possible that a simple family can only be expanded a finite amount
    // before it conflicts with a known minimal prime. If so, we should not
    // jump right to sieving, but try to eliminate it quickly.
    println!("---- INTERMEDIARY PHASE ----");

    unsolved_families
        .into_iter()
        .filter_map(|family| intermediate_process_family(base, family, primes))
        .collect()
}

/// If the family becomes prime, adds it to `primes`. If it contains another
/// prime, just discards it. If it can't do either of those things, returns
/// the family, incremented to as far as we searched.
fn intermediate_process_family(
    base: u8,
    family: SimpleFamily,
    primes: &mut CandidateSequences,
) -> Option<SimpleFamily> {
    // Figure out how many iterations we need to check (if any)
    let mut repeats_until_prime = None;
    for p in primes.iter() {
        match family.will_contain_at(p) {
            None => {
                // no luck, move to the next prime
            }
            Some(n) => {
                if n <= family.min_repeats {
                    // we can discard this immediately!
                    println!("  Discarding {}, contains prime {}", family, p);
                    return None;
                }

                // Otherwise, take the the running minimum of these
                println!("  {} will contain {} after {} more repeats", family, p, n);
                repeats_until_prime = Some(repeats_until_prime.map_or(n, |m| n.min(m)));
            }
        }
    }

    let repeats_until_prime = match repeats_until_prime {
        Some(n) => n,
        None => {
            println!("  {} will not contain any known minimal primes", family);
            return Some(family);
        }
    };

    // Either this family hits a prime quickly, or after it gets too long,
    // will contain another prime and be discarded. Let's find out which.
    let mut family = family;
    while family.min_repeats < repeats_until_prime {
        // Test if it's prime
        let value = family.value(base);

        if is_prime(&value, None).probably() {
            println!("  Saving {}, is prime", family);
            let seq = family.sequence();
            primes.insert(seq);
            return None;
        }

        // not yet, increment and try again
        family.min_repeats += 1;
    }

    // Didn't become prime, discard it
    println!("  {} did not become prime, discarding", family);
    None
}

fn second_stage(
    cmd: &SearchArgs,
    unsolved_families: Vec<SimpleFamily>,
    primes: CandidateSequences,
) -> (CandidateSequences, Vec<SimpleFamily>) {
    println!("---- SIEVING PHASE ----");
    let base = cmd.base;
    let mut primes = primes;
    let mut unsolved_branches = vec![];
    let mut remaining_branches: Vec<_> = unsolved_families
        .iter()
        .map(|simple| {
            let seq = Sequence::from_family(&simple, base);
            (simple, seq)
        })
        .collect();

    // Okay, now we have a collection of simple familes, and the sequences
    // they correspond to. Let's do some sieving.
    let mut n_lo = 0;
    let mut n_len = 16;

    while !remaining_branches.is_empty() {
        let mut sequences_to_sieve = vec![];
        let mut slices_to_sieve = vec![];

        let lo = n_lo;
        let hi = std::cmp::min(n_lo + n_len, cmd.n_hi);

        // Real quick, check if this can be eliminated via a minimal prime
        for (simple, seq) in std::mem::take(&mut remaining_branches) {
            if let Some(p) = primes
                .iter()
                .find(|p| simple.will_contain_at(p).is_some_and(|n| n < n_lo))
            {
                println!("{} can be eliminated, since it contains {}", simple, p);
                continue;
            }

            // Do we give up on this sequence?
            if n_lo >= cmd.n_hi {
                println!("Reached limit on n for {}", simple);
                unsolved_branches.push(simple.clone());
                continue;
            }

            sequences_to_sieve.push(simple);
            slices_to_sieve.push(SequenceSlice::new(seq, lo, hi))
        }

        // Now sieve all these slices at once
        println!("Sieving {} families", slices_to_sieve.len());
        sieve::sieve(base, &mut slices_to_sieve, cmd.p_max);

        for (simple, slice) in std::iter::zip(sequences_to_sieve, slices_to_sieve) {
            println!(
                "Investigating family {} for n from {} to {}",
                slice.seq,
                slice.n_lo(),
                slice.n_hi(),
            );
            match sieve::last_resort(base, &slice) {
                Some((i, p)) => {
                    let digitseq =
                        DigitSeq(p.to_radix_be(base.into()).into_iter().map(Digit).collect());
                    println!("Found prime at exponent {}: {}", i, digitseq);
                    primes.insert(digitseq);
                }
                None => {
                    println!("Unable to find prime in the given range: {}", simple);
                    remaining_branches.push((simple, slice.seq))
                }
            }
        }

        n_lo = hi;
        n_len *= 2;
    }

    primes.sort();
    (primes, unsolved_branches)
}

fn do_sieve(cmd: &SieveArgs) {
    let x =
        sieve::find_first_prime(cmd.base, cmd.k, cmd.c, 1, cmd.n_lo, cmd.n_hi, cmd.p_max).unwrap();
    println!("{}, {}", x.0, x.1);
}

#[cfg(test)]
mod tests {
    use regex::Regex;

    use super::*;

    use std::io;

    struct IncompleteBranches {
        /// There are some branches that we know are composite, but the
        /// program can't prove it yet.
        composites: Vec<&'static str>,
        /// These are branches that do eventually become prime, but take
        /// an extremely long time to reach that point, so we don't run
        /// them all the way, just verify that they're still in our search
        /// set.
        eventual_primes: Vec<&'static str>,
    }

    enum Status {
        /// Completely solved; all branches eliminated.
        Complete,
        /// Not complete, we can eliminate all non-simple families,
        /// but some simple families can't be resolved.
        IncompleteSimple(IncompleteBranches),
        /// The number of families grows dramatically, with no signs of
        /// being reducible. This means there's something about our
        /// first stage that can be improved.
        Explodes,
        /// Something else!
        Other,
    }

    macro_rules! declare_test_for_base {
        ($name:ident, $base:literal, $status:expr) => {
            #[test]
            fn $name() {
                test_for_base($base, $status);
            }
        };
    }

    declare_test_for_base!(test_base_2, 2, Status::Complete);
    declare_test_for_base!(test_base_3, 3, Status::Complete);
    declare_test_for_base!(test_base_4, 4, Status::Complete);
    declare_test_for_base!(test_base_5, 5, Status::Complete);
    declare_test_for_base!(test_base_6, 6, Status::Complete);
    declare_test_for_base!(test_base_7, 7, Status::Complete);
    declare_test_for_base!(
        test_base_8,
        8,
        Status::IncompleteSimple(IncompleteBranches {
            composites: vec!["10*1"], // 8^n + 1, sum of cubes
            eventual_primes: vec![]
        })
    );
    declare_test_for_base!(
        test_base_9,
        9,
        Status::IncompleteSimple(IncompleteBranches {
            composites: vec![
                // (9^n - 1) / 8: unusual
                // for even n, this is difference of squares (the 8 doesn't matter)
                // for odd n, this is an even number:
                //   9^(2k+1) mod 16 = 9 * 81^k = 9
                "1*"
            ],
            eventual_primes: vec![]
        })
    );
    declare_test_for_base!(test_base_10, 10, Status::Complete);
    declare_test_for_base!(test_base_11, 11, Status::Complete);
    declare_test_for_base!(test_base_12, 12, Status::Complete);
    declare_test_for_base!(
        test_base_13,
        13,
        Status::IncompleteSimple(IncompleteBranches {
            composites: vec![],
            // 32021 digits :(
            eventual_primes: vec!["80*111"]
        })
    );
    declare_test_for_base!(test_base_14, 14, Status::Complete);
    declare_test_for_base!(test_base_15, 15, Status::Complete);
    declare_test_for_base!(
        test_base_16,
        16,
        Status::IncompleteSimple(IncompleteBranches {
            composites: vec![
                // (4*16^(n+1) - 49) / 15
                // difference of squares
                "4*1",
                // (8*16^(n+1) + 15^2 - 8*16) / 15
                // = (2^(4n+7)+97) / 15
                // unknown
                "8*F", //
                // 9*16^n - 1
                // difference of squares
                "8F*", //
                // (2535*16^(n+2) - 1215) / 15
                // unknown
                "A8F*AF",
                // (2153*16^(n+1) + 97) / 15
                // unknown, but there's that 97 again
                "8F8*F",
            ],
            eventual_primes: vec!["88F*", "90*91", "F8*F"]
        })
    );
    declare_test_for_base!(test_base_17, 17, Status::Explodes);
    declare_test_for_base!(test_base_18, 18, Status::Complete);
    declare_test_for_base!(test_base_19, 19, Status::Explodes);
    // 20 fails to understand that [G]*[I]*9 is composite, and expands it forever.
    declare_test_for_base!(test_base_20, 20, Status::Other);
    // ramps up to tens of thousands of branches, but then drops to 400ish and
    // oscillates for a while. unfortunately, it starts going back up again.
    declare_test_for_base!(test_base_21, 21, Status::Explodes);
    declare_test_for_base!(test_base_22, 22, Status::Explodes);
    declare_test_for_base!(test_base_23, 23, Status::Explodes);
    declare_test_for_base!(
        test_base_24,
        24,
        Status::IncompleteSimple(IncompleteBranches {
            // (6*24^(n+1) - 121) / 23
            // (2^(3n+2) * 3^(n+2) - 11^2) / 23
            // feels like difference of squares on even n, something
            // else on odd n
            composites: vec!["6*1"],
            eventual_primes: vec![]
        })
    );
    declare_test_for_base!(test_base_25, 25, Status::Explodes);
    declare_test_for_base!(test_base_26, 26, Status::Explodes);
    declare_test_for_base!(test_base_27, 27, Status::Explodes);
    declare_test_for_base!(test_base_28, 28, Status::Explodes);
    declare_test_for_base!(test_base_29, 29, Status::Explodes);
    // 30 is solvable, but takes too long, even with --release
    declare_test_for_base!(
        test_base_30,
        30,
        Status::IncompleteSimple(IncompleteBranches {
            composites: vec![],
            eventual_primes: vec!["C0*1"] // 1024 digits
        })
    );

    fn test_for_base(base: u8, status: Status) {
        // Complete is just a shortcut for "nothing incomplete", so let's
        // reduce our casework.
        let expected_incomplete = match status {
            Status::Complete => IncompleteBranches {
                composites: vec![],
                eventual_primes: vec![],
            },
            Status::IncompleteSimple(branches) => branches,
            Status::Explodes | Status::Other => {
                // this is real bad. just run it a short amount to make sure
                // it's not panicking or anything.
                // TODO: compare the results of this to the expected ones,
                // to make sure we're not discovering fake primes or anything
                do_search(&SearchArgs {
                    base: base,
                    max_weight: Some(5),
                    max_iter: Some(10_000),
                    with_sieve: false,
                    n_hi: 0,
                    p_max: 0,
                    tree_log: false,
                });
                return;
            }
        };

        // Figure out the equivalent command
        let cmd = SearchArgs {
            base: base,
            max_weight: None,
            max_iter: None,
            with_sieve: true,
            n_hi: 500,
            // seems to work better than p = 1M, should this be backported
            // to the actual CLI command?
            p_max: 1_000,
            tree_log: false,
        };

        // First stage
        let mut results = first_stage(&cmd);

        // Remove any composite branches that are expected to be present.
        for composite in expected_incomplete.composites {
            match results
                .simple_families
                .iter()
                .position(|family| family.pattern() == composite)
            {
                Some(i) => {
                    results.simple_families.remove(i);
                }
                None => panic!(
                    "Expected to find composite family {}, but it was not present",
                    composite
                ),
            }
        }

        // Do intermediate and second stages
        let unsolved_families =
            intermediate_stage(cmd.base, results.simple_families, &mut results.primes);

        let (primes, unsolved) = second_stage(&cmd, unsolved_families, results.primes);

        // Compare the primes we got to the primes we expect, except for the ones we
        // know we're missing.
        compare_primes(
            base,
            &primes,
            expected_incomplete
                .eventual_primes
                .iter()
                .map(|s| Regex::new(&format!("^{s}$")).unwrap())
                .collect(),
        );

        // We should also check that these eventual primes show up in our unsolved
        // list. Otherwise that'd mean we forgot them somehow.
        let mut unsolved: Vec<_> = unsolved.iter().map(|f| f.pattern()).collect();
        unsolved.sort();
        assert_eq!(unsolved, expected_incomplete.eventual_primes);

        //     match status {
        //         Status::Complete => {
        //             // Simulate it for the full duration
        //             let max_weight = get_max_weight(base);
        //             let results = calculate_full(base, max_weight);
        //             compare_primes(base, &results.primes, true);
        //             assert!(
        //                 results.simple_families.is_empty() && results.other_families.is_empty(),
        //                 "Some branches were not eliminated!\n{}\n{}",
        //                 results.simple_families.iter().format("\n"),
        //                 results.other_families.iter().format("\n")
        //             );
        //         }
        //         Status::StrayBranches { unresolved } => {
        //             // Simulate it for the full duration
        //             let max_weight = get_max_weight(base) + 1;
        //             let results = calculate(base, max_weight);
        //             compare_primes(base, &results.primes, false);

        //             assert!(
        //                 results.other_families.is_empty(),
        //                 "Got some non-simple families:\n{}",
        //                 results.other_families.iter().format("\n")
        //             );

        //             let mut simple_strings: Vec<_> = results
        //                 .simple_families
        //                 .iter()
        //                 // TODO: better method here
        //                 .map(|x| x.to_string().split(" -- ").next().unwrap().to_owned())
        //                 .collect();
        //             simple_strings.sort();

        //             assert_eq!(
        //                 simple_strings,
        //                 unresolved,
        //                 "Didn't get the expected unsolved branches:\n{}\nvs\n{}",
        //                 results.simple_families.iter().format("\n"),
        //                 results.other_families.iter().format("\n")
        //             );
        //         }
        //         Status::Unsolved { max_weight } => {
        //             let results = calculate(base, max_weight);
        //             compare_primes(base, &results.primes, false);
        //             assert!(
        //                 !results.simple_families.is_empty() || !results.other_families.is_empty(),
        //                 "All branches were eliminated, this test should be marked Complete!"
        //             );
        //         }
        //     }
    }

    fn calculate(base: u8, max_weight: usize) -> SearchResults {
        search_for_simple_families(base, Some(max_weight), None, false, false)
    }

    fn calculate_full(base: u8, max_weight: usize) -> SearchResults {
        do_search(&SearchArgs {
            base,
            max_weight: Some(max_weight),
            max_iter: None,
            with_sieve: true,
            n_hi: 5000,
            p_max: 1_000_000,
            tree_log: false,
        })
    }

    fn iter_ground_truth(base: u8) -> impl Iterator<Item = String> {
        use io::BufRead;
        let file_path = format!("mepn-data/minimal.{}.txt", base);
        let file = std::fs::File::open(file_path).expect("open ground truth file");
        io::BufReader::new(file)
            .lines()
            .map(|line| line.expect("read line"))
    }

    fn get_max_weight(base: u8) -> usize {
        iter_ground_truth(base)
            .map(|s| s.len())
            .max()
            .expect("nonempty ground truth")
    }

    fn compare_primes(base: u8, primes: &CandidateSequences, exceptions: Vec<Regex>) {
        let mut truth_iter = iter_ground_truth(base).peekable();
        let mut iter = primes.iter().map(|seq| seq.to_string()).peekable();

        let mut fail = false;
        loop {
            // We want to consume both iterators until they're both gone.
            //
            // If there's a mismatch, we need to be intelligent about which
            // iterator to increase (the one with the smaller prime.)
            // It turns out to be simpler to combine th
            let cmp = match (truth_iter.peek(), iter.peek()) {
                // If we get two different primes, we only want to advance
                // the iterator with the lesser prime, so that we can re-sync.
                (Some(p_truth), Some(p_got)) => {
                    // First compare them by length, then lexicographically;
                    // otherwise 9 will compare larger than 10.
                    p_truth
                        .len()
                        .cmp(&p_got.len())
                        .then_with(|| p_truth.cmp(p_got))
                }
                // It turns out to be easier to combine this case
                // with the one above by treating None as the biggest
                // prime.
                (None, Some(_)) => std::cmp::Ordering::Greater,
                // same
                (Some(_), None) => std::cmp::Ordering::Less,
                // both exhausted! break
                (None, None) => break,
            };
            match cmp {
                // truth < actual, so we skipped something we were
                // supposed to see. this might be okay, if it's in
                // our exception list
                std::cmp::Ordering::Less => {
                    let p = truth_iter.next().unwrap();
                    if exceptions.iter().any(|re| re.is_match(&p)) {
                        // great, allowable exception
                    } else {
                        println!("Didn't see expected prime: {}", p);
                        fail = true;
                    }
                }
                // truth > actual, so there's something in actual that
                // shouldn't be there. this is always a failure.
                std::cmp::Ordering::Greater => {
                    let p = iter.next().unwrap();
                    println!("Got extra unexpected prime: {}", p);
                    fail = true;
                }
                // all is well, just increment both iterators
                std::cmp::Ordering::Equal => {
                    iter.next();
                    truth_iter.next();
                }
            }
        }

        assert!(!fail, "Some mismatches between primes");
    }
}
