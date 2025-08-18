use std::ops::ControlFlow;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;

use clap::Parser;
use itertools::Itertools;
use log::{info, LevelFilter};
use num_prime::nt_funcs::is_prime;

use crate::data_structures::CandidateSequences;
use crate::digits::{Digit, DigitSeq};
use crate::families::{Family, Sequence, SimpleFamily};
use crate::search::{SearchTree, Stats};
use crate::sieve::SequenceSlice;

mod data_structures;
mod digits;
mod families;
mod logging;
mod search;
mod sieve;

#[derive(Parser)]
struct Args {
    /// What to do
    #[command(subcommand)]
    command: Command,

    /// Log level
    #[arg(long, global = true, default_value_t = LevelFilter::Info)]
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
    /// Do not continue onwards to sieving.
    #[arg(long)]
    first_stage_only: bool,
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

    // Set up Ctrl-C handling
    let stop_signal = Arc::new(AtomicBool::new(false));
    ctrlc::set_handler({
        let stop_signal = stop_signal.clone();
        move || {
            println!("Ctrl-C received! Performing cleanup...");
            let oldval = stop_signal.swap(true, std::sync::atomic::Ordering::Relaxed);
            if oldval {
                // not the first time this was pressed, killing the process
                std::process::exit(128 + 2);
            }
        }
    })
    .expect("Error setting Ctrl-C handler");

    match args.command {
        Command::Search(cmd) => {
            do_search(&cmd, &stop_signal);
        }
        Command::Sieve(cmd) => {
            do_sieve(&cmd);
        }
    }
}

fn do_search(cmd: &SearchArgs, stop_signal: &AtomicBool) -> SearchResults {
    let mut results = first_stage(
        cmd.base,
        cmd.max_weight,
        cmd.max_iter,
        cmd.tree_log,
        stop_signal,
    );

    if cmd.first_stage_only {
        println!("Done with first stage, stopping early!");
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
        primes.clone_and_sort_and_iter().format(", ")
    );
    println!("{} branches unsolved", unsolved.len());
    for x in &unsolved {
        println!("{x}");
    }

    SearchResults {
        primes,
        simple_families: unsolved,
        other_families: vec![],
        stats: results.stats,
    }
}

pub struct SearchResults {
    pub primes: CandidateSequences,
    pub simple_families: Vec<SimpleFamily>,
    pub other_families: Vec<Family>,
    pub stats: Stats,
}

fn first_stage(
    base: u8,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    tree_log: bool,
    stop_signal: &AtomicBool,
) -> SearchResults {
    let mut tree = SearchTree::new(base, tree_log);

    tree.explore_until(|tree| {
        if stop_signal.load(std::sync::atomic::Ordering::Relaxed) {
            info!("Interrupted! Stopping now...");
            return ControlFlow::Break(());
        }

        let weight = match tree.frontier.min_weight() {
            Some(w) => w,
            // this means the frontier is empty!
            None => return ControlFlow::Break(()),
        };

        if let Some(max) = max_weight {
            if weight > max {
                info!("Reached weight cutoff; stopping...");
                return ControlFlow::Break(());
            }
        }

        if let Some(max) = max_iter {
            if tree.ctx.iter >= max {
                info!("Reached iteration cutoff; stopping...");
                return ControlFlow::Break(());
            }
        }

        if tree.all_nodes_simple_and_checked() {
            info!("All remaining families are simple; stopping...");
            return ControlFlow::Break(());
        }

        info!(
            "Iteration {} - Weight {} - {} branches",
            tree.ctx.iter,
            weight,
            tree.frontier.len()
        );

        ControlFlow::Continue(())
    });

    let results = tree.into_results();

    println!("---- BRANCHES REMAINING ----");
    for f in results.simple_families.iter() {
        println!("{f}");
    }
    for f in results.other_families.iter() {
        println!("{f}");
    }
    println!("---- MINIMAL PRIMES ----");
    println!("{}", results.primes.clone_and_sort_and_iter().format(", "));
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
    let branch_stats = &results.stats.branch_stats;
    println!(
        "{} branches eliminated with leading zeros",
        branch_stats.leading_zeros
    );
    println!(
        "{} branches eliminated for containing a prime",
        branch_stats.contains_prime
    );
    println!(
        "{} branches eliminated by discovering a new prime",
        branch_stats.is_new_prime
    );
    println!(
        "{} branches eliminated for compositeness",
        branch_stats.detected_composite
    );
    println!(
        "{} branches simplified into simple families",
        branch_stats.simplified
    );
    println!(
        "{} branches split on a limited digit",
        branch_stats.split_on_limited_digit
    );
    println!(
        "{} branches split on incompatible digits (same core)",
        branch_stats.split_on_incompatible_same_core
    );
    println!(
        "{} branches split on incompatible digits (different cores)",
        branch_stats.split_on_incompatible_different_cores
    );
    println!(
        "{} branches split on a necessary digit",
        branch_stats.split_on_necessary_digit
    );
    println!(
        "{} branches explored generically",
        branch_stats.explored_generically
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
                    println!("  Discarding {family}, contains prime {p}");
                    return None;
                }

                // Otherwise, take the the running minimum of these
                println!("  {family} will contain {p} after {n} more repeats");
                repeats_until_prime = Some(repeats_until_prime.map_or(n, |m| n.min(m)));
            }
        }
    }

    let repeats_until_prime = match repeats_until_prime {
        Some(n) => n,
        None => {
            println!("  {family} will not contain any known minimal primes");
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
            println!("  Saving {family}, is prime");
            let seq = family.sequence();
            primes.insert(seq);
            return None;
        }

        // not yet, increment and try again
        family.min_repeats += 1;
    }

    // Didn't become prime, discard it
    println!("  {family} did not become prime, discarding");
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
    let (mut remaining_branches, unsievable_branches): (Vec<_>, Vec<_>) = unsolved_families
        .into_iter()
        .map(|simple| match Sequence::try_from_family(&simple, base) {
            Ok(seq) => Ok((simple, seq)),
            Err(_) => Err(simple),
        })
        .partition_result();

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
        // TODO: shouldn't this come _after_ sieving?
        for (simple, seq) in std::mem::take(&mut remaining_branches) {
            if let Some(p) = primes
                .iter()
                .find(|p| simple.will_contain_at(p).is_some_and(|n| n < n_lo))
            {
                println!("{simple} can be eliminated, since it contains {p}");
                continue;
            }

            // Do we give up on this sequence?
            if n_lo >= cmd.n_hi {
                println!("Reached limit on n for {simple}");
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
                    println!("Found prime at exponent {i}: {digitseq}");
                    primes.insert(digitseq);
                }
                None => {
                    println!("Unable to find prime in the given range: {simple}");
                    remaining_branches.push((simple, slice.seq))
                }
            }
        }

        n_lo = hi;
        n_len *= 2;
    }

    // Before we go, check whether any of the unsievable branches can be eliminated with
    // a minimal prime.
    // TODO: we really gotta have some quick way to do this, it happens a lot
    for family in unsievable_branches {
        println!("Checking if we get a lucky break on {family}");
        if let Some(still_unsolved) = intermediate_process_family(base, family, &mut primes) {
            unsolved_branches.push(still_unsolved);
        }
    }

    (primes, unsolved_branches)
}

fn do_sieve(cmd: &SieveArgs) {
    let x =
        sieve::find_first_prime(cmd.base, cmd.k, cmd.c, 1, cmd.n_lo, cmd.n_hi, cmd.p_max).unwrap();
    println!("{}, {}", x.0, x.1);
}

#[cfg(test)]
mod tests {
    use std::io;

    use regex::Regex;

    use super::*;

    struct IncompleteBranches {
        /// There are some branches that we know are composite, but the
        /// program can't prove it yet.
        // composites: Vec<&'static str>,
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
    declare_test_for_base!(test_base_8, 8, Status::Complete);
    declare_test_for_base!(test_base_9, 9, Status::Complete);
    declare_test_for_base!(test_base_10, 10, Status::Complete);
    declare_test_for_base!(test_base_11, 11, Status::Complete);
    declare_test_for_base!(test_base_12, 12, Status::Complete);
    declare_test_for_base!(
        test_base_13,
        13,
        Status::IncompleteSimple(IncompleteBranches {
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
            eventual_primes: vec!["88F*", "90*91", "F8*F"]
        })
    );
    declare_test_for_base!(
        test_base_17,
        17,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "A0*1", // 1357 digits
                "A*GF", // 2016 digits
                "6*E9", // 4663 digits
                "1F*",  // 7093 digits
                "49*",  // 111334 digits!
                "F19*", // unsolved!
            ],
        })
    );
    declare_test_for_base!(test_base_18, 18, Status::Complete);
    declare_test_for_base!(
        test_base_19,
        19,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "6*FA",   // 507 digits
                "8C96*",  // 626 digits
                "40G*1",  // 1334 digits
                "E9*E",   // 1465 digits
                "G*1",    // 2035 digits
                "FA6*",   // 9293 digits
                "E0*111", // 16416 digits
                "90*G",   // 42996 digits
                "4F0*6",  // 49850 digits
                "FG6*",   // 110986 digits
                "EE16*",  // known to be unsolved
                // this is a tricky one! it gets eliminated
                // by FA6*, but that one's not discovered yet
                "FFFA6*",
            ]
        })
    );
    declare_test_for_base!(test_base_20, 20, Status::Complete);
    declare_test_for_base!(
        test_base_21,
        21,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "A9F*K",   // 605 digits
                "A*6FK",   // 1634 digits
                "40*9G",   // 47336 digits
                "CF*0K",   // 479150 digits
                "G0*FK",   // unsolved!
                "A*6FFFK"  // killed by A*6FK
            ],
        })
    );
    declare_test_for_base!(
        test_base_22,
        22,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "I*AF",   // 628 digits
                "K0*EC1", // 764 digits
            ]
        })
    );
    declare_test_for_base!(
        test_base_23,
        23,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "C0*MLC",   // 507 digits
                "IF*A",     // 527 digits
                "K*EL",     // 582 digits
                "FA*CC",    // 924 digits
                "F*G",      // 1092 digits
                "FFFFFK*C", // 1119 digits
                "9F*A",     // 1308 digits
                "696E*",    // 1358 digits
                "G*09",     // 1381 digits
                "E0*KLE",   // 1658 digits
                "E*L6",     // 1713 digits
                "9EE6E*",   // 2187 digits
                "IK*FFF",   // 2605 digits
                "K*LLF",    // 2808 digits
                "KLF*",     // 2874 digits
                "EL*6",     // 3261 digits
                "K*L",      // 3762 digits
                "FFFK*C",   // 4465 digits
                "F*KC",     // 5569 digits
                "IIE*L",    // 8122 digits
                "G*9",      // 9526 digits
                "AIF*",     // 21145 digits
                "K9AE*",    // 23278 digits
                "96E*",     // 25513 digits
                "80*1",     // 119216 digits
                "9E*",      // 800874 digits
                "F*GI",     // eventually gets killed by F*G
                "F*KG",     // same
                "F*KAG",    // same
                "G*69",     // eventually gets killed by G*9
                "AIIF*",    // eventually gets killed by AIF*
                "80*81",    // eventually gets killed by 80*1
            ]
        })
    );
    declare_test_for_base!(test_base_24, 24, Status::Complete);
    // It works but takes a long time
    declare_test_for_base!(test_base_25, 25, Status::Other);
    declare_test_for_base!(
        test_base_26,
        26,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "40*GL",  // 512 digits
                "K0*IP",  // 656 digits
                "G*OO9",  // 1108 digits
                "G*9",    // 1160 digits
                "A*06F",  // 1296 digits
                "KIA*F",  // 1301 digits
                "F*PCF",  // 1572 digits
                "M*P",    // 8773 digits
                "I*GL",   // unsolved!
                "A*6F",   // unsolved!
                "AM*P",   // eventually killed by M*P
                "K0*IPP", // eventually killed by K0*IP
                "KKIA*F", // eventually killed by KIA*F
                "G*O9",   // eventually killed by G*9
                "IG*9",   // same
            ]
        })
    );
    // takes a really long time (first stage ends at weight 100!)
    // but does eventually start sieving.
    // the ones that takes the longest to resolve look like
    // K[0]*K0KKK...KKKF6A
    // [0K]*K0KKK...KKKF6A
    // [0K]*K[0]*KKK...KKKF6A
    declare_test_for_base!(test_base_27, 27, Status::Other);
    declare_test_for_base!(
        test_base_28,
        28,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "4O*09", // 617 digits
                "LK*F",  // 927 digits
                "A*6F",  // 1425 digits
                "QO*69", // 4242 digits
                "O4O*9", // 94538 digits
                "OA*F",  // unsolved!
            ]
        })
    );
    // after about a half hour and 860k iterations, it's starting
    // to go down from 130k branches, but it's not done yet.
    // after an hour and 1.3M iterations, it's stalled out around
    // 6k branches at weight 29, and going back up :(
    declare_test_for_base!(test_base_29, 29, Status::Explodes);
    declare_test_for_base!(
        test_base_30,
        30,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec!["C0*1"] // 1024 digits
        })
    );

    fn test_for_base(base: u8, status: Status) {
        // Complete is just a shortcut for "nothing incomplete", so let's
        // reduce our casework.
        let expected_incomplete = match status {
            Status::Complete => IncompleteBranches {
                eventual_primes: vec![],
            },
            Status::IncompleteSimple(branches) => branches,
            Status::Explodes | Status::Other => {
                // this is real bad. just run it a short amount to make sure
                // it's not panicking or anything.
                // TODO: compare the results of this to the expected ones,
                // to make sure we're not discovering fake primes or anything
                do_search(
                    &SearchArgs {
                        base,
                        max_weight: Some(5),
                        max_iter: Some(10_000),
                        first_stage_only: true,
                        n_hi: 0,
                        p_max: 0,
                        tree_log: false,
                    },
                    &AtomicBool::new(false),
                );
                return;
            }
        };

        // Figure out the equivalent command
        let cmd = SearchArgs {
            base,
            max_weight: None,
            max_iter: None,
            first_stage_only: false,
            n_hi: 500,
            // seems to work better than p = 1M, should this be backported
            // to the actual CLI command?
            p_max: 1_000,
            tree_log: false,
        };

        // First stage
        let mut results = first_stage(base, None, None, false, &AtomicBool::new(false));

        // Remove any composite branches that are expected to be present.
        // TODO: all composites are detected right now, but re-use this for
        // unsolved families maybe?
        // for composite in expected_incomplete.composites {
        //     match results
        //         .simple_families
        //         .iter()
        //         .position(|family| family.pattern() == composite)
        //     {
        //         Some(i) => {
        //             results.simple_families.remove(i);
        //         }
        //         None => panic!(
        //             "Expected to find composite family {}, but it was not present",
        //             composite
        //         ),
        //     }
        // }

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
        let unsolved: Vec<_> = unsolved.iter().map(|f| f.pattern()).collect();
        assert_eq!(
            sort_and_dedup(unsolved),
            sort_and_dedup(expected_incomplete.eventual_primes)
        );
    }

    fn iter_ground_truth(base: u8) -> impl Iterator<Item = String> {
        use io::BufRead;
        let file_path = format!("mepn-data/minimal.{base}.txt");
        let file = std::fs::File::open(file_path).expect("open ground truth file");
        io::BufReader::new(file)
            .lines()
            .map(|line| line.expect("read line"))
    }

    fn compare_primes(base: u8, primes: &CandidateSequences, exceptions: Vec<Regex>) {
        let mut truth_iter = iter_ground_truth(base).peekable();
        let mut iter = primes
            .clone_and_sort_and_iter()
            .map(|seq| seq.to_string())
            .peekable();

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
                        println!("Didn't see expected prime: {p}");
                        fail = true;
                    }
                }
                // truth > actual, so there's something in actual that
                // shouldn't be there. this is always a failure.
                std::cmp::Ordering::Greater => {
                    let p = iter.next().unwrap();
                    println!("Got extra unexpected prime: {p}");
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

    // TODO: we shouldn't be getting duplicates in the first place!
    // so we shouldn't need to dedup
    fn sort_and_dedup<T>(mut list: Vec<T>) -> Vec<T>
    where
        T: Ord,
    {
        list.sort();
        list.dedup();
        list
    }
}
