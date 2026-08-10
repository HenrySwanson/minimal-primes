use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use clap::Parser;
use itertools::Itertools;
use log::{LevelFilter, info};
use num_prime::buffer::PrimeBufferExt;

use crate::families::{Family, SimpleFamily};
use crate::search::{DiesAt, SearchContext, SearchTree, print_stats};
use crate::sequence::Sequence;
use crate::sieve::SieveContext;

mod candidates;
mod digits;
mod families;
mod logging;
mod search;
mod sequence;
mod sieve;
mod tree;

const LOG_EVERY_N: usize = 10_000;

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
    /// Explores the search tree for minimal primes in the given base.
    ///
    /// Will stop once all familes are simple.
    Search(SearchArgs),
    /// Sieves through a sequence of the form k b^n + c.
    Sieve(SieveArgs),
    /// Finds all minimal primes in the given base.
    Solve(SolveArgs),
    /// Reconstructs and prints the search tree from a trace file produced by
    /// `--trace`.
    Tree(TreeArgs),
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
    /// Whether to skip printing the actual primes and branches
    #[arg(long)]
    stats_only: bool,
    /// Log the details of search-tree exploration in detail.
    ///
    /// This will produce a JSONL file under `results/`.
    #[arg(long)]
    trace: bool,
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

#[derive(clap::Args)]
struct SolveArgs {
    /// Base, e.g., decimal, binary, etc.
    base: u8,
    /// upper bound for n
    #[arg(long, default_value_t = 5_000)]
    n_hi: usize,
    /// When doing the first round of sieving, how many primes do we use?
    #[arg(long, default_value_t = 200)]
    p_initial: u64,
    /// What's the ultimate limit on how many primes we should sieve with?
    #[arg(long, default_value_t = 10_000_000)]
    p_max: u64,
    /// Log the details of search-tree exploration in detail.
    ///
    /// This will produce a JSONL file under `results/`.
    #[arg(long)]
    trace: bool,
}

#[derive(clap::Args)]
struct TreeArgs {
    /// Path to a trace file produced by `--trace`.
    path: PathBuf,
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
        Command::Solve(cmd) => {
            do_solve(&cmd, &stop_signal);
        }
        Command::Tree(cmd) => {
            tree::print_tree(&cmd.path).expect("failed to print tree");
        }
    }
}

fn do_search(cmd: &SearchArgs, stop_signal: &AtomicBool) {
    let mut ctx = SearchContext::new(cmd.base, cmd.trace);
    let results = first_stage(&mut ctx, cmd.max_weight, cmd.max_iter, stop_signal);

    if !cmd.stats_only {
        println!("---- BRANCHES REMAINING ----");
        for f in results.simple_families.iter() {
            println!("{f}");
        }
        for f in results.other_families.iter() {
            println!("{f}");
        }
        println!("---- MINIMAL PRIMES ----");
        println!("{}", ctx.primes.clone_and_sort_and_iter().format(", "));
        println!("------------");
        println!(
            "{} primes found, {} simple branches and {} non-simple branches remaining",
            ctx.primes.len(),
            results.simple_families.len(),
            results.other_families.len()
        );
    }

    println!("---- STATS ----");
    print_stats(&ctx.stats);

    println!(
        "Final set of primes ({}): {}",
        ctx.primes.len(),
        ctx.primes.clone_and_sort_and_iter().format(", ")
    );

    println!("{} branches unsolved", results.simple_families.len());
    for x in &results.simple_families {
        println!("{x}");
    }

    if !results.other_families.is_empty() {
        println!("Not all remaining branches are simple! Must bail out now.");
        for x in &results.other_families {
            println!("{x}");
        }
    }
}

fn do_solve(cmd: &SolveArgs, stop_signal: &AtomicBool) -> RemainingNodes {
    let mut ctx = SearchContext::new(cmd.base, cmd.trace);
    let results = first_stage(&mut ctx, None, None, stop_signal);

    println!(
        "{} primes found, {} simple branches and {} non-simple branches remaining",
        ctx.primes.len(),
        results.simple_families.len(),
        results.other_families.len()
    );

    if !results.other_families.is_empty() {
        println!("Not all remaining branches are simple! Must bail out now.");
        return results;
    }

    let mut ctx = SieveContext::from(ctx);
    let unsolved = second_stage(cmd, results.simple_families, &mut ctx);

    println!(
        "Final set of primes ({}): {}",
        ctx.primes.len(),
        ctx.primes.clone_and_sort_and_iter().format(", ")
    );
    println!("{} branches unsolved", unsolved.len());
    for x in &unsolved {
        println!("{x}");
    }

    RemainingNodes {
        simple_families: unsolved,
        other_families: vec![],
    }
}

pub struct RemainingNodes {
    pub simple_families: Vec<SimpleFamily>,
    pub other_families: Vec<Family>,
}

fn first_stage(
    ctx: &mut SearchContext,
    max_weight: Option<usize>,
    max_iter: Option<usize>,
    stop_signal: &AtomicBool,
) -> RemainingNodes {
    let mut tree = SearchTree::new(ctx);

    let mut prev_weight = 0;
    let mut counter = 0;

    loop {
        // Check if we need to stop for any reason
        if stop_signal.load(std::sync::atomic::Ordering::Relaxed) {
            info!("Interrupted! Stopping now...");
            break;
        }

        let weight = match tree.nodes.min_weight() {
            Some(w) => w,
            // this means the frontier is empty!
            None => break,
        };

        if let Some(max_weight) = max_weight
            && weight > max_weight
        {
            info!("Reached weight cutoff; stopping...");
            break;
        }

        if let Some(max_iter) = max_iter
            && ctx.iter >= max_iter
        {
            info!("Reached iteration cutoff; stopping...");
            break;
        }

        if !tree.any_nodes_to_solve() {
            info!("All remaining families are simple; stopping...");
            break;
        }

        // Don't log every single time, that's annoying to read
        let mut should_print = false;
        if weight != prev_weight {
            prev_weight = weight;
            counter = 0;
            should_print = true;
        }

        counter += 1;
        if counter == LOG_EVERY_N {
            counter = 0;
            should_print = true;
        }

        if should_print {
            let num_complex = tree.num_nodes_to_solve();
            let num_simple = tree.nodes.len() - num_complex;
            info!(
                "Weight {} - Iteration {} - {} complex branches - {} simple branches",
                weight, ctx.iter, num_complex, num_simple
            );
        }

        // Now, finally, we can explore a new node
        if tree.explore_once(ctx).is_break() {
            break;
        }
    }

    tree.into_results()
}

/// If the family becomes prime, adds it to `primes`. If it contains another
/// prime, just discards it. If it can't do either of those things, returns
/// the family, incremented to as far as we searched.
fn fast_forward_if_potentially_prime(
    family: SimpleFamily,
    ctx: &mut SieveContext,
) -> Option<SimpleFamily> {
    let dies_at = find_dies_at(&family, ctx);

    let (repeats_until_prime, killer_prime) = match dies_at {
        DiesAt::KilledBy(n, p) => (n, p),
        DiesAt::Unknown => {
            println!("  {family} will not contain any known minimal primes");
            return Some(family);
        }
    };

    // Either this family hits a prime quickly, or after it gets too long,
    // will contain another prime and be discarded. Let's find out which.
    let mut family = family;
    while family.min_repeats < repeats_until_prime {
        // Test if it's prime
        let value = family.value(ctx.base);

        if ctx.prime_buffer.is_prime(&value, None).probably() {
            println!("  Saving {family}, is prime");
            let seq = family.contract();
            ctx.primes.insert(seq);
            return None;
        }

        // not yet, increment and try again
        family.min_repeats += 1;
    }

    // Didn't become prime, discard it
    println!("  {family} killed by {killer_prime}");
    None
}

fn find_dies_at(family: &SimpleFamily, ctx: &mut SieveContext) -> DiesAt {
    let mut dies_at = DiesAt::Unknown;

    for p in ctx.primes.iter() {
        match family.will_contain_at(p) {
            None => {
                // no info, move to the next prime
            }
            Some(n) => {
                // Take the the running minimum of these
                println!("  {family} will contain {p} after {n} more repeats");
                dies_at.update(n, p);

                // We might be able to bail out instantly!
                if n <= family.min_repeats {
                    println!("  Discarding {family}, contains prime {p}");
                    break;
                }
            }
        }
    }

    dies_at
}

fn second_stage(
    cmd: &SolveArgs,
    unsolved_families: Vec<SimpleFamily>,
    ctx: &mut SieveContext,
) -> Vec<SimpleFamily> {
    // It's possible that a simple family can only be expanded a finite amount
    // before it conflicts with a known minimal prime. If so, we should not
    // jump right to sieving, but try to eliminate it quickly.
    println!("---- INTERMEDIARY PHASE ----");

    let unsolved_families: Vec<_> = unsolved_families
        .into_iter()
        .filter_map(|family| fast_forward_if_potentially_prime(family, ctx))
        .collect();

    println!("---- SIEVING PHASE ----");
    let base = ctx.base;
    let (mut remaining_branches, unsievable_branches): (Vec<_>, Vec<_>) = unsolved_families
        .into_iter()
        .map(
            |simple| match Sequence::try_from_family(&simple.bare, base) {
                Ok(seq) => Ok((simple, seq)),
                Err(_) => Err(simple),
            },
        )
        .partition_result();

    // Okay, now we have a collection of simple familes, and the sequences
    // they correspond to. Let's do some sieving.

    // Start slow with a small range
    let mut n_range = 0..16;

    // We'll track some timing stats from the previous round, which will allow
    // us to adaptively pick p_max for the next round. For now, it's None, since
    // there's no history to go on yet.
    // Also put the previous size of n_range in there since we need it.
    let mut prev_round_stats: Option<(sieve::SieveStats, usize)> = None;

    // Put reasonable bounds on p_max, just in cast the adaptive algorithm
    // decides to get silly with it. Remember we have to stay within a u32.
    const MIN_P_MAX: u64 = 100;

    while !remaining_branches.is_empty() {
        // clamp the range
        n_range.end = std::cmp::min(n_range.end, cmd.n_hi);

        // Reached the end of our sieve?
        if n_range.start >= cmd.n_hi {
            println!("Reached limit on n, stopping sieving...");
            break;
        }

        // Pick p_max for this upcoming round by asking the magic oracle.
        // Details are specific to implementation in `sieve` but it depends
        // on timing stats from the previous round.
        let p_max = prev_round_stats
            .as_ref()
            .and_then(|(prev_stats, prev_range_len)| {
                sieve::suggest_next_p_max(
                    prev_stats,
                    *prev_range_len,
                    n_range.len(),
                    remaining_branches.len(),
                )
            })
            .map(|p| p.clamp(MIN_P_MAX, cmd.p_max))
            .unwrap_or(cmd.p_initial);

        let stats = sieve::do_one_round(ctx, &mut remaining_branches, &n_range, p_max);
        prev_round_stats = Some((stats, n_range.len()));

        // Double the range for next time
        n_range = n_range.end..(n_range.end * 2);
    }

    // Before we go, check any of our remaining branches and see if they can be
    // eliminated by an existing minimal prime.
    let mut output = vec![];
    for (family, _) in remaining_branches {
        if let Some(still_unsolved) = fast_forward_if_potentially_prime(family, ctx) {
            output.push(still_unsolved);
        }
    }

    // Definitely also check the unsievable branches, since there's nothing else
    // we can do with them.
    // TODO: we really gotta have some quick way to do this, it happens a lot
    for family in unsievable_branches {
        println!("Checking if we get a lucky break on {family}");
        if let Some(still_unsolved) = fast_forward_if_potentially_prime(family, ctx) {
            output.push(still_unsolved);
        }
    }

    output
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
    use crate::candidates::CandidateSequences;

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
        #[allow(dead_code)]
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
                "FFFA6*", // killed by the above
                "E0*111", // 16416 digits
                "90*G",   // 42996 digits
                "4F0*6",  // 49850 digits
                "FG6*",   // 110986 digits
                "EE16*",  // known to be unsolved
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
                "A*6FFFK", // killed by the above
                "40*9G",   // 47336 digits
                "CF*0K",   // 479150 digits
                "G0*FK",   // unsolved!
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
                "F*GI",     // killed by the above
                "F*KG",     // same
                "F*KAG",    // same
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
                "G*69",     // killed by the above
                "AIF*",     // 21145 digits
                "AIIF*",    // killed by the above
                "K9AE*",    // 23278 digits
                "96E*",     // 25513 digits
                "80*1",     // 119216 digits
                "80*81",    // killed by the above
                "9E*",      // 800874 digits
            ]
        })
    );
    declare_test_for_base!(test_base_24, 24, Status::Complete);
    // Takes 30s, but it works!
    declare_test_for_base!(
        test_base_25,
        25,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "LC*8",     // 509 digits
                "1F0*1",    // 520 digits
                "GOF*I",    // 538 digits
                "LE8*",     // 543 digits
                "L0*91",    // 634 digits
                "LO*66KC",  // 652 digits
                "L08L*",    // 683 digits
                "8C*FFL",   // 711 digits
                "G0F*C",    // 726 digits
                "G0F*CI",   // killed by the above
                "L*O8",     // 743 digits
                "LIC*K6",   // 759 digits
                "M6*9",     // 773 digits
                "9*CM",     // 801 digits
                "8LF*L",    // 948 digits
                "L4E*",     // 1153 digits
                "80*FAG",   // 1155 digits
                "LKI*",     // 1187 digits
                "4F*OOOO",  // 1222 digits
                "60GF*I",   // 1223 digits
                "L0*AI",    // 1223 digits
                "6L*000KI", // 1419 digits
                "AOF*I",    // 1608 digits
                "6L*0KI",   // 1769 digits
                "CMF*6",    // 1912 digits
                "8CF*L",    // 2181 digits
                "8C*FFL",   // killed by the above
                "8L*O",     // 3169 digits
                "LC*KK6",   // 3311 digits
                "C*LKC",    // 4302 digits
                "E1*F1",    // 5510 digits
                "60*LK6",   // 5554 digits
                "ME1*",     // 6393 digits
                "FOK*O",    // 7039 digits
                "FKOK*O",   // killed by the above
                "9MF*9",    // 7138 digits
                "A*FO",     // 7986 digits
                "A*EFO",    // killed by the above
                "O*L8",     // 10177 digits
                "KF*66",    // 10954 digits
                "1*8",      // 11104 digits
                "C1*8",     // killed by the above
                "80*1KO",   // 16610 digits
                "MF1*",     // 16886 digits
                "L*I8",     // 18247 digits
                "G06*FC",   // 18470 digits
                "LF*KI",    // 27169 digits
                "4F*OO",    // 42786 digits
                "LO*KC",    // 66380 digits
                "E*FOO",    // 98399 digits
                "MF*0F6",   // 109992 digits
                "96*M",     // 136967 digits
                "6MF*9",    // unsolved
                "CM1*",     // unsolved
                "EE1*",     // unsolved
                "E1*E",     // unsolved
                "EFO*",     // unsolved
                "F1*F1",    // unsolved
                "F0*KO",    // unsolved
                "F0K*O",    // unsolved
                "LOL*8",    // unsolved
                "M1*F1",    // unsolved
                "M10*8",    // unsolved
                "OL*8",     // unsolved
            ]
        })
    );
    declare_test_for_base!(
        test_base_26,
        26,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "40*GL",  // 512 digits
                "K0*IP",  // 656 digits
                "K0*IPP", // killed by the above
                "G*OO9",  // 1108 digits
                "G*9",    // 1160 digits
                "G*O9",   // killed by the above
                "IG*9",   // same
                "A*06F",  // 1296 digits
                "KIA*F",  // 1301 digits
                "KKIA*F", // killed by the above
                "F*PCF",  // 1572 digits
                "M*P",    // 8773 digits
                "AM*P",   // killed by the above
                "I*GL",   // unsolved!
                "A*6F",   // unsolved!
            ]
        })
    );
    // Takes 30s but does solve!
    declare_test_for_base!(
        test_base_27,
        27,
        Status::IncompleteSimple(IncompleteBranches {
            eventual_primes: vec![
                "L8*",      // 512 digits
                "LM8*",     // killed by the above
                "K*POK",    // 514 digits
                "K*0POK",   // killed by the above
                "M0*QL8",   // 518 digits
                "9*GG",     // 527 digits
                "9*GIG",    // killed by the above
                "9*GGGG",   // same
                "C90*4",    // 538 digits
                "MA*FK",    // 697 digits
                "Q0*FFFA",  // 790 digits
                "Q0*EEFA",  // killed by the above
                "Q0*EEEFA", // same
                "FA*CA",    // 853 digits
                "Q*FFFA",   // 858 digits
                "Q*FFFFFA", // killed by the above
                "A0*LPP",   // 860 digits
                "F*I8",     // 891 digits
                "F*IFF8",   // killed by the above
                "P*4IP",    // 909 digits
                "L0L*6E",   // 920 digits
                "C*94",     // 960 digits
                "F*AIP",    // 1009 digits
                "PP4L0*E",  // 1036 digits
                "I*LE",     // 1067 digits
                "EI*LE",    // killed by the above
                "GGI*LE",   // same
                "M0*KLG",   // 1119 digits
                "P0*OE",    // 1305 digits
                "L*6E",     // 1349 digits
                "L*666E",   // killed by the above
                "QQL*G",    // 1662 digits
                "OLC*M",    // 1889 digits
                "M0*L8",    // 2317 digits
                "QL*E",     // 2561 digits
                "CQCL*E",   // killed by the above
                "Q0*FA",    // 2858 digits
                "4A*PPP",   // 3225 digits
                "FA*IA",    // 4022 digits
                "LPIP*",    // 4442 digits
                "100O*8",   // 4551 digits
                "FKI*K",    // 4825 digits
                "L*IG",     // 5567 digits
                "O16*8",    // 6237 digits
                "P0P*IP",   // 6359 digits
                "Q*FA",     // 7688 digits
                "IL*G",     // 7881 digits
                "4AP*",     // 8885 digits
                "CF*IA",    // 10080 digits
                "A0*F9P",   // 13201 digits
                "Q0*964",   // 17277 digits
                "KL*G",     // 17471 digits
                "KKL*G",    // killed by the above
                "16*8",     // 19397 digits
                "L*GLG",    // 35567 digits
                "ME*9G",    // 49643 digits
                "CA0F*A",   // 88887 digits
                "L*G",      // 101106 digits
                "L*0G",     // killed by the above
                "IL*0G",    // same
                "ICL*G",    // same
                "A0*PM",    // 109006 digits
                "80*9A",    // unsolved
                "999G*",    // unsolved
                "CL*E",     // unsolved
                "EI*F8",    // unsolved
                "F*9FM",    // unsolved
                "9G*",      // this one's composite! alternating squares and div 7
                            // why don't we catch it?
            ]
        })
    );
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
    // This is solvable, but it takes 4 minutes :(
    declare_test_for_base!(test_base_29, 29, Status::Other);
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
                        stats_only: false,
                        trace: false,
                    },
                    &AtomicBool::new(false),
                );
                return;
            }
        };

        // Figure out the equivalent command
        let cmd = SolveArgs {
            base,
            n_hi: 500,
            // seems to work better than p = 1M, should this be backported
            // to the actual CLI command?
            p_initial: 1_000,
            p_max: 1_000_000,
            trace: false,
        };

        let mut ctx = SearchContext::new(base, false);
        let results = first_stage(&mut ctx, None, None, &AtomicBool::new(false));
        let mut ctx = SieveContext::from(ctx);
        let unsolved = second_stage(&cmd, results.simple_families, &mut ctx);

        // Compare the primes we got to the primes we expect, except for the ones we
        // know we're missing.
        compare_primes(
            base,
            &ctx.primes,
            expected_incomplete
                .eventual_primes
                .iter()
                .map(|s| Regex::new(&format!("^{s}$")).unwrap())
                .collect(),
        );

        // We should also check that these eventual primes show up in our unsolved
        // list. Otherwise that'd mean we forgot them somehow.
        let unsolved: Vec<_> = unsolved.iter().map(|f| f.bare.to_string()).collect();
        pretty_assertions::assert_eq!(
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
