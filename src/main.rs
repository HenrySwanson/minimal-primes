use clap::Parser;
use data_structures::CandidateSequences;
use itertools::Itertools;
use log::LevelFilter;
use num_bigint::{BigInt, BigUint};
use search::{is_substring_of_simple, search_for_simple_families, SearchContext};
use sequences::{Digit, DigitSeq};
use sieve::{Sequence, SequenceSlice};
use std::{ops::ControlFlow, sync::atomic::AtomicBool, usize};

mod composite;
mod data_structures;
mod logging;
mod math;
mod search;
mod sequences;
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

fn do_search(cmd: &SearchArgs) -> (CandidateSequences, Vec<search::SearchNode>) {
    let mut ctx = first_stage(cmd);

    if !cmd.with_sieve {
        return (ctx.primes, ctx.frontier.iter().cloned().collect());
    }

    if !ctx.frontier.all_simple() {
        println!("Not all remaining branches are simple! Must bail out now.");
        return (ctx.primes, ctx.frontier.iter().cloned().collect());
    }

    while intermediate_stage(&mut ctx).is_continue() {}

    let (primes, unsolved) = second_stage(cmd, ctx);

    println!(
        "Final set of primes ({}): {}",
        primes.len(),
        primes.iter().format(", ")
    );
    println!("{} branches unsolved", unsolved.len());
    for x in &unsolved {
        println!("{}", x);
    }

    (primes, unsolved)
}

fn first_stage(cmd: &SearchArgs) -> search::SearchContext {
    let ctx = search_for_simple_families(cmd.base, cmd.max_weight, cmd.max_iter, cmd.with_sieve);

    println!("---- BRANCHES REMAINING ----");
    for f in ctx.frontier.iter() {
        println!("{}", f);
    }
    println!("---- MINIMAL PRIMES ----");
    println!("{}", ctx.primes.iter().format(", "));
    println!("------------");
    println!(
        "{} primes found, {} branches unresolved",
        ctx.primes.len(),
        ctx.frontier.len()
    );
    println!("---- STATS ----");
    println!(
        "{} branches explored",
        ctx.stats.borrow().num_branches_explored
    );
    println!(
        "{} primality tests ({}ms)",
        ctx.stats.borrow().num_primality_checks,
        ctx.stats.borrow().duration_primality_checks.as_millis()
    );
    println!(
        "{} substring tests ({}ms)",
        ctx.stats.borrow().num_substring_checks,
        ctx.stats.borrow().duration_substring_checks.as_millis()
    );
    println!(
        "{} simple substring tests ({}ms)",
        ctx.stats.borrow().num_simple_substring_checks,
        ctx.stats
            .borrow()
            .duration_simple_substring_checks
            .as_millis()
    );

    ctx
}

fn intermediate_stage(ctx: &mut SearchContext) -> ControlFlow<(), ()> {
    println!("---- INTERMEDIARY PHASE ----");
    // It's possible that a simple family can only be expanded a finite amount
    // before it conflicts with a known minimal prime. If so, we should not
    // jump right to sieving, but should instead continue normal searching.

    // Check if any family is potentially able to contain any prime.
    let should_search = ctx.frontier.iter().any(|family| {
        let simple = match family {
            search::SearchNode::Simple(s) => s,
            _ => unreachable!("found non-simple family after all_simple()"),
        };

        ctx.primes
            .iter()
            .any(|p| match is_substring_of_simple(p, simple) {
                search::SubstringResult::Never => false,
                search::SubstringResult::Yes => true,
                search::SubstringResult::Eventually(_) => true,
            })
    });

    if should_search {
        println!(
            "Searching one extra round: iter={}, num primes={}, num_branches={}",
            ctx.iter,
            ctx.primes.len(),
            ctx.frontier.len()
        );
        ctx.search_one_level();
        ControlFlow::Continue(())
    } else {
        ControlFlow::Break(())
    }
}

fn second_stage(
    cmd: &SearchArgs,
    ctx: SearchContext,
) -> (CandidateSequences, Vec<search::SearchNode>) {
    println!("---- SIEVING PHASE ----");
    let base = ctx.base;
    let mut primes = ctx.primes;
    let mut unsolved_branches = vec![];
    let mut remaining_branches: Vec<_> = ctx
        .frontier
        .iter()
        .map(|family| {
            let simple = match &family {
                search::SearchNode::Simple(x) => x,
                _ => unreachable!("found non-simple family after all_simple()"),
            };

            // Compute the sequence for this family: xy*z
            let x = simple.before.value(base);
            let y = simple.center.0;
            let z = simple.after.value(base);

            let b_z = BigUint::from(base).pow(simple.after.0.len() as u32);
            let d = u64::from(base) - 1;
            let k = (x * d + y) * &b_z;
            let c = BigInt::from(d * z) - BigInt::from(y * b_z);

            // Try to fit it into the appropriate ranges
            let k = u64::try_from(k)
                .unwrap_or_else(|e| panic!("Can't convert {} to u64", e.into_original()));
            let c = i64::try_from(c)
                .unwrap_or_else(|e| panic!("Can't convert {} to i64", e.into_original()));

            let seq = Sequence::new(k, c, d);

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
                .find(|p| match is_substring_of_simple(p, simple) {
                    search::SubstringResult::Never => false,
                    search::SubstringResult::Eventually(n) => n_lo >= n,
                    search::SubstringResult::Yes => true,
                })
            {
                println!("{} can be eliminated, since it contains {}", simple, p);
                continue;
            }

            // Do we give up on this sequence?
            if n_lo >= cmd.n_hi {
                println!("Reached limit on n for {}", simple);
                unsolved_branches.push(search::SearchNode::Simple(simple.clone()));
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
    use std::io;

    use super::*;

    // TODO: name exactly which branches can't be resolved
    enum Status {
        /// Completely solved; all branches eliminated.
        Complete,
        /// Can get all the minimal primes, but there's some branches
        /// we can't realize are composite.
        StrayBranches { unresolved: usize },
        /// Not solved yet. Limit the number of iterations.
        Unsolved { max_weight: usize },
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
    declare_test_for_base!(test_base_8, 8, Status::StrayBranches { unresolved: 1 });
    declare_test_for_base!(test_base_9, 9, Status::StrayBranches { unresolved: 1 });
    declare_test_for_base!(test_base_10, 10, Status::Complete);
    declare_test_for_base!(test_base_11, 11, Status::Complete);
    declare_test_for_base!(test_base_12, 12, Status::Complete);
    declare_test_for_base!(test_base_13, 13, Status::Unsolved { max_weight: 7 });
    declare_test_for_base!(test_base_14, 14, Status::Complete);
    declare_test_for_base!(test_base_15, 15, Status::Complete);
    declare_test_for_base!(test_base_16, 16, Status::Unsolved { max_weight: 10 });
    declare_test_for_base!(test_base_17, 17, Status::Unsolved { max_weight: 4 });
    declare_test_for_base!(test_base_18, 18, Status::Complete);
    declare_test_for_base!(test_base_19, 19, Status::Unsolved { max_weight: 3 });
    // 20 is solvable, but requires --release to be practical
    declare_test_for_base!(test_base_20, 20, Status::Unsolved { max_weight: 50 });
    declare_test_for_base!(test_base_21, 21, Status::Unsolved { max_weight: 4 });
    declare_test_for_base!(test_base_22, 22, Status::Unsolved { max_weight: 5 });
    declare_test_for_base!(test_base_23, 23, Status::Unsolved { max_weight: 3 });
    declare_test_for_base!(test_base_24, 24, Status::StrayBranches { unresolved: 1 });
    declare_test_for_base!(test_base_25, 25, Status::Unsolved { max_weight: 3 });
    declare_test_for_base!(test_base_26, 26, Status::Unsolved { max_weight: 3 });
    declare_test_for_base!(test_base_27, 27, Status::Unsolved { max_weight: 3 });
    declare_test_for_base!(test_base_28, 28, Status::Unsolved { max_weight: 3 });
    declare_test_for_base!(test_base_29, 29, Status::Unsolved { max_weight: 2 });
    // 30 is solvable, but takes too long, even with --release
    declare_test_for_base!(test_base_30, 30, Status::Unsolved { max_weight: 50 });

    fn test_for_base(base: u8, status: Status) {
        match status {
            Status::Complete => {
                // Simulate it for the full duration
                let max_weight = get_max_weight(base);
                let (primes, frontier) = calculate_full(base, max_weight);
                compare_primes(base, &primes, true);
                assert!(
                    frontier.is_empty(),
                    "Some branches were not eliminated!\n{}",
                    frontier.iter().format("\n")
                );
            }
            Status::StrayBranches { unresolved } => {
                // Simulate it for the full duration
                let max_weight = get_max_weight(base) + 1;
                let (primes, frontier) = calculate(base, max_weight);
                compare_primes(base, &primes, false);
                assert_eq!(
                    frontier.len(),
                    unresolved,
                    "Didn't get the expected number of unsolved branches: {}",
                    frontier.iter().format("\n")
                );
            }
            Status::Unsolved { max_weight } => {
                let (primes, frontier) = calculate(base, max_weight);
                compare_primes(base, &primes, false);
                assert!(
                    !frontier.is_empty(),
                    "All branches were eliminated, this test should be marked Complete!"
                );
            }
        }
    }

    fn calculate(base: u8, max_weight: usize) -> (CandidateSequences, Vec<search::SearchNode>) {
        let ctx = search_for_simple_families(base, Some(max_weight), None, false);
        (ctx.primes, ctx.frontier.iter().cloned().collect())
    }

    fn calculate_full(
        base: u8,
        max_weight: usize,
    ) -> (CandidateSequences, Vec<search::SearchNode>) {
        do_search(&SearchArgs {
            base,
            max_weight: Some(max_weight),
            max_iter: None,
            with_sieve: true,
            n_hi: 5000,
            p_max: 1_000_000,
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

    fn compare_primes(base: u8, primes: &CandidateSequences, require_all: bool) {
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
                // supposed to see. this might be okay.
                std::cmp::Ordering::Less => {
                    let p = truth_iter.next().unwrap();
                    if require_all {
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
