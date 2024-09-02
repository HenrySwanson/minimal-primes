use clap::Parser;
use data_structures::CandidateSequences;
use itertools::Itertools;
use num_bigint::{BigInt, BigUint};
use search::{is_substring_of_simple, search_for_simple_families};
use sequences::{Digit, DigitSeq};
use std::{sync::atomic::AtomicBool, usize};

mod composite;
mod data_structures;
mod math;
mod search;
mod sequences;
mod sieve;

static LOGGING_ENABLED: AtomicBool = AtomicBool::new(false);

macro_rules! debug_println {
    ($($arg:tt)*) => {
        if crate::LOGGING_ENABLED.load(std::sync::atomic::Ordering::Relaxed) {
            println!($($arg)*);
        }
    };
}

pub(crate) use debug_println;

#[derive(Parser)]
struct Args {
    /// What to do
    #[command(subcommand)]
    command: Command,

    /// Enable logging
    #[arg(long)]
    log: bool,
}

#[derive(clap::Subcommand)]
enum Command {
    /// Searches for minimal primes in the given base.
    Search(SearchArgs),
    /// Sieves through a sequence of the form k b^n + c.
    Sieve(SieveArgs),
    /// Do both search and sieve!
    Full(FullArgs),
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

    /// Stop exploring after all remaining families are simple.
    #[arg(long)]
    stop_when_simple: bool,
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
struct FullArgs {
    #[clap(flatten)]
    search_args: SearchArgs,

    /// upper bound for n
    #[arg(default_value_t = 5_000)]
    n_hi: usize,

    /// max p to sieve with
    #[arg(default_value_t = 1_000_000)]
    p_max: u64,
}

fn main() {
    let args = Args::parse();
    LOGGING_ENABLED.store(args.log, std::sync::atomic::Ordering::Relaxed);

    match args.command {
        Command::Search(cmd) => {
            do_search(&cmd);
        }
        Command::Sieve(cmd) => {
            do_sieve(&cmd);
        }
        Command::Full(cmd) => {
            do_full(&cmd);
        }
    }
}

fn do_full(cmd: &FullArgs) -> (CandidateSequences, Vec<search::SearchNode>) {
    let mut ctx = do_search(&cmd.search_args);

    if !ctx.frontier.all_simple() {
        println!("Not all remaining branches are simple! Must bail out now.");
        return (ctx.primes, ctx.frontier.iter().cloned().collect());
    }

    println!("---- INTERMEDIARY PHASE ----");
    // It's possible that a simple family can only be expanded a finite amount
    // before it conflicts with a known minimal prime. If so, we should not
    // jump right to sieving, but should instead continue normal searching.

    loop {
        // Check if any family is potentially able to contain any prime.
        let should_loop = ctx.frontier.iter().any(|family| {
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

        if should_loop {
            println!(
                "Searching one extra round: iter={}, num primes={}, num_branches={}",
                ctx.iter,
                ctx.primes.len(),
                ctx.frontier.len()
            );
            ctx.search_one_level();
        } else {
            break;
        }
    }

    println!("---- SIEVING PHASE ----");
    let mut primes = ctx.primes;
    let mut leftover_branches = vec![];

    // TODO: sieve them all at the same time!
    for family in ctx.frontier.iter().cloned() {
        let simple = match &family {
            search::SearchNode::Simple(x) => x,
            _ => unreachable!("found non-simple family after all_simple()"),
        };

        // Compute the sequence for this family: xy*z
        let base = cmd.search_args.base;
        let x = simple.before.value(base);
        let y = simple.center.0;
        let z = simple.after.value(base);

        let b_z = BigUint::from(base).pow(simple.after.0.len() as u32);
        let d = u64::from(base) - 1;
        let k = (x * d + y) * &b_z;
        let c = BigInt::from(d * z) - BigInt::from(y * b_z);

        // Try to fit it into the appropriate ranges
        let k = match u64::try_from(k) {
            Ok(k) => k,
            Err(e) => {
                println!("Can't convert {} to u64", e.into_original());
                leftover_branches.push(family);
                continue;
            }
        };
        let c = match i64::try_from(c) {
            Ok(c) => c,
            Err(e) => {
                println!("Can't convert {} to i64", e.into_original());
                leftover_branches.push(family);
                continue;
            }
        };

        // Great, we have a sequence! Sieve it and see what we can get :)
        let n_lo = simple.num_repeats;
        let n_hi = cmd.n_hi;
        println!(
            "Investigating family {} -> ({}*{}^n+{})/{} for n from {} to {}",
            simple, k, base, c, d, n_lo, n_hi,
        );
        match sieve::find_first_prime(base, k, c, d, n_lo, n_hi, cmd.p_max) {
            Some((i, p)) => {
                let digitseq =
                    DigitSeq(p.to_radix_be(base.into()).into_iter().map(Digit).collect());
                println!("Found prime at exponent {}: {}", i, digitseq);
                primes.insert(digitseq);
            }
            None => {
                println!("Unable to find prime in the given range: {}", simple);
                leftover_branches.push(family);
            }
        }
    }

    primes.sort();
    println!(
        "Final set of primes ({}): {}",
        primes.len(),
        primes.iter().format(", ")
    );
    (primes, leftover_branches)
}

fn do_search(cmd: &SearchArgs) -> search::SearchContext {
    let ctx =
        search_for_simple_families(cmd.base, cmd.max_weight, cmd.max_iter, cmd.stop_when_simple);

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

fn do_sieve(cmd: &SieveArgs) {
    let x =
        sieve::find_first_prime(cmd.base, cmd.k, cmd.c, 1, cmd.n_lo, cmd.n_hi, cmd.p_max).unwrap();
    println!("{}, {}", x.0, x.1);
}

#[cfg(test)]
mod tests {
    use std::io;

    use search::SearchContext;

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
                let final_ctx = calculate(base, max_weight);
                compare_primes(&final_ctx, true);
                assert!(
                    final_ctx.frontier.is_empty(),
                    "Some branches were not eliminated!\n{}",
                    final_ctx.frontier.iter().format("\n")
                );
            }
            Status::StrayBranches { unresolved } => {
                // Simulate it for the full duration
                let max_weight = get_max_weight(base) + 1;
                let final_ctx = calculate(base, max_weight);
                compare_primes(&final_ctx, false);
                assert_eq!(
                    final_ctx.frontier.len(),
                    unresolved,
                    "Didn't get the expected number of unsolved branches: {}",
                    final_ctx.frontier.iter().format("\n")
                );
            }
            Status::Unsolved { max_weight } => {
                let final_ctx = calculate(base, max_weight);
                compare_primes(&final_ctx, false);
                assert!(
                    !final_ctx.frontier.is_empty(),
                    "All branches were eliminated, this test should be marked Complete!"
                );
            }
        }
    }

    fn calculate(base: u8, max_weight: usize) -> SearchContext {
        search_for_simple_families(base, Some(max_weight), None, false)
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

    fn compare_primes(ctx: &SearchContext, require_all: bool) {
        let mut truth_iter = iter_ground_truth(ctx.base).peekable();
        let mut iter = ctx.primes.iter().map(|seq| seq.to_string()).peekable();

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
