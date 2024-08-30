use clap::Parser;
use itertools::Itertools;
use search::search_for_simple_families;
use std::sync::atomic::AtomicBool;

mod composite;
mod data;
mod math;
mod search;

static LOGGING_ENABLED: AtomicBool = AtomicBool::new(false);

const DEFAULT_MAX_WEIGHT: usize = 10;

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
    /// Base, e.g., decimal, binary, etc.
    #[arg(default_value_t = 10)]
    base: u8,

    /// Stop exploring when families get above this weight.
    /// If neither --max-weight or --max-iter is specified,
    /// --max-weight defaults to 10.
    #[arg(long)]
    max_weight: Option<usize>,

    /// Stop exploring after a specific number of iterations.
    /// If neither --max-weight or --max-iter is specified,
    /// --max-weight defaults to 10.
    #[arg(long)]
    max_iter: Option<usize>,

    /// Enable logging
    #[arg(long)]
    log: bool,
}

fn main() {
    let mut args = Args::parse();
    if args.max_iter.is_none() && args.max_weight.is_none() {
        args.max_weight = Some(DEFAULT_MAX_WEIGHT);
    }

    LOGGING_ENABLED.store(args.log, std::sync::atomic::Ordering::Relaxed);

    let ctx = search_for_simple_families(args.base, args.max_weight, args.max_iter);

    println!("---- BRANCHES REMAINING ----");
    for f in ctx.frontier.iter() {
        println!("{}", f);
    }
    println!("---- MINIMAL PRIMES ----");
    println!("{}", ctx.primes.iter().map(|(seq, _)| seq).format(", "));
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
        search_for_simple_families(base, Some(max_weight), None)
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
        let mut iter = ctx.primes.iter().map(|(seq, _)| seq.to_string()).peekable();

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
