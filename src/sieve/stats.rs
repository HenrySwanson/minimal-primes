use std::ops::Range;
use std::time::Duration;

/// Stats for a single round of sieving (one n_range). A fresh one is created
/// (and thus implicitly reset) at the start of every [crate::sieve::do_one_round]
/// call, and printed at the end of it, so these numbers are always specific
/// to one round rather than cumulative across the whole run.
#[derive(Debug, Default)]
pub struct SieveStats {
    /// number of sequences being sieved this round
    pub num_slices: usize,

    /// number of primes iterated over while running BSGS this round
    pub num_bsgs_primes: usize,
    /// total time spent inside baby_step_giant_step this round
    pub duration_bsgs: Duration,

    /// number of (slice, n) candidates the BSGS sieve eliminated this round
    pub num_eliminated_by_sieve: usize,
    /// number of (slice, n) candidates that survived sieving and had to be
    /// checked individually
    pub num_candidates_after_sieve: usize,

    /// number of individual candidates actually run through a primality test
    /// this round (<= num_candidates_after_sieve, since we stop checking a
    /// slice as soon as we find a prime in it)
    pub num_primality_tests: usize,
    /// total time spent on those primality tests
    pub duration_primality_tests: Duration,

    /// number of new primes found this round
    pub num_primes_found: usize,
}

impl SieveStats {
    pub fn record_bsgs(&mut self, count: usize, elapsed: Duration) {
        self.num_bsgs_primes += count;
        self.duration_bsgs += elapsed;
    }

    pub fn record_primality_test(&mut self, elapsed: Duration) {
        self.num_primality_tests += 1;
        self.duration_primality_tests += elapsed;
    }

    pub fn print(&self, n_range: &Range<usize>, p_max: u64) {
        println!(
            "---- SIEVE STATS (n = {}..{}, p_max = {}) ----",
            n_range.start, n_range.end, p_max
        );
        println!(
            "{} sequences sieved over {} values of n",
            self.num_slices,
            n_range.len()
        );
        println!(
            "{} primes used in BSGS ({}ms total, {}/prime, {}/prime/seq)",
            self.num_bsgs_primes,
            self.duration_bsgs.as_millis(),
            format_avg(self.duration_bsgs, self.num_bsgs_primes),
            format_avg(self.duration_bsgs, self.num_bsgs_primes * self.num_slices),
        );
        println!(
            "{} candidates eliminated by sieving, {} left over for primality testing",
            self.num_eliminated_by_sieve, self.num_candidates_after_sieve,
        );
        println!(
            "{} primality tests performed ({}ms total, {}/test)",
            self.num_primality_tests,
            self.duration_primality_tests.as_millis(),
            format_avg(self.duration_primality_tests, self.num_primality_tests),
        );
        println!("{} new primes found", self.num_primes_found);
    }
}

/// Formats `d / n` as a human-readable per-unit duration, or "n/a" if n is 0.
fn format_avg(d: Duration, n: usize) -> String {
    if n == 0 {
        return "n/a".to_string();
    }
    let micros_per = d.as_micros() as f64 / n as f64;
    if micros_per >= 1000.0 {
        format!("{:.3}ms", micros_per / 1000.0)
    } else {
        format!("{:.3}us", micros_per)
    }
}

/// e^-gamma, from Mertens' third theorem: the fraction of integers surviving
/// sieving by every prime up to P is asymptotically e^-gamma / ln(P).
/// e^-γ, where γ is the Euler-Mascheroni constant. Computed with WolframAlpha.
const E_NEG_GAMMA: f64 = 0.561_459_483_566_885_1;

/// Suggests the right p_max for the next round of sieving, based on the
/// timing statistics from the previous round.
///
/// The idea is that the more we sieve, the more numbers we can eliminate,
/// which means the fewer (expensive) primality tests we have to run, but
/// the sieving itself takes time. The tricky part is that it's not as easy
/// as "one more prime in the sieve => M fewer primality tests".
///
/// Mertens' third theorem tells us the probability of a number surviving a
/// sieve of primes up to P -- it is approximately e^-γ / ln(P), where γ is
/// the Euler-Mascheroni constant.
///
/// Meanwhile, the cost of doing that sieve is the cost of one round of BSGS,
/// times the number of primes less than P (approximately P/ln(P)).
///
/// So, if N is the number of candidates, cost_bsgs is the cost-per-prime, and
/// cost_test is the cost of a primality test, the total cost is:
///
///   cost_bsgs * P / ln(P) + cost_test * N * e^-gamma / ln(P)
///
/// Taking the derivative with respect to P, and setting it equal to zero, we get:
///
///   cost_bsgs * (ln(P) - 1) / ln(P)^2 = cost_test * N * e^-gamma / P / ln(P)^2
///
///   P * (ln(P) - 1) = cost_test * N * e^-gamma / cost_bsgs
///
/// If this calculation can't be done for some reason (empty stats?), returns
/// `None` instead.
pub fn suggest_next_p_max(
    prev: &SieveStats,
    prev_range_len: usize,
    next_range_len: usize,
    next_num_slices: usize,
) -> Option<u64> {
    if prev.num_bsgs_primes == 0 || prev.num_primality_tests == 0 || prev.num_slices == 0 {
        return None;
    }

    // The cost of BSGS scales up with sqrt(range * num_slices) (see the BSGS
    // code above), so take the previous cost and scale it accordingly.
    let scale = ((next_range_len * next_num_slices) as f64
        / (prev_range_len * prev.num_slices) as f64)
        .sqrt();
    let cost_bsgs = scale * (prev.duration_bsgs.as_secs_f64() / prev.num_bsgs_primes as f64);

    // The cost of primality testing grows as the number of digits gets larger.
    // I don't have a hard-and-fast rule here, but it seems to grow about 3x
    // each time the number of digits doubles.
    let scale = 3.0;
    let cost_test =
        scale * prev.duration_primality_tests.as_secs_f64() / prev.num_primality_tests as f64;

    // Remember, we want to solve P * (ln(P) - 1) = cost_test * N * e^-gamma / cost_bsgs.
    let total_candidates = (next_range_len * next_num_slices) as f64;
    let rhs = total_candidates * E_NEG_GAMMA * cost_test / cost_bsgs;

    // We approximate a solution by doing some fixed-point iteration on it.
    let mut p = rhs;
    for _ in 0..10 {
        // If p < e, then the denominator below is negative. I don't see
        // how this could happen, but bail out and return something small.
        if p < std::f64::consts::E {
            return Some(2);
        }
        p = rhs / (p.ln() - 1.0);
    }

    Some(p.round() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suggest_next_p_max_no_history() {
        let stats = SieveStats::default();
        assert_eq!(suggest_next_p_max(&stats, 16, 32, 100), None);
    }

    #[test]
    fn test_suggest_next_p_max_solves_the_balance_condition() {
        // Made-up but self-consistent stats: 1000 primes cost 1ms total in
        // BSGS (1us/prime), and 500 primality tests cost 100ms total
        // (200us/test).
        let stats = SieveStats {
            num_slices: 10,
            num_bsgs_primes: 1000,
            duration_bsgs: Duration::from_millis(1),
            num_primality_tests: 500,
            duration_primality_tests: Duration::from_millis(100),
            ..Default::default()
        };

        let p = suggest_next_p_max(&stats, 16, 16, 10).expect("should have a suggestion");

        // Check it actually (approximately) solves p*(ln(p)-1) = rhs for the
        // rhs this scenario implies, rather than just trusting the formula
        // was transcribed correctly. Mirrors suggest_next_p_max's own math,
        // including the flat 3x fudge factor for primality-test cost growth.
        let t_bsgs = stats.duration_bsgs.as_secs_f64() / stats.num_bsgs_primes as f64;
        let t_test =
            3.0 * stats.duration_primality_tests.as_secs_f64() / stats.num_primality_tests as f64;
        let n = 16.0 * 10.0;
        let rhs = n * E_NEG_GAMMA * t_test / t_bsgs;

        let p = p as f64;
        let residual = (p * (p.ln() - 1.0) - rhs).abs() / rhs;
        assert!(
            residual < 0.01,
            "p={p} doesn't satisfy p*(ln(p)-1)={rhs} closely enough (residual {residual})"
        );
    }

    #[test]
    fn test_suggest_next_p_max_grows_with_more_candidates() {
        // Same per-unit costs, but the next round covers more ground (more
        // candidates to potentially save primality tests on) -- should want
        // to sieve deeper.
        let stats = SieveStats {
            num_slices: 10,
            num_bsgs_primes: 1000,
            duration_bsgs: Duration::from_millis(1),
            num_primality_tests: 500,
            duration_primality_tests: Duration::from_millis(100),
            ..Default::default()
        };

        let p_small = suggest_next_p_max(&stats, 16, 16, 10).unwrap();
        let p_large = suggest_next_p_max(&stats, 16, 16, 1000).unwrap();
        assert!(
            p_large > p_small,
            "expected deeper sieving for a round with more candidates: {p_large} <= {p_small}"
        );
    }
}
