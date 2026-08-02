use num_prime::buffer::NaiveBuffer;

use crate::candidates::CandidateSequences;
use crate::search::SearchContext;

/// Context needed for sieving. This makes the most sense if it takes place
/// after searching, where we have a list of minimal primes discovered so far,
/// but I think it might make sense without it too.
pub struct SieveContext {
    pub base: u8,
    pub primes: CandidateSequences,
    pub prime_buffer: NaiveBuffer,
}

impl From<SearchContext> for SieveContext {
    fn from(ctx: SearchContext) -> Self {
        Self {
            base: ctx.base,
            primes: ctx.primes,
            prime_buffer: ctx.prime_buffer,
        }
    }
}
