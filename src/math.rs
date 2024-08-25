use num_bigint::BigUint;

pub fn big_one() -> BigUint {
    BigUint::from(1_u8)
}

pub fn gcd(mut a: BigUint, mut b: BigUint) -> BigUint {
    // Since we're working with a binary representation, it's more efficient
    // to use this algorithm than the standard Euclidean one.

    // First, pull out all the factors of 2.
    // gcd(2^i a, 2^j b) = 2^k gcd(a, b) when a, b odd and k = min(i, j);
    // But since x.trailing_zeros() is None when x is 0, take this opportunity
    // to check for the edge case of a = 0 or b = 0.
    let i = match a.trailing_zeros() {
        Some(i) => i,
        None => return b,
    };
    let j = match b.trailing_zeros() {
        Some(j) => j,
        None => return a,
    };

    // How many of those 2s do we want to keep?
    let k = i.min(j);
    a >>= i;
    b >>= j;

    loop {
        // Now they're both odd
        debug_assert!(a.bit(0));
        debug_assert!(b.bit(0));

        // Swap so a is larger
        if a < b {
            std::mem::swap(&mut a, &mut b);
        }

        // Subtract, just like with Euclid
        // gcd(u, v) = gcd(u, v-u)
        a -= &b;

        // Now a is even; remove its 2s (again checking for the case of a = 0).
        // gcd(2^i a, b) = gcd(a, b) when b is odd
        match a.trailing_zeros() {
            Some(i) => {
                a >>= i;
            }
            None => {
                // gcd(0, b) = b, then add in the 2^k we removed earlier
                return b << k;
            }
        }
    }
}

pub fn gcd_reduce(numbers: impl IntoIterator<Item = BigUint>) -> BigUint {
    gcd_reduce_with(BigUint::ZERO, numbers)
}

pub fn gcd_reduce_with(init: BigUint, numbers: impl IntoIterator<Item = BigUint>) -> BigUint {
    let mut g = init;
    for n in numbers {
        g = gcd(g, n);
        if g == big_one() {
            // g will always be 1; break now and save some time
            break;
        }
    }
    g
}
