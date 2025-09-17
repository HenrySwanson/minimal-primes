use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::Zero;

use crate::families::SimpleFamily;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sequence {
    pub k: u64,
    pub c: i64,
    pub d: u64,
}

// TODO: merge with sequence somehow?
pub struct BigSequence {
    pub k: BigUint,
    pub c: BigInt,
    pub d: u8,
}

impl Sequence {
    pub fn new(k: u64, c: i64, d: u64) -> Self {
        assert_ne!(k, 0);
        assert_ne!(c, 0);
        assert_ne!(d, 0);
        // k and c have to be opposites mod d
        assert_eq!(k.checked_add_signed(c).unwrap() % d, 0);

        // Do some quick reduction to put it in lowest terms
        let gcd = k.gcd(&c.unsigned_abs()).gcd(&d);

        Self {
            k: k / gcd,
            // casting is okay because 0 < gcd <= |c|
            c: c / (gcd as i64),
            d: d / gcd,
        }
    }

    // TODO: better error type
    pub fn try_from_family(simple: &SimpleFamily, base: u8) -> Result<Self, String> {
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
            .map_err(|e| format!("Can't convert {} to u64 for {}", e.into_original(), simple))?;
        let c = i64::try_from(c)
            .map_err(|e| format!("Can't convert {} to i64 for {}", e.into_original(), simple))?;

        Ok(Sequence::new(k, c, d))
    }

    pub fn compute_term(&self, n: u32, base: u64) -> BigUint {
        let bn = BigUint::from(base).pow(n);
        let kbnc = if self.c > 0 {
            self.k * bn + self.c.unsigned_abs()
        } else {
            self.k * bn - self.c.unsigned_abs()
        };
        let (q, r) = kbnc.div_rem(&self.d.into());
        debug_assert_eq!(r, BigUint::ZERO);
        q
    }

    pub fn check_term_equal(&self, base: u64, p: u64, n: usize) -> bool {
        let mut x = u128::from(p);
        x *= u128::from(self.d);
        if self.c > 0 {
            let c = u128::from(self.c.unsigned_abs());

            x = match x.checked_sub(c) {
                Some(x) => x,
                None => return false,
            };
        } else {
            x += u128::from(self.c.unsigned_abs());
        }

        if x % u128::from(self.k) != 0 {
            return false;
        }
        x /= u128::from(self.k);

        for _ in 0..n {
            if x % u128::from(base) != 0 {
                return false;
            }
            x /= u128::from(base);
        }

        x == 1
    }
}

impl BigSequence {
    pub fn new(k: BigUint, c: BigInt, d: u8) -> Self {
        assert_ne!(k, BigUint::ZERO);
        assert_ne!(c, BigInt::ZERO);
        assert_ne!(d, 0);
        // k and c have to be opposites mod d
        assert!(((BigInt::from(k.clone()) + &c) % d).is_zero());

        // Do some quick reduction to put it in lowest terms
        let gcd = BigUint::from(d).gcd(&k).gcd(c.magnitude());
        let gcd: u8 = gcd.try_into().expect("gcd is <= d");

        Self {
            k: k / gcd,
            // casting is okay because 0 < gcd <= |c|
            c: c / gcd,
            d: d / gcd,
        }
    }

    pub fn from_family(simple: &SimpleFamily, base: u8) -> Self {
        // Compute the sequence for this family: xy*z
        let x = simple.before.value(base);
        let y = simple.center.0;
        let z = simple.after.value(base);

        let b_z = BigUint::from(base).pow(simple.after.0.len() as u32);
        let d = base - 1;
        let k = (x * d + y) * &b_z;
        let c = BigInt::from(d * z) - BigInt::from(y * b_z);

        BigSequence::new(k, c, d)
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}*b^n+{})/{}", self.k, self.c, self.d)
    }
}

impl std::fmt::Display for BigSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}*b^n+{})/{}", self.k, self.c, self.d)
    }
}
