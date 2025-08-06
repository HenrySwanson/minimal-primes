use num_integer::Integer;
use num_traits::One;

/// TODO: the Integer trait expects these by reference, but we
/// usually want to discard the input! should we reimplement GCD
/// but taking by value / mut ref?
pub fn nontrivial_gcd<T: Integer + One>(a: &T, b: &T) -> Option<T> {
    let g = a.gcd(b);
    if g.is_one() {
        None
    } else {
        Some(g)
    }
}
