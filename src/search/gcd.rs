use num_integer::Integer;
use num_traits::One;

/// Returns the GCD of the two inputs, unless it's 1, in which case
/// it returns None.
pub fn nontrivial_gcd<T: Integer + One>(a: &T, b: &T) -> Option<T> {
    // TODO: the Integer trait expects these by reference, but we
    // usually want to discard the input! should we reimplement GCD
    // but taking by value / mut ref?
    let g = a.gcd(b);
    if g.is_one() {
        None
    } else {
        Some(g)
    }
}
