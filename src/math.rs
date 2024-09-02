use num_integer::Integer;

pub fn gcd_reduce<T: Integer>(numbers: impl IntoIterator<Item = T>) -> T {
    gcd_reduce_with(T::zero(), numbers)
}

pub fn gcd_reduce_with<T: Integer>(init: T, numbers: impl IntoIterator<Item = T>) -> T {
    let mut g = init;
    for n in numbers {
        g = g.gcd(&n);
        if g == T::one() {
            // g will always be 1; break now and save some time
            break;
        }
    }
    g
}
