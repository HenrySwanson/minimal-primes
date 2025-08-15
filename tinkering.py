import dataclasses
from typing import Tuple


def factor(n):
    factors = []
    while n % 2 == 0:
        factors.append(2)
        n //= 2
    for k in range(3, n, 2):
        if k * k > n:
            break
        while n % k == 0:
            factors.append(k)
            n //= k
    if n > 1:
        factors.append(n)
    return factors


@dataclasses.dataclass
class SimpleFamily:
    x: str
    y: str
    z: str


def parse_simple_family(s) -> SimpleFamily:
    before, after = s.split("*")
    return SimpleFamily(x=before[:-1], y=before[-1], z=after)


def to_sequence(f: SimpleFamily, base: int) -> Tuple[int, int, int]:
    x = int("0" + f.x, base=base)
    y = int("0" + f.y, base=base)
    z = int("0" + f.z, base=base)
    d = base - 1
    k = (x * d + y) * base ** len(f.z)
    c = d * z - y * base ** len(f.z)
    return k, c, d


def factor_simple_family(f, n, base):
    for i in range(n):
        s = f.x + f.y * i + f.z
        value = int(s, base=base)
        print(s, value, factor(value))

def foo(s, base, n=10):
    f = parse_simple_family(s)
    print(to_sequence(f, base))
    factor_simple_family(f, n, base)