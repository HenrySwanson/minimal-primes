Analysis of mepn/minimal.c: https://github.com/curtisbright/mepn/blob/master/minimal.c

has_divisor(family f) -> bool
-----------------------------
- do a bunch of gcds and check which ones are >1
- use i to mean "x_1 ... x_i y ... x_n" for any y in L_i
- one divisor
  - g = gcd(∅, i) for all i
  - everything in the family is divisible by g
  - this is just Corollary 4
- two divisors
  - fix some m
  - g1 = gcd(∅, i, mm) for all i≠m
  - g2 = gcd(m, im, mmm) for all i≠m
  - this is similar to Corollary 6
- two divisors (different approach)
  - focuses on two (nonempty) cores, i and j (i < j)
  - g1 = gcd(∅, ii)
  - g2 = gcd(i, ikk) for all k (including i)
  - g3 = gcd(j, jkk) for all k (temp2)
  - g4 = gcd(ij, iiij, ijjj) (temp3)
  - idk what this is
- 3/4/5/... divisors
  - only works with one nonempty core, so our pattern is xLz
  - example with 3: check gcd(xz, xL^3z), gcd(xLz, xL^4z), and gcd(xL^2z, xL^5z)
- TODO: check sq and cube divisors
- TODO: something with residues?

split(family f, char insplit, out list unsolved) -> bool
--------------------------------------------------------
- if insplit=0, add to list (split=0) and return false
- for i in 2...5
  - check if x_1 ... x_i y y ... x_n contains a prime
  - if so, split the pattern, add the children to the list (split=2), and return true
- TODO: there's something about splitting and gcds?
- if all that fails, add to list with split=1, return false

split2
------
???

examine(family f) -> bool
-------------------------
should we keep this family? true=keep, false=discard
- if contract(f) contains a prime, return false
- if it is prime, add to list and return false
- reduce the cores; if they're all empty, discard
- special case: simplify y*y^ny* to y*y^n (i think)
- return !has_divisor(f)

explore(family f, bool side, int pos, out list unsolved)
--------------------------------------------------------
explores the pos-th (nonempty) core in f
- pos is interpreted mod number of cores in f
- for each y in L_i
  - produces child patterns where x_i L_i is replaced with x_i y L_i or x_i L_i y, depending on side
  - 1 is split left, 0 is split right (not that it matters)
- those are all added to list (split=1)
- so is the child with L_i deleted entirely (split=1)  

main
----
- loop until only simple families left
  - repeat
    - split() everything in the list
    - keep everything that either a) couldn't split or b) examine returns true
    - same thing with split2(), with the same retention
    - repeat until no more splitting happens
  - explore everything we have left


open questions
--------------
- what's split for? i think it tracks if a family was splittable, but then why 1 vs 2?
