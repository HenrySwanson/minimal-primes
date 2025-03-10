A [*minimal prime*](https://en.wikipedia.org/wiki/Minimal_prime_(recreational_mathematics)) (in base B) is a prime that, when written out in base B, doesn't contain a substring that is also prime.

For example, in base 10, 409 is a minimal prime, but 353 is not, because it contains 53 as a substring. The substring does not need to be contiguous, so, for example, 127 is not a minimal prime (contains 17).

In base 10, there are exactly 26 minimal primes: 2, 3, 5, 7, 11, 19, 41, 61, 89, 409, 449, 499, 881, 991, 6469, 6949, 9001, 9049, 9649, 9949, 60649, 666649, 946669, 60000049, 66000049, and 66600049.

As a consequence of [Higman's Lemma](https://en.wikipedia.org/wiki/Higman%27s_lemma), the set of minimal primes in any base is always finite. Can we find all of them?

---

I got nerdsniped by this problem from a coworker, and did base 10 by hand. Then I got further into the rabbit hole and and am trying to see how far I can get by computer. Turns out for some bases it's still an open problem.

Some of these methods I found myself, but others I'm taking from a paper I found by Curtis Bright: https://cs.uwaterloo.ca/~cbright/reports/mepn.pdf

He's also got a GitHub repository here, which has been invaluable for testing my output: https://github.com/curtisbright/mepn-data/tree/master/data 
