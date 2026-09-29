# Separating opposite pairs without losing the centering

**The full Paley conjecture and the Proximity Prize remain unproved.**
This argument separates the contribution of opposite pairs from the
distinct-coordinate moment in pass33. Both the forward transform and
its inverse have controlled Gaussian coefficient norms. Consequently
the missing moment hierarchy can be posed as a one-sided estimate on
centered, distinct, opposite-free relations. That estimate is not proved
here. A separate lower bound shows why its principal term cannot be
omitted, already at degree ten.

These are ordinary uniform arguments supported by exact finite checks.
They have not received separate-author review or Lean verification, and
no literature novelty is claimed.

## 1. Counts and normalization

Let p be an odd prime and A=-A be a subset of F_p^* of even size
n=2N. The combinatorial transforms below do not require a multiplicative
subgroup. Define

    I_r = number of ordered, distinct zero-sum r-tuples in A,
    O_r = number of ordered, distinct zero-sum r-tuples in A
          with no two entries opposite,
    a_r = 2^r (N)_r.

Use I_0=O_0=1. There are (n)_r ordered distinct tuples in all, and
a_r ordered distinct opposite-free tuples in all. The latter formula
chooses r of the N opposite classes, one sign in each, and an ordering.

For even orders put

    Q_s = I_(2s) - (n)_(2s)/p,
    B_s = O_(2s) - 2^(2s)(N)_(2s)/p,
    G_s = (2s-1)!! n^s,      G_0=1.

In particular Q_0=B_0=1-1/p, O_2=0, and

    B_1 = -n(n-2)/p.                                      (1)

Do not substitute n^(2s)/p for the principal term of either restricted
count. At six terms O_6 includes both product-balanced and unbalanced
words; it is not the previously isolated remainder D6.

## 2. An exact positive transform

Assume 2s<=N. A distinct word has a unique set of full opposite pairs.
After removing its j pairs, its remaining k=2s-2j entries are distinct
and opposite-free. A fixed ordered remainder occupies k opposite
classes. There are binomial(N-k,j) choices of additional classes whose
full pairs will be inserted. There are (2s)!/k! ways to place and order
the added entries while keeping the order of the remainder. Thus

    I_(2s) = sum_(t=0)^s C_(s,t) O_(2t),
    C_(s,t) = (2s)!/(2t)! binomial(N-2t,s-t).              (2)

The same decomposition without the zero-sum requirement gives

    (n)_(2s) = sum_(t=0)^s C_(s,t) 2^(2t)(N)_(2t).

Subtracting its p-th part from (2) proves the centered identity

    Q_s = sum_(t=0)^s C_(s,t) B_t.                        (3)

Every C_(s,t) is nonnegative. The restriction 2s<=N is sufficient for
all displayed binomial coefficients and for the inverse below; it is
not claimed to be the maximal domain of the forward counting identity.

For example,

    I_2 = O_2 + 2N O_0,
    I_4 = O_4 + 12(N-2)O_2 + 12N(N-1)O_0,
    I_6 = O_6 + 30(N-4)O_4
          + 360 binomial(N-2,2)O_2 + 720 binomial(N,3)O_0.

The all-pair term alone has size comparable to G_s at depths
s^2/n tending to zero. Opposite pairs are therefore not a negligible
Gaussian-scale error.

## 3. Generating functions and the inverse transform

Work first in the rational group algebra of the additive group F_p,
with basis X^a and multiplication X^a X^b=X^(a+b). For each pair
{a,-a}, the generating factor for arbitrary distinct selections is

    1 + z(X^a+X^(-a)) + z^2.

The corresponding opposite-free factor is 1+w(X^a+X^(-a)). Taking
the coefficient of X^0 after multiplication over the N pairs gives

    sum_r I_r z^r/r! = (1+z^2)^N sum_r O_r w^r/r!,
    w = z/(1+z^2).                                       (4)

The same substitution relates the total-count series (1+z)^n and
(1+2w)^N. Therefore (4) remains true for the centered series

    Q(z) = sum_r [I_r-(n)_r/p] z^r/r!,
    B(w) = sum_r [O_r-2^r(N)_r/p] w^r/r!.

Since w=z+O(z^3), it has a unique inverse z=z(w) over rational formal
series, and B(w)=(1+z(w)^2)^(-N)Q(z(w)). Formal residue substitution,
using dw/dz=(1-z^2)/(1+z^2)^2, yields

    [w^r] B(w) = [z^r] Q(z)(1-z^2)(1+z^2)^(r-N-1).       (5)

This coefficient calculation is purely formal and uses no analytic
convergence assumption. For r=2s, k=2t, j=s-t, and r<=N, the
coefficient multiplying [z^k]Q is

    binomial(r-N-1,j) - binomial(r-N-1,j-1)
      = (-1)^j [binomial(N-r+j,j)+binomial(N-r+j-1,j-1)].

The equality follows by expanding a negative integral power with
binomial(-d,j)=(-1)^j binomial(d+j-1,j). With M=N-2t, define

    c(M,0)=1,
    c(M,j)=binomial(M-j,j)+binomial(M-j-1,j-1), j>=1.

Here M>=2j whenever j>=1. Restoring factorials in (5) proves

    B_s = sum_(t=0)^s (-1)^(s-t) (2s)!/(2t)!
                      c(N-2t,s-t) Q_t.                  (6)

For M>=3, c(M,j) counts j-element independent vertex sets of a cycle
of length M. Excluding a specified vertex gives the first binomial
coefficient; including it and excluding its neighbors gives the
second. The exceptional M=2,j=1 formula is directly c(2,1)=2.
In either case

    c(M,j) <= binomial(M,j) <= N^j/j!,     M<=N.           (7)

For example, the inverse at degree four is

    O_4 = I_4 - 12(N-2)I_2 + 12N(N-3)I_0.

The finite verifier checks the inverse matrices directly, as well as
counting the cycle independent sets by exhaustive binary masks. Those
checks support, but do not replace, the formal derivation above.

## 4. The Gaussian hierarchy is preserved up to constants

For j=s-t, both forward and absolute inverse coefficients are at most

    (2s)!/(2t)! * (n/2)^j/j!
       = binomial(s,t) G_s/G_t.                          (8)

Consequently if K>=1 and the one-sided bounds

    B_t <= K^t G_t,          0<=t<=s,

hold, positivity of (3) and the binomial theorem imply

    Q_s <= (1+K)^s G_s.                                  (9)

There is no need to assume B_t>=0 or a lower bound on B_t in this
direction. The t=0 hypothesis holds automatically because B_0<=1.

For s>=2, use the [pass33 estimate](parallel33-centered-distinct-moments-2026-09-05.md).
Write

    eta(a) = sum_(x in A) exp(2 pi i a x/p),
    T_t = (1/p) sum_(a != 0) |eta(a)|^(2t),
    alpha = n(p-n)/(p-1),
    epsilon_(2t)(alpha)=product_(j=1)^(2t-1)(1+j/sqrt(alpha))-1.

For 2s<p and epsilon_(2s)<=1/2, pass33 gives T_s<=2Q_s. Combining
this with (9) gives

    T_s <= 2(1+K)^s G_s.                                 (10)

Conversely, suppose T_t<=K^t G_t for 1<=t<=s and the same epsilon
condition holds. Since epsilon increases with t, pass33 gives
0<Q_t<=2T_t for 2<=t<=s. At t=1 one has directly

    Q_1 = n-n(n-1)/p = T_1+n/p <= 2T_1,

since p-n>=1; Q_0<=1. Taking absolute values in (6) and using (8)
therefore gives

    |B_s| <= 2(1+K)^s G_s.                               (11)

Thus Gaussian bounds for the full centered moment hierarchy and for
the opposite-free centered hierarchy are equivalent up to absolute
changes in constants, on these parameter ranges. An estimate at just
one degree is not being promoted to all smaller degrees.

In the quartic window n^4/4<=p<=n^4, the conditions 2s<=N, 2s<p,
and epsilon_(2s)<=1/2 hold eventually when s=O(log n), by pass33.
The finite constants in the epsilon bound can be large. If A=H is a
multiplicative subgroup, coset constancy gives

    M^(2s) <= (p/n) T_s,
    M = max_(a != 0) |eta(a)|.

Since G_s<=(2sn)^s, (10) would imply

    M <= (2p/n)^(1/(2s)) sqrt(2(1+K)sn).

At s>=log(2p/n) this is at most sqrt(2e(1+K)sn). The missing input
is the uniform upper estimate on B_t. The transformation does not
provide it, improve the currently recorded period exponent 71/72,
or prove a reduction from this subgroup target to every part of the
full classical Paley conjecture or the official prize statement.

## 5. Distinct, opposite-free ten-term relations are forced

This section works for every symmetric A subset F_p^*, without the
subgroup assumption or the earlier condition 2s<=N. Let

    E_t = (1/p) sum_a |eta(a)|^(2t),  t>=1.

At integer s, E_s counts all zero-sum 2s-words. Parseval gives E_1=n.
For a fixed positional pair required to be equal, Fourier inversion
counts its words as

    (1/p) sum_a eta(a)^(2s-2) eta(2a)
      <= (1/p) sum_a |eta(a)|^(2s-2)|eta(2a)|
      <= E_(s-1/2).

The last inequality is Holder and the fact that a->2a permutes F_p.
This is the repeated-word argument of pass33 with its unnecessary
subgroup restriction removed: the Fourier formula counts the sum over
all possible common values directly, without scaling to a fixed value.

Interpolating between E_1 and E_s and using Jensen E_s>=n^s shows,
for s>=2, that the number R of words with any repeated coordinate obeys

    R/E_s <= binomial(2s,2)/sqrt(n).                       (12)

For a fixed positional opposite pair there are exactly n E_(s-1)
words: choose its first entry and then a zero-sum remaining word.
Normalized Lyapunov on all p frequencies gives

    E_(s-1) <= E_s^((s-1)/s).

The principal frequency also gives E_s>=n^(2s)/p. Consequently the
number P of words having some opposite pair satisfies

    P/E_s <= binomial(2s,2) n E_s^(-1/s)
           <= binomial(2s,2) p^(1/s)/n.                  (13)

The complement of the union of these two bad classes is precisely
the distinct, opposite-free count. Combining (12) and (13) proves

    O_(2s) >= [1-binomial(2s,2)(n^(-1/2)+p^(1/s)/n)] E_s. (14)

It is valid even when the bracket is negative, but only useful when
positive. For fixed s>4 and p<=n^4, the bracket tends to one.
In particular at s=5,

    O_10 >= [1-45(n^(-1/2)+n^(-1/5))] n^6               (15)

when the bracket is nonnegative and p<=n^4. A simple explicit
sufficient threshold is n>=2^35: then n^(-1/5)<=1/128 and
n^(-1/2)<=1/2^17, so

    45(1/128+1/131072) < 1/2,
    O_10 >= n^6/2.                                       (16)

This conclusion applies in particular to every eligible symmetric
subgroup. For fixed C, the raw upper bound O_10<=C n^5 cannot hold
in an eligible case with n>=2^35 and n>2C. No existence of a quartic
prime/subgroup pair at every such order, or infinite sequence of such
pairs, is asserted. This is not a counterexample to a centered Gaussian
bound, the Paley conjecture, or the prize. The term subtracted in B_5
is 2^10(N)_10/p, which itself has order n^10/p for large n; its size
explains why the raw count is the wrong Gaussian target.

## 6. Verification scope and next missing input

The [exact verifier](../experiments/parallel34_verify_2026_09_05.py)
computes I and O using different dynamic programs: one processes
individual elements; the other processes opposite classes with three
choices. A third convolution counts unrestricted words. It checks raw
and centered transforms, total cardinalities, both union bounds,
Lyapunov and Jensen inequalities, the inverse matrices, Gaussian
coefficient bounds, and the explicit degree-ten threshold.

The [results](../results/parallel34_verification_2026_09_05.json)
include eight multiplicative subgroups, a symmetric non-subgroup,
orders through twelve, and an actual quartic case p=33713,n=16.
Binary-mask cycle counts and coefficient-matrix checks use independent
finite constructions. All comparisons are rational or integral.

What remains is a uniform quantitative upper estimate after the
principal subtraction, for example B_t<=K^t G_t uniformly through
t=O(log n). A bound on the raw number of short relations does not
address the higher-degree principal term. The transform identifies
the cancellation that must be proved; it does not establish it.
