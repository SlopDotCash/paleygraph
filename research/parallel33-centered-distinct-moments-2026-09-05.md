# Removing repeated coordinates from centered subgroup moments

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass gives an explicit reduction to distinct-coordinate counts at
growing moment depth, with the principal-frequency terms subtracted
exactly. It also improves the recorded upper bound on repeated six-word
relations. No upper bound on the remaining distinct-coordinate aggregate
is proved. These are ordinary arguments with exact finite checks;
separate-author review and Lean verification remain outstanding.

## 1. The centered distinct-coordinate estimate

Let p be an odd prime, A=-A a subset of F_p with 0<n=|A|<p, and
r=2s an even integer with 4<=r<p. A need not be a subgroup in this
section. Define

    eta(a) = sum_(x in A) exp(2 pi i a x/p),
    T_s = (1/p) sum_(a != 0) |eta(a)|^(2s),
    I_(2s) = #{(x_1,...,x_(2s)) in A^(2s):
                sum x_i=0 and all x_i are distinct},
    Q_s = I_(2s) - (n)_(2s)/p,
    alpha = n(p-n)/(p-1),

where (n)_r=n(n-1)...(n-r+1), with value zero when r>n.
Symmetry makes eta real. Consequently T_s is the ordinary zero-sum
2s-word count with n^(2s)/p subtracted, and T_s>0.

Put

    epsilon_r(alpha) = product_(j=1)^(r-1)(1+j/sqrt(alpha)) - 1.

Then

    |Q_s-T_s| <= epsilon_(2s)(alpha) T_s.                  (1)

In particular, when epsilon_(2s)(alpha)<1,

    (1-epsilon) T_s <= Q_s <= (1+epsilon) T_s.             (2)

The centered distinct count is positive under this last condition.
Positivity is not asserted when the error factor is at least one.

## 2. Center each collision partition before estimating it

Let pi be a set partition of the r labeled positions, with k blocks
of sizes m_1,...,m_k. Let W_pi count solutions of

    m_1 y_1 + ... + m_k y_k = 0,       y_i in A,

with no distinctness requirement between the y_i. Fourier inversion
gives

    V_pi := W_pi - n^k/p
           = (1/p) sum_(a != 0) product_(j=1)^k eta(m_j a).       (3)

Every m_j is nonzero modulo p because 1<=m_j<=r<p. Thus multiplication
by m_j permutes the nonzero frequencies. Write

    M_t = (1/(p-1)) sum_(a != 0) |eta(a)|^(2t).

Parseval gives M_1=alpha, and T_s=((p-1)/p)M_s. Holder followed by
the monotonicity of normalized L^q norms gives

    |V_pi| <= ((p-1)/p) M_(k/2)
            <= ((p-1)/p) M_s^(k/r).

Jensen also gives M_s>=alpha^s. Since k<=r, these inequalities imply

    |V_pi|/T_s <= M_s^(-(r-k)/r)
                 <= alpha^(-(r-k)/2).                     (4)

The symmetry of A is used to identify the full even zero-sum count
with the nonnegative Fourier moment. The scalar-permutation argument
requires the characteristic condition r<p. It is not legitimate to
apply (4) to a partition having a block size divisible by p.

The exact inclusion-exclusion coefficient of pi is

    mu(pi) = product_(B in pi) (-1)^(|B|-1)(|B|-1)!.

It satisfies

    I_r = sum_pi mu(pi) W_pi,
    (n)_r = sum_pi mu(pi) n^(number of blocks of pi).

One way to prove the first identity is to fix a word and sum mu(pi)
over partitions refined by its equality partition. On an equality
class of size b, the coefficient sum is the signed sum over permutations
of b letters by cycle count. It equals 1 at b=1 and 0 at b>1, by the
falling-factorial identity evaluated at 1. The word survives precisely
when all its coordinates are distinct. Applying the same reasoning
without the zero-sum restriction proves the second identity.

Subtracting the second identity divided by p from the first yields

    Q_s = sum_pi mu(pi) V_pi.                              (5)

The discrete partition contributes exactly T_s. The sum of |mu(pi)|
over partitions with k blocks is the unsigned Stirling number c(r,k):
each factor (|B|-1)! counts cyclic orderings of its block. The cycle
insertion recurrence gives

    sum_k c(r,k) x^k = x(x+1)...(x+r-1).

Apply (4) to the other partitions in (5), then use this generating
polynomial with x=sqrt(alpha). The result is exactly (1).

This proof does not treat the number of all repeated words as a small
absolute error after centering. Their own principal terms are removed
in (3) before Holder is applied. Using n^r/p for both counts would
lose the identity (5).

## 3. Uniformity at logarithmic depth

Since log(1+u)<=u for u>=0,

    epsilon_r(alpha) <= exp(r(r-1)/(2 sqrt(alpha))) - 1.    (6)

In the quartic window n^4/4<=p<=n^4 with n>=4,

    alpha >= n(1-n/p) >= n-4/n^2 >= n-1.

Therefore epsilon_(2s)(alpha)=O(s^2/sqrt(n)) whenever
s^2/sqrt(n) tends to zero. In particular (2) gives

    Q_s = (1+o(1)) T_s uniformly for s=O(log n),           (7)

along admissible quartic fields. More generally this holds for
s=o(n^(1/4)). The constants and the characteristic condition are
explicit; no theorem at fixed s is being promoted to growing s.

The estimates can be coarse at modest orders. For example, taking
s=3 log_2 n and replacing sqrt(alpha) by floor(sqrt(n-1)) gives
an exact rational error upper bound below one at n=2^30, but not at
the sampled n=2^28. The recorded upper bounds are about 0.32249 at
n=2^32, s=96, and 0.02773 at n=2^40, s=120. These are parameter
bounds for all eligible fields, not assertions that particular primes
were constructed or moment values computed at these orders.

For a multiplicative subgroup H of size n, eta is constant on each
nonzero H-coset. Its maximum M therefore satisfies

    M^(2s) <= (p/n) T_s.

As a conditional consequence, if epsilon<=1/2 and one proves

    Q_s <= (K s n)^s,

then M <= (2p/n)^(1/(2s)) sqrt(K s n). Taking
s>=log(2p/n) makes the first factor at most sqrt(e).
Thus a Gaussian-scale upper bound for the centered distinct count
would suffice for the subgroup square-root scale up to logarithms.
The missing Gaussian upper bound is not proved here. Neither this
subgroup criterion nor its proof is asserted to settle the full
classical arbitrary-set Paley conjecture or the official code prize.

## 4. A stronger bound for repeated six-word relations

This section assumes H<=F_p^*, -1 in H, n=|H|, and p odd. Put

    E_t = (1/p) sum_a |eta(a)|^(2t)

for real t>=1, so E_1=n and at integer s the quantity E_s counts
zero-sum 2s-words. Let R_(2s) count all such words with at least
one pair of equal coordinates. In particular it includes the
opposite-free repeated class from pass25.

Fix an equal positional pair and call its common value h. Scaling by
h^(-1) shows that the remaining 2s-2 entries sum to -2. Consequently
the number for this marked pair is n r_(2s-2)(-2). Since r_k(t) is
constant on nonzero H-cosets,

    n r_(2s-2)(-2)
      = (1/p) sum_a eta(a)^(2s-2) eta(2a)
      <= (1/p) sum_a |eta(a)|^(2s-2) |eta(2a)|
      <= E_(s-1/2).                                      (8)

The sign in eta(2a) is immaterial because -H=H. The final inequality
is Holder with exponents (2s-1)/(2s-2) and 2s-1; a -> 2a is a
permutation because p is odd. A union bound now gives

    R_(2s) <= binomial(2s,2) E_(s-1/2)
      <= binomial(2s,2) n^(1/(2s-2)) E_s^(1-1/(2s-2)).      (9)

The last step interpolates between E_1=n and E_s. Jensen gives
E_s>=n^s, and hence

    R_(2s)/E_s <= binomial(2s,2)/sqrt(n).                   (10)

This is an uncentered statement. It alone does not imply (1); the
partition-by-partition principal subtraction in Section 2 is needed
for that conclusion.

At s=3, interpolation between E_2 and E_3 also gives

    R_6 <= 15 sqrt(E_2 E_3).                              (11)

The already used [MRSS source, Theorem 3 and Corollary 7](https://arxiv.org/html/1712.00410v1)
provides E_2<<n^(49/20)(log n)^(1/5) and
E_3<<n^4 log n for multiplicative subgroups of size at most sqrt(p).
These hypotheses hold in the quartic window. Substitution into (11)
proves

    R_6 << n^(129/40) (log n)^(3/5).                      (12)

This improves the project's earlier n^(69/20)(log n)^(1/5) estimate
for its opposite-free repeated subset, and applies to the larger class
of all repeated six-words. The power saving relative to the previous
bound is 9/40. It is not an improvement of the full E_3 bound, since
the distinct-coordinate part remains uncontrolled.

## 5. Exact verification and limits

The [verifier](../experiments/parallel33_verify_2026_09_05.py) computes
zero-sum counts for all integer partition types at orders 4,6,8. A
separate subset dynamic program computes distinct zero-sum counts
directly, as coefficients of product_(x in A)(1+z X^x) modulo X^p-1.
It does not use inclusion-exclusion. The two methods agree exactly.

The 33 set/order cases cover nine multiplicative subgroups and two
additional symmetric sets, including a symmetric interval containing
zero. All principal terms, mixed-partition bounds, Jensen bounds,
Stirling identities, and norm comparisons use rational/integer
arithmetic. The parameter table likewise stores exact rational upper
bounds; its decimal column is only for display.

The [record](../results/parallel33_verification_2026_09_05.json)
contains 740 checks, finishing in about 0.63 seconds in the recorded
run. It also checks a failure outside the characteristic hypothesis:
at p=3, A={1,2}, r=6, the partition (3,3) has centered count 8/3,
while T_3=2/3. Its zero block coefficients invalidate (4).

The uniform arguments remain ordinary proofs, not a conclusion from
these finite checks. No Lean job, hosted proof submission, or separate
author review was performed. The remaining centered distinct-coordinate
upper bound, the full signed spectral estimate, and all requested full
targets remain open in this project.
