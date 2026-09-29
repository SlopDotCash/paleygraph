# Reflection smoothing: a positive moment bound with an exact recovery loss

This tests a concrete affine alternative to the uniform-translation
smoothing in pass41. It averages over square scalings about random centers,
smooths the test set at the same time, and retains an exactly recoverable
pairing in expectation. The resulting moment bound is positive, but the
specified recovery certificate is always at least as large as direct
sixth-moment Holder. Thus this route does not improve the original Paley
discrepancy bound.

The affine lane returned the key identities before reaching its usage
limit. Root independently derived the recurrences, constant, and recovery
comparison below, and wrote the exact checker. No Lean build was started.

## 1. Reflections as orthogonal projections

Let p be prime with p=1 mod4. Let chi be the quadratic character,
let A be a nonempty subset of F_p of size m, and let
V subset F_p have size n>=2. Put

    F(x)=sum_(v in V) chi(x-v), M_j(F)=sum_x F(x)^j.

For a center c define

    P_c G(x)=[G(x)+G(2c-x)]/2.

Since -1 is a square, reflection of this character convolution corresponds
to reflection of V with the same sign. Each P_c is a self-adjoint
orthogonal projection in counting L2, preserves total mass, and contracts
every Lq norm for q>=1. Starting with F_0=F and A_0=1_A, apply the same
independent uniformly chosen centers to both:

    F_j=P_(c_j)F_(j-1), A_j=P_(c_j)A_(j-1), t=2^r.

Always sum F_j=0, sum A_j=m=|A|, and 0<=A_j<=1. The shifted branches
are not treated as independent samples; the projection recurrence accounts
for their correlations.

## 2. Exact sixth-moment recurrence

For a fixed zero-mean G, averaging over c makes 2c-x uniform for each x.
Binomial expansion gives

    E_c M6(P_c G)
      =M6(G)/32+15 M4(G)M2(G)/(32p)+5 M3(G)^2/(16p).           (1)

Every term is nonnegative. Therefore B_r=E M6(F_r) satisfies

    B_r>=M6(F)/t^5.                                           (2)

Also E M2(F_j)=M2(F)/2^j. Under n^4<=p, the same individual Weil
inputs checked in pass41 give

    M2(F)<=pn, M4(F)<=6pn^2, M6(F)<=pn^3(15+5n).

Pointwise L4 contraction gives M4(F_j)<=M4(F), and Cauchy gives
M3(F_j)^2<=M4(F_j)M2(F_j). Substituting these into (1) yields

    B_j<=B_(j-1)/32+(150/32)pn^3/2^(j-1).

Summing the geometric recurrence gives the useful explicit bound

    B_r<=M6(F)/t^5+(10pn^3/t)(1-t^(-4)).                      (3)

If t>=n^(1/5) and t>=2, then

    B_r/(pn^3)<=5+10/t+5/t^5<=325/32.                        (4)

Choosing r=ceil(log_2(n^(1/5))) supplies such a t with
t<2n^(1/5). This proves a positive uniform smoothed sixth-moment
bound without asserting independence of the reflection branches.

## 3. Exactly recovering the original pairing in expectation

Write S=sum_(a in A) F(a). Since P_c is self-adjoint and idempotent,
and E_c P_c G=G/2 for zero-mean G,

    E <A_r,F_r>=S/t.                                         (5)

The second moment of the smoothed indicator is also exact:

    E ||A_r||2^2=m^2/p+(m-m^2/p)/t.                          (6)

Indeed its one-step recurrence is E_c||P_c A||2^2=
||A||2^2/2+m^2/(2p), and its initial squared norm is m.

Holder in x and interpolation give

    |<A_r,F_r>|<=m^(2/3)||A_r||2^(1/3)||F_r||6.

Holder in the probability space, then (5)-(6), gives

    (|S|/(mn))^6
      <= t^5 [1/m+(t-1)/p] B_r/n^6.                          (7)

This recovery avoids an unspecified unsmoothing error: its cost is an
explicit factor. But (2) shows that the right-hand side of (7) is at least

    [1+m(t-1)/p] M6(F)/(m n^6)
      >= M6(F)/(m n^6).                                      (8)

The final expression is precisely direct sixth-moment Holder for 1_A.
Thus this particular recovery certificate can never improve that direct
certificate, even if B_r is known exactly. This is a comparison of upper
bound expressions, not a lower bound on the true discrepancy.

On the critical slice n=p^(1/4+o(1)), using (4) in (7) recovers the
direct Weil threshold |A|>p^(1/2+epsilon): the normalized sixth-power
bound has a main term O(p/(|A|n^2)). It does not reach the desired
original Gaussian moment scale. The obstruction concerns this joint
reflection/projection/Holder calculation. It does not exclude every
affine method, signed average, or a more character-sensitive recovery.

## 4. Exact validation and source scope

The [checker](../experiments/parallel42_classical_affine.py) enumerates
reflection-center tuples in six fixed cases, retaining nonzero third
moments, and checks (1)-(8) by integer arithmetic and rational fractions.
The [results](../results/parallel42_classical_affine_2026_09_06.json)
separate the exact certificate comparison from any claim about the original
Paley bound. The only imported character estimate is the individual Weil
bound already scoped and reviewed in
[pass41](parallel41-classical-smoothing-2026-09-06.md).
