# Recovering endpoint energy from first-precision triple mass

The Paley conjecture and Proximity Prize remain unproved. This note gives
a more direct quantitative target for the existing energy route. It does
not prove the required uniform bound on that target. The recovery bounds below
are ordinary derivations by root; the imported energy theorem is
identified separately. No novelty, independent review or Lean
formalization is claimed.

## 1. Setup and the two-sided estimate

Fix a dyadic N>=4 and an odd prime p=1 mod N. For each dyadic s dividing
N, let H_s be the subgroup of order s and E_s its ordered additive energy.
At a step s=2k, put K=H_k and L=H_s\K. Let c_b be the multiplicities of
the primitive polynomial R_s modulo p, and retain the notation from
[pass 48](parallel48-collision-eliminants-2026-09-06.md):

    Y_s = sum_(c_b>=3) c_b(c_b-2),
    B_s = ||1_K*1_L||_2^2 = s sum_b c_b^2,
    T_s = <1_K*1_K, 1_K*1_L>,
    E_s = 2E_k + 6B_s + 8T_s.                           (1)

These are actual first-precision field multiplicities. No p-adic
valuation or higher-precision correction is included in Y_s. Define

    Ytot_N = sum_(s=4,8,...,N) Y_s,
    S_N = sum_(s=4,8,...,N) sqrt(Y_s),
    ell_N = log_2(N)-1.

Then

    3N + 6N Ytot_N <= E_N,
    sqrt(E_N) <= sqrt(14N^2-25N) + 2sqrt(2N) S_N,        (2)
    E_N <= 28N^2 + 16N ell_N Ytot_N.                       (3)

The lower bound and the upper bound in (3) compare the same unweighted
total triple mass, with a logarithmic gap. Thus the intermediate energy
target E_N<<N^(7/3) necessarily implies Ytot_N<<N^(4/3). Conversely,
Ytot_N<<N^(4/3) implies E_N<<N^(7/3) log N. For the bound without that extra
logarithm it suffices to prove S_N<<N^(2/3), or the stronger condition
Ytot_N<<N^(4/3)/ell_N. None of these uniform inputs is proved here.

## 2. Subtract the two-representation baseline before Cauchy

Write a(x)=(1_K*1_K)(x), r(x)=(1_K*1_L)(x), and
q(x)=max(r(x)-2,0). Pointwise r(x)<=2+q(x), while sum_x a(x)=k^2.
Therefore Cauchy gives

    T_s <= 2k^2 + sqrt(E_k sum_x q(x)^2).

The primitive-fiber identity says that each c_b is repeated at exactly
s different nonzero sums. Consequently

    sum_x q(x)^2 = s Z_s,
    Z_s = sum_(c_b>=3) (c_b-2)^2 <= Y_s.

We obtain the useful replacement for the unsplit Cauchy estimate:

    T_s <= 2k^2 + sqrt(E_k sY_s).                       (4)

Also B_s<=s^2/2+sY_s, as already proved in pass 48. Combining this with
(1) and (4) gives

    E_s <= 2E_k + 7s^2 + 6sY_s + 8sqrt(E_k sY_s).       (5)

The baseline contribution is quadratic. Only the tail beyond two
representations enters the square-root term.

## 3. Recovery from any intermediate subgroup

In fact, for every dyadic 2<=M<=N,

    (N/M)E_M + 6N sum_(M<s<=N) Y_s <= E_N,              (6)

    sqrt(E_N)
      <= sqrt((N/M)E_M + 14N(N-M))
         + 2sqrt(2N) sum_(M<s<=N) sqrt(Y_s).             (7)

All sums over s in this note are dyadic. For (6), use B_s>=sY_s and
T_s>=0 in (1), and iterate E_s>=2E_k+6sY_s. Each contribution at level
s is multiplied by N/s. No lower bound on the smaller levels is assumed.

Here is a proof of the radical estimate that keeps its constants. Put
e_s=E_s/s^2 and y_s=Y_s/s. Equation (5) becomes

    e_s <= e_k/2 + 7 + 6y_s + 4sqrt(e_k y_s).

With Phi(u)=sqrt(u^2/2+7), it follows that

    sqrt(e_s) <= Phi(sqrt(e_k)) + 2sqrt(2y_s).           (8)

Indeed, squaring the right side gives at least the previous right side:
Phi(u)>=u/sqrt(2), and the y_s square term is 8y_s rather than 6y_s.
The function Phi is increasing and has Lipschitz constant 1/sqrt(2)
on the nonnegative real axis, since

    Phi'(u)=u/[2sqrt(u^2/2+7)] <= 1/sqrt(2).

Start the baseline recursion b_M=sqrt(e_M), b_s=Phi(b_k).
Its square at the endpoint is

    b_N^2 = (M/N)e_M + 14(1-M/N).

Inductively bounding the positive difference sqrt(e_s)-b_s in (8) gives

    sqrt(e_N) <= b_N
       + 2sqrt(2) sum_(M<s<=N) (s/N)^(1/2) sqrt(Y_s/s).

Multiplying by N proves (7). At M=2, E_2=6, so (6)-(7) give (2).
Finally use (a+b)^2<=2a^2+2b^2 and S_N^2<=ell_N Ytot_N to obtain (3).
The number of levels with positive Y_s can replace ell_N in (3).

## 4. Why this improves the previous sufficient criterion

Pass 48 used

    V_N = sum_s (3/4)^(log_2(N/s)) Y_s/s.

Writing j=log_2(N/s), Cauchy with weights (3/2)^j gives

    S_N^2 <= [sum_j (2/3)^j] [sum_j (3/2)^j Y_s]
           <= 3N V_N.                                 (9)

Thus the earlier V_N<<N^(1/3) hypothesis implies the new sufficient
S_N<<N^(2/3) hypothesis. The new condition does not give exponentially
larger weights to collisions at smaller levels. For arbitrary
nonnegative mass profiles, the reverse comparison has no uniform
constant: a single positive mass at depth j has N V_N/S_N^2=(3/2)^j.
This observation compares the recovery estimates; it is not a claim
that such an arbitrary profile comes from a field.

The logarithmic step S_N^2<=ell_N Ytot_N cannot be improved using only
nonnegative masses and the degree budget. For example, let N=2^(3m),
c=2^m, and at each level 4c<=s<=N assign one cluster of size c and fill
the remaining degree s/4 with singletons. Then Y_s=c(c-2) at every
one of the 2m-1 assigned levels, and

    S_N^2 = (2m-1) Ytot_N.

These are individually admissible multiplicity lists, not a compatible
subgroup tower. They establish only the limitation of recovering S_N
from its total Ytot_N by the compared scalar constraints.

## 5. Smaller subgroups can be handled with the existing bound

Even the trivial E_M<=M^3 in (7) handles the base term when M<=N^(2/3).
The known uniform energy estimate allows a substantially higher cutoff.
The previously audited MRSS source, [Corollary 12, equation (32)](https://arxiv.org/abs/1712.00410), states

    E(H_M) << M^(49/20) log^(1/5) M,  when M<=sqrt(p).

The complete statement was reread in the archived primary text in
this pass. No newer or stronger estimate is imported. In the quartic
range N^4/4<=p<=N^4 and N>=4, every M<=N satisfies its size hypothesis.
Write L_N=1+log N and, for sufficiently large N, choose the largest
dyadic M no greater than

    N^(80/87) / L_N^(4/29).                             (10)

The finite small-N cases can be absorbed in the implicit constants.
Since (49/20)-1=29/20, the propagated base energy in (7) satisfies

    (N/M) E_M << N M^(29/20) L_N^(1/5) << N^(7/3).

The exponent and logarithm balances are exact:

    1 + (80/87)(29/20) = 7/3,
    -(4/29)(29/20) + 1/5 = 0.

Consequently it now suffices to prove just

    sum_(M<s<=N) sqrt(Y_s(p)) << N^(2/3),                (11)

uniformly over the target quartic primes, with M given by (10).
The smaller levels have been absorbed into the existing source theorem.
The same conclusion would follow from a sufficiently small total mass
over these upper levels, using their count in Cauchy. Equation (11) is
still unproved. In terms of powers of N, the unresolved interval starts
at 80/87, with the displayed logarithmic correction. The number of
remaining levels is

    (7/87) log_2 N + (4/29) log_2 L_N + O(1).

This is a shorter growing interval, not a bounded number of levels.

The [previous energy-to-triangle implication](parallel43-energy-prime-exceptions-2026-09-06.md)
then gives the intermediate W<<N^(17/3) if (11) holds. This would still
not prove the full Paley conjecture or the prize goal. Merely substituting
the current 49/20 bound for each upper-level mass recovers the existing
energy exponent; it does not establish (11). Explicitly, (1) gives
Y_s<=E_s/(6s), so the known theorem only yields a square-root mass sum
of order N^(29/40)L_N^(1/10). The gap from 29/40 to the required 2/3
is 7/120, and the recovery returns energy power 1+2(29/40)=49/20.

The critical content from [pass 49](parallel49-critical-content-2026-09-06.md)
still gives Y_s<=v_p(Ccrit_s), so it supplies a stronger arithmetic
version of (11). It is sufficient, not necessary: higher-precision mass
and reciprocal cancellation can overcount the actual first-precision
quantity. The direct field version (11) avoids those costs.

## 6. Exact finite verification and scope

The [standard-library checker](../experiments/parallel50_tower_mass.py)
uses all 70 certified splitting-prime cases from pass 48. At every
dyadic step from order four upward it independently counts parent, child
and mixed pairs, recovers the primitive multiplicities, and checks the
tail inequality, recurrence and both endpoint bounds. It also checks
both bounds starting at every intermediate cutoff.

All 390 level checks, 460 cutoff checks and 1,689,552 literal pair
enumerations passed. Radical inequalities are checked against rational
lower bounds for their proposed upper-bound expressions; passing those
checks therefore does not depend on floating-point rounding. Rational
upper bounds separately validate the comparison in (9). The exact
cutoff exponents and abstract equal-mass examples are also checked.

The three order-128 quartic cases still have total tower mass Ytot_128=3
and energies 54912, 57984 and 54912. These repeat the earlier finite
evidence; no new prime classification is claimed. Full output is in
the [results](../results/parallel50_tower_mass_2026_09_06.json).

No new uniform energy, triangle, period or absolute exception exponent
is established. No Lean process or Prove2Me submission was changed.
The next mathematical task is (11), or a stronger sufficient estimate,
in the upper part of the actual splitting-prime tower.
