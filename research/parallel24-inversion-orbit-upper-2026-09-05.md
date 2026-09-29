# A bounded-moment representative in every inversion family

**Status: a uniform existence theorem for an unsigned B_h inverse of
every input set. The moment of the original input remains unbounded
at the desired scale.** This is stronger placement than a statement
about almost all sets in the whole slice, but its quantifiers do not
supply the uniform hypothesis of SS-B*. No novelty claim is made.

Let p be an odd prime, χ(0)=0, and C⊂F_p have size k≥1. Write

    F_C(x)=Σ_(c∈C)χ(x−c), M_(2r)(C)=Σ_x F_C(x)^(2r),
    D_z={1/(c−z):c∈C}, z∉C,
    Q_C(x,z)=Σ_(c∈C)χ(x−c)χ(z−c),
    A_r=(2r−1)!!, r≥1.

A B_h set has only trivial equal h-term sums, with repetition allowed.
For h≥2 and p>h define, as in [pass 20](parallel20-inversion-moments-2026-09-05.md),

    N_h(k)=binom(k+h−1,h),
    Q_h(k)=k+(2h−2)binom(N_h(k),2),
    U_r(p,k)=A_r p²k^r+(2r−1)²p k^(2r).

The letter Q with two field arguments denotes the character kernel;
the subscript-h expression is the pole-count bound.

## 1. Finite theorem and quantifiers

If p>Q_h(k), then for **every** set C of size k there is a pole z∉C
such that D_z is B_h and

    M_(2r)(D_z) ≤ U_r(p,k)/(p−Q_h(k)).                 (1)

More generally, for every real λ>1, at least

    (1−1/λ)(p−Q_h(k))

poles have both the B_h property and the bound λ times the right side
of (1). A real lower bound for an integer count has its usual meaning.
This is an unsigned moment estimate. The chosen pole may depend on C,
r, h, and λ. It need not work simultaneously at other moment orders.

## 2. Exact inversion identity

For x≠z, set t=1/(x−z). Then

    t−1/(c−z)=(c−x)/[(x−z)(c−z)]

and multiplicativity, including both factors χ(−1), gives

    F_(D_z)(t)=χ(x−z)Q_C(x,z).                        (2)

At t=0 we have F_(D_z)(0)=F_C(z), while Q_C(z,z)=k because z∉C.
Thus

    M_(2r)(D_z)=Σ_x Q_C(x,z)^(2r)−k^(2r)+F_C(z)^(2r)
               ≤Σ_x Q_C(x,z)^(2r).                  (3)

The boundary term is included exactly. The argument works for both
congruence classes of odd primes; it does not assume χ(−1)=1.

## 3. Squared complete sums and good poles

Expanding before summing in the two row variables gives the exact
nonnegative identity

    Σ_(x,z)Q_C(x,z)^(2r)
      =Σ_(c_1,...,c_(2r)∈C) [Σ_x χ(Π_i(x−c_i))]².   (4)

For an ordered tuple in which every element has even multiplicity,
the inner sum has absolute value at most p. There are at most A_r k^r
such tuples: pair the positions and assign a value to every pair.
Every such tuple admits at least one pairing, and overcounting is safe.

For any other tuple, Weil's quadratic-character bound gives absolute
value at most (2r−1)√p. The imported input is
[McDonald–Sahay–Wyman, Lemma 2.1](https://arxiv.org/html/2210.03789v2),
specialized to a monic polynomial and the quadratic character. To see
explicitly that repeated roots cause no issue, remove the even-positive
multiplicity roots. If there are s odd-multiplicity roots and e removed
roots, then s≥2; the sum differs from its squarefree degree-s sum at
at most e points. Its absolute value is therefore at most
(s−1)√p+e≤(s+e−1)√p≤(2r−1)√p.

It follows that

    Σ_(x,z)Q_C(x,z)^(2r) ≤ U_r(p,k).                  (5)

This is the standard squared-complete-sum ingredient, already derived
for an interval in [the sigma sum-product note, Lemma D](sigma-sumproduct-2026-09-05.md).
Its proof only uses cardinality, so it applies to arbitrary C. The
new combination here is with the inversion identity and the explicit
count of good poles, not a new form of Weil's theorem.

For completeness, the pole count does not assume independence of
additive relations and character moments. For two distinct h-element
multisets of C, cancel their common summands. Their reciprocal-sum
difference is Σ_c ν_c/(c−z), with nonzero integers |ν_c|≤h and Σν_c=0.
The numerator after clearing denominators is nonzero: evaluate it at
one of its distinct poles. Its degree is at most 2h−2. Since p>h,
none of its nonzero coefficients ν_c vanishes modulo p. Summing this
root bound over unordered pairs of multisets, and excluding C, bounds
the total number of bad poles by Q_h(k).

Let G be the good poles. All terms in (3)–(5) are nonnegative, so

    Σ_(z∈G)M_(2r)(D_z)≤U_r(p,k), |G|≥p−Q_h(k)>0.

Averaging proves (1). At most (p−Q_h(k))/λ good poles can exceed
λU_r/(p−Q_h(k)); subtracting this number from |G| proves the stated
positive-fraction version. No unproved distribution assertion enters.

## 4. The single-size slice

Fix r≥3 and 2≤h≤floor((r+1)/2), before letting p grow. At
k=floor(p^(1/(r+1))), the theorem applies eventually, since

    Q_h(k)/p → 0                    if 2h<r+1,
    Q_h(k)/p → (h−1)/(h!)²          if 2h=r+1.

Indeed N_h(k)∼k^h/h!, and k^(r+1)/p→1. The boundary constant is
at most 1/4. Also k^r/p→0. Dividing (1) by p k^r proves that every
C has a B_h inverse satisfying

    M_(2r)(D_z) ≤ [A_r/(1−b_h)+o(1)]p k^r,          (6)

where b_h=0 in the strict interior and b_h=(h−1)/(h!)² at the boundary.
The error is uniform in C, for the fixed r,h. In particular, at r=3,
h=2, k=floor(p^(1/4)), some Sidon inverse of every C satisfies

    M_6(D_z) ≤ (20+o(1))p k³.                       (7)

This is a proved existence statement at the actual character kernel.
The coefficient 20 comes from 15/(1−1/4), not an empirical fit.

## 5. Why this has not proved the original bound

The inverse that transports the original moment has the specific
coefficients w_d=χ(c−z)=χ(d). The earlier exact identity is

    M_(2r)(D_z,w)=M_(2r)(C)−|F_C(z)|^(2r)+k^(2r)
                  ≥M_(2r)(C).                       (8)

Equation (1) controls M_(2r)(D_z,1). It supplies no upper bound for
M_(2r)(D_z,w). The uniform-over-all-B_h-sets hypothesis in pass 20
allowed removal of signs by averaging completions; having one good
inverse in each family does not bound those particular completions.

For example, let E={d∈D_z:w_d=−1}. Then F_(D_z,w)=F_(D_z)−2F_E.
The triangle inequality introduces the moment of this specific E;
the present existence theorem does not bound E itself. Inverting E
again supplies a different unsigned set and another transported sign
vector, with the same missing step. This is the precise remaining
quantifier gap, rather than a claim that sign removal is impossible.

Consequently no uniform SS-B* estimate, arbitrary-two-set Paley
cancellation exponent, full spectral edge, or prize theorem follows
from (1) alone. The full goal remains open.

## 6. Verification

The [exact verifier](../experiments/parallel24_inversion_orbit_upper_2026_09_05.py)
and [result ledger](../results/parallel24_inversion_orbit_upper_2026_09_05.json)
check actual finite-field inversion coordinates, unsigned and signed
boundary terms, the two-row moment bound, B_h pole counts, and the
good-pole existence and fraction estimates. A separate small-field
ordered-tuple expansion checks (4) and each individual Weil bound.
The script checks p²k⁶≤2⁶³−1 before its fixed-width integer array
computations, so all sums and powers in the stated coverage are exact;
the scalar bounds use Python integers. Finite tests support implementation
and boundary checks but do not prove the asymptotic theorem. Primary
HTML was read live in this pass and its existing archive is pinned.
No PDF, human referee report, or Lean proof is claimed for this lane.
