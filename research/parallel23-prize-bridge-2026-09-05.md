# What a subgroup-period theorem would still need to reach the pinned prize

Status: a conditional list formula for a restricted polynomial family,
an obstruction to one proposed strengthening, and exact finite counting
obligations for the pinned official certificate. No prime-field period
bound, general list bound, MCA upper bound, or prize certificate is proved.
All official parameters below concern the archived September-4 contract,
not a claim about the current website, leaderboard or latest repository.

The new point is quantitative. The first a coefficient conditions of
a root polynomial correspond to polynomial phases of degrees through a.
For a proportional to the domain size, a uniform square-root bound for
all those phases is impossible by an elementary pigeonhole argument on
the actual prime-field domain. A successful bridge must use more than
the maximum of these phases. This does not contradict the linear-period
target and does not assert that the associated codes have large lists.

## 1. The pinned target and the integer budget

The archived [profile](../sources/official-prize-2026-09-04/ProximityPrize/Benchmark/IRSProfile.lean)
at contract b34c0131cfa36b51111521541d7d3e35c8791082 has

    p=2130706433,       q=p^6,
    H=μ_(262144)⊂F_p*⊂F_q,
    n=262144,          k=131072,
    C_r=RS[F_q,H,k] with r interleaved scalar rows.

Its base code is C_8 and its twice-interleaved list code is C_16.
Column agreement requires every scalar row to agree. The archived
[IRS equality, lines 669 onward](../sources/official-prize-2026-09-04/dependencies/ArkLib/ProofSystem/ToyProblem/Impl/IRS.lean)
identifies the combination-round certificate as

    Γ(δ)=MCA(C_8,δ)+Λ(C_16,δ)/q.                       (1)

The [lower contract](../sources/official-prize-2026-09-04/ProximityPrize/Benchmark/TargetLower.lean)
requires Γ(δ)≤2^(−128), an admissible radius, and a separate
spot-check score inequality. Fix an integer 1≤j≤131071 and put

    δ=j/n,       s=n−j≥k+1.

This radius lies strictly below rate-half capacity and is admissible.
If B_s is the maximum number of MCA-bad scalars and L_s the maximum
C_16 list cardinality, (1) gives the exact integer condition

    Γ(j/n)≤2^(−128)
      iff B_s+L_s≤floor(q/2^128)
      iff B_s+L_s≤274980728111395087.                  (2)

No numerical bound on B_s or L_s is obtained here. Equation (2)
records their joint budget; separate bounds using this whole budget
for each term would not suffice.

The working subgroup theorem concerns p in a quartic window around
the subgroup size. The displayed p is less than n^4/4, while q is
greater than n^4. Thus granting that theorem in its stated window
does not instantiate it at this profile. The [trace audit](official-profile-and-trace.md)
also rules out changing p to q while retaining H in the base field:
nonzero trace-zero frequencies have full period n. The following
analysis grants any needed linear-period estimate as an additional
hypothesis, so the remaining problem is not hidden by that mismatch.

## 2. An exact conditional bridge for monic received polynomials

Let p be prime, F_p⊂F_q any finite extension, D⊂F_p have n<p
distinct points, and 1≤k<s=k+a≤n with a≥1. Let W∈F_q[X]
be monic of degree s, and let c_j be its coefficient of X^(s−j),
for 1≤j≤a. At radius δ_a=1−s/n,

    |List_(RS[F_q,D,k])(W,δ_a)|
      =#{T⊂D: |T|=s, (-1)^j e_j(T)=c_j for 1≤j≤a}.
                                                               (3)

If any c_j is outside F_p, this list is empty. Otherwise define
v_1,…,v_a∈F_p recursively by Newton's identities:

    v_j=−j c_j−Σ_(i=1)^(j−1)c_i v_(j−i).

Because a<p, the coefficient conditions in (3) are equivalent to

    Σ_(x∈T)x^j=v_j,       1≤j≤a.                      (4)

Proof of (3): for a nearby polynomial F of degree less than k,
W−F is monic of degree s and has at least s roots in D. It therefore
equals the unique polynomial Q_T=∏_(x∈T)(X−x) for a size-s set T.
The coefficients of degrees s−1 through k must match W. Conversely
these matches make W−Q_T have degree less than k. Distinct root
sets give distinct codewords, and no extra root is possible. The
lower coefficients of W may be arbitrary extension-field elements;
only the displayed top coefficients must lie in the base field.

Define the polynomial-phase maximum

    M_a(D)=max_(b∈F_p^a, b≠0)
             |Σ_(x∈D)e_p(b_1x+…+b_a x^a)|.

For an integer U≥max(1,M_a(D)), put t=min(s,n−s). Fourier inversion
on F_p^a and the same generating-function argument as in the
[one-coefficient subset note](subset-sums-and-lists.md) give

    |N_s(v)−binom(n,s)/p^a|
       ≤(p^a−1)/p^a · binom(U+t−1,t),                 (5)

where N_s(v) counts (4). For completeness, for a nonzero b the
generating product is

    ∏_(x∈D)(1+z e_p(P_b(x))).

The jth coefficient of its formal logarithm uses the phase jP_b.
For j≤s<p, its coefficient vector remains nonzero, so M_a bounds
every power sum in this logarithm. Taking absolute values in its
exponential bounds the s-coefficient by binom(U+s−1,s). Complementing
subsets gives t. Summing the p^a−1 nonprincipal frequencies proves
(5). This is a valid conditional upper and lower list estimate for
(3), over every extension F_q of the same base field.

For a=1, M_a is exactly the linear additive-character maximum in
the working subgroup target when D=H. For a>1, (5) requires a new
input. Even for a subgroup, pure monomials only give

    Σ_(h∈H)e_p(bh^j)=gcd(j,n)·η_(H^j)(b).

Mixed phases such as b_1h+b_2h² are not covered by that identity.
Moreover H^j may lie outside the original theorem's quartic window.

An exact small-field warning against reusing the same numerical
linear bound is D=F_17*, where M_1(D)=1. Among its eight-element
subsets there are 54 with both first and second power sums zero.
Their uniform mean is binom(16,8)/17²=12870/289. Inserting U=1
into the two-dimensional version of (5) would allow error at most
288/289, whereas the actual excess is 2736/289. The proposed
substitution is false. This example concerns that numerical
substitution, not the conjectured linear square-root bound.

## 3. Why uniform higher polynomial-phase square-root bounds fail

The following statement holds for every actual D⊂F_p with n distinct
points, without replacing it by an abstract kernel.

**Proposition.** Let 1≤a<n and let R≥1 be an integer satisfying
R^n<p^a. Then a nonzero polynomial P=Σ_(j=1)^a b_jX^j exists with

    Re Σ_(x∈D)e_p(P(x)) > n(1−20/R²).                (6)

Proof: the p^a coefficient vectors give p^a distinct evaluation
vectors, since degree at most a<n makes evaluation injective.
Represent each coordinate in {0,…,p−1} and divide this interval
into R real intervals [lp/R,(l+1)p/R). There are at most R^n
boxes. Two distinct evaluation vectors fall in one box. Their
difference is the evaluation of a nonzero polynomial with zero
constant coefficient and degree at most a; in every coordinate it
has an integer representative d_x with |d_x|<p/R. Therefore

    Re e_p(P(x))=cos(2πd_x/p)
      ≥1−2π²(d_x/p)²>1−20/R²,

using π²<10. Summing proves (6).

For any sequence with p→∞ and a/n→ε>0, choose
R=floor(p^(a/(2n))). Eventually R→∞ and R^n<p^a, so

    M_a(D)/n → 1.                                    (7)

In particular this applies to actual dyadic subgroups in the working
quartic prime family. It precludes M_a=O(√(n log p)) with a constant
independent of n when a is proportional to n. Granting the linear
period conjecture does not remove this unconditional obstruction.
At a fixed positive gap below capacity, a/n is positive in (3).
Thus a maximum bound for all the higher polynomial phases cannot
be the proposed fixed-gap bridge.

There is also an exact finite instance in the pinned official field:

    a=26215,       s=k+a=157287,
    δ_a=104857/262144,       R=8.

Since p>2^30 and 30a−3n=18>0, p^a>8^n. Proposition (6) gives
a nonzero degree-at-most-26215 polynomial over F_p with

    Re Σ_(h∈H)e_p(P(h)) > (11/16)n=180224.            (8)

This is an existence proof, not an explicitly constructed polynomial
at the production size. It does not assert a large list or MCA error
at this radius: Fourier inversion still sums many coefficients with
phases, and a large individual polynomial-phase maximum need not
prevent cancellation of that aggregate. It also does not contradict
a linear-period estimate, since the polynomial can have higher degree.
The particular finite bound (8) is not asserted at every proposed
official radius, including radii nearer 0.464.

## 4. Exact remaining counts for arbitrary official words

Here is a concrete statement that would actually close (2), rather
than a replacement by monomial words. Fix j,s as in Section 1.
Every r-row received word has a unique coordinatewise polynomial
interpolant W=(W_1,…,W_r) with degrees less than n over F_q.
For every size-s subset T⊂H, write

    rem_T(W_i)=W_i mod Q_T,       deg rem_T(W_i)<s,
    ρ_T(W)=(coefficients k,…,s−1 of each rem_T(W_i)).

Thus ρ_T(W) lies in F_q^(r(s−k)), and W has a codeword explanation
on T exactly when ρ_T(W)=0. If it does, the explanation is unique:
its rows are rem_T(W_i), each of degree less than k.

Define the set of **distinct** explanation vectors

    F_s(W)={ (rem_T(W_i))_(i=1)^r : |T|=s, ρ_T(W)=0 }.

Then |F_s(W)| is exactly the list cardinality at radius j/n.
Counting subsets T instead would overcount a codeword once for every
size-s subset of its agreement set; if W is itself a codeword, every
T qualifies but the list can still be a singleton.

For two r-row words W,Z define

    R_s(W,Z)={γ∈F_q: some |T|=s satisfies
                  ρ_T(W)+γρ_T(Z)=0 and ρ_T(Z)≠0}.

This is exactly the set of MCA-bad scalars at radius j/n. Linearity
of polynomial remainder gives the folded explanation criterion. If
ρ_T(Z)=0, a folded explanation also forces ρ_T(W)=0, so the pair
is jointly explainable and is not an MCA witness. When ρ_T(Z)≠0,
at most one γ can solve the displayed vector equation.

To justify using size exactly s, start with any larger MCA witness.
Its direction Z cannot have a joint codeword explanation there.
Choose k points and interpolate every row. A violating row supplies
one additional point, giving a nonexplainable (k+1)-subset. Extend
it to size s within the witness. Nonexplainability persists and the
folded explanation restricts. This uses s≥k+1, which is why the
radius range was fixed explicitly above.

Consequently the precise missing counting assertion at j/n is

    max_(W,Z∈(F_q[X]_<n)^8) |R_s(W,Z)|
      + max_(V∈(F_q[X]_<n)^16) |F_s(V)|
        ≤274980728111395087.                          (9)

Equation (9) is an exact polynomial form of the needed certificate,
not a claimed easier theorem or a proved estimate. It states all
field, domain, row-width and radius quantifiers. A useful analytic
bridge must prove bounds for these arbitrary remainder fibres and
ratios, or supply a justified domination theorem. The restricted
top-coefficient family in (3) is only one special case. Earlier
[monomial extremality counterexamples](prize-reduction-audit.md)
already exclude one unrestricted shortcut to that domination.

Thus even granting square-root cancellation in every relevant linear
prime-field period does not, through the established argument chain,
supply either uniform maximum in (9). This is a finding about the
missing proof, not a theorem that no logical implication can exist.

## 5. Verification and limits

The [companion verifier](../experiments/parallel23_prize_bridge_2026_09_05.py)
uses exact polynomial arithmetic and subset enumeration to check
the coefficient/power-sum correspondence and restricted list formula,
including extension-valued centers in F_25. It compares the general
remainder description against direct interpolation and codeword
enumeration in small prime fields, checks the false U=1 substitution,
constructs a small actual polynomial-phase box collision, and verifies
the official integer parameters and budget. The [results](../results/parallel23_prize_bridge_2026_09_05.json)
pin the files used. No production-size list or MCA maximum is computed.

The official archived source bytes and their existing manifest are
checked locally; no claim about current prize status is made and no
new source version is substituted. No submission, old-file edit, Lean
build, human referee review or formal proof verification is claimed.
The full Paley, subgroup and official-prize goals remain unproved.
