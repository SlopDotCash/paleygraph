# An actual sixth-moment upper bound for almost all Sidon inputs

**Status: an almost-all theorem for actual prime-field character kernels
is proved below. No uniform bound for the exceptional Sidon sets is
proved.** This extends the local typical-set work from the fourth to the
sixth moment and conditions the result on the Sidon property at the
requested size. No literature novelty or formal proof certification is
claimed. The uniform Paley, subgroup, spectral, and prize targets remain
open.

The [quartic-star interface](parallel21-positive-upper-review-2026-09-05.md)
and its [independent review](parallel22-quartic-star-independent-review-2026-09-05.md)
identify the missing uniform upper bound. The result here bounds the
actual quantity on a quantified family; it is not another equivalent
formulation of the uniform hypothesis or a relaxed sign-model example.
It does not remove the exceptional family from that hypothesis.

## 1. The restricted theorem

Let p be an odd prime, χ its quadratic character with χ(0)=0, and
C⊂F_p have size n. Write

    F_C(x)=Σ_(c∈C)χ(x−c),     M_6(C)=Σ_x F_C(x)^6.

Choose C uniformly among n-element subsets. For every n≥6 with n⁴≤p,
define U_s(j)=E(ξ_1+…+ξ_j)^(2s) for independent uniform signs ξ_i,
where j=n−1, and put

    A=36U_5(j)+240U_4(j)+472U_3(j)+240U_2(j)+36U_1(j),
    B=U_4(j)+2U_3(j)+U_2(j),
    V(p,n)=n(√(pA)+15√B)²,
    G(p,n)=p(30n²−16n).                                (1)

These are explicit quantities; U_s(j) can be evaluated by the finite
binomial distribution. A proved finite bound is

    P[M_6(C)>15pn³ | C is Sidon]
        ≤ min{1, (8/7)V(p,n)/G(p,n)²}.                 (2)

Here Sidon means that equal two-term sums, including repetitions, have
identical summand multisets. In particular, on the exact slice
n=⌊p^(1/4)⌋,

    P[M_6(C)>15pn³ | C is Sidon]
        ≤ (216/5+o(1)) n²/p
         = (216/5+o(1))p^(−1/2).                       (3)

The estimate holds for every sufficiently large prime, without averaging
over primes. The constants do not depend on the chosen prime congruence
class. For every member satisfying M_6(C)≤15pn³, its actual off-set
quartic-star energy obeys

    L_out(C)≤(2/3)p²n³.                                (4)

The assertion is about a proportion tending to one of Sidon sets.
Neither (2) nor (3) asserts that all Sidon sets satisfy this bound.

## 2. Moment domination for a sampled character row

A complete row χ(x−c), as c varies, contains one zero and (p−1)/2
copies of each sign. A random n-subset samples this population without
replacement. Its even moments are bounded by those of n independent
signs:

    E_C F_C(x)^(2s)≤U_s(n),
    E_C M_(2s)(C)≤pU_s(n).                              (5)

Here is a direct proof that avoids importing a sampling comparison
theorem. Randomly place the zero among p positions. Uniformly pair the
remaining p−1 positions and independently assign one + and one − to
each pair, with either orientation equally likely. This produces a
uniform permutation of the original balanced population: every balanced
sign arrangement is compatible with exactly ((p−1)/2)! matchings, so
all arrangements have equal probability.

Restrict this random population to any fixed n positions. A pair with
both endpoints selected contributes zero; a pair with exactly one
selected endpoint contributes an independent uniform sign. Conditional
on the pairing and zero position, the sum is therefore a sum of ℓ≤n
independent signs. Adding independent mean-zero signs cannot decrease
the expected convex function t↦t^(2s), so U_s(ℓ)≤U_s(n). This proves
(5), including the exceptional zero column.

The exact Rademacher sixth moment is

    U_3(n)=15n³−30n²+16n.

Consequently, if μ=E_C M_6(C), then

    15pn³−μ≥G(p,n)>0.                                  (6)

We also use U_s(j)≤(2s−1)!!j^s. Expanding the moment and pairing the
2s positions proves this: only words with even multiplicities survive,
and each is counted at least once by the pairings. For fixed s,

    U_s(j)=(2s−1)!!j^s+O_s(j^(s−1)).                   (7)

The leading words have s distinct indices, each occurring twice; every
other surviving word uses fewer indices. Thus (7) is a polynomial
count, not a limiting independence assumption about character rows.

## 3. A self-contained variance inequality on the subset slice

For any real function f of n-element subsets of a p-element set,
1≤n<p, let A be a uniform (n−1)-subset and b a uniform element of
its complement. Then

    Var_C f(C) ≤ [n(p−n+1)/p]
                       E_A Var_(b∉A) f(A∪{b}).        (8)

For completeness, define D_n(f) as half the expected squared difference
across a uniform oriented edge of the Johnson graph: the two sets differ
by one deletion and one insertion. With q=p−n, conditional sampling of
two elements outside A gives

    E_A Var_b f(A∪{b})=[q/(q+1)]D_n(f).                 (9)

Let g(A)=E_b f(A∪{b}). If A,A′ are adjacent (n−1)-sets, the common
term f(A∪A′) cancels in g(A)−g(A′). The remaining q terms have common
denominator q+1. Cauchy–Schwarz, followed by averaging the edge and its
extra element, gives

    D_(n−1)(g)≤[q/(q+1)]²D_n(f).                       (10)

The resulting n-set edge is uniform by permutation symmetry. Starting
from Var(f)=((p−1)/p)D_1(f) at n=1, total variance, (9), and induction
therefore give

    Var_n(f)
      ≤ [q/(q+1)+(n−1)(q+1)/p·q²/(q+1)²]D_n(f)
      = [nq/p]D_n(f).

Combining with (9) proves (8). The argument is for arbitrary slice
functions; it does not assume independence of the elements of C.

## 4. Full character structure bounds the sixth-moment variance

Fix A of size n−1, abbreviate F=F_A, and define the actual character
transform

    (Qh)(b)=Σ_x χ(x−b)h(x).

Its full-field identity is

    Σ_b(Qh)(b)²=pΣ_xh(x)²−(Σ_xh(x))²≤pΣ_xh(x)².        (11)

This follows by summing χ(x−b)χ(y−b) over b: the value is p−1
on the diagonal and −1 otherwise. It holds in both odd prime
congruence classes, with the orientation shown. This use of the actual
complete translate kernel is essential to the variance estimate.

Set

    H=6F⁵+20F³+6F,       Z=F⁴+F²,
    h_A(b)=(QH)(b)−15Z(b).

The insertion formula, with b∉A, is exactly

    M_6(A∪{b})=M_6(A)+15M_4(A)+15M_2(A)+p−1+h_A(b).    (12)

Indeed χ(x−b)^(2j)=1−1_(x=b) for j≥1. Expanding (F+χ(x−b))⁶
therefore produces the two subtractions −15F(b)⁴−15F(b)². These
exceptional-row terms must be kept; the odd powers give QH.

Put m=p−n+1. Variance over b outside A is unchanged by the constant
in (12) and is at most m^(-1)Σ_bh_A(b)². Inequality (8) consequently
yields

    Var_C M_6(C)≤(n/p)E_AΣ_bh_A(b)².                   (13)

The polynomial identities

    H²=36F^10+240F⁸+472F⁶+240F⁴+36F²,
    Z²=F⁸+2F⁶+F⁴

have positive coefficients. Applying (5) with |A|=j gives
E_AΣ_xH(x)²≤pA and E_AΣ_xZ(x)²≤pB, with A,B from (1).
By (11), E_AΣ_b(QH)(b)²≤p²A. The triangle inequality in the
Hilbert space of functions of (A,b) now gives

    [E_AΣ_bh_A(b)²]^(1/2)≤p√A+15√(pB).

Together with (13), this proves the explicit variance bound

    Var_C M_6(C)≤n(√(pA)+15√B)²=V(p,n).                (14)

No genus-two estimate or assumed near orthogonality of quartic traces
was used. The favorable negative rank-one term in (11) was simply
omitted for this upper bound.

Chebyshev and (6) prove the all-set probability bound V/G². Since
U_5(j)=945j⁵+O(j⁴) and U_4(j)=105j⁴+O(j³), as n→∞ with n⁴≤p,

    A=34020n⁵(1+O(1/n)),       B=105n⁴(1+O(1/n)),
    V(p,n)=34020pn⁶(1+o(1)),
    V(p,n)/G(p,n)²=(189/5+o(1))n²/p.                  (15)

This is a variance upper bound, not an exact variance formula. The
constant is not claimed optimal; the earlier local fourth-moment
variance calculation does not supply the sixth-moment variance here.

## 5. Conditioning on Sidon sets and controlling the quartic-star energy

Every nontrivial two-sum relation has either three or four distinct
labels. The number of relations 2a=b+c with b,c distinct and unordered
is p(p−1)/2. The number of equalities between unordered, disjoint pairs,
where the two pairs are themselves unordered, is

    p(p−1)(p−3)/8.

To see the latter, choose ordered a≠b, then choose c different from
a,b and (a+b)/2; d=a+b−c is forced and all four labels are distinct.
Each relation has eight such orderings. A uniform n-set contains any
specified j distinct labels with probability (n)_j/(p)_j. Thus the
union bound gives

    P[C is not Sidon]
      ≤ (n)_3/[2(p−2)]+(n)_4/[8(p−2)]
      = n(n−1)(n−2)(n+1)/[8(p−2)]≤1/8,               (16)

where the final inequality uses n≥3 and n⁴≤p. This counting does not
need a Poisson approximation or independence between relation events.
Dividing the unconditional exceptional probability by the Sidon
probability, at least 7/8, proves (2). Equations (15) and (16) give
(3), since (8/7)(189/5)=216/5.

Finally, let e_3(x) be the cubic elementary symmetric sum of the
character row on C. If N(x) is its number of nonzero entries, N is n
or n−1, and the exact identity is

    6e_3(x)=F_C(x)³−(3N(x)−2)F_C(x).

For n≥2, the coefficient 3N−2 is nonnegative. Squaring and dropping
the nonpositive fourth-power term therefore gives

    Σ_xe_3(x)²≤M_6(C)/36+(n²/4)M_2(C).

The exact second moment is M_2(C)=pn−n². If M_6(C)≤15pn³,
the displayed bound is at most (2/3)pn³. The full transform (11)
gives

    L_out(C)≤Σ_t U_C(t)²
      =pΣ_xe_3(x)²−T_3(C)²≤(2/3)p²n³,

which proves (4) with all zero and on-set rows accounted for.

## 6. Exact checks and what remains unresolved

The [verifier](../experiments/parallel23_classical_upper_2026_09_05.py)
uses exact integers and rational comparisons. Its exhaustive small-field
part checks moment domination, the insertion formula, the full transform,
the slice variance inequality for both M_6 and an unrelated function,
the radical variance certificate, and the counts of relation types.
The actual larger-field part samples genuine Sidon sets on the requested
slice and computes their character moments, together with a full integer
quartic-star transform for a selected set at each prime. Larger analytic
probability certificates use actual primes but do not evaluate those
large character kernels. These scopes are separate in the
[results](../results/parallel23_classical_upper_2026_09_05.json).

The script uses the integer upper enclosure

    V≤n[pA+225B+30 ceil(√(pAB))]

for its probability certificates. Its direct comparisons with the
unrounded radical variance bound square only after checking that the
other side is positive. No floating-point acceptance threshold is used.
The scientific claim rests on the all-parameter proofs above; the sampled
Sidon success counts are not statistical evidence for a worst-case claim.

Compared with [pass 6](parallel6-classical-2026-09-04.md), this gives a
sixth-moment coefficient 15 and a quantitative upper bound for actual
quartic-star energy after conditioning on Sidon sets. Elementary random
set concentration already gives other typical-set cancellation results;
no claim that the general random-set method is new is made.

The obstruction to a uniform conclusion is quantitative and explicit:
(3) leaves up to O(p^(−1/2)) of the Sidon slice uncontrolled. That is a
vanishing proportion but still allows an enormous number of individual
sets. The [pass-20 transfer](parallel20-inversion-moments-2026-09-05.md)
requires the input bound for **every** Sidon set of the prescribed size,
and begins, under its convenient threshold, at n≥25. Neither the
proportion estimate nor the five small-field samples satisfy that
universal quantifier. They also do not show that an arbitrary larger
set contains enough good subsets for the sampling transfer.

The uniform positive T₆/quartic-star estimate remains open, as does any
new worst-case Paley exponent from this lane. The result is restricted
to the stated almost-all M₆ theorem; no broader higher-order or B_h
extension is asserted in this pass.
