# Repeated-incidence row capacity and the surviving concentration pattern

The full Paley conjecture and Proximity Prize remain unproved. This pass
derives a necessary inequality from actual incidence rows. It excludes
the earlier interval model with Sidon filler under the available cubic
incidence bound. A different filler preserves the correlation obstruction
and satisfies its total-X consequence, so no uniform exponent improves.

## 1. A necessary inequality for every subset of rows

Let C be a symmetric matrix of nonnegative integers with row sums at most
n. Define

    h_i=sum_j C_ij(C_ij-1),
    X_D=sum_(i,j in D) C_ij(C_ij-1)(C_ij-2),
    X=sum_(i,j) C_ij(C_ij-1)(C_ij-2).

For a nonempty set D of d rows, put

    H_D=sum_(i in D) h_i,
    R_D=max_(j outside D) h_j,

where the maximum is0 if the complement is empty. Then

    H_D <= nd max(1,sqrt(R_D)) + 6^(1/3)d^(2/3)X_D^(2/3).  (1)

Consequently

    X >= X_D >= [H_D-nd max(1,sqrt(R_D))]_+^(3/2)
                      / (sqrt(6)d).                        (2)

To prove (1), split the ordered entries with row in D into entries whose
column is outside D, internal entries of size at most2, and internal
entries of size at least3. Symmetry gives

    C_ij(C_ij-1)<=h_j.

Thus an outside entry c contributes c(c-1)<=c sqrt(R_D). An internal
entry c<=2 contributes at most c. These two groups together contribute
at most nd max(1,sqrt(R_D)), using their shared row-mass budget.
For c>=3,

    [c(c-1)]^3 <= 6[c(c-1)(c-2)]^2,

because 6(c-2)^2-c(c-1)=(c-3)(5c-8)>=0. Holder over at most d^2
internal entries bounds their contribution by
6^(1/3)d^(2/3)X_D^(2/3), proving (1). No matrix realization or asymptotic
incidence theorem is assumed in this lemma.

## 2. The row data are determined by the actual difference function

For the subgroup incidence matrix, let a_i=C_0i and
r_i=sum_j a_j a_(j+i). The previously derived square identity and exact
row sums give

    h_i=r_i-(n-1)1_(i=0),
    sum_i h_i=(n-1)(n-2).                                  (3)

The total X is precisely the nontrivial shifted-product excess. The
[scoped source input](parallel41-level-collision-bound-2026-09-06.md)
gives X<<n^2(1+log n) when n^2<p. Equations (1)-(3) therefore impose a
restriction on the quotient autocorrelation that an arbitrary function
need not satisfy. This uses the diagonal consequences of the full
[weighted multiplication law](parallel45-weighted-matrix-algebra-2026-09-06.md).

## 3. Exclusion of the earlier Sidon-filled interval family

Use the exact family from
[pass44](parallel44-interval-and-matrix-review-2026-09-06.md):

    n=t^5, A=t^3, d=t, I={1,...,d},
    l=(n-2-Ad)/2, B=2l^2+1,
    S={Bj+j^2+10d+1:1<=j<=l},
    a=A1_I+2 1_S+1_{d+1}.

Here t>=2 is dyadic and the cyclic quotient has prime order m between
n^3 and2n^3, large enough to avoid wrap in all compared differences.
No claim is made that nm+1 is prime. Define h from (3), and take
D=I-I, of size2d-1. Directly summing the short differences gives

    H_D=A^2d^2+n-2-2A >= t^8.                              (4)

Outside D, the high interval's self-correlation vanishes. The nonzero
differences of S have multiplicity at most1, and its spacing ensures a
translate of I meets S in at most one point. Bounding the remaining
mixed and singleton terms gives

    R_D<=6A+8.                                             (5)

For all dyadic t>=1024, (4)-(5) imply

    H_D-n(2d-1)max(1,sqrt(R_D)) >= t^8/2.

For example ceil(sqrt(6A+8))<=4sqrt(A), so the subtracted term is at
most8t^(15/2)<=t^8/4 in this range. Applying (2), with2d-1<=2t, yields

    X >= t^11/(8sqrt(3)) = n^(11/5)/(8sqrt(3)).              (6)

This contradicts X=O(n^2 log n) as t grows. Thus, for all sufficiently
large t, the particular Sidon-filled family cannot extend to actual subgroup
incidence matrices obeying that uniform source bound. This does not refute its earlier role
as an abstract obstruction for the smaller list of scalar constraints.

The checker verifies (4)-(5) directly at t2 and t4. A closed-form integer
check at t4096 forces X>192n^2. It does not enumerate a group or construct
an actual subgroup at that enormous parameter. The general contradiction
with any fixed implied constant follows from (6), not this finite check.

## 4. A different filler survives every total-X row inequality

Keep n,A,d,I,l and the singleton as above, but replace S by the long
interval J={10n+1,...,10n+l}. Work in a prime cyclic group with
n^3<m<2n^3; all relevant sums and differences have no wrap ambiguity.
Then

    a=A1_I+2 1_J+1_{d+1}

still has mass n-1, avoids the identity, and has exactly one odd positive
value. Also A2<=3n^(7/5), and its dyadic weak cubic constant is Q=n^2.
The singleton and whole-group coset mass caps both hold with constant1.
Those are all subgroup-coset caps because m is prime.

The weighted quotient energy K remains O(n^3). Expanding the convolution
into nine component terms and applying Cauchy gives

    K <= 9[A^4d^3+16l^3+1+8A^2d^2l+2A^2d+8l]
      <=126n^3.                                            (7)

Here the mixed interval energy is at most d^2l. The high interval alone
still supplies K>=(2/3)n^3 and, for the same short shifts as in pass44,

    ||sum_x a(x)^2a(x+u)^2a(x+v)^2||_(3/2)
      >=128^(-1/3)n^(61/15).                               (8)

The new family satisfies the necessary total-X consequence

    H_D<=nd max(1,sqrt(R_D))+6^(1/3)d^(2/3)X_model^(2/3)    (9)

for every D with X_model=2n^2, for all t>=32. No local quantities X_D
are assigned to it. To see (9), h_max<=A2<=3n^(7/5). If the positive
part in (2) is nonzero, necessarily

    R_D < (h_max/n)^2 <=9n^(4/5)<n/2.

For all nonzero shifts |s|<=floor(l/2), the filler alone supplies
h_s>=2l>=n/2. Consequently D must contain at least l-1>=n/4 indices.
But H_D<=sum h<n^2, and then

    H_D^3 < n^6 <=6d^2(2n^2)^2.

This proves (9) in the remaining case. If the positive part is zero,
(9) holds immediately.

One can set X_model=192n^2 to satisfy also the previously compared
scalar conditions: divisibility by6, X_model>=3sum(a_i)_3, and
K<=n[2(n-1)^2-(n-1)+X_model], using (7). This is only a formal budget,
not the excess of a supplied incidence matrix. At t2 and t4, an exact
threshold check verifies every subset inequality (9) with this budget.
At fixed R_D, the strongest subset is {i:h_i>R_D}: adding indices with
h_i<=R_D cannot increase the positive part divided by d^(2/3), since
R_D<n^2 and h_i-n max(1,sqrt(R_D))<=0. For R_D>=n^2 all bounds are
trivial. A floor square root in the checker makes its test conservative.

No C or prime-field subgroup realizing the modified family is claimed.
The family shows that the total-X row-capacity conditions, together with
the compared scalar constraints, still does not improve the standalone
correlation power. The local X_D allocations and off-diagonal multiplication
identities remain unimposed and are the next constraints to investigate.

## 5. Verification scope

The [checker](../experiments/parallel46_row_capacity.py) verifies6,656
subset inequalities on five actual field matrices and the exact identity
(3). One positive lower certificate occurs in a dense example used only
to test the general matrix lemma; that example does not satisfy n^2<p
and is not used with the sparse source theorem. Four abstract examples,
two symbolic large-parameter checks, and399 conservative threshold tests
for the modified filler also pass. Results are in the
[certificate](../results/parallel46_row_capacity_2026_09_06.json).

Root derived and reviewed the ordinary proofs and exact checks. The
parallel agents remain explicitly errored at their usage limit; no
separate-agent, Lean, external review, or novelty claim is made. No
uniform triangle, energy, period, or full-goal estimate improves here.
