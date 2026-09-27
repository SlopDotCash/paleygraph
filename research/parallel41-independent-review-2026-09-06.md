# Pass 41 independent review and final weighted-triangle bound

The full Paley conjecture and Proximity Prize remain unproved. Root has
checked the ordinary arguments below and independently reproduced their
finite certificates. This is review within the task, not external peer
review or Lean certification. No new Lean build or submission was needed.

## 1. Final uniform triangle theorem

Let p be an odd prime, H=-H a multiplicative subgroup of size n with
n^2<p, and L=1+log n. Define

    f(x)=r_(H-H)(x), E_*=sum_(x!=0) f(x)^2,
    W=sum_(x,y!=0,x!=y) f(x)^2 f(y)^2 f(x-y)^2.

The final bound from this pass is

    W << n^(86/15) L^(8/15).                                    (1)

The power improves our previous generic n^6 estimate. It is still above
17/3 by 1/15. This is not a claim of a new best literature bound. The
[initial independent audit](parallel41-level-collision-bound-2026-09-06.md)
obtained L^(53/15); the geometric summation below removes L^3.

Use the quotient group G=F_p^*/H, and write

    a(C)=f(C)=|(H-1) intersect C|, b=a^2,
    T(u,v)=sum_C b(C)b(uC)b(vC),
    S(q)=sum_C a(C)a(qC), K=sum_q S(q)^2,
    B=E^times((H-1) minus {0}),
    X=B-[2(n-1)^2-(n-1)].

Field-ratio aggregation followed by Cauchy gives K<=nB. The exact
normalization, retaining every multiplicity, is

    W=n sum_(u,v) rho(u,v)T(u,v),
    rho(u,v)=#{x in uH,y in vH:y-x=1},
    X=sum_(u,v) (rho(u,v))_3.

The rho<=2 contribution is at most 2E_*^3/n^2. On rho>=3,
rho^3<=(9/2)(rho)_3. Thus Holder gives

    W_high <= (9/2)^(1/3) n X^(1/3) ||T||_(3/2).                (2)

All norms here use counting measure. Nonzero edges are required; including
zero edges or summing rho^3 without this restriction changes the estimate.

### A distribution-sensitive correlation lemma

Split positive a into levels D_i with T_i<=a<2T_i, where T_i=2^i.
Let d_i=|D_i| and Q=max_i T_i^3 d_i. Empty levels contribute nothing.
The following functional inequality is valid for nonnegative a on any
finite abelian group:

    ||T||_(3/2) << K^(1/5) Q^(26/15).                           (3)

For completeness, here is the full summation argument proposed by root
and checked independently by the correlation lane. Put

    t_i(q)=|D_i intersect qD_i|,
    m_ijk(u,v)=#{C in D_k:uC in D_i,vC in D_j},
    P_T=T_i T_j T_k, P_d=d_i d_j d_k.

Since S(q)>=T_i^2 t_i(q), sum_q t_i(q)^2<=K/T_i^4.
Also t_i(q)<=d_i and

    sum_(u,v) m_ijk(u,v)^2 = sum_q t_i(q)t_j(q)t_k(q).

Holder therefore bounds this last sum by K P_d^(1/3)P_T^(-4/3).
The exact l1 mass of m is P_d, and its maximum is at most min(d_i,d_j,d_k).
Interpolation gives two bounds for the same weighted block:

    ||P_T^2 m_ijk||_(3/2)
      <= min[ P_T^2 P_d^(2/3) min(d_i,d_j,d_k)^(1/3),
              K^(1/3) P_T^(14/9) P_d^(4/9) ].                 (4)

Order the three scales and write their maximum as z, with the other two
equal to z*2^(-r), z*2^(-s), r>=s>=0. From d_i<=Q/T_i^3,
the two arms in (4) are bounded respectively by

    Q^(7/3)/z,
    K^(1/3)Q^(4/3) z^(2/3) 2^(-2(r+s)/9).                    (5)

For positive A,D, summing over all dyadic z gives

    sum_z min(A/z,D z^(2/3)) << A^(2/5)D^(3/5),

by splitting at z=(A/D)^(3/5) and summing the two geometric series.
For each r,s, (5) thus sums to at most

    C K^(1/5)Q^(26/15) 2^(-2(r+s)/15).

The sum over r>=s>=0 is another convergent geometric series. At most six
orders of the scales and the factor 64 between actual block weights and
P_T^2 m_ijk affect only the absolute constant. This proves (3), with no
factor depending on the number of levels.

The actual subgroup tail a_j<<n^(2/3)j^(-1/3) gives Q<<n^2. Equations
(2) and (3) yield the more informative bound

    W <= 2E_*^3/n^2 + C n^(67/15) X^(1/3) K^(1/5)
      <= 2E_*^3/n^2 + C n^(14/3) X^(1/3) B^(1/5).             (6)

The shifted-energy theorem gives B<<n^2 L and X<=B. Its insertion in
(6) gives the high contribution n^(86/15)L^(8/15). The elementary scoped
energy bound E(H)<<n^(5/2) already bounds the low contribution by n^(11/2),
which is smaller. Thus (1) does not depend on the stronger MRSS energy
input; MRSS is used only when comparing outcomes.

Closing the earlier retained-triangle energy inequality with (1) gives
only E(H)<<n^(221/90)L^(23/90). Its power exceeds the imported 49/20 by
1/180. Accordingly, the energy input 49/20 and period input 71/72 are
unchanged. The repeated-six-word and full conjecture targets are also
unchanged. A smaller residual gap is not a measure of distance to a full
Paley proof.

### Removing diagonal incidence from the next target

The [actual-subgroup refinement](parallel41-level-refinement-2026-09-06.md)
records the diagonal formulas and the separate hereditary estimate.

The exact symmetries rho(u,v)=rho(v,u)=rho(u^(-1),v/u), together with
rho(1,C)=a(C)=rho(C^(-1),C^(-1)), give

    X_dist = sum_(1,u,v pairwise distinct) (rho(u,v))_3
           = X-3 sum_C(a(C))_3+2(a(1))_3 >= 0.                (7)

For W with three distinct edge cosets, X_dist can replace X in (2).
The previously checked coincident-coset contribution is <<n^(17/3)L.
Hence X_dist<<n^(9/5) up to logarithms would reach the target triangle
power via (6). Equation (7) is proved; that proposed upper bound is not.
The S3 symmetries preserve T as well as rho and give no sign cancellation.

A stronger consequence for the **total** excess uses its dependence on a.
The row rho(1,C)=a(C) is part of X, so sum_C(a(C))_3<=X. Since a is
integer valued and sum_C a(C)=n-1,

    sum_C a(C)^3 <= 4(n-1)+(9/2)X,
    Q << min(n^2,n+X).                                         (8)

Use (8) in (2)-(3) before replacing Q by its generic tail bound. This gives

    W_high << n X^(1/3)(nB)^(1/5) min(n^2,n+X)^(26/15).        (9)

Both algebra lanes independently checked this root deduction. With the
exact B=2(n-1)^2-(n-1)+X, a bound X<<n^(2-delta), 0<=delta<=1,
would make the exponent in (9) at most 86/15-31delta/15. Therefore

    X << n^(61/31) up to logarithms                            (10)

would reach the triangle power 17/3. A fixed power saving dominates any
fixed logarithmic factor in keeping B=O(n^2). The low contribution is
already smaller. Equation (10) is a sufficient **unproved input**, not a
consequence of the current published estimate X<<n^2 L. Under such a
subquadratic input, (8) also excludes the formerly extremal level
a~n^(3/5), d~n^(1/5), because that level alone has cubic mass of order
n^2. Under the currently available quadratic input it remains a possible
scale for these inequalities. This narrows the next question without
claiming an actual subgroup extremizer or a full-conjecture reduction.

The conditional criterion can in fact be proved without importing either
an energy bound or a weak-tail bound. Cauchy gives A2^2<=A1*A3, and
A1=n-1. Using (8) without its n^2 alternative gives the explicit bound

    W << n^(5/2)(n+X)^(3/2)
         + n^(6/5)X^(1/3)(n^2+X)^(1/5)(n+X)^(26/15).        (11)

For X<<n^(2-delta), its two powers are 11/2-3delta/2 and
86/15-31delta/15. The same delta=1/31 suffices. This formulation makes
clear that the conditional reduction requires the stated collision input;
it does not hide the needed improvement in a separate energy premise.

Root also checked the refinement's hereditary calculation: split each
level's multiplicative energy into 2|R_i|^2-|R_i| and its nontrivial
quadruples, which are a subset of those counted by X. Its mixed-level
interpolation replaces the blanket B by X+E_*^4/n^8, at the cost of
three logarithms. The diagonal criterion above gives a better sufficient
total-X power. The finite example X_dist=504>216=X-X_dist was independently
recomputed from shifted ratios at p=67403009,n=128; it refutes that
specific diagonal-domination inequality, not every constant-factor version.

## 2. Shell recovery review

The [real-subfield argument](parallel41-shell-arithmetic-recovery-2026-09-06.md)
is valid under its stated dyadic conductor and prime-splitting hypotheses.
If h is real, |Norm_+(h)|=lambda*p with p not dividing lambda, and h
vanishes at u, then it vanishes at exactly u and u^(-1) modulo p. The
real adjugate makes c=lambda*p/h integral. Its norm valuation shows that
c has nonzero evaluations at those two roots and vanishes elsewhere.
These facts justify division by p in f=h*bar(F)/p, and fc=lambda*bar(F)
then proves f(u)=0. Merely observing h(u)=0 before division would not
have been enough. The norm is exactly p*lambda^2*tau.

For lambda=1 the proposed maps are mutually inverse on the full integral
lattices, preserve the normalized norms, and identify the indicated
quotients. The restriction to centered actual coset representatives does
not preserve norm minimization. Root independently checked this with
32-by-32 integer determinants: the recovered cofactor is 1217; centering
the reverse image of the known cofactor-641 element changes the cofactor
to 178771841. No minimum-cofactor or uniform-real-generator claim follows.

The [unit-distortion follow-up](parallel41-shell-unit-distortion-2026-09-06.md)
also checks out. A real integer element of odd prime norm has odd constant
coefficient. Its squared coefficient norm S is therefore odd; S=1 gives
a trivial unit, while S=3 gives one of +/-1 +/- (X^j+X^(-j)), whose
cyclotomic norm is one by the X->X^3 automorphism. Thus S>=5. Fourier
diagonalization gives condition number >=sqrt(5)p^(-2/N), after every
real unit choice. In the quartic window this rules out a full-space
relative 1+o(1) isometry, and the exact n>=1024 threshold forces relative
error greater than 1/3. It does not rule out useful constant distortion,
estimates on actual cosets, or vector-dependent units. Root separately
recomputed the three specified unit multiples by integer determinants.

## 3. Classical smoothing review

The [smoothing argument](parallel41-classical-smoothing-2026-09-06.md)
correctly distinguishes deliberately independent shifts from the original
character summands. The partitions of six indices without singletons are
6, 4+2, 3+3, and 2+2+2; retaining all of them gives the stated exact
formula, including the nonzero third-moment-square term. Under n^4<=p,
t=ceil(n^(1/5)) gives expected smoothed sixth moment <=(325/32)pn^3.

This is a bound on a weighted convolution. It is not an original Sidon
moment estimate. The L2 error has expected square (1+1/t)(pn-n^2), and
its paired expectation is the original discrepancy. Jensen further shows
that the expected sixth moment of the error is at least the original
sixth moment. Thus the proposed unsmoothing shortcuts still require the
missing cancellation; the positive smoothed estimate alone cannot enter
the full Paley reduction.

## 4. Primary inputs and independent finite evidence

Root visually checked Theorem 4, d=2, on PDF page 10 (printed 197) of
the archived published Shkredov 2013 paper. Its n<sqrt(p) hypothesis
matches this pass. The adjacent three-set incidence lemma contains the
injectivity hypothesis identified in pass38; it is not applied here.
Root also checked [Shkredov Theorem 6](https://arxiv.org/html/1504.04522v1#Thmsatz6)
for the shifted multiplicative energy, [MRSS Corollary 12](https://arxiv.org/pdf/1712.00410)
for the comparison exponent, and [Volostnov Theorem 5](https://arxiv.org/pdf/1712.09355v1)
for the individual Weil character estimate. All uses retain their prime-field
and size hypotheses. The literature supplies the inputs, not the new
unproved X_dist or unsmoothed-moment estimates.

The [independent verifier](../experiments/parallel41_independent_verify.py)
checks exact full-field W normalization, shifted-energy excess, diagonal
subtraction, coset symmetries, quotient energy, and mixed-level collision
identities in six small fields, including two cases with positive high-rho
contribution. It also recomputes three larger, fixed subgroup excesses
using O(n^2) shifted-ratio counts without scanning their fields. It separately
checks the shell determinants and four further
lattice vectors, and enumerates 4,167 shift tuples in three smoothing cases,
including signed nonzero third moments. The [results](../results/parallel41_independent_verification_2026_09_06.json)
record success. The agent's two original verifiers were also rerun.
Finite checks support these identities and implementations; the displayed
ordinary proofs and their source hypotheses support the uniform conclusions.
