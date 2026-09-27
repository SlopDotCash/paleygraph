# Actual subgroup constraints on the weighted correlation bound

The strongest verified sufficient criterion in this lane is

\[
 \boxed{X(H)\ll n^{61/31}\ \text{up to logarithms}
 \quad\Longrightarrow\quad
 W\ll n^{17/3}\ \text{up to logarithms}.}                         \tag{1}
\]

Here `X` is the nontrivial multiplicative-energy excess of
`(H-1)\{0}`, defined below. The required bound on `X` is not proved.
The implication uses the actual relation between the diagonal cyclotomic
numbers and difference multiplicities. It does not require the imported
MRSS energy bound. This note also gives a sharper off-diagonal excess
parameter and a valid hereditary refinement of the intermediate level
energy estimate.

The root lane supplied the `X`-to-tail-constant observation during this
audit. Its normalization, constants, and exponents are independently
verified here. The logarithm-free correlation estimate is from
[the concurrent correlation optimization](parallel41-correlation-optimization-2026-09-06.md),
whose two block bounds and geometric summation were read and checked.
No new external theorem is invoked in the refinements below.

## Exact quantities and the additional diagonal constraint

Let `p` be odd, `H=-H` have order `n>=2`, and write `G=F_p^*/H`.
Set

\[
 R_0=(H-1)\setminus\{0\},\quad a(C)=r_{H-H}(C)=|R_0\cap C|,
 \quad A_j=\sum_C a(C)^j.
\]

Thus `A1=n-1` and `E_*=nA2`. Define

\[
 B=E^\times(R_0),\qquad
 X=B-\bigl(2(n-1)^2-(n-1)\bigr)\ge0,
\]

and

\[
 \rho(u,v)=\#\{x\in uH:1+x\in vH\}.
\]

The previously verified exact identities give

\[
 X=\sum_{u,v}(\rho(u,v))_3,
 \qquad \rho(1,C)=a(C)=\rho(C^{-1},C^{-1}).                       \tag{2}
\]

For every nonnegative integer `c`,
`c^3<=4c+(9/2)(c)_3`: check `c=0,1,2` directly and use
`c^3<=(9/2)(c)_3` for `c>=3`. Hence

\[
 A_3\le4(n-1)+\frac92\sum_C(a(C))_3
 \le4(n-1)+\frac92 X.                                           \tag{3}
\]

For actual dyadic levels `D_i={C:T_i<a(C)<=2T_i}`, let

\[
 Q=\max_i T_i^3|D_i|,
 \qquad K_{\rm quot}=\sum_q\left(\sum_Ca(C)a(qC)\right)^2.
\]

Every one of the level sums is bounded by `A3`, so

\[
 Q\le A_3\ll n+X,
 \qquad K_{\rm quot}\le nB\ll n(n^2+X).                          \tag{4}
\]

The second inequality is Cauchy on the `n` field ratios within each
quotient ratio. It retains the factor `n`; no quotient fibers are
discarded. The exact expression for `B` shows that the last bound in
(4) does not need a published estimate for shifted energy.

## The sufficient power saving in the excess

Write `b=a^2` and
`T(u,v)=sum_C b(C)b(uC)b(vC)`. The audited correlation bound is

\[
 \|T\|_{3/2}\ll K_{\rm quot}^{1/5}Q^{26/15},
\]

and the exact normalization gives

\[
 W=n\sum_{u,v}\rho(u,v)T(u,v)
 \le2nA_2^3+C nX^{1/3}K_{\rm quot}^{1/5}Q^{26/15}.              \tag{5}
\]

The high-incidence estimate uses
`sum_(rho>=3)rho^3<=(9/2)X`; it does not apply that estimate to the
many cells of size one or two.

Cauchy gives `A2^2<=A1 A3`. Equations (3)-(5) therefore prove the
fully quantitative, combinatorial bound

\[
 \boxed{W\ll n^{5/2}(n+X)^{3/2}
 +n^{6/5}X^{1/3}(n^2+X)^{1/5}(n+X)^{26/15}.}                    \tag{6}
\]

In particular, if `X<=n^(2-delta)` with `0<=delta<=1`, the two powers
in (6) are

\[
 \frac{11}{2}-\frac{3\delta}{2},
 \qquad \frac{86}{15}-\frac{31\delta}{15}.                       \tag{7}
\]

Taking `delta=1/31` makes the second power exactly `17/3`; the first
is `169/31<17/3`. This proves (1). Fixed powers of logarithms in the
hypothesis introduce only logarithmic factors in the conclusion.
The implication itself does not require a quartic relation between
`p` and `n`, or any imported additive-energy estimate.

The available scoped theorem gives only `X<<n^2(1+log n)` when
`n^2<p`. It supplies no fixed positive `delta` in (7). Subtracting the
quadratic trivial energy from this coarse bound does not supply the
missing saving. Thus (1) is a sharper sufficient condition, not a
completed estimate for the actual subgroup family.

## Removing the diagonal rich cells exactly

Let `kappa=a(1)`, and put

\[
 X_{\rm off}=\sum_{\substack{u,v:\1,u,v\ {m pairwise\ distinct}}}
 (\rho(u,v))_3.
\]

The three coincidence regions are `u=1`, `v=1`, and `u=v`.
Their pairwise intersections are all `(1,1)`. By (2) and
inclusion-exclusion,

\[
 \boxed{X_{\rm off}
 =X-3\sum_C(a(C))_3+2(\kappa)_3\ge0.}                            \tag{8}
\]

Consequently, the high-incidence estimate for `W_distinct` uses
`X_off^(1/3)` in place of `X^(1/3)`. Its exact low-incidence mass also
has the sharper expression

\[
 \sum_{1,u,v\ {m distinct}}T(u,v)
 =A_2^3-3A_2A_4+2A_6.                                           \tag{9}
\]

This is the ordered distinct-coordinate identity applied to `b=a^2`.
Thus a refined version of (5) is

\[
 W_{\rm distinct}\le
 2n(A_2^3-3A_2A_4+2A_6)
 +C nX_{\rm off}^{1/3}K_{\rm quot}^{1/5}Q^{26/15}.               \tag{10}
\]

The separate coincident-edge contribution remains controlled by the
pass 39 estimate. Formula (10) retains a potentially useful saving
when most rich cells lie in the three coincidence regions.

Symmetry alone gives no further saving in the pairing. Both `rho`
and `T` are invariant under the generators
`(u,v)->(v,u)` and `(u,v)->(u^(-1),v/u)`. For `T`, this follows by
permuting its three factors and changing `C` to `uC`. Averaging either
factor over this six-element symmetry group therefore leaves it
unchanged. This identifies precisely why simply symmetrizing the
weighted sum supplies no cancellation.

There is also a finite obstruction to the concrete candidate
`X_off<=X-X_off`: off-diagonal rich mass need not be bounded by the
diagonal rich mass. In the already archived quartic example

\[
 p=67403009,\quad H=\langle64701253\rangle,\quad n=128,
\]

the exact rich-cell records in `results/mixed_period_collisions.json`
give `X=720`, diagonal rich mass `216`, and `X_off=504`.
This is a counterexample to that stated inequality, not to every
comparison with an unspecified constant. Two additional bounded
checks constructed `H`, its shifted products, and its coset counts
directly: `(p,n)=(6700417,64)` gives `(X,diagonal,X_off)=(114,78,36)`;
`(215535361,128)` gives `(1740,1056,684)`.
No search for a zero-diagonal witness was continued.

## A hereditary refinement that preserves the actual excess

There is another valid refinement of the step `E^times(R_i)<=B` in
the original level proof. Let `R_i` be the part of `R0` in `D_i`,
`m_i=|R_i|`, and let `X_i` count its nontrivial product equalities.
Then

\[
 B_i=E^\times(R_i)=2m_i^2-m_i+X_i,
 \qquad0\le X_i\le X,
 \qquad\sum_iX_i\le X.                                         \tag{11}
\]

The last inequality holds because the `R_i` are disjoint and each
counted nontrivial quadruple lies entirely in one of them. This is
an exact hereditary statement; it does not assume a new theorem
bounding energies proportionally to subset density.

Here is one resulting bound, retaining a dyadic logarithm for clarity.
Let `V` be a weak-cubic constant such that
`|D_i|<=V/T_i^3`, and write `J=A2`. Under `n^2<p`, the scoped subgroup
tail permits `V=Cn^2`. Then

\[
 \boxed{\|T\|_{3/2}
 \ll n^{1/5}V^{26/15}
       (X+J^4/V^2)^{1/5}(1+\log n)^3.}                          \tag{12}
\]

For verification, classify a level as trivial type when
`X_i<=2m_i^2`, and as excess type otherwise. In the first case
`B_i<<T_i^2d_i^2`; in the second `B_i<=2X`. The ordinary and energy
arms for a mixed triple are, respectively,

\[
 \prod_l T_l^2d_l^{7/9},\qquad
 n^{1/3}\prod_l T_l^{14/9}d_l^{4/9}B_l^{1/9}.                  \tag{13}
\]

For a trivial-type level, `d_i<=min(J/T_i^2,V/T_i^3)` gives

\[
 T_i^2d_i^{7/9}\le V^{4/9}J^{1/3},\qquad
 T_i^{16/9}d_i^{2/3}\le V^{4/9}J^{2/9}.
\]

For an excess-type level, the corresponding bounds are
`V^(7/9)T_i^(-1/3)` and
`V^(4/9)T_i^(2/9)X^(1/9)`. Taking the ordinary arm to power `2/5`
and the energy arm to power `3/5` cancels these latter powers of `T_i`.
If there are `k` excess-type indices in the triple, the result is

\[
 n^{1/5}V^{4/3+2k/15}J^{4(3-k)/15}X^{k/15}.
\]

With `Y=X+J^4/V^2`, both `X<=Y` and `J^4<=V^2Y` hold, so this is
at most `n^(1/5)V^(26/15)Y^(1/5)` for every `k=0,1,2,3`.
Summing the bounded number of dyadic triples proves (12).

Taking `V=Cn^2` in (12) yields the additional estimate

\[
 W_{\rho\ge3}\ll
 n^{14/3}X^{1/3}(X+E_*^4/n^8)^{1/5}(1+\log n)^3.               \tag{14}
\]

For `W_distinct`, the first `X^(1/3)` may be replaced by
`X_off^(1/3)`. This sharpens the dependence on the actual excess in
the old level-energy arm. Its sufficient criterion using MRSS alone
would be `X<<n^(15/8)` up to logarithms, but the independently checked
diagonal-to-tail argument (1) is stronger: it needs only `X<<n^(61/31)`.

All refinements are ordinary proofs of identities, upper bounds, or
conditional implications. The existing uniform power gap is not closed.
Only this note was written; no Lean build, heavy subgroup scan, or
central-file change was performed.
