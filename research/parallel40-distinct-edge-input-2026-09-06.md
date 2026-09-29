# Actual difference levels do not restore distinct-edge injectivity

The proposed repair “exclude coincident edge cosets, then apply the old
unweighted normalization to the actual difference levels” fails in a
quartic-regime example. The failure occurs at a normalized pair with
positive additive incidence. This note also gives an exact identity
expressing the remaining normalization collisions through multiplicative
intersections of the actual levels. It does not establish an improved
uniform bound for `W_distinct` or disprove the old numerical energy bound.

The starting point is
[pass 39's edge reduction](parallel39-edge-coset-reduction-2026-09-06.md).
The published incidence hypothesis and the existing weighted layer repair
are recorded in
[the pass 38 audit](parallel38-independent-incidence-2026-09-06.md).
Neither the weighted layer repair nor generic Young is asserted as new here.

## Exact normalization after excluding coincident cosets

Let `H=-H` have order `n`, put `f(x)=r_(H-H)(x)`, and set

\[
 E=\sum_x f(x)^2,\qquad d=E/(16n^2),\qquad
 S_i=\{x\ne0:2^{i-1}d<f(x)\le2^i d\}.
\]

Write `D_i` for the set of `H` cosets contained in `S_i`, in
`G=F_p^*/H`. For a fixed level triple `(i,j,k)`, define

\[
 m(u,v)=\#\{(A,B,C)\in D_i\times D_j\times D_k:
 A,B,C\text{ pairwise distinct},\ A/C=u,\ B/C=v\},
\]

and

\[
 \rho(u,v)=\#\{(x,y)\in uH\times vH:y-x=1\}.
\]

Thus `m` is supported on pairs for which `1,u,v` are pairwise distinct.
The number of `y-x=z` solutions in `S_i x S_j x S_k` with distinct edge
cosets is exactly

\[
 I_{ijk}^{\mathrm{dist}}=n\sum_{u,v}m(u,v)\rho(u,v).                 \tag{1}
\]

There are `n` choices of `z` in each denominator coset. Setting
`alpha=y`, `beta=x` shows that this sign convention agrees with `W` in
pass 39. Its weighted contribution for this level triple is exactly

\[
 n\sum_{u,v}\rho(u,v)
 \sum_{\substack{(A,B,C)\text{ counted by }m(u,v)}}
 f(A)^2 f(B)^2 f(C)^2.                                             \tag{2}
\]

In particular, even the weakened injectivity candidate
`m(u,v)<=1 whenever rho(u,v)>0` would need a separate proof.

## A quartic-regime counterexample to that candidate

Take

\[
 p=215535361,\qquad H=\langle25525303\rangle,\qquad n=128.
\]

The exact checker verifies primality, order `128`, and `-1 in H`.
Here `n^4/2 < p < n^4`. The nonzero difference multiplicities on the
quotient have histogram

| Value of `f` | Number of `H` cosets |
|---:|---:|
| 1 | 1 |
| 2 | 48 |
| 4 | 3 |
| 6 | 3 |

Consequently `E=61056`, `d=477/2048`, and the old levels are:
`S3` has the one coset of value `1`; `S4` has the 48 cosets of value `2`;
`S5` has the six cosets of values `4` or `6`. In particular,
`|S4|=6144` and `|S5|=768`. Numerically, the two size expressions in the
original incidence lemma obey

\[
 |S_4||S_5|^2=3623878656=(27/256)n^5,
 \qquad n|S_4||S_5|^2<p^3.
\]

These comparisons do not specify the source's implicit constants. The
counterexample below concerns injectivity itself and does not assert a
violation of the published incidence conclusion under its full hypotheses.

For `(i,j,k)=(4,5,5)`, take `(u,v)=(674430H,2246831H)`. The two
preimage triples are

| `A,B,C` | `f(A),f(B),f(C)` |
|---|---|
| `817794H, H, 1492225H` | `2,6,6` |
| `674430H, 2246831H, H` | `2,4,6` |

Each triple has three distinct cosets. The first normalizes to `(u,v)`
because, with `g=25525303`, the following congruences hold modulo `p`:

\[
 \frac{817794}{1492225\cdot674430}=g^{44},\qquad
 \frac{1}{1492225\cdot2246831}=g^{89}.
\]

The second normalization is immediate. There are exactly three normalized
solutions, so `m(u,v)=2` and `rho(u,v)=3`:

\[
 (x,y)=(19584002,19584003),\quad
 (108122689,108122690),\quad
 (162237922,162237923).
\]

As short membership certificates, the pairs of exponents for
`(x/u,y/v)` in powers of `g` are `(53,117)`, `(126,17)`, and `(81,73)`.
Thus this is a collision that contributes to `W_distinct`; it is not
merely a repeated pair with zero additive incidence.

Over the whole level triple `(4,5,5)`, the domain of distinct coset
triples has size `48*6*5=1440`; its normalized image has size `1438`.
Exactly two image points have multiplicity `2`, both with `rho=3`.
Independent enumeration in the field gives

\[
 I_{455}^{\mathrm{dist}}=52736,
 \qquad n\sum_{\{m>0\}}\rho=51968.
\]

Discarding multiplicity loses `768`. The same field enumeration gives
the actual weighted contribution (2) as `144850944`, agreeing with
the weighted quotient calculation. The single displayed repeated fiber
contributes
`128*3*(2^2*6^2*6^2+2^2*4^2*6^2)=2875392` to (2).

There is no positive repeated fiber with all three levels equal in this
particular quartic example. Mixed level triples already occur in the old
argument, so this does not rescue it. For comparison, the bounded check
also finds a same-level failure at `p=193`, `H=<64>`, `n=16`: all seven
support cosets lie in `S4`, and `I444_dist=4992` versus `2304` after
discarding multiplicity. That second example is outside the quartic window.

## A reduction to multiplicative intersections of actual levels

For each distinct level index `l` appearing in `(i,j,k)`, let `nu_l`
be its number of appearances, and put

\[
 t_l(q)=|D_l\cap qD_l|,\qquad (t)_a=t(t-1)\cdots(t-a+1).
\]

There is an exact collision identity

\[
 \mathcal C_{ijk}:=\sum_{u,v}m(u,v)(m(u,v)-1)
 =\sum_{q\ne1}\prod_l (t_l(q))_{\nu_l}.                           \tag{3}
\]

To prove it, count ordered pairs of different distinct-coset triples
with equal normalized image. There is a unique common multiplier
`q!=1` taking the first triple to the second. For each level, its chosen
coordinates must be distinct elements whose multiples by `q` remain
in that level, giving the falling factorial. Different levels are
disjoint, so no further cross-level exclusions are needed. Using
`D_l cap q^{-1}D_l` in this counting gives the same cardinality as
`D_l cap qD_l`.

For example, `(3)` becomes
`sum_(q!=1) t_i(q)t_j(q)(t_j(q)-1)` for `(i,j,j)`, and
`sum_(q!=1) t_i(q)(t_i(q)-1)(t_i(q)-2)` for `(i,i,i)`.

Let `R_col=max_{m>=2} rho`, with value zero if there are no repeated
fibers. Equation (1) gives the valid upper estimate

\[
 I_{ijk}^{\mathrm{dist}}
 \le n\sum_{\{m>0\}}\rho+\frac{nR_{\mathrm{col}}}{2}\mathcal C_{ijk}, \tag{4}
\]

since `m-1 <= m(m-1)/2` for positive integer `m`. This separates the
distinct-pair term from an explicit collision term of the actual
difference levels. Applying any incidence theorem to the first term
still requires that theorem's own size hypotheses.

In the quartic `(4,5,5)` example, the only nonzero terms of (3) are
`q=1492225H` and `q=2246831H`. Each has `t4(q)=1`, `t5(q)=2`, so
`C455=4`. Here `R_col=3` and (4) is an equality: its correction is
`128*3*4/2=768`. The obstruction therefore has a concrete structural
description: the high level contains the three successive quotient
powers `q^{-1},1,q`, while the lower level contains a compatible
successive pair.

Controlling these intersection sums uniformly, with their additive
weights and across all retained levels, remains unresolved. This note
does not turn (4) into the desired `n^(17/3)` upper estimate for
`W_distinct` or an energy exponent improvement.

## Reproduction and limits

Run `python3 experiments/parallel40_distinct_levels.py`. The checker
uses exact integer arithmetic, constructs differences from `H x H`,
and examines only eight specified subgroups and their support cosets.
It does not scan all elements of the largest field or run Lean.
The complete results are in
[the exact output](../results/parallel40_distinct_levels_2026_09_06.json).
The first three examples' normalized `W_distinct` values also agree with
the independently enumerated pass 39 results. The two new positive
witnesses have direct field checks of their whole indicated level triple,
both unweighted and weighted, plus independent checks of (3).

The concrete failed candidate and identities (1)-(4) are the bounded
outcome. There is no claimed literature novelty, unrestricted incidence
counterexample, general energy improvement, or full-conjecture proof.
