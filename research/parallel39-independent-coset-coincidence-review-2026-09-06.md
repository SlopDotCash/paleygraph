# Pass 39: independent review of the edge-coset coincidence bound

**Verdict:** the proposed coincidence estimate is correct. For a symmetric multiplicative subgroup `H=-H` of order `n`, the part of the nondegenerate triangle sum in which at least two edge cosets coincide is at most `3 F3* F4*/n`. The complementary region with three distinct edge cosets remains unestimated. In particular, distinctness does not itself restore the missing injectivity hypothesis from the published incidence lemma.

No central source file, accepted frontier bound, or Lean proof was changed by this review.

## Exact identity and Holder bound

Let `f(x)=(H circ H)(x)` and `Fk*=sum_(x!=0) f(x)^k`. For `G=F_p^*/H`, write `f(A)` for the common value on the coset `A`. Consider

\[
W=\sum_{\substack{\alpha,\beta\ne0\\\alpha\ne\beta}}
f(\alpha)^2f(\beta)^2f(\alpha-\beta)^2.
\]

Let `P` denote the sub-sum with `beta/alpha in H`. Setting `beta=h alpha` gives exactly

\[
P=\sum_{h\in H\setminus\{1\}}\sum_{\alpha\ne0}
f(\alpha)^4f((1-h)\alpha)^2. \tag{1}
\]

For every quotient coset `D`,

\[
|\{h\in H\setminus\{1\}:1-h\in D\}|=f(D). \tag{2}
\]

To verify (2), fix `d in D`. A solution `a-b=d`, with `a,b in H`, determines `h=b/a`, and `1-h=d/a in D`. Conversely, given such an `h`, the unique pair is

\[
a=d/(1-h),\qquad b=hd/(1-h),
\]

which belongs to `H x H`. This is a bijection, not an averaging identity.

There are `n` choices of `alpha` in each coset `A`, so (1) and (2) give

\[
P=n\sum_{A,D\in G}f(A)^4f(D)f(DA)^2. \tag{3}
\]

For each fixed `A`, Holder's inequality with exponents `3` and `3/2`, followed by the bijection `D -> DA` of `G`, yields

\[
\sum_D f(D)f(DA)^2
\le\left(\sum_Df(D)^3\right)^{1/3}
\left(\sum_Df(DA)^3\right)^{2/3}
=\sum_Df(D)^3.
\]

Since `Fk*=n sum_A f(A)^k`, this proves

\[
P\le n\left(\frac{F_4^*}{n}\right)
\left(\frac{F_3^*}{n}\right)=\frac{F_3^*F_4^*}{n}. \tag{4}
\]

All powers and factors of `n` in the proposed identity are therefore correct.

## The three coincidence regions

The three nonzero edges of the triangle with vertices `0,alpha,beta` have cosets represented by `alpha`, `beta`, and `alpha-beta`. Permuting its vertices gives bijections of the nondegenerate summation domain and permutes these edges up to sign. The weights are preserved because `f(-x)=f(x)`; the cosets are preserved under a change of sign because `-1 in H`.

It follows that each of the three pairwise-coset-equality regions has total weight `P`. Their intersections need not be empty, but the union bound is sufficient because the summands are nonnegative. Writing `Wcoin` for the weight of their union,

\[
W_{\mathrm{coin}}\le3P\le\frac{3F_3^*F_4^*}{n}. \tag{5}
\]

The symmetry hypothesis is used explicitly here. No disjointness of the three regions is claimed.

## Range of the moment input and energy implication

The published Shkredov paper, `sources/parallel38-energy-source/shkredov-2013-published.pdf`, PDF page 10, printed page 197, Corollary 2, states `E3(H) << n^3 log n` and `E_l(H)=n^l+O(n^((2l+3)/3))` for `l>=4`, under `n << p^(2/3)`. Its notation `E_l` is the moment of difference multiplicities, not the `2l`-variable additive energy `T_l`.

The `x=0` contribution is `n^l`. Consequently, in that source range,

\[
F_3^*\ll n^3\log n,\qquad F_4^*\ll n^{11/3},
\quad\text{and hence}\quad
W_{\mathrm{coin}}\ll n^{17/3}\log n. \tag{6}
\]

The quartic regime satisfies the range requirement. The estimate also bounds any cutoff-restricted sub-sum of the coincidence region used in the old energy argument.

For clarity, the Cauchy-Schwarz denominator in that old argument uses the **full** third moment `F3=n^3+F3*`. With this distinction preserved, the already derived lower estimate takes the form

\[
\frac{E(H)^6}{n^6F_3}\ll W_{\mathrm{retained}}.
\]

Apart from the separately handled degenerate-edge case, (6) therefore reduces the remaining task to the part of `Wretained` with three pairwise distinct edge cosets. If the coincidence contribution alone accounts for a fixed positive fraction of this lower bound, then

\[
E(H)^6\ll n^5F_3F_3^*F_4^*
\ll n^{44/3}(\log n)^2,
\]

which gives `E(H) << n^(22/9) log^(1/3)n` in that branch. This is a conditional branch conclusion; there is no estimate here for the distinct-coset branch.

## Distinct edge cosets do not remove normalization multiplicities

The generic injectivity obstruction from pass 38 persists even after excluding all edge-coset equalities. Let `A=B=C=D` be a subgroup of the quotient `G` of cardinality at least three. For normalized pairs `(u,v)` with `1,u,v` pairwise distinct, the multiplicity

\[
m(u,v)=|\{c\in D:cu,cv\in D\}|
\]

is still `|D|` whenever `u,v in D`. Distinctness of `1,u,v` does not change this fact.

There is a concrete positive-incidence example. Take `p=97`, `H=<5^12>` of order `8`, and `S=<5^4>` of order `24`. The three `H` cosets in `S` have least representatives `1,4,16`. The only normalized coset pairs with positive counts in the equation `x-y=1` are `(4,16)` and `(16,4)`, each with count `1`. Both have `1,u,v` pairwise distinct, and each has quotient multiplicity `3`. Thus all the direct incidence count `48` lies in this distinct-coset region, while discarding multiplicity gives only `16`.

This example concerns general invariant sets `S`; it is not asserted to be a dyadic level set of `f_H`. It demonstrates that distinctness **alone** is insufficient to recover the published incidence hypothesis. A proof tailored to the actual difference-multiplicity levels could still use additional structure, which has not been established here.

## Independent finite checks

Three exact integer enumerations checked (2), the equality between direct `P` and (3), (4), and `Wcoin<=3P`:

| p | n | P | F3* F4*/n | Wcoin | W with distinct edge cosets |
|---:|---:|---:|---:|---:|---:|
| 97 | 8 | 2880 | 9800 | 8640 | 8448 |
| 1153 | 8 | 2304 | 9800 | 6912 | 7680 |
| 257 | 16 | 479408 | 789168 | 950448 | 689664 |

The groups were checked to have the indicated orders and contain `-1`. These computations support the implementation of the identities; the proof above establishes the universal coincidence estimate. No off-region or full-conjecture claim follows from the finite checks.
