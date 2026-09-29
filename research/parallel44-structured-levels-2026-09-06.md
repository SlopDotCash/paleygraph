# Actual subgroup mass on cosets and a structured triangle case

The full Paley conjecture and prize remain open. This note proves a
restriction on actual difference multiplicities inside a coset of any
larger multiplicative subgroup. It gives the desired triangle power when
the large-multiplicity levels admit a specified nested coset partition.
That structural hypothesis is not proved for general actual levels.
No novelty or full uniform improvement is claimed.

Write H=-H, n=|H|, G=F_p^*/H, and a(C)=r_(H-H)(C). Throughout the
quartic application p lies in [c n^4,C n^4], for fixed positive c,C and
sufficiently large n. Let E_*=sum_(x!=0) r_(H-H)(x)^2 and let W be
the same nonzero weighted-triangle sum as in passes41–43.

## 1. A stronger mass bound inside a subgroup coset

Let K<=G have order d, let S be its preimage in F_p^*, and let D=lambda K
be any coset. Then |S|=nd and the actual shifted-set identity gives

    sum_(C in D) a(C)=|(H-1) intersect lambda S|.

Apply the two-subgroup intersection estimate to H and S. It yields

    sum_(C in D) a(C) << (n|S|)^(1/3)=n^(2/3)d^(1/3).       (1)

In particular, if a(C)>=T on every C in D, then

    T^3 d^2 << n^2.                                         (2)

This is stronger than the general weak-cubic-tail bound T^3 d<<n^2.
It applies because D is a whole subgroup coset, not merely a set of d
cosets of H.

The precise source is Mit'kin's lemma as stated in
[Shkredov, arXiv:1504.04522v1, Lemma2, equation5](https://arxiv.org/html/1504.04522v1#Thmsatz2).
For a singleton pair of coset representatives it assumes
(|Gamma||Pi|)^2<p^3 and |Gamma||Pi|>=33^3 and bounds the intersection
by O((|Gamma||Pi|)^(1/3)). Root read the primary HTML statement and
its mathematical source text. For Gamma=H,Pi=S, the first condition
holds for every S<=F_p^* when n^2<p. If n^2d<33^3, the trivial bound
by n absorbs the finitely bounded n into an absolute constant. Thus (1)
is uniform in S. No comparability of the two subgroup sizes is assumed.

## 2. Three nested subgroup cosets

Suppose three pieces D_i=lambda_i K_i consist of quotient subgroup cosets
on which T_i<=a<2T_i. Suppose their underlying subgroups form a chain;
after permuting the three edge roles, write K_1<=K_2<=K_3 and
d_1<=d_2<=d_3. Such a permutation preserves the count since H=-H.
Let S_i be their preimages and A_i=lambda_i S_i.

Consider the retained incidence

    I=#{(z,x,y) in A_1 x A_2 x A_3:y-x=z}.

For a fixed z0 in A_1, multiplication by S_1 maps its normalized
solutions bijectively to those for every other z in A_1. Hence

    I=nd_1 * #{x in A_2:x+z0 in A_3}
      << n^(5/3)d_1(d_2d_3)^(1/3).                         (3)

The second step is the same two-subgroup lemma, now for S_2,S_3.
Because a>=T_i>=1 on each piece, d_i<=n-1. Thus its field-size condition
(n^2d_2d_3)^2<p^3 holds for sufficiently large n in the quartic range.
The 33^3 condition is again absorbed for bounded n.

The contribution of this piece triple to W is at most
64 T_1^2 T_2^2 T_3^2 I. Applying (2) separately gives

    W_(D1,D2,D3)
      << n^(17/3) d_1^(-1/3)d_2^(-1)d_3^(-1)
      <= O(n^(17/3)).                                      (4)

All field incidence multiplicities are retained. Neither normalization
injectivity nor a bound on the global shifted-product excess is assumed.

## 3. Only the large levels need this hypothesis

If any edge has f(x)<=T0, positivity and a bijective change between the
other two edge coordinates give

    W_(some edge small)<=3T0^2 E_*^2.                        (5)

Take T0=n^(23/60). The existing scoped bound
E(H)<<n^(49/20)(1+log n)^(1/5) makes (5) at most
n^(17/3)(1+log n)^(2/5).

Partition a>=T0 into dyadic levels. Assume each such level is a disjoint
union of at most r whole cosets of quotient subgroups, and that every
underlying subgroup in the partition belongs to a single chain. There
are O(log n) levels. Summing (4) over every ordered piece triple proves

    W << n^(17/3)[(1+log n)^(2/5)+r^3(1+log n)^3].          (6)

Thus a polylogarithmic value of r would close this intermediate triangle
power up to logarithms. This is a sufficient conditional statement. It is
not a claimed description of actual large-multiplicity levels. In
particular, a low-doubling set need not be a subgroup coset, and no such
replacement is made here.

## 4. What the restriction excludes, and what it leaves

The abstract extremizer a=A*1_K for the previous functional bound can
have A=n^(3/5), |K|=n^(1/5). It meets Q=n^2 and K_quot=n^3, but (2)
would require n^(11/5)<<n^2. Consequently it cannot model actual
subgroup-coset concentration in an unbounded family.

The [interval construction](parallel44-interval-and-matrix-review-2026-09-06.md)
shows why this observation alone does not improve the general functional
exponent: a short interval in a large prime-order quotient group retains
the same power while satisfying every subgroup-coset mass bound. It is
an abstract function, not actual subgroup difference data. Controlling
that more general concentration, or using additional arithmetic identities,
remains necessary. The nested-coset theorem does not prove the full Paley
or official prize conclusions even if its structural hypothesis is met.

The [checker](../experiments/parallel44_structured_levels.py) verifies 24
actual mass identities and three nested incidence normalizations at a
quartic prime. It records whether the source size hypotheses hold in each
mass example instead of treating an unspecified implied constant as 1.
The [certificate](../results/parallel44_structured_levels_2026_09_06.json)
is finite validation. Root completed the ordinary proof review; no
separate-agent, Lean, or external peer review is claimed.
