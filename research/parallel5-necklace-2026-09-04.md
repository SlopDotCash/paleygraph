# The full three-gap necklace family

**Status: a uniform growing-family bound is proved.** Let p≡1 mod 4
be prime, χ(0)=0, and A={0}, B={1}, C={0,1}. For a,b,c≥1 put

\[
 w_{a,b,c}=B A^{a-1}B A^{b-1}C A^{c-1},\qquad k=a+b+c.
\]

Then

\[
 \boxed{|N(w_{a,b,c})|\le
 (3a+1)\min(b,c)\,p^{(k+1)/2}+5p^{k/2}+2.}
 \tag{1}
\]

This includes every A^r B A^s B A^t C with r,s,t≥0, and all ranks
a,b,c, including ranks at least p. The proof retains all zero coordinates,
both other boundary fibers, and the constant term in the quartic kernel.
There is no exceptional equal-rank tensor left untreated.

Together with the preceding results, it proves the individual degree-two
necklace estimate for **all 243 words of length five**. At length six the
covered finite catalog grows to **639 of 729 words**; the 90 words with
multiplicities (2,2,2) remain outside these proofs. These are coverage
counts, not a percentage of a Paley conjecture. No cancellation estimate
for the full weighted necklace aggregate or spectral edge is proved here.

## 1. Full-field rows and multiplicative kernels are different at zero

Set S_xy=χ(x−y), D=diag(χ(x)), C_0=DS, and

\[
 K_j=S(DS)^{j-1},\qquad
 \kappa_j(x)=(K_j)_{1x}\quad(x\in F_p),\qquad
 k_j=\kappa_j|_{F_p^*},\qquad t_j=\kappa_j(1).
\]

The symbol κ_j denotes the **full-field** row. In particular,

\[
 \kappa_j(0)=(-1)^{j-1},\qquad
 (K_j)_{0y}=(-1)^{j-1}\chi(y),\qquad (K_j)_{00}=0.
 \tag{2}
\]

It must not be replaced by a convolution function extended by zero at
zero. The [earlier kernel proof](parallel3-necklace-2026-09-04.md) gives
the multiplicative convolution identity and Mellin eigenvalues
J(ψ,χ)^j for k_j on F_p*. It also gives

\[
 \begin{split}
 &K_j(x,y)=\chi(x)k_j(y/x)\quad(x,y\ne0),\\
 &k_j(x^{-1})=\chi(x)k_j(x),\qquad
 \sum_{x\ne0}k_j(x)=\sum_{x\ne0}\chi(x)k_j(x)=(-1)^j.
 \end{split}
 \tag{3}
\]

The elementary estimates needed below are

\[
 \|K_j\|\le p^{j/2},\quad
 \|\kappa_j\|_2\le p^{j/2},\quad |t_j|\le p^{j/2},
 \quad
 \sum_{x\ne0}|k_j(x)|^2
 =\frac{(p-3)p^j+2}{p-1}\le p^j.
 \tag{4}
\]

The matrix bounds follow from S²=pI−J and ||D||≤1. The last identity
is multiplicative Parseval: the two exceptional Jacobi eigenvalues have
absolute value one, and each of the remaining p−3 has absolute value √p.
No pointwise rank-dependent trace bound is used to sum the outside variable.

## 2. The exact identity, with every boundary term

Define the quartic anchor correlation

\[
 H(x,y,z)=\sum_v\chi(v)\chi(v-x)\chi(v-y)\chi(v-z).
\]

The [preceding pass](parallel4-necklace-2026-09-04.md) proves by weighted
anchor averaging and scaling the C coordinate to 1 that

\[
 N(w_{a,b,c})=\sum_{x,y\in F_p}
 (K_a)_{xy}\kappa_b(y)\kappa_c(x)H(x,y,1).
 \tag{5}
\]

For x,y≠0 it also gives

\[
 H(x,y,1)+1=\chi(xy)K_2(x^{-1}-1,y^{-1}-1).
 \tag{6}
\]

For x∈F_p\{0,1}, let

\[
 q_x(t)=\frac{1-xt}{t(1-x)},\qquad
 T_{a,b}(x)=
 \sum_{\substack{t\ne0\\t\ne x^{-1}}}
 \chi(t)k_a(t)k_b(xt)k_2(q_x(t)).
 \tag{7}
\]

The variable t=1 is included. The excluded t=x^(-1) corresponds to
y=1 and is retained separately below. Put

\[
 \begin{split}
 M^\circ_{a,b,c}&=\sum_{x\ne0,1}\chi(1-x)k_c(x)T_{a,b}(x),\\
 W_{r,s}&=\sum_{x\ne0}\chi(1-x)k_r(x)k_s(x),\\
 U_{a,b,c}&=\kappa_c^T K_a\kappa_b.
 \end{split}
 \tag{8}
\]

Here U is a full-field matrix expression. The exact master identity is

\[
 \boxed{\begin{aligned}
 N(w_{a,b,c})={}&M^\circ_{a,b,c}-U_{a,b,c}
      -t_cW_{a,b}-t_bW_{a,c}\\
 &+p\bigl[(-1)^{a+c}t_b+(-1)^{a+b}t_c\bigr]-2(-1)^k.
 \end{aligned}}
 \tag{9}
\]

For completeness, the boundary accounting is as follows. Restricting (5)
to x,y≠0 contributes the correction

\[
 p[(-1)^{a+c}t_b+(-1)^{a+b}t_c]-4(-1)^k,
\]

because H(0,y,1)=p·1_(y=1)−1−χ(y), (2)–(3) hold, and
K_a(0,0)=0. On the restricted domain the −1 in (6) subtracts

\[
 U_G=\sum_{x,y\ne0}k_c(x)(K_a)_{xy}k_b(y).
\]

Equation (2) gives **U=U_G+2(−1)^k**, with one contribution from
each zero coordinate. This changes the constant correction from −4(−1)^k
to −2(−1)^k in (9).

For H+1, the x=1 fiber contributes −t_c W_(a,b), and the y=1 fiber
contributes −t_b W_(a,c); their common point contributes zero. On the
remaining domain, y=xt and (3), (6) give exactly

\[
 (K_a)_{xy}k_b(y)[H(x,y,1)+1]
 =\chi(1-x)\chi(t)k_a(t)k_b(xt)k_2(q_x(t)).
 \tag{10}
\]

Alternatively, restoring t=x^(-1) to (7) with the full-field value
κ_2(0)=−1 adds precisely −t_b k_a(x), using (3). Its outside sum is
−t_b W_(a,c). Thus neither the zero-extension convention nor a boundary
fiber is hidden in the notation of (9).

## 3. Uniform cancellation in the inner Möbius correlation

We prove

\[
 \boxed{|T_{a,b}(x)|\le(3a+1)b\,p^{(a+b)/2}
       \quad(x\ne0,1).}
 \tag{11}
\]

The imported input is the hypergeometric theorem reviewed in
[Katz, *G₂ and hypergeometric sheaves*, Section 2, printed pp. 3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf),
and the curve trace formula, Euler–Poincaré formula, and Deligne weight
bound used in [Katz, *Gauss sums, Kloosterman sums, and monodromy groups*,
§2.3 and §3.6](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
Their precise application and normalization to the rank-j sheaves F_j
were proved in the [third pass](parallel3-necklace-2026-09-04.md).
The [archived hypergeometric PDF](../sources/katz-g2-hypergeometric.pdf)
has SHA256 `0bf485bcc9dde2afebc8268af679566bcff65aa6d6ed8c9f4b0387d7e9e954c1`.

Recall that F_j has F_p trace −k_j, is geometrically irreducible of
rank j, and is pure of weight j−1 away from 0,1,∞. It is tame.
Its monodromy at infinity has only the quadratic eigencharacter χ,
with one unipotent block of size j. Its actual middle-extension stalk at
1 has dimension j−1 and trace −t_j. The normalizing constant twist is
fixed over F_p; the extension-field signs from that twist are not chosen
anew. These local and weight statements hold without j<p.

Fix x≠0,1 and consider

\[
 V_x=P^1\setminus\{0,1,x^{-1},\infty\},\qquad
 X_x=G_m\setminus\{x^{-1}\}.
\]

On X_x form the tensor of the actual pullbacks/middle extensions

\[
 \mathcal T=F_a\otimes[t\mapsto xt]^*F_b
               \otimes q_x^*F_2\otimes L_\chi(t).
 \tag{12}
\]

The map q_x is a morphism X_x→G_m. Its special values are

\[
 q_x(0)=\infty,\quad q_x(1)=1,\quad q_x(x^{-1})=0,\quad
 q_x(\infty)=\frac{x}{x-1}\notin\{0,1,\infty\}.
 \tag{13}
\]

On V_x, the tensor has rank 2ab and weight a+b−1, and all four
missing points are tame. Its F_p trace on X_x is minus the summand
of (7), including t=1. The minus sign comes from its three factors F_j.

**The global invariant vanishes at infinity for every pair of ranks.**
At t=∞, F_a and [t↦xt]^*F_b each have scalar quadratic
semisimplification, while q_x^*F_2 is unramified by (13). The Kummer
factor L_χ(t) contributes a third quadratic character. Thus every
inertia eigenvalue of the tensor at infinity is −1. A unipotent factor
does not change those eigenvalues. There is no inertia invariant or
coinvariant there, hence H_c²(V_x,𝒯)=0. This directly excludes global
invariants, even when ranks are equal or some tensor constituent might
otherwise match another factor; no independence hypothesis is used.

Since H_c⁰(V_x,𝒯)=0 and χ_c(V_x)=−2, Euler–Poincaré gives

\[
 \dim H_c^1(V_x,\mathcal T)=4ab.
\]

The point t=1 is included in X_x. Its **actual tensor stalk** has
dimension (a−1)·b·1=(a−1)b, and trace t_a k_b(x), which is minus
the t=1 summand of (7). The factor at xt=x is lisse because x≠0,1.
The tensor of these stalks injects into the invariant subspace of the
generic tensor. Consequently 𝒯 has no punctual subsheaf. The open–closed
exact sequence therefore gives

\[
 \boxed{\dim H_c^1(X_x,\mathcal T)=4ab-(a-1)b=(3a+1)b.}
 \tag{14}
\]

We have not replaced this tensor stalk by the potentially larger stalk
of the middle extension of its generic tensor. The point t=x^(-1)
remains excluded, so no unspecified Frobenius trace at q_x(t)=0 is
imported; its literal finite-sum value was restored in (9).

Deligne bounds the weights of H_c¹(V_x,𝒯) by a+b. Its quotient
H_c¹(X_x,𝒯) has the same upper weight bound. The trace formula and
(14) prove (11). The argument includes a=1 and b=1 and all ranks
at least p; its four points remain distinct for every allowed x.

## 4. From the inner bound to the necklace bound

By (4) and Cauchy–Schwarz,

\[
 \sum_{x\ne0,1}|k_c(x)|\le p^{(c+1)/2}.
\]

Hence (11) gives |M°|≤(3a+1)b p^((k+1)/2). The other quantities
in (9) satisfy

\[
 \begin{split}
 |U|&\le p^{k/2},\\
 |t_cW_{a,b}|, |t_bW_{a,c}|&\le p^{k/2},\\
 p|t_b|, p|t_c|&\le p^{k/2}.
 \end{split}
\]

For W this uses only the last bound in (4), not a conjectured signed
aggregate estimate. The last line uses a+c≥2 and a+b≥2. Equation (9)
now proves (1) with b in place of min(b,c). Swapping x,y in (5)
interchanges b,c exactly, so choose the smaller of them.

The [previous projective transfer](parallel-necklace-2026-09-04.md)
adds at most kp^(k/2) under any permutation of the three labels. Thus
the normalized permuted-family bound is

\[
 \frac{|N|}{p^{k/2+1}}
 \le\frac{(3a+1)\min(b,c)}{\sqrt p}
       +\frac{k+5}{p}+2p^{-k/2-1}.
 \tag{15}
\]

The coefficient is at most (3k+1)²/24, hence O(k²). This yields
uniform decay for k=o(p^(1/4)), including logarithmic depth. It does not
supply cancellation between the exponentially many other words in the
spectral aggregate.

## 5. Fixed-length consequences and verification boundary

At length five, the (2,2,1) words are label permutations of this family,
and (3a+1)min(b,c)≤10. Only four positions at most contribute to the
required inversion-transfer error. Therefore every such word obeys

\[
 |N|\le10p^3+9p^{5/2}+2\le15p^3\qquad(p\ge5).
\]

The other length-five words have multiplicities (3,1,1), already covered
by the two-singleton-exception theorem, or use at most two labels. For
two labels, every nonconstant length-five word has minority count at
most two, so the prior bound 30p^(5/2) applies; the constant-word bound
is also below 15p³. Thus the same constant covers all 243 words.

The exact finite coverage counts are:

| Length | At most two labels | Two singleton exceptions | Three-gap class | Remaining |
|---|---:|---:|---:|---:|
| 5 | 93 | 60 of type (3,1,1) | 90 of type (2,2,1) | 0 |
| 6 | 189 | 90 of type (4,1,1) | 360 of type (3,2,1) | 90 of type (2,2,2) |

At length six, the at-most-two-label claim uses the previously proved
complete binary length-six result and its projective transfer, rather
than a minority-count assertion that would miss the balanced binary words.

The [exact verifier](../experiments/parallel5_necklace_2026_09_04.py)
checks the Möbius maps, every derived inner factorization and bound in
its bounded fields/ranks, the included t=1 term, the excluded t=x^(-1)
term, the full U versus U_G correction, the master identity, the final
bound and reversal, and all 243 length-five words. Original necklaces
are also summed literally in small cases. Integer comparisons handle
half-integral exponents without floating-point acceptance.

[Results](../results/parallel5_necklace_2026_09_04.json) retain counts,
values, missing words at length six, and source hashes. These checks
supplement the uniform proof and do not validate cohomological theorems
by enumeration. The full degree-two all-word estimate, the required
weighted spectral aggregate, both Paley targets, and the prize bridge
remain unproved.
