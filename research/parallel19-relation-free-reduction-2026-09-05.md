# A relation-free reduction for the full classical Paley target

Status: an unconditional combinatorial reduction and a restricted
single-size sufficient moment criterion. The required character estimate
on the restricted sets is unproved. No new cancellation exponent for
arbitrary sets, subgroup period bound, clique bound, or prize implication
is claimed. No novelty claim is made for the elementary extraction method.

Let p be an odd prime, χ(0)=0, and
S(A,B)=Σ_{a∈A,b∈B}χ(a−b). A B_h set is one whose equal h-term
sums, with repetitions allowed, have identical multisets of summands.
This property is hereditary and implies B_j for every 1≤j≤h, by
padding the shorter sums with a fixed element. In particular B_2 is
the ordinary additive Sidon property, including exclusion of nontrivial
three-term arithmetic progressions. This is not just a restriction
on sums of distinct elements.

## 1. An explicit extension certificate

Fix h≥2 and assume p>h. If S is a nonempty B_h set of size s, define

\[
\mathcal F_h(S)=S\cup\bigcup_{j=1}^h j^{-1}
                 \bigl(hS-(h-j)S\bigr),                 \tag{1}
\]

where tS is the set of t-term sums, with repetitions, and 0S={0}.
For S=∅ set \(\mathcal F_h(S)=∅\).
Then, for x∉S,

\[
S\cup\{x\}\text{ is }B_h
\quad\Longleftrightarrow\quad x\notin\mathcal F_h(S).
                                                               \tag{2}
\]

To prove the forward implication, membership gives
jx+u₁+⋯+u_{h−j}=v₁+⋯+v_h with all u_i,v_i∈S.
This is a nontrivial h-term equality since x∉S. Conversely, suppose
adding x creates a nontrivial equality. Cancel all copies of x common
to its two sides, leaving j≥1 extra copies on one side. If no copies
remain on either side, padding would contradict the B_h property of
S. Otherwise pad both sides with the same copies of an element of S
until the side without x again has length h. The equality then has
exactly the displayed form, with 1≤j≤h. Division by j is valid since
p>h. The empty-set case follows because a singleton is B_h.

Counting ordered sums, even though some may coincide, gives

\[
|\mathcal F_h(S)|\le f_h(s):=s+\sum_{j=1}^h s^{2h-j}.
                                                               \tag{3}
\]

Set \(L_h(k)=f_h(k-1)\), for k≥1. If a set D has size greater
than L_h(k), the greedy extension rule (2) constructs a B_h subset
of size exactly k: at every previous size s<k, at most
f_h(s)≤L_h(k) possible elements of D are forbidden. This proof
does not import a sharp Sidon-extraction theorem.

Consequently every B⊂F_p admits a disjoint partition

\[
B=C_1\sqcup\cdots\sqcup C_q\sqcup R,\quad
|C_i|=k,\quad C_i\text{ is }B_h,\quad |R|\le L_h(k).      \tag{4}
\]

Repeatedly extract a block while the remaining size exceeds L_h(k).
For k=1 the remainder is empty. For k≥2, L_h(k)≥k−1, so
every extraction has enough elements. The process terminates because
each step removes k elements. In particular

\[
L_h(k)\le(h+1)k^{2h-1}.                                 \tag{5}
\]

At h=2 the explicit remainder bound is
L₂(k)=(k−1)+(k−1)²+(k−1)³. It is deliberately conservative.

## 2. The complete two-set target can be restricted to B_h sets

Let K(a,b) be any real or complex kernel with |K(a,b)|≤1.
Partition A and B as in (4), with remainders of sizes r_A,r_B.
If every pair of complete blocks satisfies
\(|\sum_{a\in C_i,b\in D_j}K(a,b)|\le\eta k^2\), then

\[
\begin{split}
\frac{|\sum_{A\times B}K|}{mn}
&\le\eta\frac{(m-r_A)(n-r_B)}{mn}
  +\frac{r_A}{m}+\frac{r_B}{n}-\frac{r_Ar_B}{mn}\\
&\le\eta+\frac{L_h(k)}m+\frac{L_h(k)}n.                 \tag{6}
\end{split}
\]

The exceptional pairs are precisely those with at least one endpoint
in a remainder; their number is mr_B+nr_A−r_Ar_B. No sign or
positivity assumption on K is used.

Fix any h≥2. Suppose the classical Paley cancellation statement holds
when **both** input sets are B_h: for every σ>0 there is δ_h(σ)>0
such that, uniformly for B_h sets C,D of sizes greater than p^σ,

\[
|S(C,D)|\ll_{h,\sigma}|C||D|p^{-\delta_h(\sigma)}.       \tag{7}
\]

Then it holds for arbitrary input sets. Indeed, given ε>0, choose
θ=ε/[2(2h−1)] and k=⌊p^θ⌋. For sufficiently large p,
k>p^{θ/2}, so (7) applies to every block. For m,n>p^ε,
(5)–(6) give

\[
\frac{|S(A,B)|}{mn}
\ll_{h,\varepsilon}p^{-\delta_h(\theta/2)}+p^{-\varepsilon/2}.
                                                               \tag{8}
\]

The reverse implication is immediate by restriction. Thus for every
fixed h, the full all-exponents classical conjecture is **equivalent**
to its restriction to B_h input pairs. In particular, a proof on
Sidon pairs at every positive size exponent would suffice.
This does not assert cancellation on those pairs; it proves that
requiring the absence of finitely many orders of additive relations
does not remove the full difficulty in the all-exponents formulation.

## 3. Restricting the sufficient single-size moment criterion

Write \(F_C(x)=\sum_{c\in C}\chi(x-c)\) and
\(M_{2r}(C)=\sum_x|F_C(x)|^{2r}\).
Suppose a bound \(M_{2r}(C)\le U\) is known only for B_h sets
of one size k. Partition B as in (4). Hölder for each block and
the trivial remainder bound give, for every nonempty A and B,

\[
\boxed{\frac{|S(A,B)|}{mn}
\le\left(\frac{U}{m k^{2r}}\right)^{1/(2r)}
    +\frac{L_h(k)}n.}                                  \tag{9}
\]

Indeed |S(A,C_i)|≤m^{1−1/(2r)}U^{1/(2r)}, the sum has q≤n/k
complete blocks, and |S(A,R)|≤m|R|. The assertion remains valid
when q=0. More precisely, if U is the maximum of the actual block
moments, with U=0 for no blocks, the entirely integer inequality is

\[
(|S(A,B)|-m|R|)_+^{2r}\le q^{2r}m^{2r-1}U.             \tag{10}
\]

Now take an unbounded sequence of fixed integers r_j≥2, integers
h_j≥2 with h_j/r_j→0, and β_j≥0 with β_j→0. Put
k_j(p)=⌊p^{1/(r_j+1)}⌋. The following restricted hypothesis suffices
for the full classical Paley conjecture:

\[
\boxed{M_{2r_j}(C)\le C_j p^{1+\beta_j}k_j(p)^{r_j}
\quad\text{for every }B_{h_j}\text{ set }C\text{ of size }k_j(p).}
\tag{SS-B}
\]

For each j, C_j and the sufficiently-large-prime threshold may depend
on j; all sets in that field and slice must satisfy the same bound.
No estimate at a moment order growing with p is assumed.

For a fixed j write θ=1/(r_j+1), h=h_j and r=r_j. Equation (9),
k≥p^θ/2 and m,n>p^ε imply

\[
\frac{|S(A,B)|}{mn}
\le(C_j2^r)^{1/(2r)}p^{(\theta+\beta_j-\varepsilon)/(2r)}
 +(h+1)p^{(2h-1)\theta-\varepsilon}.                    \tag{11}
\]

Choose one j large enough that θ+β_j<ε/2 and
(2h_j−1)θ<ε/2. This is possible precisely from the stated limits.
Both terms have a fixed power saving; for example, after absorbing
the fixed constants into the prime threshold, δ=ε/(8r_j) works.
The case h_j=2 for every j is allowed. One can also choose h_j
unbounded, for example h_j=⌊√r_j⌋ once r_j≥4.

Unlike the [earlier unrestricted SS criterion](parallel5-classical-2026-09-04.md),
SS-B quantifies only over sets with all additive energies through h_j
equal to their permutation-only minima. It is a restricted
domain for the proposed estimate, but its conclusion is still the full
all-exponents classical conjecture. No converse from Paley to SS-B
is proved. Merely assuming low additive energy does not establish SS-B.

## 4. The extraction loss matters

For fixed r,h and β, the displayed argument gives a saving when

\[
\varepsilon>\max\{(1+ (r+1)\beta)/(r+1),\ (2h-1)/(r+1)\}.
                                                               \tag{12}
\]

Thus a fourth-moment bound only for Sidon sets on k=⌊p^{1/3}⌋
does **not**, through this argument, prove Paley above ε=1/3:
the remainder condition would require ε>1. High moment orders
are essential to the sufficient sequence SS-B.

Nor can a B_h block always be extracted at a size independent of h.
If C⊂{0,…,N−1} is B_h modulo p and h(N−1)<p, then its
\(\binom{|C|+h-1}{h}\) multiset sums are distinct integers in
{0,…,h(N−1)}. Hence

\[
\binom{|C|+h-1}{h}\le h(N-1)+1,\qquad
|C|^h\le h!\,[h(N-1)+1].                               \tag{13}
\]

This elementary example explains why extraction incurs a size loss.
It does not prove optimality of (3)–(5), nor rule out better extraction
or transfer methods. In particular the condition h_j=o(r_j) is a
sufficient condition for this proof, not a claimed necessary condition
for every possible approach.

The earlier forced-row construction refutes Gaussian high moments on
B_h sets at size p^{1/r} for r>h. It does not refute SS-B: here the
size is p^{1/(r+1)}, and its forced-row lower contribution divided
by p k^r is only O(p^{−1/(r+1)}log p). This distinction is about
that particular obstruction, not a proof that SS-B is true.

## 5. Relation to the transformation route and validation

For completeness, return to p≡1 mod4 and the Weil normalization of
pass 18. Its simplest mixed product does not
add a missing estimate. With g(s,t)=u_swu_t, d=s'−s and e=t'−t,
its quotient has trace 2−de. When d,e≠0 the Weil character is
χ(de); when exactly one is zero it is √p times the corresponding
nonzero difference character, and when both are zero it is p.
Summing these four cases reproduces exactly the factorized
Hilbert–Schmidt identity already proved in pass 18: if
s_U=S(U,U), s_V=S(V,V), the unnormalized sum is
p mn+√p(ns_U+ms_V)+s_Us_V. Division by m²n² gives the squared
Hilbert–Schmidt norm. This is why the
present pass pursues the restricted moment reduction instead.

The [exact verifier](../experiments/parallel19_relation_free_2026_09_05.py)
tests the forbidden-extension identity against independent multiset-sum
enumeration, constructs partitions, checks minimal energies, verifies
(6) and (10) in exact arithmetic, and checks the exponent ledger.
Finite tests complement the proofs; they do not prove SS-B or (7).

The definition and broader context of additive Sidon extraction were
checked in [Shkredov, arXiv:2103.14670v1, Introduction](https://arxiv.org/html/2103.14670v1).
The sharp extraction results discussed there are not imported: (1)–(13)
are proved here. No new external theorem is a proof input. The full
classical conjecture, subgroup target and official prize bridge remain
open; independent mathematical review is outstanding.
