# A two-anchor Schur kernel and one balanced necklace class

**Status: a new operator bound and one balanced length-six class are
proved.** The other four balanced classes remain open in this pass.
Let p≡1 mod 4 be prime, χ(0)=0, and write

\[
 S_{xy}=\chi(x-y),\quad D_A=\operatorname{diag}(\chi(x)),\quad
 D_B=\operatorname{diag}(\chi(x-1)),\quad D_C=D_A D_B.
\]

For the real symmetric matrices

\[
 Q=S D_C S,\qquad W_C=S\circ Q,
\]

where the circle denotes entrywise multiplication, we prove

\[
 \boxed{\|W_C+S\|\le2p,\qquad
        \|W_C\|\le2p+\sqrt p.}
 \tag{1}
\]

The same bounds hold, after translation/scaling, for any two distinct
anchors. This is a bound for the displayed Schur kernel, not for the
two-anchor localized Paley adjacency matrix or its spectral edge.

Put f_0=S e_0, f_1=S e_1, and t_4=tr((D_A S)^4)/(p−1). There is
an exact necklace identity

\[
 \boxed{N(AABBCC)=f_0^T W_C S W_C f_1-p t_4+2.}
 \tag{2}
\]

Consequently

\[
 \boxed{|N(AABBCC)|\le7p^{7/2}.}
 \tag{3}
\]

Its twelve-word dihedral orbit is also its entire label-permuted orbit,
so all twelve words have exactly the same value; no transfer error is
needed for this class. Combined with [pass five](parallel5-necklace-2026-09-04.md),
this gives finite coverage of **651 of 729** degree-two words of length
six. The remaining **78** are listed in the results. These counts are
not a measure of progress toward the full Paley quantifiers.

## 1. The single-anchor input and an exact extra coordinate

Let K_2=S D_A S and W_A=S∘K_2. On F_p*, the
[second pass](parallel2-necklace-2026-09-04.md) proves

\[
 \|W_A^*\|\le2p,\qquad W_A^*\mathbf1=2\mathbf1,
 \qquad (W_A)_{0y}=-1\quad(y\ne0).
 \tag{4}
\]

The analytic input for the first inequality is
[Katz, *Sato-Tate theorems for finite-field Mellin transforms*,
Theorem 15.1 and Corollary 4.2](https://web.math.princeton.edu/~nmk/mellin186.pdf).
The precise TwLeg trace normalization, all good-character conditions,
and the exceptional trivial eigenvalue 2 were checked in that earlier
proof. No additional sheaf identification, rank assumption, or purity
claim is imported here.

Adjoin a coordinate ∞ to W_A and set

\[
 \widetilde W_{xy}=(W_A)_{xy}\quad(x,y\in F_p),\qquad
 \widetilde W_{\infty,y}=\widetilde W_{y,\infty}
      =p\,1_{y=0}-1,\qquad \widetilde W_{\infty,\infty}=0.
 \tag{5}
\]

The entries in this extra row have not been approximated by constants.
Let n=p−1 and v=1_(F_p*) on the finite coordinates. The span of
e_∞,e_0,v is invariant, and the matrix in this (not normalized) basis is

\[
 \begin{pmatrix}0&n&-n\\n&0&-n\\-1&-1&2\end{pmatrix}.
 \tag{6}
\]

Its three eigenvectors/eigenvalues are

\[
 \begin{array}{c|c}
 e_\infty-e_0&-(p-1)\\
 e_\infty+e_0+v&0\\
 (p-1)(e_\infty+e_0)-2v&p+1.
 \end{array}
\]

The orthogonal complement is the zero-sum subspace of F_p*, where the
operator is W_A*. Thus (4) and p+1≤2p imply

\[
 \boxed{\|\widetilde W\|\le2p.}
 \tag{7}
\]

This argument handles every added-coordinate mode explicitly. A
rank-only estimate of the border is unnecessary.

## 2. Projective conjugation, including the pole

For x≠0 write T(x)=x^(-1)−1 and put T(0)=∞. This is a bijection
from F_p onto P¹(F_p)\{−1}. The quartic identity already derived in
[pass four](parallel4-necklace-2026-09-04.md) says

\[
 Q_{xy}+1=\chi(xy)K_2(T(x),T(y))\quad(x,y\ne0).
\]

Also S_(T(x),T(y))=χ(xy)S_(xy). Therefore

\[
 (W_C+S)_{xy}=(W_A)_{T(x),T(y)}\quad(x,y\ne0).
 \tag{8}
\]

At the pole use the exact quadratic correlation

\[
 Q_{0y}=p\,1_{y=1}-1-\chi(y).
\]

For y≠0 this gives

\[
 (W_C+S)_{0y}=p\,1_{y=1}-1
       =p\,1_{T(y)=0}-1=\widetilde W_{\infty,T(y)}.
\]

Both diagonal entries indexed by 0/∞ are zero. Hence (8) extends to
the whole finite matrix:

\[
 \boxed{W_C+S\text{ is permutation-conjugate to }
 \widetilde W|_{P^1(F_p)\setminus\{-1\}}.}
 \tag{9}
\]

The first bound in (1) follows by principal compression of (7), and the
second by ||S||=√p. For anchors a≠b, scaling x=a+(b−a)u changes
W_(a,b) and S by the common sign χ(b−a), while Q is unchanged. Thus
both norms in (1) are uniform in the anchor pair.

## 3. Normalize the adjacent doubleton vertices

Form the graph G of the necklace AABBCC: six cycle vertices, two
anchor vertices, and the eight label spokes. Its total edge count is
14. Let Z(G) be its full affine character sum, with both anchor values
also summed. Simultaneous affine substitutions preserve every summand's
overall scaling sign because 14 is even. If the anchors are distinct,
their contribution is p(p−1)N(AABBCC). If they coincide, translation
reduces their contribution to p N_0, where

\[
 N_0=\operatorname{tr}((D_A S)^4(P_0S)^2),
 \qquad P_0=I-e_0e_0^T.
 \tag{10}
\]

The diagonal at a collapsed C label is χ(x)²=P_0(x,x), not the
constant function one. Thus

\[
 Z(G)=p(p-1)N(AABBCC)+pN_0.
 \tag{11}
\]

Instead normalize the last two, adjacent C vertices to values 1 and 0.
Their equal-value contribution is zero because their cycle edge has
character zero. With these two values fixed, summing each old anchor
produces Q on its two remaining neighbors. The four remaining cycle
vertices have path-label order A,A,B,B. The path sum is exactly

\[
 F=f_0^T W_C S W_C f_1.
\]

There are p(p−1) distinct ordered values for the two C vertices, and
affine invariance gives Z(G)=p(p−1)F. Comparison with (11) yields

\[
 N(AABBCC)=F-\frac{N_0}{p-1}.
 \tag{12}
\]

## 4. Evaluate the coincident-anchor contribution

Let C*=D_A S restricted to F_p* and let u=(χ(x))_(x≠0). Both
D_A S and P_0S have zero row at 0, so their product in (10) has trace

\[
 N_0=\operatorname{tr}((C^*)^4(S^*)^2).
\]

The full identity S²=pI−J gives

\[
 (S^*)^2=pI-J-uu^T.
\]

Since C*1=−1 and C*u=−u, while both vectors have squared norm p−1,

\[
 \boxed{N_0=(p-1)(p t_4-2).}
 \tag{13}
\]

This proves (2), including its +2 correction. Nothing at the collapsed
anchor or at the projective pole has been discarded.

Finally ||f_0||=||f_1||=√(p−1), so

\[
 |N(AABBCC)|\le
 (p-1)\sqrt p\,(2p+\sqrt p)^2+p|t_4|+2.
 \tag{14}
\]

Use |t_4|≤p² and √p≤p/2 for p≥5. Dividing the resulting bound
by p^(7/2) gives at most 25/4+1/√p+2/p^(7/2)<7. This proves (3).

## 5. What this does not settle

The 90 balanced words have the following disjoint orbits under rotation,
reversal, and label permutation:

| Representative | Orbit size | Status in this pass |
|---|---:|---|
| AABBCC | 12 | proved |
| AABCBC | 36 | open |
| AABCCB | 18 | open |
| ABACBC | 18 | open |
| ABCABC | 6 | open |

The first orbit already consists solely of rotations/reversals of
AABBCC, which proves the exact-value assertion made above.

Normalizing an adjacent equal-label pair in the next two classes gives
path pairings ABAB and ABBA, rather than AABB. In terms of the same
Q and W_C their fixed-pair contractions are respectively

\[
 \begin{split}
 F_{ABAB}&=\sum_{x,y}S_{xy}
                 (S D_A Q)_{xy}(Q D_B S)_{xy},\\
 F_{ABBA}&=\sum_{x,y}(f_0)_x(f_1)_y Q_{xy}(S W_C S)_{xy}.
 \end{split}
\]

These exact formulas are checked, but **are not additional word bounds**.
The operator inequalities proved here give only the order-p⁴ scale
for these contractions, which does not establish their required
o(p⁴) estimate. Applying the AABB chain bound to either contraction
would be invalid. The last two classes have no adjacent equal labels;
the fixed-pair argument above cannot remove their coincident-C fiber
using a zero cycle edge. They remain open as well.

No growing all-word or signed aggregate theorem follows from this pass.
The uniform bound (1) is a reusable operator estimate; its use for other
contractions still needs a proof.

## Verification

The [verifier](../experiments/parallel6_necklace_2026_09_04.py) checks
the full pole-inclusive conjugation (9), every added-coordinate
eigenvector, the norm bounds using exact integer positivity certificates,
the two evaluations of Z(G), the exact coincident-anchor correction,
the necklace identity and bound, and the twelve-word exact-value orbit.
It also verifies the two unresolved contraction formulas without marking
them as proved bounds, and independently enumerates all five balanced
orbits. Small-field graph and necklace sums use their original definitions.

[Results](../results/parallel6_necklace_2026_09_04.json) record exact
values, coverage limits, source hashes, and terminal check counts. The
finite tests supplement the proof; they do not prove the imported
TwLeg theorem or any assertion about the remaining 78 words. Both Paley
targets, the full weighted spectral aggregate, and the prize bridge
remain unproved.
