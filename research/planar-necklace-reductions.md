# Planar reductions for the remaining binary necklaces of length six

**Status: Paley, the full localization estimates, and the prize remain
unproved.** This note resolves the three length-six singleton patterns
left outside the previous one-/two-occurrence argument. It also evaluates
every two-block singleton necklace, with arbitrary block lengths.

Let p≡1 mod 4 be prime. Retain
`S_xy=χ(x−y)`, `D_b=diag(χ(x−b))`, `C_b=D_b S`, and
`t_j=tr(C_0^j)/(p−1)` for j≥1, with t₀=1. A binary word w denotes
the necklace `N(w)=tr(C_(w_1)⋯C_(w_k))`. Its anchors 0,1 may be
replaced by any two distinct elements, by an affine substitution.

The exact identities are

\[
 \boxed{N(0^r1^s)=p\,t_r t_s-t_{r+s}\quad(r,s\ge1),}
 \tag{1}
\]

\[
 \boxed{N(001011)=p(t_5+2t_3)-t_6,}
 \tag{2}
\]

\[
 \boxed{N(010101)=p(p-5)t_4+p^2(p-2)-t_6.}
 \tag{3}
\]

In particular N(000111)=p t₃²−t₆. The proofs use affine averaging,
finite Fourier duality of planar graphs, and projective changes of variables.
They are ordinary mathematical proofs, not numerical pattern fits or Lean
certificates. No novelty claim is made.

## Consequences from the existing analytic input

The [preceding note](localized-necklace-identities.md) derives

\[
 |t_j|\le(j-1)p^{(j-1)/2}
 \tag{4}
\]

from [Lu–Zheng–Zheng, Lemma 2.1 (2.3), with R=1](https://arxiv.org/html/1305.3405v3#S2).
No new analytic theorem is assumed here. Formula (1) immediately gives

\[
 |N(0^r1^s)|\le(r+s)^2p^{(r+s)/2}.
 \tag{5}
\]

Its normalized bound is `(r+s)²/p`, so this subfamily is also controlled
at logarithmic length. No claim is made for the exponentially many other
words needed at that depth.

For the balanced words of length six, (1)–(4) give, respectively,

\[
 \begin{aligned}
 |N(000111)|&\le4p^3+5p^{5/2},\\
 |N(001011)|&\le4p^3+4p^2+5p^{5/2},\\
 |N(010101)|&\le3p^{7/2}+p^3+5p^{5/2}.
 \end{aligned}
 \tag{6}
\]

Rotation, reversal, and exchange of the two anchors leave N unchanged.
The three cyclic gaps between three marked vertices sum to six, so the
balanced cases have representatives 000111, 001011, and 010101 under
these symmetries. Every unbalanced word has a minority anchor occurring
at most twice and is covered by the preceding note. Constant words have
value `(p−1)t₆`. Thus, with a deliberately loose explicit constant,

\[
 \boxed{|N(w)|\le36p^{7/2}\quad\text{for every }w\in\{0,1\}^6.}
 \tag{7}
\]

After normalization by p⁴ this tends to zero. This covers all 64 binary
singleton words at length six. It does not cover all lengths or words
using the doubleton label {0,1}.

## Affine averaging proves the two-block identity

For j≥1 put `Q_j=Σ_(b∈F_p) C_b^j`. Translation makes its entries
depend only on x−y. Simultaneous nonzero scaling fixes C_b entries
because each contains two quadratic-character factors. Thus Q_j has one
diagonal value and one off-diagonal value. Its diagonal value is
`tr(C_0^j)=(p−1)t_j`, by summing the translated diagonal over b.
Every C_b kills the constant vector, hence Q_j does too. Therefore

\[
 \boxed{\sum_b C_b^j=t_j(pI-J),}
 \tag{8}
\]

where J is the all-ones matrix. It follows that

\[
 \sum_b\operatorname{tr}(C_0^r C_b^s)
 =t_s\operatorname{tr}(C_0^r(pI-J))
 =p(p-1)t_r t_s.
\]

The b=0 term is `(p−1)t_(r+s)`, and the other p−1 terms are equal
by scaling. This proves (1), without a restriction on either block length.

## Two general graph identities

For a finite graph G, define the affine character partition function

\[
 Z_p(G)=\sum_{x:V(G)\to F_p}\prod_{uv\in E(G)}\chi(x_u-x_v).
 \tag{9}
\]

Repeated edges are retained. A loop makes the product zero. As χ(−1)=1,
edge orientation does not change the sum.

### Planar Fourier duality

For a connected plane graph G with plane dual G*, v vertices and e edges,

\[
 \boxed{Z_p(G)=p^{v-e/2-1}Z_p(G^*).}
 \tag{10}
\]

Use the quadratic Gauss expansion
`χ(t)=p^(−1/2)Σ_y χ(y)e_p(yt)` on each edge. Summing all vertex
variables yields p^v times the indicator that the edge variables form a
divergence-free flow. Such flows are exactly differences of dual face
potentials. The potential-to-flow map has kernel the p constant potentials,
so the flow sum is Z_p(G*)/p. This proves (10), with every factor of p
included. The argument also accounts for dual loops and parallel edges.

### Projective changes of variables and infinity terms

Use homogeneous representatives `[x:1]` for finite x and `[1:0]` for ∞.
Define the projective edge weight by χ of their determinant. It agrees
with χ(x−y) at finite points, equals 1 between ∞ and a finite point,
and is zero between two copies of ∞.

If every vertex of G has even degree and G has an even number of edges,
the product is invariant under PGL₂(F_p): changes of representatives
cancel by even degrees, and the determinant of the transformation contributes
`χ(det M)^e=1`. Let the sum over all projective assignments be Ẑ_p(G).
In particular, the sum with any specified vertex fixed to ∞ is
Ẑ_p(G)/(p+1), independently of the chosen vertex.

For any graph, the assignments with precisely I at ∞ contribute
Z_p(G−I) if I is independent, and zero otherwise. Hence

\[
 \widehat Z_p(G)=\sum_{I\text{ independent}}Z_p(G-I),
 \qquad
 \widehat Z_p(G;x_v=\infty)
   =\sum_{\substack{I\text{ independent}\\v\in I}}Z_p(G-I).
 \tag{11}
\]

The second expression has the same value for all v under the preceding
even-degree/even-edge hypotheses. The affine summands must not be silently
replaced by projective sums; their boundary terms drive the corrections below.

## Turning a binary necklace into a planar graph

Attach two anchor vertices a,b to a cycle of length k, connecting each
cycle vertex to the anchor named by its binary label. Call the graph G_w.
Place a inside the cycle and b outside. This gives a plane embedding with
v=k+2, e=2k, and k faces whenever both anchors occur. Summing the anchors
first gives

\[
 Z_p(G_w)=p(p-1)N(w)+p(p-1)t_k.
\]

Formula (10) has exponent one here. Thus

\[
 \boxed{N(w)=\frac{Z_p(G_w^*)}{p-1}-t_k.}
 \tag{12}
\]

The dual can be constructed explicitly: an a-face lies between each pair
of consecutive a-spokes, and likewise for b. The a-faces form a cycle,
as do the b-faces; every original cycle edge adds a cross-edge between
its two incident faces. This prescription retains multiplicities when
there are only two spokes or two faces share multiple edges.

## The alternating word: cube to octahedron

For w=010101, G_w is a cube and G_w* is an octahedron O=K_(2,2,2).
Write its three nonedge pairs as `(x₁,x₂)`, `(y₁,y₂)`, `(z₁,z₂)`.
It has degree four at every vertex and twelve edges, so projective
invariance applies.

First suppose x₁≠x₂. Any nonzero summand has y₁ distinct from both.
PGL₂ acts sharply transitively on ordered triples of distinct points;
normalize `(x₁,x₂,y₁)=(0,∞,1)`. There are p(p²−1) such triples.
Writing y₂=t, z₁=u, z₂=v, the remaining sum is

\[
 W=\sum_{t\ne0}\chi(t)E(t)^2,\qquad
 E(t)=\sum_u\chi(u(u-1)(u-t)).
 \tag{13}
\]

All zero factors and the allowed coincidences within nonedge pairs are
included. In particular, t=1 is allowed. If x₁=x₂, move the common
value to ∞. The remaining four vertices are finite and contribute the
four-cycle partition `tr(S⁴)=p²(p−1)`. There are p+1 choices for the
common point. Therefore

\[
 \widehat Z_p(O)=p(p^2-1)(W+p).
 \tag{14}
\]

The earlier kernel K₂=S D₀ S satisfies `K₂(1,t)=E(t)`.
Homogeneity then gives

\[
 (p-1)t_4=\operatorname{tr}(D_0K_2D_0K_2)
          =(p-1)\sum_{t\ne0}\chi(t)E(t)^2.
\]

Thus W=t₄ exactly; we do not separately estimate the absolute elliptic
traces. Finally, O has six singleton independent sets and three independent
pairs. A single ∞ leaves a wheel with four rim vertices, contributing
p(p−1)t₄. A pair at ∞ leaves a four-cycle. Equation (11) gives

\[
 \widehat Z_p(O)-Z_p(O)=6p(p-1)t_4+3p^2(p-1).
\]

Combining with (14),

\[
 Z_p(O)=(p-1)\bigl[p(p-5)t_4+p^2(p-2)\bigr].
\]

Equation (12) proves (3). The cancellation in the weighted cubic expression
from the preceding note has been obtained here by duality and by retaining
the signed second elliptic moment W.

## The mixed word: completion reduces the dual to a wheel

For w=001011, label its three a-faces a,b,c and its three b-faces d,e,f,
in cyclic order. The dual H has edges

\[
 ab,ac,bc,\ de,df,ef,\ af,bf,bd,cd,ce,cf.
 \tag{15}
\]

Equivalently it is K₆ with the path edges da,ae,eb removed. Add a new
vertex g connected to a,c,e,f, and call the resulting graph H⁺. Every
degree is even, and it has sixteen edges. Vertex c is universal.

Fixing c=∞ forces every other vertex finite. Removing c leaves a wheel
with hub f and rim `a,b,d,e,g`, so this conditional sum equals
`p(p−1)t₅`.

Fixing g=∞ instead permits b and d to be at ∞, but not both because
bd is an edge. Thus (11) and projective invariance give

\[
 Z_p(H)+Z_p(H-b)+Z_p(H-d)=p(p-1)t_5.
 \tag{16}
\]

In H−b, vertex a has degree two, with neighbors c,f that are adjacent.
Summing a gives `p 1_(x_c=x_f)−1`; the equality term vanishes on the
edge cf. The remaining graph is K₄, so `Z_p(H−b)=−Z_p(K₄)`.
The same argument at e works in H−d. Fixing one K₄ vertex by translation
leaves a monochromatic triangle, giving `Z_p(K₄)=p(p−1)t₃`.
Hence

\[
 Z_p(H)=p(p-1)(t_5+2t_3).
\]

Equation (12) proves (2). ∎

## Verification and the unresolved scope

[The exact script](../experiments/planar_necklace_reductions.py) checks
arbitrary two-block lengths with total at most twelve, all 64 binary
length-six words, and the affine matrix average. It constructs the planar
duals from oriented face boundaries and checks Euler's formula and opposite
edge incidences. In small fields it independently sums the original graphs,
their duals, the completed graph conditioned at two different vertices,
both deletion terms, and the affine/projective octahedron. It also directly
computes the signed elliptic moment (13). These checks include boundary
assignments and do not substitute the claimed formulas when summing graphs.

[Results](../results/planar_necklace_reductions.json) record exact values,
faces, edges, trace sequences, and source hashes. Arithmetic is in Python
integers, including NumPy object matrices. The finite checks supplement the
uniform proofs; they do not establish any additional asymptotic assertion.

The general three-occurrence problem at larger lengths is still open here,
as are necklaces using all three labels {a},{b},{a,b}. Even a proof of all
fixed-length necklace estimates would leave the spectral-edge/growing-depth
requirements and the full Paley quantifiers. No estimate for arbitrary
small sets, no uniform thin-subgroup bound, and no prize reduction follows
from the present identities.
