# Two middle convolutions complete the length-six necklace estimate

**Status: all 729 degree-two words of length six have a uniform
O(p^(7/2)) bound.** The remaining 42 words from the seventh pass
are covered here. This is a fixed-length result. The full
growing-depth signed aggregate, localized spectral edge, classical
Paley conjecture and prize reduction remain unproved.

The new bounds, for every prime p≡1 mod 4, are

\[
 |N(ABBACC)|\le46p^{7/2},\qquad
 |N(ABCBAC)|,\ |N(ABCABC)|\le45p^{7/2}.                 \tag{1}
\]

Their full label-permuted dihedral classes have sizes 18,18,6 and
bounds 48,47,45 times p^(7/2), respectively. Combined with the
earlier 687 words, this completes length six. Ordinary mathematical
proofs and finite checks are distinguished below; no novelty or
formal-verification claim is made.

## 1. Definitions and a common contraction

Retain χ(0)=0 and the full F_p-coordinate matrices

\[
 S_{xy}=\chi(x-y),\quad D_a=\operatorname{diag}\chi(x-a)\ (a=0,1),
 \quad D_C=D_0D_1,\quad Q=S D_C S.
\]

The labels A,B,C mean D_0,D_1,D_C, and
N(w)=tr(D_(w1)S⋯D_(wk)S). Define, for a,b∈{0,1},

\[
 L_{ab}=S D_a Q D_b S
       =S D_a S D_C S D_b S,\qquad
 J^{a,b,c}_y=\sum_x\chi(x-c)Q_{xy}(L_{ab})_{xy}.        \tag{2}
\]

For y∉{0,1} and c∈{0,1,y}, we prove

\[
 \boxed{|J^{a,b,c}_y|
       \le32p^{5/2}+26p^2+2p^{3/2}.}                 \tag{3}
\]

For y=0,1 the same sum, with any c, has the simpler bound p³:
||Q||≤p and ||L_ab||≤p², so column Cauchy–Schwarz gives it directly.

## 2. The rank-four intermediate object

All middle convolutions use the quadratic Kummer sheaf and the
arithmetic convention from the
[seventh-pass proof](parallel7-necklace-2026-09-04.md). For fixed
y∉{0,1}, let E_y be the rank-two sheaf with actual trace

\[
 e_y(x)=-Q_{xy}-1.
\]

Its local types at 0,1,y are J₂, and at infinity χI₂. It is
geometrically irreducible, tame and pure of weight one. In particular,
|e_y(x)|≤2√p at every finite point. The exact exceptional values are

\[
 e_y(0)=\chi(y),\quad e_y(1)=\chi(1-y),\quad
 e_y(y)=\chi(y(y-1)).                                  \tag{4}
\]

Write K_b=S D_b S. The affine Legendre sheaf B_(b,y) has trace
−(K_b)_(xy), local type J₂ at b,y, and χJ₂ at infinity; it is
lisse at the other anchor 1−b.

Multiply this sheaf by the Kummer character χ(x(x−1)) and take
middle extension. Denote the resulting rank-two object G_(b,y).
Its types are

| Point | G_(b,y) |
|---|---|
| b | χJ₂ |
| 1−b | χI₂ |
| y | J₂ |
| infinity | χJ₂ |

The quadratic factors have no invariants at either finite anchor,
so the zero masks give the actual middle-extension traces.

Let B'_(b,y)=MCχ(G_(b,y)). It is irreducible, tame and pure of
weight two. The finite invariant codimensions sum to 2+2+1=5,
and G_(b,y)(infinity)⊗χ has one invariant, so its rank is four.
Its complete local types are

| Point | B'_(b,y) |
|---|---|
| b | J₃⊕1 |
| 1−b | J₂⊕J₂ |
| y | χ⊕1³ |
| infinity | 1⊕χI₃ |

These types follow from the finite quotient rule and the infinity
auxiliary representation M=χJ₂⊕I₃. After twisting M by χ and
quotienting by invariants one obtains 1⊕χI₃. This records the full
types needed for another convolution, not just the type at y.

Put V_b=Q D_b S. The actual trace of B'_(b,y) is

\[
 v_{b,y}(x)=(V_b)_{xy}-\tau_{b,y},\qquad
 |\tau_{b,y}|\le\sqrt p.                              \tag{5}
\]

The infinity correction is constant in x. Its exact arithmetic
value is unnecessary for this proof.

## 3. The second convolution always has rank six

For any a∈{0,1}, twist B'_(b,y) by χ(x−a) and take middle
extension, obtaining T_(a,b,y). At x=a its local type is quadratic
times a unipotent representation, so its invariant space is zero.
Thus χ(x−a)v_(b,y)(x) is its actual trace at every finite point:
there is no missing punctual term.

The finite invariant codimensions of T_(a,b,y) are as follows.

| Point | a=b | a≠b |
|---|---|---|
| a | 4, from χJ₃⊕χ | 4, from χJ₂⊕χJ₂ |
| 1−a | 2, from J₂⊕J₂ | 2, from J₃⊕1 |
| y | 1, from χ⊕1³ | 1, from χ⊕1³ |

The sum is seven in both cases. At infinity the type is χ⊕1³.
After the convolution twist it is 1⊕χI₃, with exactly one
invariant. Therefore

\[
 C_{a,b,y}=MC_\chi(T_{a,b,y})
 \quad\text{has rank }7-1=6.                          \tag{6}
\]

Every object used as input is geometrically irreducible of rank
at least two. None is a punctual, rank-one Kummer or
Artin–Schreier exception. The preservation theorem therefore makes
C_(a,b,y) geometrically irreducible. It is tame and pure of weight
three, with possible singularities only at 0,1,y,infinity.
No restriction rank<p is required, including when p=5.

For completeness, its infinity auxiliary representation is
χ⊕J₂⊕J₂⊕J₂. After twisting and removing invariants the output
type is χJ₂⊕χJ₂⊕χJ₂. At y the output type is J₂⊕1⁴.
These are consistent with rank six.

Let w_(a,b,y)(x) be its actual middle-extension trace. The compact
versus middle correction at infinity is a scalar κ_(a,b,y),
constant in x, with |κ_(a,b,y)|≤p because T has weight two.
Using S(χ(x−a))=p e_a−1 gives exactly

\[
 \boxed{w_{a,b,y}(x)=-(L_{ab})_{xy}
       +\tau_{b,y}(p\,1_{x=a}-1)-\kappa_{a,b,y}.}      \tag{7}
\]

At all finite x, |w_(a,b,y)(x)|≤6p^(3/2). Formula (7) keeps both
infinity corrections and the finite spike at x=a.

## 4. Rank mismatch excludes the invariant

On U_y=P¹\\{0,1,y,infinity}, consider

\[
 E_y\otimes C_{a,b,y}\otimes\mathcal L_{\chi(x-c)},
 \qquad c\in\{0,1,y\}.
\]

It is tame, of rank twelve and pure weight four. An invariant
would give a nonzero map from the irreducible rank-two E_y dual
to the irreducible rank-six quadratic twist of C_(a,b,y); a
trivial quotient would give a map in the reverse direction.
Schur's lemma forbids both maps. Hence H_c²=0;
also H_c⁰=0. The tame Euler characteristic gives dim H_c¹=24.
Consequently

\[
 \left|\sum_{x\notin\{0,1,y\}}
   \chi(x-c)e_y(x)w_{a,b,y}(x)\right|\le24p^{5/2}.     \tag{8}
\]

All ramification of the final Kummer factor is among the same four
punctures. There is no extra singularity depending on c.

## 5. Restore the omitted coordinates and constant terms

Suppress the subscripts on e,w,τ,κ, and put h_x=χ(x−c).
Equations (4) and (7), Q=−e−1, and Σ_x h_x=0 give

\[
 \begin{split}
 J^{a,b,c}_y={}&\sum_x h_x e(x)w(x)
       +(\tau+\kappa)\sum_xh_xe(x)+\sum_xh_xw(x)\\
       &-\tau p\,h_a(e(a)+1).                        \tag{9}
 \end{split}
\]

Only two of the three omitted points need restoring in (8), since
h_c=0. The bounds |e|≤2√p and |w|≤6p^(3/2) cost at most 24p².
The two linear terms cost at most

\[
 (p+\sqrt p)\,2p^{3/2}=2p^{5/2}+2p^2,\qquad
 6p^{5/2},
\]

respectively. At the pole a∈{0,1}, (4) gives |e(a)+1|≤2, so
the last term costs at most 2p^(3/2). Combining these estimates
with (8) proves (3). In particular the spike is not placed under
an uncorrected rank-six pointwise bound.

## 6. The three remaining classes

Set W_C=S∘Q and t₄=tr((D_0S)^4)/(p−1), with |t₄|≤p².
The adjacent-C graph normalization from
[pass six](parallel6-necklace-2026-09-04.md) gives

\[
 \begin{split}
 N(ABBACC)
 &=\sum_{x,y}\chi(x)\chi(y-1)Q_{xy}(S W_C S)_{xy}
       -p t_4+2\\
 &=\sum_y J^{0,1,y}_y-p t_4+2.                       \tag{10}
 \end{split}
\]

For the second equality, move the symmetric S matrices through
the trace: the first contraction is
tr(D_0 Q D_1 S W_C S)=tr((S D_0 Q D_1 S)W_C).
No extra graph normalization is used in the next two cases.
Directly expanding their six matrix factors gives

\[
 \begin{split}
 N(ABCBAC)
 &=\sum_{x,y}\chi(x)\chi(y)Q_{xy}(L_{11})_{xy}
   =\sum_y\chi(y)J^{1,1,0}_y,\\
 N(ABCABC)
 &=\sum_{x,y}\chi(x)\chi(y-1)Q_{xy}(L_{10})_{xy}
   =\sum_y\chi(y-1)J^{1,0,0}_y.                      \tag{11}
 \end{split}
\]

The absolute sum of the p−2 generic fibers in each formula is at
most (p−2)(32p^(5/2)+26p²+2p^(3/2)). The two exceptional fibers
cost at most 2p³, including when an outer character makes one zero.
Thus the direct contractions in (11) are at most

\[
 32p^{7/2}+28p^3+2p^{5/2}<45p^{7/2}\quad(p\ge5).
\]

Adding p³+2 for (10) gives less than 46p^(7/2). This proves (1).
The comparisons hold already at p=5 and each normalized error
term decreases thereafter.

Rotations, reversal and exchange of A,B are exact symmetries.
The [projective exchange of B,C](parallel-necklace-2026-09-04.md)
costs at most 4p³ on a balanced word. Every label permutation uses
at most one such exchange with exact reflections on either side.
Since 4/√p<2 for p≥5, the 18-word classes represented by AABCCB
and ABACBC have bounds 48p^(7/2) and 47p^(7/2). The six words
in the ABCABC class are already its dihedral orbit and retain
45p^(7/2).

## 7. Sources, exact checks and scope

The only cohomological inputs are those already audited in pass seven:
[Katz, Rigid Local Systems](https://web.math.princeton.edu/~nmk/wholebookRLScorr.pdf),
Lemma 2.9.4, Theorems 2.9.7 and 3.3.3, Corollaries 3.3.6 and 3.3.7,
and the weight argument in Section 5.5.5.10. The archived PDF SHA256 is
ca6eb8d5d9e21076dda4e154e83dfa2821f586d6ccbe8c8f1b371b4352f753d1.
Its complete relevant pages were inspected in the preceding audit.
The rank-twelve cohomology estimate uses the same curve trace formula,
tame Euler characteristic and weight bound as the earlier proofs.

The [verifier](../experiments/parallel8_necklace_2026_09_05.py)
computes the original matrix sums and all 42 new orbit values.
It checks the finite and infinity Jordan-block transformations for
all four (a,b) pairs, formal coefficients in (9), original graph
normalization and small literal necklace sums. It does not guess κ:
its finite trace-bound checks allow every τ,κ satisfying the modulus
bounds |τ|≤√p and |κ|≤p, using exact comparisons in Q(√p).
The [results](../results/parallel8_necklace_2026_09_05.json)
record source hashes and exact counts.

These computations check formulas and finite examples, not geometric
irreducibility, purity, or the uniform theorem. Those are supplied by
the written argument and its stated sources. The finite word catalog
does not substitute for the growing-depth aggregate required by the
Paley problem.
