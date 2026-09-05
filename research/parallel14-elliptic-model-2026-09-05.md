# Three-anchor localization on an elliptic group

Status: exact representation and a source-dependent uniform scalar Fourier
bound. No new asymptotic operator bound, clique bound, or Paley proof.
The elementary identities below have been checked locally. The application
of curve cohomology needs independent mathematical review. No novelty
claim is made for the elliptic or character-sum ingredients.

## 1. Setup and exact covering

Let p be prime, p ≡ 1 mod 4, and let χ be its quadratic character,
extended by χ(0)=0. A three-vertex clique can be normalized by an affine
map with square multiplier to A={0,1,r}, where χ(r)=χ(r−1)=1. Put

\[
C=\{x: \chi(x)=\chi(x-1)=\chi(x-r)=1\},\qquad m=|C|,
\]
\[
E_r:y^2=x(x-1)(x-r),\qquad G=E_r(\mathbb F_p).
\]

Write O for the identity and G[2]={O,(0,0),(1,0),(r,0)}. Define

\[
\rho(P)=x(2P)\in\mathbb P^1(\mathbb F_p),\quad x(O)=\infty,
\qquad f(O)=0,\quad f(x,y)=\chi(y).
\tag{1}
\]

For P=(t,y), y≠0, the duplication formulas give

\[
x(2P)=\frac{(t^2-r)^2}{4y^2},\quad
x(2P)-1=\frac{(t^2-2t+r)^2}{4y^2},\quad
x(2P)-r=\frac{(t^2-2rt+r)^2}{4y^2}.                 \tag{2}
\]

Thus a value outside {0,1,r,∞} belongs to C. Conversely, every point of
E with x-coordinate in C is divisible by 2 in G. This uses the usual
rational-halving criterion: all three differences from the roots of the
cubic are squares. The precise statement, valid in any characteristic
other than 2 and including zero squares, is Theorem 2.1 of
[Bekker and Zarhin](https://arxiv.org/html/1702.02255v2#S2.Thmthm1).
It also makes every point of G[2] divisible by 2 under our hypotheses.
Since the kernel of doubling has four elements, all 16 geometric
4-torsion points are rational. Therefore

\[
\rho(G)=C\cup\{0,1,r,\infty\},\qquad
|\rho^{-1}(x)|=\begin{cases}8&x\in C,\\4&x\in\{0,1,r,\infty\},\end{cases}
\quad |G|=8m+16.                                      \tag{3}
\]

The preimage of the four exceptional values is exactly G[4]. For an
additional explicit inverse check, choose any of the eight triples
a²=x, b²=x−1, c²=x−r for x∈C, and set

\[
D=(1-r)a+rb-c,\qquad y=\frac{r(r-1)}D,\qquad t=r+y(a-b). \tag{4}
\]

In fact D=(a−b)(a−c)(b−c), so D is nonzero since the three squares
are distinct. Multiplying by the corresponding sums in (4) gives
t=(a+b)(a+c) and y=(a+b)(a+c)(b+c). Consequently
t−1=(a+b)(b+c), t−r=(a+c)(b+c), and y²=t(t−1)(t−r).
The three square roots on the right sides of (2), with denominator
2y, are a,b,c respectively. For example, their first two values a′,b′
satisfy a′−b′=(t−r)/y=a−b and a′²−b′²=1; hence a′+b′=a+b.
The analogous first/third difference gives c′=c.
This recovers all eight points above x; the ratios recover the triple,
so different triples give different points. The verifier checks this
inverse separately from the group-law calculation.

## 2. A complete sum-and-difference kernel identity

Define the projective character matrix K by

\[
K(u,v)=\chi(u-v)\quad(u,v\text{ finite}),\qquad
K(u,\infty)=K(\infty,u)=1\quad(u\text{ finite}),\qquad
K(\infty,\infty)=0.
\]

For all P,Q∈G, with no exceptions omitted,

\[
M_{P,Q}:=K(\rho(P),\rho(Q))=f(P+Q)f(P-Q).             \tag{5}
\]

To prove it, first note that f is even because χ(−1)=1. It is also
invariant under translation by G[2]. Indeed, for T_e=(e,0) and
{e,e_j,e_k}={0,1,r}, the addition formulas give

\[
x(P+T_e)=e+\frac{(e-e_j)(e-e_k)}{x(P)-e},\quad
y(P+T_e)=-\frac{(e-e_j)(e-e_k)y(P)}{(x(P)-e)^2}.       \tag{6}
\]

The multiplier of y is a square whenever it is defined; all pairwise
anchor differences and −1 are squares. The remaining points lie in
G[2], where f vanishes, or are covered by the same nonsingular formula.

If ρ(P)=ρ(Q), then 2P=±2Q, so P−Q or P+Q lies in G[2]. Both sides
of (5) vanish. If precisely one image is ∞, its point lies in G[2]
and the other does not. Translation invariance and evenness reduce the
right side to f(Q)²=1, or f(P)²=1. If both images are ∞, it is zero.

For distinct finite images put U=P+Q and V=P−Q. Neither U nor V
lies in G[2], and x(U)≠x(V): equality would give U=±V, forcing
2P=O or 2Q=O. The usual addition formulas now give

\[
x(U+V)-x(U-V)=\frac{-4y(U)y(V)}{(x(U)-x(V))^2}.
\]

Since U+V=2P, U−V=2Q, and χ(−4)=1, taking characters proves (5).

## 3. Exact exceptional block and centered compression

Let S_C=(χ(x−y)) indexed by C. Use the normalized indicators of all
fibers of ρ: entries 1/√8 on ordinary fibers and 1/2 on exceptional
fibers. On their span the matrix M is

\[
\begin{pmatrix}
8S_C&4\sqrt2\,\mathbf1_C\mathbf1_4^T\\
4\sqrt2\,\mathbf1_4\mathbf1_C^T&4(J_4-I_4)
\end{pmatrix}.
\]

Changing the four exceptional coordinates to their normalized constant
direction and its orthogonal complement gives exactly

\[
\begin{pmatrix}
8S_C&8\sqrt2\,\mathbf1_C\\
8\sqrt2\,\mathbf1_C^T&12
\end{pmatrix}\ \oplus\ (-4I_3).                     \tag{7}
\]

M vanishes on the orthogonal complement of all fiber indicators. The
coupling in (7) has norm 8√(2m), of the same square-root order as the
spectral target when m is proportional to p. It cannot be discarded.

For the actual localized operator remove G[4], writing G₀=G\G[4].
Let U:ℝ^C→ℝ^G₀ put 1/√8 on each ordinary fiber. Then

\[
U^TU=I,\quad U^TM_{G_0}U=8S_C,\quad U^TJ_{G_0}U=8J_C. \tag{8}
\]

With the prior Paley projection P=(I−J/p+S/√p)/2, its principal
compression to C is therefore

\[
P_C=\tfrac12 I+
 U^T\left(\frac{M_{G_0}}{16\sqrt p}-\frac{J_{G_0}}{16p}\right)U.
\tag{9}
\]

For three anchors, the normalized deviation
Z=(P_C−I/2)/(√7/8) is unitarily equivalent to the restriction to the
fiber-constant subspace of

\[
\frac{M_{G_0}-J_{G_0}/\sqrt p}{2\sqrt{7p}}.          \tag{10}
\]

The operator in (10) vanishes on the complementary subspace. Thus its
norm is exactly ‖Z‖ (also with the empty-space convention). This is
an exact restatement of the missing edge estimate for three anchors;
the central J term has not been replaced by an estimate or dropped.

## 4. Four invariant character blocks

The maps

\[
\tau_0(x)=r/x,\qquad \tau_1(x)=(x-r)/(x-1),\qquad
\tau_r(x)=r(x-1)/(x-r)                               \tag{11}
\]

permute C and, with the identity, form a Klein four group. They are
the x-coordinate actions of translation by the three nonzero
2-torsion points on 2P. Alternatively the fractional-linear difference
formula proves directly that each preserves χ(x−y): its determinant
is a square and its denominator is a square at every x∈C. They
therefore preserve S_C and J_C.

For a character ε=(1,ε₀,ε₁,ε_r), ε_r=ε₀ε₁, its orthogonal projector
is (I+ε₀T₀+ε₁T₁+ε_rT_r)/4. If t_i counts fixed points of τ_i in C,
the dimension is

\[
d_\varepsilon=\frac{m+\varepsilon_0t_0+\varepsilon_1t_1+
\varepsilon_rt_r}{4}.                               \tag{12}
\]

Each t_i≤2, by its quadratic fixed-point equation. No freeness
assumption is required. The J term acts only in the trivial-character
block. An extreme eigenvalue can occur in any of the four blocks.

For the existing finite certificate p=257, r=62, the action is free.
The character (1,−1,1,−1) has eight orbit representatives

\[
(2,16,18,26,30,60,123,141).
\]

On normalized orbit vectors with values ε(g)/2, the corresponding
block B is

```
 1 -2 -2  0 -2  0 -2  2
-2  1  2 -2  2  0  2  0
-2  2  1  2 -4  4  2  2
 0 -2  2 -3  2  0 -4  0
-2  2 -4  2 -1 -2  0  4
 0  0  4  0 -2 -1 -4  2
-2  2  2 -4  0 -4 -1 -2
 2  0  2  0  4  2 -2 -1
```

Its entries are B_ij=Σ_g ε(g)χ(r_i−g(r_j)). For
w=(−1,1,−3,3,−3,1,2,2), wᵀw=38 and wᵀBw=−406. The corresponding
Rayleigh quotient for Z is −812/(19√1799)<−1, certified by
812²−19²·1799=9905>0. This reproduces the earlier 32-coordinate
certificate exactly. It is a finite outlier, not an asymptotic
counterexample. See the [previous obstruction note](parallel12-bootstrap-obstruction-2026-09-05.md).

## 5. Uniform scalar Fourier bound

Translation invariance in (6) lets f descend to an even real function
f̄ on the finite abelian group H=G/G[2]. It vanishes only at 0 and is
±1 elsewhere. Put L=|H|=2m+4. Since all G[4] is rational, |H[2]|=4.
For each character η of H define

\[
F_\eta=\sum_{h\in H}\bar f(h)\overline{\eta(h)}\in\mathbb R.
\]

Then the following bound holds uniformly in η, including the trivial
character:

\[
|F_\eta|\le\sqrt p.                                 \tag{13}
\]

Proof using standard curve cohomology. Lift η to a character of G.
The full sum Σ_{P∈G} f(P)η̄(P) equals 4F_η. On the open curve
V=E\E[2], f is the trace of the rank-one Kummer sheaf associated
with y and the quadratic character. The divisor of y is

\[
\operatorname{div}(y)=(0,0)+(1,0)+(r,0)-3O.
\]

Thus it has nontrivial tame quadratic local monodromy at all four
removed points. The Lang isogeny 1−Frob_p:E→E is an étale G-torsor;
pushing out by η̄ gives a rank-one sheaf, lisse on the complete curve
E, whose trace at a rational point is η̄(P). It has finite-order
monodromy and weight zero. Tensoring with it does not change any
of the four nontrivial local monodromies.

The tensor is consequently geometrically nontrivial, lisse of rank
one on V, tame at the four punctures, and pure of weight zero.
Its H_c^0 vanishes since V is an affine connected curve; H_c^2
vanishes by absence of geometric coinvariants. The tame
Euler characteristic is 2−2·1−4=−4, so dim H_c^1=4. The trace
formula and the weight-at-most-one bound give |4F_η|≤4√p.
All constants are independent of the order of η: the Lang twist is
unramified on E. This proves (13).

The exact primary inputs are in
[Katz, Gauss Sums, Kloosterman Sums, and Monodromy Groups](../sources/katz-gauss-kloosterman-monodromy.pdf):
§§2.3.1–2.3.3, printed pp.32–33 (Euler characteristic and trace),
§3.6, printed pp.40–41 (H_c^1 weights), and §4.3, printed pp.59–60
(Lang-character sheaves on any smooth connected commutative group,
their trace convention, and tame Kummer sheaves). Complete PDF
spreads 21,25,34,35 were rendered and visually inspected in this pass.
The proof above applies those statements to this elliptic curve; the
source itself is not being represented as stating our final formula.

## 6. Fourier matrix: four products, still not diagonal

On H put M̄_{s,t}=f̄(s+t)f̄(s−t), and e_α(s)=α(s)/√L. The exact
Fourier-basis formula is

\[
\langle e_\alpha,\bar M e_\beta\rangle
=\frac1L\sum_{\gamma^2=\alpha\beta^{-1}}
              F_\gamma F_{\gamma\beta}.              \tag{14}
\]

There are exactly four terms if αβ⁻¹ is a square character, and zero
otherwise. To derive (14), insert f̄(h)=L⁻¹Σ_γ F_γγ(h) twice.
Summation over s,t imposes γδ=α and γδ⁻¹β=1; equivalently
δ=γβ and γ²=αβ⁻¹. Character orthogonality gives the stated factor.
There are four square roots because |H/2H|=|H[2]|=4.
The four cosets of the square-character subgroup give another
description of the four symmetry blocks.

Parseval also gives Σ_η F_η²=L(L−1). These identities and (13) do
not by themselves establish the needed norm bound for (10). In
particular, (14) has nonzero off-diagonal entries. Reading its
products as eigenvalues would be an invalid diagonalization.
Any future use must control their correlations and the restriction
away from G[4], with the centered term retained.

The ordinary Weil-representation gates involving additive quadratic
phases have not been identified with our quadratic-character
multipliers. No theorem about those gates has been imported as an
operator-norm result here.

## 7. Verification and remaining obligation

The [exact verifier](../experiments/parallel14_elliptic_model_2026_09_05.py)
enumerates points and the group law independently of the matrix
formula. It checks every kernel entry in its retained cases, the
explicit square-root inverse, all exceptional fibers, centered
compression, fixed-point character projectors, the reduced finite
certificate, and the Fourier transform over exact cyclotomic rings.
For each finite case it also certifies (13) simultaneously for all
characters via positive leading minors of pI−A_f², where A_f is
the real symmetric convolution matrix of f̄. Fourier orthogonality
identifies its eigenvalues with F_η. Finite positivity is additional
evidence; the general proof of (13) is the cohomological one above.
Results and source hashes are in the
[machine-readable report](../results/parallel14_elliptic_model_2026_09_05.json).

The next obligation in this model is an asymptotically sharp norm
bound for the explicitly centered, punctured kernel (10), uniformly
over eligible r as required by the intended spectral application.
This pass has not established it, or generalized an edge bound to
arbitrary anchor depth. The earlier logarithmic energy range in
[pass thirteen](parallel13-principal-budget-2026-09-05.md) is unchanged.
The subgroup target, classical arbitrary-set cancellation, and the
precise official Reed–Solomon prize bridge remain separate open
obligations.
