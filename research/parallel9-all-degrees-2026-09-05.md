# Growing convolution chains for every fixed necklace degree

**Status: source-dependent proof completed and checked locally.**
Independent mathematical review and formal verification have not been obtained.
The claim below concerns individual necklaces, not the signed spectral
aggregate, a Paley clique bound, the two-set conjecture or the prize.

## Theorem and scope

Let p≡1 mod 4 be prime. Let A⊂F_p contain a≥1 distinct anchors,
with a<p, and let Z₁,…,Z_k be nonempty subsets of A. Put
d_Z(x)=∏_(z∈Z)χ(x−z), S_xy=χ(x−y), and
N=tr(D_Z₁ S⋯D_Z_k S).

For k≥3 the uniform estimate is

\[
 \boxed{|N|\le3a(2a+2)^{k-2}p^{(k+1)/2}.}             \tag{T}
\]

For k=1, N=0. For k=2, |N|≤a²p by S_xy²=1_(x≠y)
and the one-variable squarefree-polynomial character bound.
Thus (T) proves Kunisky's Conjecture 1.14 at every fixed degree:
after dividing by p^(k/2+1), the bound tends to zero at each fixed a,k.
It also gives a quantitative growing-length range, but its exponential
constant does not supply the signed aggregate needed for extreme
eigenvalues. No novelty or formal-verification claim is made.

Here and below, lisse sheaves have coefficients in an algebraic closure
of Q_l for a prime l≠p. Fix an embedding of their coefficient field into
C when taking trace magnitudes. Arithmetic Kummer sheaves use the
normalization with trace χ; all middle convolutions inherit that
normalization. Involution statements are used geometrically, so their
arithmetic Tate twists do not enter rank calculations. The k=2 estimate
uses |Σ_x d_Z(x)|≤(|Z|−1)√p and a diagonal contribution at most p:
((a−1)²+1)p≤a²p.

## 1. The principal middle-convolution chain

Fix y∉A. Use quadratic additive middle convolution MC, and let F₀
be the rank-one Kummer sheaf with trace χ(x−y). For any sequence
of nonempty subsets T₁,…,T_j⊂A, define

\[
 G_i=\operatorname{ME}(F_{i-1}\otimes\mathcal L_{d_{T_i}}),
 \qquad F_i=\operatorname{MC}(G_i).
\]

ME means middle extension after twisting on the lisse open.
All possible finite singularities are A∪{y}; infinity is also allowed.
The local eigencharacters are always trivial or quadratic.

The initial twist G₁ has |T₁|+1 finite quadratic singularities.
It is neither constant, an Artin–Schreier sheaf, nor a Kummer
translate with just one finite singularity. The middle-convolution
irreducibility theorem applies. At y, every F_i has a pseudoreflection:
χ⊕1^(d_i−1) for even i, and J₂⊕1^(d_i−2) for odd i.
In either case the invariant codimension is one.

We next prove d_i strictly increases, so each later twist is an
irreducible object of rank at least two, outside every excluded
rank-one or punctual case. Each F_i is tame and pure of usual
weight i, and F_i[1] is pure of perverse weight i+1.

## 2. A rank-growth lemma

For a tame local representation with only the eigencharacters ±1,
write i⁺,i⁻ for the numbers of trivial and quadratic Jordan blocks,
and j⁺,j⁻ for the numbers of those blocks of size one. Invariants
of the corresponding character twist have dimensions i⁺ or i⁻.

Suppose F=MC(G), with ranks d and e respectively, and let δ=d−e.
The finite quotient rule in Katz gives, at any finite anchor v,

\[
 i_v^+(F)=\delta+i_v^+(G),\qquad
 i_v^-(F)=i_v^+(G)-j_v^+(G).
\]

Therefore t_v(F):=i_v⁺(F)−i_v⁻(F)=δ+j_v⁺(G)≥δ.
At infinity the auxiliary-representation rule similarly gives

\[
 i_\infty^-(F)=\delta+i_\infty^-(G),\qquad
 i_\infty^+(F)=i_\infty^-(G)-j_\infty^-(G).
\]

Thus t_∞(F):=i_∞⁻(F)−i_∞⁺(F)=δ+j_∞⁻(G)≥δ.
These formulas count blocks, not dimensions of the ±1 eigenspaces.

For the first step, with u=|T₁|,

\[
 d_1=u+1-1_{u\ {\rm even}}=2\lceil u/2\rceil,\qquad
 \delta_1=d_1-1\ge1.
\]

Now suppose δ_i=d_i−d_(i−1)>0. Since MC is geometrically an
involution, applying MC without a twist to F_i has rank d_(i−1).
Its rank formula is

\[
 d_{i-1}=a d_i+1-\sum_{v\in A}i_v^+(F_i)-i_\infty^-(F_i).
\]

A twist at a nonempty T⊂A exchanges the two block counts at every
v∈T, and also at infinity exactly when |T| is odd. Define
F(T)=T, together with infinity in the odd case. Comparison with
the untwisted formula gives

\[
 d_{i+1}-d_i=-\delta_i+\sum_{v\in F(T)}t_v(F_i)
       \ge(|F(T)|-1)\delta_i\ge\delta_i.               \tag{R}
\]

The last step uses |F(T)|=|T|+(|T| mod 2)≥2. Induction proves

\[
 d_i\ge i+1,\qquad
 \delta_{i+1}\ge\delta_i\ge1,\qquad
 d_{i+1}\le a d_i+1.                                 \tag{R'}
\]

The inverse and rank formulas are used only on irreducible,
non-punctual objects. The positive computed rank and the local
pseudoreflection rule close the induction; there is no assumption
that a rank is smaller than p.

## 3. Retaining the raw convolution, including zero masks

The matrix chain integrates over every finite field coordinate.
It is not obtained by repeatedly discarding the corrections
between compact and middle convolution.

Let K₀=F₀[1], and let Cχ=Lχ[1] on A¹. For a nonempty T⊂A,
define the functor

\[
 {\cal T}_T(K)=j_{T!}\bigl(j_T^*K\otimes\mathcal L_{d_T}\bigr),
 \quad j_T:\mathbb A^1\setminus T\hookrightarrow\mathbb A^1,
 \qquad \Phi_T(K)={\cal T}_T(K)*_!C_\chi.
\]

Both functors are exact on perverse sheaves. For the first, the
open immersion is affine and quasifinite; for the second, the
Kummer perverse sheaf has Katz's property P. They preserve the
fixed finite singularity set A∪{y} and tameness. The first
preserves an upper weight bound; the second raises it by at most one.

Put K_i=Φ_(T_i)(K_(i−1)). If
R₀(x,y)=χ(x−y) and R_i=S D_(T_i)R_(i−1), the trace formula gives

\[
 \operatorname{Tr}(\operatorname{Frob}_x\mid K_i)
       =(-1)^{i+1}R_i(x,y).                           \tag{F}
\]

This includes x in the finite singularity set. Extension by zero
implements exactly the zero-valued character mask.

There is an exact sequence of mixed perverse sheaves

\[
 0\longrightarrow E_i\longrightarrow K_i
       \longrightarrow F_i[1]\longrightarrow0,
 \qquad E_i\text{ has weights }\le i.                 \tag{W}
\]

For i=0 the kernel is zero. For the induction, twist and extend by
zero the previous sequence. On the top constituent, the surjection
from extension by zero to middle extension has punctual kernel at
T_i. Its stalk is the new invariant space after twisting F_(i−1),
and has weights at most i−1. The already lower-weight kernel also
has weights at most i−1. Compact convolution raises these bounds
to at most i. Finally, the compact-to-middle exact sequence has
constant kernel V[1], where V is the surviving infinity-invariant
space of G_i. It has usual weights at most i−1, hence perverse
weights at most i. The quotient is the pure F_i[1], of weight i+1.
Composing the surjections proves (W).

This is the step that protects the strict weight gap through all
subsequent operations. Finite and infinity corrections need not
have either zero trace or a chosen arithmetic value.

For a concrete nonzero finite correction, take two anchors α,β and
the twists {α},{β},{α}. After the first two convolutions the rank is
four and the local representation at α is χ⊕1³. The third twist
changes this to 1⊕χ³, so extension by zero differs from middle
extension by a one-dimensional punctual kernel at α. Frobenius on
that invariant line is invertible, so its trace is nonzero. Its
arithmetic value is not guessed or used. The third middle-convolution
rank is six; the relevant infinity-invariant dimension is zero.
This example prevents an incorrect identification of raw chains with
the principal pure middle extensions.

## 4. A uniform size bound for all correction constituents

For a perverse sheaf supported on this stratification, let c(K)
be the sum, over geometric simple constituents with multiplicity,
of their generic ranks; assign mass one to a one-dimensional
punctual constituent. It is additive in exact sequences.

For a simple middle extension of rank r, extension by zero at T
after a Kummer twist adds at most |T|r≤ar punctual constituents.
Compact convolution of the new rank-r middle extension has middle
rank at most (a+1)r and constant-kernel rank at most r. A punctual
constituent becomes a rank-one Kummer translate under convolution.
Consequently

\[
 c(\Phi_TK)\le(2a+2)c(K),\qquad c(K_i)\le(2a+2)^i.     \tag{C}
\]

The exceptional simple objects cause no larger contribution:
constant perverse sheaves convolve to zero; a Kummer translate can
give a punctual middle convolution plus a rank-one constant kernel;
punctual objects either are deleted by the mask or become Kummer
translates. Nontrivial Artin–Schreier objects are absent by tameness.
All other simple middle convolutions are irreducible by Katz.
The same case distinction proves preservation of the stratification
and tameness: middle convolution introduces no new finite singular
point, its kernel is constant, and convolution of a punctual piece
at v has its sole finite singularity at v. Exact extensions remain
lisse on U_y and tame. Frobenius may permute geometric constituents;
the mass still bounds the total geometric stalk dimension.

Write c_i=c(K_i). Since F_i has rank d_i, (W) implies c(E_i)=c_i−d_i.
On U_y=P¹\\(A∪{y,infinity}), E_i is a lisse sheaf shifted by one,
whose usual weights are at most i−1. Therefore its trace has
absolute value at most (c_i−d_i)p^((i−1)/2) there.

At finite singularities, the ordinary cohomology of E_i lies in
degrees −1,0 and has weights at most i−1,i respectively. Total
stalk dimension is at most its constituent mass. The top middle
extension has no degree-zero cohomology and stalk weights at
most i. Thus at every finite x

\[
 |R_i(x,y)|\le c_i p^{i/2}.                            \tag{P}
\]

In (P) the generic lower part is better by p^(-1/2); that stronger
bound is used in the final sum.

## 5. Close the necklace

Take j=k−2≥1 and choose the interior labels in the order that
makes R_j(x,y) the k−1-edge path between the endpoint variables.
Then

\[
 N=\sum_y d_{Z_k}(y)\sum_x
       d_{Z_1}(x)\chi(x-y)R_j(x,y).
\]

For y∉A, the main middle-extension contribution on U_y is a
rank-d_j sheaf twisted by a rank-one Kummer sheaf. Since d_j≥2,
it has neither an invariant nor a trivial quotient. It is tame
at the a+2 punctures and of usual weight j. Thus H_c¹ has
dimension a d_j, giving an upper bound
a d_j p^((j+1)/2) for the open main sum.

The open error from E_j costs at most
(c_j−d_j)p^((j+1)/2). At x=y the character vanishes. At least
one anchor is killed by the nonempty endpoint mask Z₁, so at
most a−1 finite omitted points remain; by (P) they cost at most
(a−1)c_j p^(j/2). The whole inner sum is therefore bounded by

\[
 \bigl(c_j+(a-1)d_j+(a-1)c_j/\sqrt p\bigr)
       p^{(j+1)/2}
 \le(2a-1)c_j p^{(j+1)/2}.
\]

There are at most p generic y. At most a−1 exceptional y∈A
survive the nonempty Z_k mask. Each of their inner sums is at
most p^((j+2)/2), by the norm of the original matrix chain and
the norm of the endpoint character vector. Hence

\[
 |N|\le\bigl((2a-1)(2a+2)^j+(a-1)/\sqrt p\bigr)
       p^{(j+3)/2}
 \le3a(2a+2)^{k-2}p^{(k+1)/2}.
\]

This proves (T) using the primary inputs recorded below.

## 6. Consequences and the remaining spectral obstruction

For each fixed a and k, the normalized maximum in
[Kunisky's Conjecture 1.14](https://arxiv.org/html/2303.16475v1#S1.SS5)
is at most 3a(2a+2)^(k−2)/√p for k≥3, uniformly in the anchor
positions. If the union has size less than a, either use its actual
size or enlarge it; both give the displayed bound. The cases k=1,2
were handled above. The published Theorem 1.17 therefore gives the
corresponding weak empirical spectral convergence at each fixed
degree a. This is a consequence through that imported theorem, not
a proof of the minimum-eigenvalue conjecture.

For fixed a and 0<ε<1/2, (T) also gives normalized decay throughout
k≤(1/2−ε)log(p)/log(2a+2), with bound at most
3a(2a+2)^(−2)p^(−ε) when k≥3. The theorem itself allows a,k to vary,
but its coefficient must always be retained.

The [aggregate criterion](parallel2-spectral-transfer-2026-09-04.md)
needs much more. Even for the S-only portion of its trace expansion,
summing absolute values of the (b−1)^k nonempty-label words, b=2^a,
and using (T) gives only

\[
 \frac{3a}{(2a+2)^2}\sqrt p\,
       \bigl((2a+2)\sqrt{b-1}\bigr)^k.
\]

This expression is an upper bound supplied by this method, not a
lower bound on the actual signed aggregate. It does not have the
required subexponential growth in k. Rank-one J terms and the anchor
mask corrections of the full trace also need their own accounting.
Consequently the argument supplies neither that aggregate nor an
improved extreme eigenvalue, clique bound, arbitrary-two-set bound,
uniform subgroup square-root estimate or official prize certificate.

## 7. Primary inputs and verification

The following are imported theorems, rather than numerical claims.
Page indices for PDFs are one-based.

- Katz, [Rigid Local Systems](../sources/katz-rigid-local-systems.pdf),
  corrected manuscript: Section 2.3.3 (PDF 45) gives affine-open
  perverse exactness; Section 2.6.1 (PDF 50) gives exact compact
  convolution with a property-P object. Section 2.9.1 and Lemma
  2.9.4 (PDF 63) give the Kummer property and the exact constant
  infinity kernel. Theorem 2.9.7 (PDF 67), Theorem 3.3.3 (PDF 104),
  Corollaries 3.3.6–3.3.7 (PDF 107,109), and Section 5.5.5.10
  (PDF 139) supply geometric inversion, irreducibility, tame local
  monodromy, rank, and purity of the middle image. All these complete
  pages were visually inspected in this or the preceding audits.
- Beilinson–Bernstein–Deligne,
  [Faisceaux pervers](../sources/bbd-faisceaux-pervers.pdf),
  Sections 5.1.8–5.1.9 (printed 126–127; PDF 66),
  5.1.13–5.1.14 (printed 128–129; PDF 67), and 5.3.1–5.3.2
  (printed 134–135; PDF 70) give the stalk convention for upper
  weights, compact-image and external-product weight bounds,
  stability under perverse subquotients, and purity of intermediate
  extensions. Complete PDF spreads 65,66,67,69,70,71 were visually
  inspected. In particular, a perverse object of weights ≤w has
  ordinary H^q stalk weights ≤w+q; overlooking q=−1 would lose the
  strict saving used here.
- Katz, [Gauss Sums, Kloosterman Sums, and Monodromy Groups](../sources/katz-gauss-kloosterman-monodromy.pdf), printed
  pages 32–33 and 38–41, previously inspected, supplies the tame
  Euler characteristic, trace formula, and H_c¹ weight bound.
- [Kunisky's versioned primary text](../sources/kunisky-2303.16475v1.html),
  Definition 1.13, Conjecture 1.14 and Theorem 1.17, fixes the
  normalization and the weak-convergence implication. The online
  primary text was rechecked in this pass. No literature-priority
  claim is made.

The [exact verifier](../experiments/parallel9_all_degrees_2026_09_05.py)
enumerates every nonempty-label sequence through depth 40 at a=1,
8 at a=2, 5 at a=3, 4 at a=4 and 3 at a=5. It checks 114,510 rank
transitions and inverse local types, 1,691,778 local block identities,
and the stated nonzero finite-invariant example. The all-degree proof
is the symbolic induction (R), not this bounded enumeration.

Original character matrices over five field/anchor configurations
check 1,170 endpoint identities, 980,530 trace-sign entries,
340 literal path entries, 6,821 masked rows, 420 order-one/two
identities and 410 literal necklaces. All arithmetic is integral.
The small-field upper bounds with these exponential constants can be
vacuous and are explicitly not used as evidence of asymptotic decay.

The [results](../results/parallel9_all_degrees_2026_09_05.json)
pin source and proof hashes. Root audited the induction, all correction
weights, exceptional constituents, stalk dimensions and endpoint
normalizations. This is a locally checked written proof dependent on
the cited theorems, without an independent worker review, a Lean
formalization or a claimed solution of the full goal.
