# A crossing contraction with no invariant summand

**Status: the 36-word balanced crossing class is proved.** For every prime
p≡1 mod 4, with χ(0)=0,

\[
 |N(ABABCC)|\le56p^{7/2}.
 \tag{1}
\]

Every word obtained by cyclic rotation, reversal, or permutation of the
three labels satisfies |N(w)|≤58p^(7/2). This is the class represented by
AABCBC in [pass six](parallel6-necklace-2026-09-04.md). Together with the
earlier results, finite length-six coverage is **687/729**. The remaining
42 words form the AABCCB, ABACBC, and ABCABC classes, of sizes 18,18,6.
No estimate for their contractions, growing-depth signed aggregate, or
localized adjacency spectral edge is claimed.

The new mechanism is an irreducible rank-four middle convolution on each
side of a crossing. Their local monodromies cannot match after the final
quadratic twist. This removes the invariant term that separate operator
norms leave uncontrolled.

## 1. Exact contraction and all coincident fibers

Use the full F_p-coordinate matrices

\[
 S_{xy}=\chi(x-y),\quad D_A=\operatorname{diag}\chi(x),\quad
 D_B=\operatorname{diag}\chi(x-1),\quad D_C=D_A D_B,
 \quad Q=S D_C S.
\]

Set

\[
 U=S D_A Q,\qquad V=Q D_B S,\qquad
 I_y=\sum_x\chi(x-y)U_{xy}V_{xy}.
\]

The adjacent-C normalization from pass six gives the exact identity

\[
 N(ABABCC)=\sum_y I_y-p t_4+2,
 \qquad t_4=\operatorname{tr}((D_A S)^4)/(p-1).
 \tag{2}
\]

For completeness, the underlying graph has six cycle vertices, two
anchors, and fourteen edges. Averaging distinct anchors contributes
p(p−1)N; coincident anchors contribute pN₀, where the collapsed C mask is
P₀=I−e₀e₀ᵀ, and

\[
 N_0=\operatorname{tr}((D_A S)^4(P_0S)^2)
     =(p-1)(p t_4-2).
\]

Normalizing the adjacent C vertices to 1,0 has no equal-value fiber,
because their joining edge has character zero. Its remaining four-vertex
path has pairing ABAB, and its sum is Σ_y I_y. This proves (2), including
the +2 term. In particular, no zero-coordinate or coincident-anchor
condition is silently imposed on U,V.

## 2. Imported middle-convolution input

We use quadratic additive middle convolution MCχ over F_p, with the
arithmetic normalization in which its trace is minus the additive
convolution of the input traces, minus the trace of the surviving
infinity-invariant space. Middle extension is taken at every finite
singularity. The source is
[Katz, *Rigid Local Systems*, corrected manuscript](https://web.math.princeton.edu/~nmk/wholebookRLScorr.pdf):

* Lemma 2.9.4, PDF page 63 (chapter 2, page 26): the exact constant
  infinity correction between compact and middle convolution.
* Theorem 2.9.7, PDF page 67 (chapter 2, page 30), and Theorem 3.3.3,
  PDF page 104: preservation of the relevant irreducible objects.
* Corollary 3.3.6, PDF pages 107–109 (chapter 3, pages 14–16), and
  Corollary 3.3.7 immediately after it: tame local monodromy and rank rules.
* Section 5.5.5.10, PDF page 139: the middle image has weight w+1
  when the input is pure of weight w.

The PDF fetched for this audit has 1,115,958 bytes and SHA256
`ca6eb8d5d9e21076dda4e154e83dfa2821f586d6ccbe8c8f1b371b4352f753d1`.
Page numbers above are one-based PDF indices; the manuscript has its own
chapter-local numbering.

Only rank-two, geometrically irreducible, everywhere tame inputs are
used. They have finite singularities and are neither punctual objects,
rank-one Kummer translates, nor translation-invariant Artin–Schreier
objects. Thus the irreducibility theorem applies, including p=5.
No hypothesis rank<p is being inferred or added.

## 3. The two rank-two inputs, with actual stalks

Fix y∈F_p\{0,1}. Write J₂ for a nontrivial unipotent Jordan block and
χJ₂ for its quadratic twist. Let E be the Legendre rank-two sheaf with
trace

\[
 -\sum_z\chi(z(z-1)(z-t)).
\]

It is pure of weight one, geometrically irreducible, and has local types
J₂,J₂,χJ₂ at t=0,1,∞, respectively. This is the same Legendre input
audited in [pass two](parallel2-necklace-2026-09-04.md).

First let E_y have trace −Q(x,y)−1. The exact cross-ratio formula is

\[
 Q(x,y)+1=\chi(x(1-y))
  \sum_z\chi\left(z(z-1)\left(z-
       \frac{y(1-x)}{x(1-y)}\right)\right)\quad(x\ne0).
 \tag{3}
\]

The map in (3) is a degree-one fractional linear map. Its preimages of
0,1,∞ are 1,y,0. The quadratic factor cancels the quadratic monodromy
at x=0 and introduces scalar quadratic monodromy at x=∞. Consequently
E_y has types J₂ at 0,1,y, and **χI₂ at infinity**, not χJ₂ there.

Let B_y have trace

\[
 -K_B(x,y),\qquad K_B=S D_B S.
\]

The affine identity K_B(x,y)=χ(y−1)E_sum((x−1)/(y−1)) identifies it
with an affine pullback and constant twist of Legendre. Its types are
J₂ at 1,y and χJ₂ at infinity; it is lisse at 0.

Now form the middle-extension inputs

\[
 \mathcal F_y=\mathcal L_{\chi(x)}\otimes E_y,
 \qquad
 \mathcal G_y=\mathcal L_{\chi(x(x-1))}\otimes B_y.
\]

Their local types are

| Point | F_y | G_y |
|---|---|---|
| 0 | χJ₂ | χI₂ |
| 1 | J₂ | χJ₂ |
| y | J₂ | J₂ |
| ∞ | I₂ | χJ₂ |

At 0 for F_y, and at 0,1 for G_y, the displayed quadratic local types
have no invariants. Thus multiplying the finite trace functions by the
zero-valued character masks gives their **actual middle-extension
stalks**. No missing punctual correction is hidden in these inputs.

## 4. Ranks, purity, and the infinity correction

Let A_y=MCχ(F_y), B'_y=MCχ(G_y). They are geometrically irreducible
and pure of weight two on their common lisse open set

\[
 \mathcal U_y=\mathbf P^1\setminus\{0,1,y,\infty\}.
\]

Here is a direct use of Katz's rank rule, avoiding a guessed rank.
For F_y the finite codimensions of invariants are 2,1,1, totaling 4.
The auxiliary infinity representation M has dimension 4 and
M/M^I=I₂, so M=J₂⊕J₂. After its quadratic twist it has no invariants;
the output rank is 4.

For G_y the finite codimensions are 2,2,1, totaling 5. Its auxiliary
infinity representation is χJ₂⊕I₃. After twisting it becomes
J₂⊕χI₃, with one invariant; the output rank is again 4.

At x=y, each input has J₂. The finite local rule identifies the output
modulo its invariant subspace with χ. Hence each rank-four output has
the full local type

\[
 A_y(y)\cong B'_y(y)\cong\chi\oplus1\oplus1\oplus1.
 \tag{4}
\]

Denote their actual middle-extension traces by u_y(x),v_y(x). The
compact/middle exact sequence and Σ_z χ(x−z)χ(z)=pδ₀(x)−1 give

\[
 \boxed{u_y(x)=U_{xy}+p\delta_0(x)-1,\qquad
        v_y(x)=V_{xy}-\tau_y.}
 \tag{5}
\]

For F_y the infinity-invariant space after the convolution twist is
zero. For G_y it is one-dimensional; τ_y denotes its Frobenius trace,
and |τ_y|≤√p by the weight-one bound. It is constant in x. At all finite
points, including the singularities, |u_y(x)|,|v_y(x)|≤4p.

One can determine τ_y exactly here: **τ_y=1**. In the coordinate
z=1/x, after the convolution's quadratic twist the relevant elliptic
family is

\[
 Y^2=(1-zt)(t-1)(t-y),
\]

up to the factor 1−z, whose character at z=0 is 1. The plane cubic model is
Y²Z=(Z−zT)(T−Z)(T−yZ). At z=0 it is the union of the line Z=0
and the smooth conic Y²=(T−Z)(T−yZ). They meet transversely at the
two rational points [T:Y:Z]=[1:±1:0], and the total surface is smooth
there. This is a split semistable fiber. Its dual graph has a single
cycle, fixed by Frobenius; the invariant part of H¹ therefore has
Frobenius trace 1. This holds over every finite extension as well.
The numerical bound
below deliberately uses only |τ_y|≤√p, so it does not depend on the
sharper arithmetic evaluation.

## 5. Excluding the equal-rank invariant

Consider on U_y the rank-sixteen sheaf

\[
 \mathcal H_y=A_y\otimes B'_y\otimes\mathcal L_{\chi(x-y)}.
\]

It is pure of weight four and tame at the four omitted points. If H_c²
were nonzero, irreducibility of both factors would imply

\[
 A_y^\vee\cong B'_y\otimes\mathcal L_{\chi(x-y)}
\]

geometrically. At y the left side has one quadratic and three trivial
eigenvalues. The right side has three quadratic and one trivial
eigenvalue. These local types disagree. Thus **H_c²=0**, despite the
two factors having equal rank. H_c⁰=0 on the open curve, and the tame
Euler characteristic gives dim H_c¹=32. Therefore

\[
 \left|\sum_{x\notin\{0,1,y\}}
       \chi(x-y)u_y(x)v_y(x)\right|\le32p^{5/2}.
 \tag{6}
\]

This argument gives a reusable criterion: for irreducible equal-rank
factors, unequal quadratic multiplicities after a specified local twist
exclude a global invariant. For local type χ⊕1^(r−1), the same argument
works whenever r≠2. The rank-two exception must not be discarded.

## 6. Restore every remaining term

Put a_x=χ(x−y). Formula (5), with Σ_x a_x=0, gives the exact expansion

\[
 I_y=\sum_x a_xu_y(x)v_y(x)
     +\tau_y\sum_x a_xu_y(x)+\sum_xa_xv_y(x)
     -p\chi(-y)(v_y(0)+\tau_y).
 \tag{7}
\]

Restore x=0,1 to (6), at cost at most 32p²; x=y contributes zero.
The two linear sums cost at most 4p^(5/2) and 4p², respectively.
The final pole term costs at most 4p²+p^(3/2). Hence, for y≠0,1,

\[
 |I_y|\le36p^{5/2}+40p^2+p^{3/2}.
 \tag{8}
\]

For the exceptional y=0,1, use the full matrices directly:
||U||,||V||≤p^(3/2), so |I_y|≤p³ by Cauchy–Schwarz. Also |t₄|≤p².
Combining (2) and (8) yields

\[
 |N(ABABCC)|\le
 (p-2)(36p^{5/2}+40p^2+p^{3/2})+3p^3+2
 \le36p^{7/2}+43p^3+p^{5/2}+2
 <56p^{7/2},
\]

where the last inequality holds for all p≥5. This proves (1).

Rotations and reversal preserve N exactly. The reflection exchanging
A,B is exact, while the projective exchange of B,C changes a balanced
word's full-field value by at most 4p³, as proved in
[pass one](parallel-necklace-2026-09-04.md). Every label permutation
uses at most one such exchange, with exact A,B reflections on either
side. Thus the full 36-word orbit has bound 58p^(7/2).

## Verification and limits

The [verifier](../experiments/parallel7_necklace_2026_09_04.py) retains the
original matrices and graph sums. It checks the quartic cross-ratio
identity, the actual input masks, compact-convolution normalization,
the full expansion (7), omitted-point corrections, pointwise trace and
inner-sum bounds, exceptional y=0,1, the graph identity (2), and all
36 orbit values on five small primes. It independently computes the
two corrected traces over F_(5^n), n=1,2,3,4, at six generic base-field
points. Newton identities recover degree-four Frobenius polynomials,
and exact reciprocal-polynomial inequalities check that every root has
absolute value 5. This audits signs and the infinity constant; it is
not a numerical proof of sheaf irreducibility or of the uniform theorem.

[Results](../results/parallel7_necklace_2026_09_04.json) preserve the exact
counts, values, source hashes, and the 42 unproved words. The nested
ABBA contraction and the two classes without adjacent equal labels
remain unresolved. No growing-family or full Paley conclusion follows
from this one additional crossing class.
