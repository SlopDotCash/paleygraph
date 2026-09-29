# A projection obstruction to extending the short-depth estimate

**Status: an elementary obstruction, with exact finite checks. No improved
asymptotic Paley bound is proved in this pass.** The short-depth estimate from
the [eleventh pass](parallel11-full-energy-2026-09-05.md) remains dependent on
the earlier source-dependent arguments, whose independent review is outstanding.

The question here is whether the projection identities, exact anchor data,
and the established short-depth bounds suffice to control arbitrarily long
powers. They do not. A rank-two modification preserves those inputs while
planting an extreme localized eigenvalue. It fails the exact character-entry
condition away from the anchors, so this is not a Paley counterexample.

There is a second, separate result: an actual Paley localization at p=257
has an exact rational Rayleigh certificate outside the proposed limiting
interval. This rules out an interval bound at every finite prime, while
remaining compatible with asymptotic spectral-edge convergence.

## 1. Preserve the anchors and plant an eigenvalue

Let p≡1 mod 4 be prime, S_xy=χ(x−y), J=11ᵀ, and

\[
 P=\tfrac12(I-J/p+S/\sqrt p).
\]

Then P is an orthogonal projection of rank (p−1)/2 and P1=0. Let A₀
be an a-anchor clique and C its common neighborhood, excluding the anchors.
Assume |C|≥2. Choose distinct x,y∈C and write

\[
 v=e_x-e_y,\quad s=v^Tv=2,\quad w=Pv,\quad
 q=v^TPv=1-\chi(x-y)/\sqrt p.
\]

Thus 0<q<s, wᵀw=wᵀv=q, and both v,w are perpendicular to 1.
For every z∈A₀,

\[
 \langle v,e_z\rangle=0,\qquad
 \langle v,Pe_z\rangle=P_{xz}-P_{yz}=0.
\]

The second equality follows from S_xz=S_yz=1. Since P²=P,
w is also perpendicular to e_z and Pe_z.

To plant a zero eigenvalue on v, put

\[
 u=sw-qv,\qquad
 P^-=P-\frac{ww^T}{q}
             +\frac{uu^T}{s q(s-q)}.                 \tag{1}
\]

Indeed, uᵀv=0 and uᵀu=sq(s−q). Every r∈range(P)∩w⊥ satisfies
⟨r,v⟩=⟨r,Pv⟩=0, and hence ⟨r,u⟩=0. Thus (1) removes the
unit direction w/√q from range(P) and replaces it with the unit direction
u/√[sq(s−q)]. It follows that

\[
 (P^-)^2=P^-,\quad (P^-)^T=P^-,\quad
 \operatorname{rank}P^-=(p-1)/2,\quad P^-v=0.
\]

The update has rank at most two, P⁻1=0, and P⁻e_z=Pe_z at every
anchor. A one-eigenvalue variant is

\[
 P^+=P-ww^T/q+vv^T/s,\qquad P^+v=v,                 \tag{2}
\]

with the same rank, symmetry, null constant direction, and anchor columns.

Define, for either sign,

\[
 S'=\sqrt p\,(2P'-I+J/p).
\]

It is real symmetric and obeys

\[
 S'1=0,\quad (S')^2=pI-J,\quad \|S'\|=\sqrt p,
 \quad\operatorname{rank}(S'-S)\le2,\quad S'e_z=Se_z.
                                                               \tag{3}
\]

For the square identity, use P'J=JP'=0 and J²=pJ. The exact
anchor columns preserve the anchor clique and the same common-neighborhood
mask E=diag(1_C). They also preserve every diagonal character-product mask
formed using only A₀, including the exact soft-mask correction at A₀.

The defect that prevents (3) from being a Paley example is explicit. In the
minus construction, S'v=−√p v, so

\[
 S'_{xx}-S'_{xy}=-\sqrt p.                            \tag{4}
\]

A zero diagonal and off-diagonal entry in {−1,1} cannot have this
difference. The plus construction has difference +√p instead. The
construction preserves selected character columns, not the character
matrix away from those columns. In particular, the sheaf descriptions
of the original words are not asserted for S'.

## 2. Why the existing trace estimates cannot detect this modification

Let D₁,…,D_k be the diagonal masks from the preserved anchors; each
has operator norm at most one. Since S and S' each have norm √p,
their difference has nuclear norm at most 4√p: its rank is at most
two and its operator norm is at most 2√p. Telescoping the word and
using |tr(AXB)|≤||A|| ||X||_* ||B|| gives

\[
 \left|\operatorname{tr}(D_1S'\cdots D_kS')
       -\operatorname{tr}(D_1S\cdots D_kS)\right|
       \le4k p^{k/2}.                                \tag{5}
\]

This is 4k/p after the necklace normalization p^(k/2+1), smaller
than the preceding O_(a,k)(p^(−1/2)) fixed-word error. This statement
preserves the order of those estimates, with the explicit added error
(5); it does not assert equality of all earlier numerical constants.

Now put b=2^a, N=b−1, D=bE−I and

\[
 T=(2P-I)D/\sqrt N,\quad T'=(2P'-I)D/\sqrt N,\qquad
 H_j=\|T^j\|_F^2,\quad H'_j=\|(T')^j\|_F^2.
\]

Both T and T' have operator norm at most √N. The update T'−T
has rank at most two, operator norm at most 2√N, nuclear norm
at most 4√N, and Frobenius norm at most 2√(2N). A second
telescoping argument yields

\[
 |\operatorname{tr}((T')^k)-\operatorname{tr}(T^k)|
       \le4kN^{k/2},                                 \tag{6}
\]

\[
 \|(T')^j-T^j\|_F\le2\sqrt2\,jN^{j/2},\qquad
 H'_j\le\bigl(\sqrt{H_j}+2\sqrt2\,jN^{j/2}\bigr)^2.  \tag{7}
\]

Consequently any original bound H_j≤(1+o(1))p survives the update
on a range where j²N^j=o(p). In particular, for fixed a≥2,

\[
 j\le(1/2-\epsilon)\log p/\log\Theta_a,\qquad
 \Theta_a=\frac{[2^{a-1}(a+1)-1]^2}{2^a-1}>N,
\]

the bound from pass eleven gives H'_j≤(1+o(1))p too. The inference
uses that preceding estimate as an input; the projection perturbation
argument itself is elementary and independent of its cohomological proof.
The exact nonasymptotic bound is augmented by the extra term in (7).

## 3. The same construction fails maximally at long depth

The vector v is supported on C, hence Dv=Nv. In the minus case,
(2P'−I)v=−v. Therefore

\[
 T'v=-\sqrt N\,v,\qquad H'_j\ge N^j.                 \tag{8}
\]

For fixed a≥2 and j/log p→∞, this implies

\[
 \frac{\log(H'_j+2p)}j\ge\log N>0.                   \tag{9}
\]

Thus the proposed vanishing long-depth growth condition cannot follow
from the projection identities, exact anchor data, and the established
short-depth estimates alone. It needs further input not preserved by
(1), for example exact quadratic-character identities away from A₀.

This supplies families of countermodels, rather than just a single
small-matrix example. For a=2 and A₀={0,1}, the actual Paley cell has
size (p−5)/4, so it has at least two points for every p≥13 under
consideration. Its size follows directly by expanding
(1+χ(x))(1+χ(x−1))/4, using the quadratic correlation −1, and
subtracting the two anchor contributions of 1/2. Apply (1) at each
prime. The resulting family retains the short-depth bound whenever
the original family does and violates (9) at every such prime.

This is not an obstruction to every arithmetic proof: it isolates the
insufficiency of a specified collection of inputs. It is also not a
counterexample to any asymptotic Paley conjecture.

## 4. An exact outlier in an actual Paley graph at p=257

Here S is the original character matrix, with no modification. Take
A₀={0,1,62}, a clique in F_257, and its complete common neighborhood

```
C = [2,16,18,26,30,31,32,36,58,60,61,73,114,121,122,123,
     124,129,134,135,141,185,190,196,197,198,199,208,227,235,240,249].
```

In that coordinate order, let

```
v = [-1,1,-3,3,-3,1,3,-1,-1,1,-1,-3,-3,-3,3,2,
     -2,2,1,1,2,-3,-2,1,-1,3,-2,2,-1,3,-2,3].
```

All of the following are integer identities:

\[
 v^Tv=152,\qquad \sum v_i=0,\qquad v^TS_{C,C}v=-1624.
\]

For the real symmetric normalized localization

\[
 Z=\frac{P_{C,C}-I_C/2}{\sqrt7/8}
   =\frac4{\sqrt7}\left(S_{C,C}/\sqrt{257}-J_C/257\right),
\]

the Rayleigh quotient is exactly

\[
 \frac{v^TZv}{v^Tv}=-\frac{812}{19\sqrt{1799}}<-1,
 \qquad 812^2-19^2\cdot1799=9905>0.                  \tag{10}
\]

Hence λ_min(Z)<−1 by the real symmetric variational principle.
The integer certificate is slightly stronger:

\[
 812^2\cdot144^2-19^2\cdot1799\cdot145^2
       =17702209>0.
\]

Writing r=812/(19√1799), this proves
r>145/144=(9/8+8/9)/2. The exact two-projection transfer from the
[spectral criterion](parallel2-spectral-transfer-2026-09-04.md) associates
to a localized normalized eigenvalue z a determinant-one block of T
with characteristic polynomial t²−2zt+1. If z≤−r, that block
has spectral radius at least r+√(r²−1)>9/8. Therefore the actual
p=257 matrix satisfies

\[
 \rho(T)>9/8,\qquad H_j>(81/64)^j\quad(j\ge1).        \tag{11}
\]

This excludes any polynomial-in-j energy bound valid at every fixed
prime with no exponential allowance. It does **not** exclude an
asymptotic estimate in which the allowed rate tends to zero as p→∞,
or the limiting edge conjecture. No infinite family of actual Paley
outliers with a fixed nonzero gap is established here.

The relevant published target is asymptotic; see Kunisky,
[*Spectral pseudorandomness and the road to improved clique number bounds
for Paley graphs*](https://arxiv.org/html/2303.16475v1), especially the
spectral conjectures and Section 7. Neither finite certificate changes
that target into an exact assertion at every prime.

## 5. Verification and the next obligation

The [verifier](../experiments/parallel12_bootstrap_obstruction_2026_09_05.py)
constructs both projection endpoints at p=13,17,41, including a three-anchor
case. It checks idempotence, symmetry, unchanged rank and anchor columns,
the constant null direction, the character-square identity, the explicit
entry defect, energy recurrences, and trace/power perturbation bounds.
It separately reconstructs the full actual p=257 cell and checks (10)
and (11) through their integer certificates. Acceptance uses no floating
point arithmetic. Results are in
[the exact output](../results/parallel12_bootstrap_obstruction_2026_09_05.json).

The next proof obligation remains a uniform long-depth bound for the
**actual** character matrix, with an asymptotically vanishing allowance
for edge errors. Reusing only projection identities and the current
short-depth estimates cannot close it. Independent review of passes
nine through eleven, the subgroup uniformity and exceptional-prime gap,
the classical arbitrary-two-set estimate, and the exact official
Reed–Solomon prize bridge remain outstanding.
