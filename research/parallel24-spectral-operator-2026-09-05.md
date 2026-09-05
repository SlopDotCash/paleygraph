# The next actual Krylov coefficient, with both couplings retained

**Status: the next diagonal and off-diagonal Lanczos coefficients are
evaluated, and a uniform upper estimate is proved on an actual
two-dimensional Krylov subspace. The full spectral edge is unproved.**
The second coupling is also nonvanishing. No assumption that these
directions capture an extremal eigenvector is made.

## 1. Statement

Let p≡1 mod4 be prime, p≥13, and use the actual neighborhood and
compression from [pass 23](parallel23-spectral-coupling-2026-09-05.md):

    C={t:χ(t)=χ(t−1)=1}, m=(p−5)/4,
    L(t)=Σ_y χ[y(y−1)(y−t)], S_C=(χ(x−y))_(x,y∈C),
    B=(I−J/p+S_C/√p)/2, u=1_C/√m.

Here χ(0)=0 and J is the all-ones matrix. Put

    ℓ=m^(−1)Σ_C L, g=L|_C−ℓ1_C, d=Σ_C g²,
    v=g/√d, Π=uuᵀ+vvᵀ,
    q=uᵀBu, b=uᵀBv, a=vᵀBv,
    c=||(I−Π)Bv||.

The quantities involving v are defined whenever d>0. In particular,
d>0 for every such prime p≥1024, by an explicit bound below. Then,
with absolute implied constants, as p→∞ along these primes,

    q=3/8+O(p^(−1/2)), b=1/8+O(p^(−1/2)),
    a=1/2+O(p^(−1/2)), c²=3/64+O(p^(−1/2)).       (1)

The first two coefficients refine/reuse pass 23; the evaluated next
coefficients are a and c. Since Bu=qu+bv, the space
K₂=span{u,v}=span{u,Bu} is the actual second Krylov space. Uniformly
over all unit vectors f∈K₂, the sharp subspace bounds are

    sup fᵀBf = (7+√5)/16+O(p^(−1/2)),               (2)
    ||B|_(K₂)||² = (15+√74)/64+O(p^(−1/2)).         (3)

Equation (3) retains the component of Bf outside K₂. In particular,
these are actual upper estimates, with the usual quantifier that for
every ε>0 the displayed constants plus ε apply to every sufficiently
large prime and every f in this specified subspace.

Deleting the two blocks linking K₂ and its orthogonal complement has
exact error c, and therefore error tending to √3/8:

    ||B−[ΠBΠ+(I−Π)B(I−Π)]||=c→√3/8.              (4)

Every assertion transfers to every anchor edge by residue affine
dilation. None asserts a bound for arbitrary vectors in R^C.

## 2. Exact recurrence with the anchors retained

Work first on all F_p, and set f₀(t)=χ(t(t−1)), D₀=diag(χ(t)),
D₁=diag(χ(t−1)). Let

    k(t)=k₃(t)=(S D₀ L)(t),
    K(t)=k(t)+k(1−t)+k(1−1/t), t∈C.

This k₃ is the earlier kernel S(D₀S)² evaluated at row 1.
The full-field relations S²=pI−J, L=Sf₀, Σf₀=−1 give
SL=pf₀+1. Reflection gives SD₁L(t)=k(1−t).

For t≠0, substituting y=1/z, using L(1/z)=χ(z)L(z), and restoring
the omitted z=0 term gives

    S(D₀D₁L)(t)=χ(t)k(1−1/t)+1.

The added 1 is from L(0)=−1; omitting it changes the recurrence.
Apply the exact mask

    1_C=(1+χ(t)+χ(t−1)+f₀(t))/4−(δ₀+δ₁)/2.

Since L(0)=L(1)=−1, its two anchors add another 1 at every t∈C.
Consequently

    (S_C L|_C)(t)=[p+6+K(t)]/4.                   (5)

Both C and L|_C are invariant under t↦1−t and t↦1/t. The sign
matrix S_C commutes with their permutation actions. Moreover
k(1/t)=χ(t)k(t), as follows from symmetry and homogeneity of the
rank-three kernel. Thus K is the sum over three representatives of
the S₃ orbit; it is invariant on C. If

    s₀=Σ_C k, s=Σ_C gk,

then Σ_C K=3s₀ and Σ_C gK=3s. Short S₃ orbits are included by
these bijections and require no division by orbit size.

The established row identity S_C1=(L−6)/4 and (5) imply

    b²=d/(64mp),
    a=1/2+[3s/d−ℓ]/(8√p),                        (6)
    h=K−(3s₀/m)1−(3s/d)g,
    (I−Π)Bv=h/(8√(pd)),
    c²=[||K||²−9s₀²/m−9s²/d]/(64pd).             (7)

These are exact finite identities, including the normalization of the
two orthogonal directions. In (7), h is perpendicular to both 1 and g.

## 3. The arithmetic estimates and their primary inputs

The old [seeded-kernel result](parallel6-seeded-kernels-2026-09-04.md)
controls specified kernel spans and explicitly does not capture every
spectral sector. Here the new operation is the actual recurrence (5)
and correlations between three Möbius pullbacks of k₃. Their separation
comes from different local Jordan forms, not from an assumption of
independence of sheaves or seed locations.

The primary inputs are [Katz, *G₂ and hypergeometric sheaves*, §2,
printed pp. 3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and [Katz, *Gauss sums, Kloosterman sums, and monodromy groups*,
§2.3.1–2.3.3 and the weight argument in §3.6, printed pp. 32–33,
40](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
Their primary PDF text was inspected in this lane, including the
disjoint-list hypotheses, endpoint monodromy, Euler characteristic,
trace formula, and H_c¹ weight bound. The PDFs were also opened through
the primary publisher/author URLs. No new theorem about their general
cohomological machinery is claimed here.

For r=2,3 take r trivial upstairs and r quadratic downstairs characters.
The lists are disjoint for odd p. Katz supplies a geometrically
irreducible rank-r sheaf. Twisting by the constant Gauss factor G_p^(−r)
gives the sheaf F_r with trace −k_r, pure of weight r−1 away from
0,1,∞. This normalization follows by setting x_i=u_i y_i in the raw
hypergeometric sum: RawHyp_r(t)=G_p^r k_r(t). All ramification is tame.
The constant twist is fixed over F_p, not chosen separately over each
extension. Its absolute value normalizes the weights as stated.

Write U=P¹\{0,1,∞}. The local forms of F₃ are

    at 0: J₃(1); at ∞: J₃(χ); at 1: diag(−1,1,1).

The endpoint Jordan blocks follow from the repeated character lists;
the form at 1 follows from the pseudoreflection and determinant χ³=χ.
Use the four rank-one Kummer twists with traces
χ(t)^e χ(t−1)^f, e,f∈{0,1}, on U.

The following consequences use only these fixed ranks:

* For F₂⊗F₃, even after any of those twists, an invariant or
  coinvariant would give an isomorphism between irreducible factors of
  different ranks. Hence none exists.
* For F_r⊗F_r with a nontrivial displayed twist, an isomorphism is
  excluded by inertia at 0 if e=1; if e=0,f=1, it is excluded at ∞,
  where the quadratic eigencharacter changes to the trivial one.
* For two distinct pullbacks of F₃ by t, 1−t, 1−1/t, the reflection
  occurs respectively at 1, 0, ∞. At that point the other factor has
  a single Jordan block of length three. Duality and scalar Kummer
  twists preserve these different Jordan block sizes, so no twisted
  isomorphism is possible.
* A single F₃ with any rank-one twist has no invariant or coinvariant
  because it remains irreducible of rank three.

For each resulting lisse tensor on U, H_c⁰=H_c²=0. Since χ_c(U)=−1
and Swan conductors vanish, dim H_c¹ is its rank R. If its weight is w,
the trace formula and the upper H_c¹ weight bound give R p^((w+1)/2).
The sums are on U(F_p)=F_p\{0,1}; no boundary stalk is silently added.
Expanding the exact four-term mask on this domain therefore proves

    |Σ_C k|≤3p^(3/2),  |Σ_C Lk|≤6p²,              (8)
    |Σ_C k(γ_i(t))k(γ_j(t))|≤9p^(5/2), i≠j,      (9)

for γ=(t,1−t,1−1/t). It also bounds every nontrivial mask-weighted
L² sum on U by 4p^(3/2), and every such k² sum by 9p^(5/2).

The untwisted diagonal main terms need no trace-identity assumption.
The full-field quadratic correlation gives Σ_U L²=p²−2p−3.
For k, the multiplicative-convolution description and Parseval give

    Σ_(t≠0)k(t)²=[(p−3)p³+2]/(p−1)
                 =p³−2(p²+p+1).                  (10)

Indeed the Mellin eigenvalues are J(ψ,χ)³, with |J|=√p except at
ψ=1,χ, where |J|=1. Also k(1)=Σ_t f₀(t)L(t)=U_J and the earlier
Jacobi evaluation gives |U_J|≤2p. Subtract k(1)² to restrict (10)
to U. Thus (8)–(10) and the mask prove

    Σ_C L²=p²/4+O(p^(3/2)),
    Σ_C k²=p³/4+O(p^(5/2)),
    ||K||²=3p³/4+O(p^(5/2)).                      (11)

Each diagonal pullback has the same C norm because its parameter map
permutes C. The last line keeps all six cross terms and uses (9).

## 4. Coefficients, the subspace upper estimate, and the next gap

The earlier exact mean gives ℓ=(U_J+6)/(p−5), so |ℓ|≤4 for p≥13.
In particular

    |d−p²/4|≤3p^(3/2)+(18p+3)/4.

For p≥1024 this is less than p²/8, giving d≥p²/8>0. This explicit
threshold ensures that the normalized next direction is defined; it
is not a claim that all asymptotic error constants are already small
at that prime size.

Now d=p²/4+O(p^(3/2)), s₀=O(p^(3/2)) and
s=Σ_C Lk−ℓs₀=O(p²). Equations (6) and (7) yield (1): the two
subtracted squared projections in (7) are O(p²), while ||K||² has
main term 3p³/4. For example, the displayed bounds give
|a−1/2|≤20/√p for p≥1024. The remaining O constants are also
absolute; none depends on the anchor edge or on a coefficient vector.

In the orthonormal basis (u,v), ΠBΠ has matrix [[q,b],[b,a]], tending
to [[3,1],[1,4]]/8. Its top eigenvalue is (7+√5)/16. This proves
(2), with its error rate by the norm bound for a two-by-two perturbation.
The squared norm of B acting on the same domain is the top eigenvalue
of its exact Gram matrix

    [[q²+b², b(q+a)],
     [b(q+a), b²+a²+c²]],

which tends to [[10,7],[7,20]]/64. This proves (3). The component
outside K₂ vanishes on u and has norm c on v. Thus the deleted
cross-block has eigenvalues ±c and zero elsewhere, proving (4).

The new residual h is an explicit actual function, not an abstract
countermodel. The next diagonal coefficient, if h≠0, is

    a₂=1/2+[Σ_(x,y∈C) h(x)χ(x−y)h(y)]
                 /[2√p Σ_C h²].                 (12)

For example, an explicit sufficient arithmetic input for
a₂=1/2+o(1) would be the bound

    Σ_(x,y∈C) h(x)χ(x−y)h(y)=o(p^(7/2))          (13)

for this explicit h, whose squared norm is 3p³/4+O(p^(5/2)).
Equation (13) is not proved here. It retains an additional actual
adjacency kernel connecting two variables, which is absent from the
one-variable correlations (8)–(9). Even proving (13) would control
only another fixed Krylov coefficient. Growing-depth control and the
operator on all remaining directions would still be needed.

## 5. Verification and limits

The uniquely named [exact verifier](../experiments/parallel24_spectral_operator_2026_09_05.py)
and [results](../results/parallel24_spectral_operator_2026_09_05.json)
record direct actual-field recurrences, all four Kummer-mask sums,
pairwise pullback bounds, variance identities, and exact rational/surd
Gram coefficients. The completed run passed on 88 primes from 13
through 1097, checking 11,407 recurrence coordinates, 1,056 twisted
cross-pullback bounds, and 85 nondegenerate next-coefficient identities.
Six small fields additionally check multiplicative convolution against
the full S D₀ definition. The three cases d=0 are p=13,17,29 and are
never divided by d. Six fields meet the explicit p≥1024 threshold.

The result file pins the note, verifier, previous inputs, and both
primary PDFs and their archived text. A metadata-only refresh after
the numerical run is identified separately. The analytic limits come
from the proof and its imported primary results, not finite trends.

All of K₂ and its residual lie in the S₃-invariant sector. This does
not fill the other symmetry sectors identified in the earlier seeded
work, and does not assert that an extremal eigenvector lies in K₂.
The new upper estimate is for the stated subspace, with its genuine
nonzero coupling retained. Human refereeing and formal verification
remain open; no literature novelty claim is made.
