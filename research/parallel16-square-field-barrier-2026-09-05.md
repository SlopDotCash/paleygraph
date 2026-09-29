# Sharp correlations and the Paley identity still require a field restriction

Status: a quantified obstruction over square finite fields, using actual
Paley matrices and the unchanged sharp elliptic bounds. This does not
refute the prime-field conjecture. It proves that the currently listed
field-uniform inputs alone do not yield the desired spectral edge.
No improved prime-field clique bound or full proof is claimed.

The subfield clique is classical, not a new construction. This note
tracks it through our centered projection, elliptic kernel, sharp
correlations and long-power normalization. It also identifies exactly
which hypotheses distinguish it from the abstract examples of passes
twelve and fifteen.

## 1. Actual Paley matrices over a quadratic extension

Let ℓ be an odd prime, q=ℓ², K=F_ℓ and F=F_q. Use the quadratic
character χ of F, with χ(0)=0. For t∈K*,

\[
t^{(q-1)/2}=(t^{\ell-1})^{(\ell+1)/2}=1.
\]

Consequently every pair of distinct elements of K is adjacent in
the Paley graph over F. In particular K is a clique of size ℓ=√q.
For background on this established phenomenon, see the primary
[Asgarli–Yip paper](https://arxiv.org/abs/2110.07176); its classification
of maximum cliques is not needed for the elementary argument here.

For S_{x,y}=χ(x−y) and J=11ᵀ, the same exact character identities
as in a prime field give

\[
S1=0,\quad S^2=qI-J,\quad
P=\tfrac12(I-J/q+S/\ell),\quad P^2=P=P^T,
\quad \operatorname{rank}P=(q-1)/2.                  \tag{1}
\]

For x≠y, the off-diagonal square identity is the elementary sum
Σ_z χ((x−z)(y−z))=−1. On the diagonal the sum is q−1.
These formulas require q≡1 mod4, which holds for every odd square q.

There is also a useful exact eigenvector. For v=1_K,

\[
Sv=\ell v-1,\qquad w=v-1/\ell\;1,
\qquad Sw=\ell w,\qquad Pw=w.                       \tag{2}
\]

Inside K the first sum is ℓ−1. Outside K write F=K(√ν), with
ν a nonsquare in K and x=a+b√ν, b≠0. Then
χ_F(x−t)=χ_K((a−t)²−νb²). The quadratic character sum over t∈K
is −1 because νb² is nonzero. This proves (2). Thus the full
projection has a direction with squared-mass fraction 1−1/ℓ on
the subfield. It is an actual eigenvector, not a matrix perturbation.

## 2. Centered localized outliers at every fixed anchor depth

Choose any fixed a≥1 distinct anchors A⊂K, and let C be their
common neighborhood in F, excluding A. The set B=K\A lies in C
and has n=ℓ−a points. For the normalized indicator of B, the
Rayleigh quotient of X=P_{C,C} is exactly

\[
\lambda_B=\frac12+\frac{n-1}{2\ell}-\frac{n}{2q}
=1-\frac{a+2}{2\ell}+\frac{a}{2\ell^2}.              \tag{3}
\]

Every term is retained, including J/q. Set b=2^a, N=b−1 and

\[
Z_a=\frac{X-I/2}{\sqrt N/b},\qquad
\beta_{\ell,a}=1-\frac{a+2}{\ell}+\frac a{\ell^2}.
\]

Then

\[
\lambda_{\max}(Z_a)\ge\frac b{2\sqrt N}\beta_{\ell,a}.
                                                               \tag{4}
\]

Since 0≤X≤I, the universal upper bound is ‖Z_a‖≤b/(2√N).
For every fixed a, (4) therefore implies

\[
\|Z_a\|\longrightarrow\frac{2^a}{2\sqrt{2^a-1}}
\quad(\ell\to\infty\text{ through odd primes}).       \tag{5}
\]

For a≥2 this limit is greater than 1. For a=3 it is 4/√7.
In particular these are persistent outliers, not an error at finitely
many small fields. The characteristic ℓ also tends to infinity;
excluding finitely many small characteristics does not remove them.

## 3. The complete elliptic model and its sharp constants survive

For a=3 normalize the anchors to {0,1,r} inside K. Every nonzero
element of K is a square in F, so χ(r)=χ(r−1)=1. On
E:y²=x(x−1)(x−r) over F, put G=E(F), H=G/G[2], and retain the
map ρ(Q)=x(2Q) and the lifted function f(Q)=χ(y(Q)), f(O)=0.

All identities in the [elliptic derivation](parallel14-elliptic-model-2026-09-05.md)
hold with q in place of p:

- G[4] consists of 16 rational points; ordinary fibers have eight
  points and the four exceptional fibers have four each.
- K(ρ(P),ρ(Q))=f(P+Q)f(P−Q), including all exceptional pairs.
- The full normalized fiber quotient is
  [[8S_C,8√2·1],[8√2·1ᵀ,12]] ⊕ (−4I₃).
- On H₀=H\H[2], the centered normalized kernel is
  2(M_{H₀}−J_{H₀}/√q)/√(7q), with the same nonzero spectrum as Z₃.
- The four symmetry blocks and the exact, generally off-diagonal
  Fourier-basis formula are unchanged.

These are characteristic-different-from-two group-law and character
identities. Primality of the field cardinality was not used in them.

The sharp scalar and higher-correlation constants also survive. If
O is a set of s>0 distinct shifts in H, then for every character η,

\[
\left|\sum_{h\in H}\prod_{t\in O}f(h+t)\overline{\eta(h)}\right|
\le s\sqrt q=s\ell.                                  \tag{6}
\]

For a multiset with e additional positive-even-multiplicity shifts,
the complete bound remains s√q+e, with exactly the same correction
at each excluded zero. When s=0 the exact character-orthogonality
formula remains valid. Taking s=1 gives |F_η|≤√q with constant 1.
Nothing is inflated to the constant 50 used by the modified sign
functions of pass fifteen.

To verify the field of definition, apply the
[translated-divisor proof](parallel15-elliptic-correlations-2026-09-05.md)
on E over F_q. The Lang isogeny is now 1−Frob_q. Its twist is still
unramified on E. The Kummer sheaf has exactly 4s disjoint nontrivial
tame punctures, H_c^2=0 and dim H_c^1=4s. The eigenvalue bound is
√q, and the full sum is still four times the quotient sum. The
primary [Katz GKM inputs](../sources/katz-gauss-kloosterman-monodromy.pdf)
are stated over finite fields F_q; they do not impose extension
degree one. This extension of the previous source-dependent argument
still needs independent mathematical review.

Thus the *same actual examples* satisfy the sign-kernel structure,
its original Kummer realization, sharp correlation constants and the
ambient Paley projection identity while violating a field-uniform
version of the desired edge bound. Their field cardinalities are
composite squares, which is precisely the missing restriction.

## 4. Failure of the long-depth aggregate in these actual examples

Let E_C denote the common-neighborhood diagonal projection and use
the prior transfer matrix

\[
T=(2P-I)(bE_C-I)/\sqrt N.
\]

The exact [two-projection calculation](parallel2-spectral-transfer-2026-09-04.md)
applies without any field restriction. An eigenvalue z>1 of Z_a
gives a reciprocal pair of positive eigenvalues
z±√(z²−1) of T. Equations (4)–(5), together with ‖T‖≤√N, give

\[
\rho(T)\longrightarrow\sqrt N\qquad(a\ge2\text{ fixed}). \tag{7}
\]

Consequently at depths j/log q→∞, the energy ‖T^j‖_F² cannot
have logarithm o(j), or remain bounded by a fixed power of q.
The same obstruction affects the even signed traces, not only the
energy: the other two-dimensional blocks each contribute at least
−2, and the remaining even-power rank term is nonnegative.

For a concrete uniform version, a=3 and ℓ≥29 imply

\[
\rho(T)>3\sqrt7/4,\qquad
\|T^j\|_F^2>(63/16)^j,\qquad
\operatorname{tr}(T^{2j})>(63/16)^j-q\quad(j\ge1).    \tag{8}
\]

Indeed it suffices that 96β_{ℓ,3}>79. At ℓ=29, clearing
denominators gives 96(29²−5·29+3)−79·29²=665>0.
The function 1−5/ℓ+3/ℓ² is increasing for ℓ≥29, proving the
inequality thereafter. This is equivalent to the lower bound on z
exceeding half the sum of 3√7/4 and its reciprocal. The energy
bound follows from spectral radius; the trace bound follows from
the canonical blocks just described and 2|C|≤q−1.

## 5. What this rules out, and what it leaves open

This does not disprove the prime-field statement. In the primary
[Kunisky formulation](https://arxiv.org/html/2303.16475v1#S2), the
asymptotic parameter runs through primes p≡1 mod4. Our examples
q=ℓ² do not meet that hypothesis. Nor do they establish any
equivalence with, or settle, the official Reed–Solomon prize target.

They rule out obtaining the required norm estimate solely from the
listed field-uniform identities and sharp correlation inequalities.
A successful deduction needs an additional property that excludes
this family, such as quantitative information using the prime
additive group or the lack of a proper subfield. Merely combining
the two constraints lost by the earlier abstract examples is not
enough. No such decisive additional estimate has been proved here.

The [exact verifier](../experiments/parallel16_square_fields_2026_09_05.py)
constructs F_{ℓ²} explicitly, checks its character by Euler's
criterion, and uses the actual unmodified Paley matrix. It checks
both ambient identities, the subfield eigenvector, centered Rayleigh
quotients at five anchor depths, and canonical power-trace identities.
For three fields it checks every elliptic kernel entry, all exceptional
fibers and centering, and selected sharp translated-product bounds
simultaneously for every character. The latter use exact positive
semidefinite certificates, allowing zero pivots with a required zero
row; equality can occur when √q is an integer. Finite checks are
additional evidence, not the proof of the uniform statements above.
See the [results](../results/parallel16_square_fields_2026_09_05.json).
