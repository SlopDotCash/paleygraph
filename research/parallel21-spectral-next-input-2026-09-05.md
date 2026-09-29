# Actual two-anchor Fourier mass: the constant direction is controlled

Status: an exact elementary calculation for every prime p≡1 mod4
and every two-anchor clique, with finite exact checks. It gives a
constant Fourier separation for the uniform vector on the actual
common-neighbor set. It does not give a uniform separation for all
vectors, improve the established logarithmic depth range, or prove
the spectral edge or Paley conjecture. No novelty claim is made for
the classical Jacobi-sum evaluation used here; a self-contained
derivation is included.

This investigation follows the [pass-17 determinant limitation](parallel17-prime-uncertainty-2026-09-05.md): its certificate decays
exponentially when the coordinate set has fixed positive density.
The estimate below improves that certificate on one specified
direction. It does not replace the certificate for the whole matrix.

## 1. Exact prime-field statement

Let χ(0)=0 be the quadratic character, p≡1 mod4 prime, and

    C={x∈F_p: χ(x)=χ(x−1)=1},       m=(p−5)/4.

For p≥13 this set is nonempty. Let S_(x,y)=χ(x−y), and choose a
quartic character η with η²=χ, extended by zero. Set

    J=Σ_x η(x)χ(1−x),       R_C=1_C^T S 1_C.

Then

    |J|²=p,
    R_C=(J²+conj(J)²−6p+36)/16.                    (1)

Consequently

    (9−2p)/4 ≤ R_C ≤ (9−p)/4,
    −2−1/(p−5) ≤ R_C/m ≤ −1+4/(p−5).              (2)

The normalization is unchanged for any actual edge {a,b}: the affine
map t↦a+(b−a)t preserves quadratic characters because χ(b−a)=1.
Thus (1) and (2) apply uniformly to every two-anchor clique, with no
averaging over anchors or primes.

Use the unitary additive Fourier transform and the projection

    P=(I−J_all/p+S/√p)/2

onto nonzero quadratic-residue frequencies. This frequency labeling uses
the classical positive quadratic Gauss evaluation. Its missing sign
normalization is supplied, from NIST DLMF 20.11.1–2, in
[the independent review, Section 4](parallel22-spectral-independent-review-2026-09-05.md).
The Gauss-magnitude argument below suffices for the Jacobi norm but
cannot determine that frequency label by itself.
For u=1_C/√m, define
the fraction of its total Fourier mass in those frequencies by

    q_C=||Pu||²=u^T P u.

Here J_all denotes the all-ones matrix, distinct from the Jacobi sum.
The exact identity is

    q_C=(3p+5)/(8p)+R_C/(2m√p).                    (3)

In particular q_C=3/8+O(p^(−1/2)) with an absolute constant, and
for every prime p≡1 mod4 with p≥73,

    1/4 ≤ q_C < 1/2.                              (4)

To check the lower constant without an asymptotic convention, use
R_C/m≥−17/8 for p≥13 and √p≥17/2 for p≥73:
q_C≥3/8−17/(16√p)≥1/4. The upper bound follows from
R_C<0 and (3p+5)/(8p)<1/2. The zero-frequency mass is exactly
m/p; after discarding that mass, the residue/nonresidue split tends
to one half each. The fraction 3/8 in (3) concerns total mass.

## 2. A quadratic-root identity

For any A∈F_p and B∈F_p^×,

    Σ_t χ[t(t²+At+B)]
      =Σ_z χ[z(z²−2Az+A²−4B)].                   (5)

Indeed z=t+A+B/t has, for every z, exactly
1+χ((z−A)²−4B) preimages t∈F_p^×. The product of the two
roots is B≠0, so there is no missing zero root. Therefore

    Σ_(t≠0) χ(t+A+B/t)
      =Σ_z χ(z)[1+χ((z−A)²−4B)].

The sum of χ(z) vanishes; multiplying the argument on the left by
the square t² gives (5). This proves the identity by root counting,
without importing elliptic-curve point-count invariance under isogeny.

## 3. The double character sum is a squared Jacobi sum

Put

    U=Σ_(x,y) χ[x(x−1)y(y−1)(x−y)].

I claim U=J²+conj(J)². In the expansion of the latter, group an
ordered pair of field elements by its sum s and product t. There are
1+χ(s²−4t) such pairs. Thus

    J²+conj(J)²
      =Σ_(s,t) [η(t)+conj(η(t))]χ(1−s+t)
                                 [1+χ(s²−4t)].

The term without χ(s²−4t) sums to zero over s. Also, for every
function H of t,

    Σ_t [η(t)+conj(η(t))]H(t)=Σ_z χ(z)H(z²).

For nonsquare t the coefficient on the left vanishes; for nonzero
square t its two square roots have the same quadratic character,
since χ(−1)=1. Both sides vanish at zero. It follows that

    J²+conj(J)²
      =Σ_(s,z) χ[z(1−s+z²)(s²−4z²)].

The z=0 terms vanish. Substitute s=2z(2x−1), a bijection in s
for each z≠0. The sum becomes

    Σ_x χ[x(x−1)] Σ_z χ[z(z²+(2−4x)z+1)].

For x∉{0,1}, apply (5) with A=2x−1, B=x(x−1). The inner
sum equals Σ_t χ[t(t+x)(t+x−1)]. Set y=t+x and use χ(−1)=1
to obtain U. The excluded x have zero outer weight, so no correction
has been discarded.

For completeness, the Jacobi norm follows from the elementary Gauss
sum identities. For a nontrivial multiplicative character A,

    G(A)=Σ_x A(x)exp(2πix/p),       |G(A)|²=p.

To prove the norm, set x=ty in its squared modulus. The inner sum
over y≠0 is p−1 when t=1 and −1 otherwise; the nontrivial
character sum over t is zero. If A,B,AB are all nontrivial, splitting
G(A)G(B) by z=x+y gives

    G(A)G(B)=J(A,B)G(AB).

The z=0 term vanishes because AB is nontrivial; for z≠0 substitute
x=zu. Take A=η, B=χ; all three characters are nontrivial. Hence
|J|²=p and |U|=|J²+conj(J)²|≤2p.

## 4. Keep the two exceptional coordinate rows

Let v_0(x)=χ(x), v_1(x)=χ(x−1), f=v_0v_1,
g=1+v_0+v_1+f, and d=δ_0+δ_1. The exact mask is

    1_C=g/4−d/2.

The subtraction is necessary: g/4 equals 1/2 at each anchor.
The standard quadratic correlation identity gives Σf=−1, hence
|C|=(p−1)/4−1=m. From S1=0 and S²=pI−J_all,

    g^T Sg=2p+4+U,
    g^T Sd=2p−6,
    d^T Sd=2.

For example, v_0^T Sv_1=p, v_i^T Sf=1, and (Sf)(0)=(Sf)(1)=−1.
Expanding the mask therefore gives

    R_C=(2p+4+U)/16−(2p−6)/4+2/4
       =(U−6p+36)/16,

which proves (1). These anchor corrections have constant effect on
R_C/m and cannot be omitted from an exact formula.

## 5. What this does and does not resolve

For this direction, (4) supplies a fixed lower bound on both
quadratic-residue mass and its complement. The pass-17 determinant
certificate is exponentially small at m=(p−5)/4. Thus the Fourier
calculation removes that loss for the specified uniform vector on
the actual prime-field common-neighbor set.

It is a Rayleigh-quotient estimate on a one-dimensional subspace.
It does not control u^T(P_(C,C))²u, make u an
eigenvector of P_(C,C), or control an adversarial vector with mixed
signs. The set C is not asserted to be invariant under S. A useful
exact diagnostic of the missing coupling is

    (S_C 1)(x)=(L(x)−6)/4,       x∈C,
    L(x)=Σ_y χ[y(y−1)(y−x)],

and therefore

    ||S_C u−(R_C/m)u||²
      = (1/(16m))Σ_(x∈C)(L(x)−6)²−(R_C/m)².       (6)

The squared character-trace average in (6), as well as the operator
on the subspace orthogonal to u, is not estimated here. A small
first Rayleigh quotient does not bound either quantity. This identity
is a diagnostic, not a claimed new sufficient theorem for the edge.
Likewise, a three-anchor common-neighbor set is a subset of C, but
the uniform-vector estimate does not automatically pass to subsets.

No all-vector prime-field Fourier estimate or stronger complete-power
bound was found in this bounded lane. The established [pass-13
logarithmic depth range](parallel13-principal-budget-2026-09-05.md)
remains unchanged. The computation identifies a structured direction
that is already controlled and leaves the nonconstant directions and
their coupling as the explicit missing work.

## 6. Verification and source scope

The [verifier](../experiments/parallel21_spectral_flat_direction_2026_09_05.py)
checks (5) for all parameters in selected small fields, computes
quartic Jacobi sums in exact Gaussian integers, evaluates U and the
actual common-neighbor matrix independently, and checks all the
formulas above. It also checks every affine edge at selected primes
and the exact row-variance identity (6). Its [results](../results/parallel21_spectral_flat_direction_2026_09_05.json) contain
the parameter ranges, counts and input hashes. Fourier fractions are
represented in Q(√p), with exact comparisons; no floating-point
threshold establishes a bound.

A literature search located earlier work on Jacobi-sum evaluations
for Paley four-clique counts, including the authors' [Bhowmik–Barman
abstract](https://arxiv.org/abs/2301.07021). That abstract is background
for the no-novelty claim, not an input to any equation above. A versioned
HTML fetch failed; no full text or PDF from that work was reviewed or
used. The Jacobi and anchor identities are derived explicitly here. The
Fourier-frequency identification additionally uses the classical positive
quadratic Gauss evaluation, as supplied in the linked separate-agent
review. Its review found no formula error and identified this omitted
source dependency. Separate human and formal proof review remain open. Finite checks
support the derivation but do not prove the unresolved spectral target.
