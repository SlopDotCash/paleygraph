# Independent review of the two-anchor constant-direction estimate

**Verdict: the root-count transformation, Jacobi-square identity, anchor
corrections, Fourier-mass formula and p>=73 constants are correct.** The
result concerns one specified vector and does not give the missing
all-vector spectral estimate.

One proof-scope correction is needed: identifying the positive Paley
projection with **quadratic-residue frequencies** requires the positive
quadratic Gauss-sum evaluation. The reviewed note proves Gauss-sum
magnitudes, which do not determine that sign, and pass 17 asserts the
same frequency labeling without deriving the sign. Section 4 below supplies
the missing normalization from a checked source. No displayed formula
needs an algebraic correction, but the claim that every needed identity
is explicitly derived in the reviewed note should acknowledge this
standard additional input.

This is separate-agent mathematical and code review, not human refereeing
or formal certification. The verifier was read and its recorded input
hashes were checked; it was not rerun. Repeating it would not resolve
the omitted Gauss-sign argument, since it computes the displayed surd
formula without independently evaluating the additive Fourier transform.
No prior or central artifact was edited.

## 1. Pinned inputs

| Reviewed input | SHA-256 |
|---|---|
| research/parallel21-spectral-next-input-2026-09-05.md | 60d7b38d0f1a66e5db4689cc8bc0329d5a8f9bd5c9ee68b664d98a67df11c345 |
| experiments/parallel21_spectral_flat_direction_2026_09_05.py | 32167075dbc2f8022b1345fae32795b61c6857a583395bac23492bf042783343 |
| results/parallel21_spectral_flat_direction_2026_09_05.json | 86d02185567e44d05b5fd4831bb3f2f34c6ad366e22308f4496693f891075f9c |
| research/parallel17-prime-uncertainty-2026-09-05.md | 2eaca714a8c2a139e0be257f9c2196824cbf720e566ce4ef8f0024af87f9ad83 |
| research/parallel13-principal-budget-2026-09-05.md | 50b4934169a32581859712d60269e8f8e06d829e92bfe90bd042ad68fb8684a9 |

The [spectral note](parallel21-spectral-next-input-2026-09-05.md),
[verifier](../experiments/parallel21_spectral_flat_direction_2026_09_05.py)
and [stored result](../results/parallel21_spectral_flat_direction_2026_09_05.json)
were read directly. The stored result records a subsequent Markdown-link
formatting change and a refreshed note hash; all four of its pinned
input hashes match the present files. Python syntax parsing passed.

## 2. Root counting and the squared Jacobi sum

For odd p and B!=0, the equation

\[
 z=t+A+B/t
\]

is equivalent to t^2+(A-z)t+B=0. Its discriminant is (z-A)^2-4B,
so its number of roots is 1+chi((z-A)^2-4B). Every root is nonzero
because their product is B. This includes the double-root case, where
exactly one preimage is counted.

Multiplying z=t+A+B/t by the square t^2 shows that
chi(z)=chi(t(t^2+At+B)). Summing over fibers and using sum_z chi(z)=0
proves the reviewed equation (5). The terms t=0 and z=0 vanish exactly;
there is no missing coordinate correction. This identity does not require
p=1 modulo 4, consistently with the verifier's additional test fields.

For J=J(eta,chi), grouping the two variables in J^2+conjugate(J)^2 by
s=u+v and t=uv gives exactly 1+chi(s^2-4t) ordered pairs. The two
characters combine as

\[
 [\eta(t)+\overline\eta(t)]\chi(1-s+t).
\]

The term without the discriminant character vanishes after summing over
s. The replacement

\[
 \sum_t[\eta(t)+\overline\eta(t)]H(t)=\sum_z\chi(z)H(z^2)
\]

is also exact. At nonsquares the left coefficient vanishes; at a
nonzero square the two roots have equal quadratic character since
chi(-1)=1. Both sides give zero at t=z=0. No assumption that
eta(-1)=1 is made or needed.

For z!=0, the substitution s=2z(2x-1) is bijective. The discriminant
becomes 16z^2 x(x-1), whose square factor has character one. This gives
the cubic polynomial z(z^2+(2-4x)z+1) in the inner sum. Applying the
root-count identity with

\[
 A=2x-1,\qquad B=x(x-1),\qquad A^2-4B=1
\]

is valid for x outside {0,1}; the excluded x have zero outer weight.
The transformed polynomial is t(t+x)(t+x-1). Setting y=t+x and
using chi(y-x)=chi(x-y) completes the stated identity

\[
 U=J^2+\overline J^{,2}.
\]

The Gauss-product identity used to prove abs(J)^2=p requires eta,
chi and eta*chi to be nontrivial. They have orders 4, 2 and 4,
respectively, so the condition holds. The z=0 part of the Gauss-product
expansion vanishes because eta*chi is nontrivial. The stated norm
argument is sufficient for this Jacobi norm; its missing phase only
matters later for naming the Fourier frequency set.

## 3. The anchors and exact quadratic form

Write v_i=S delta_i for i=0,1, f=v_0 v_1, g=1+v_0+v_1+f,
and d=delta_0+delta_1. At each anchor g/4 is 1/2, so

\[
 1_C=g/4-d/2
\]

is the correct mask. The pair-correlation sum is sum f=-1, giving
m=(p-5)/4. Ratios involving m and the normalized vector require
p>=13, as the reviewed text specifies when introducing nonemptiness.

The standard identities S1=0 and S^2=pI-J_all imply

\[
 Sv_i=p\delta_i-1,\quad
 v_i^TSv_i=0,\quad v_0^TSv_1=p,\quad v_i^TSf=1.
\]

Also (Sf)(0)=(Sf)(1)=-1: in either complete sum the squared
character removes just the corresponding anchor term. Therefore

\[
 g^TSg=2p+4+U,\qquad g^TSd=2p-6,\qquad d^TSd=2.
\]

The cross term from expanding the mask has coefficient -1/4, so

\[
 R_C=(2p+4+U)/16-(2p-6)/4+2/4
     =(U-6p+36)/16.
\]

All constants and signs check out. In particular, the number 36
contains the anchor contributions and must not be dropped.
Using abs(U)<=2p yields

\[
 (9-2p)/4\le R_C\le(9-p)/4,
\]

and division by m gives exactly the two bounds in the reviewed
formula (2). For p>=13 the upper bound is negative. Transporting C
by x=a+(b-a)t preserves its sign matrix whenever {a,b} is an edge,
since chi(b-a)=1. Thus the uniformity over actual two-anchor cliques
is justified without averaging.

## 4. Projection normalization: the missing sign input

From S1=0 and S^2=pI-J_all alone,

\[
 P=\tfrac12(I-J_{all}/p+S/\sqrt p)
\]

is a real orthogonal projection. These identities identify it as the
positive-S eigenspace. They do not by themselves determine whether that
eigenspace corresponds to residue or nonresidue Fourier frequencies.
The Gauss-magnitude proof in the reviewed note likewise cannot choose
between +sqrt(p) and -sqrt(p).

Here is a source-backed normalization. With the definition in
[NIST DLMF 20.11.1](https://dlmf.nist.gov/20.11.E1), its
[inversion formula 20.11.2](https://dlmf.nist.gov/20.11.E2), specialized
to (m,n)=(2,p), gives

\[
 \frac1{\sqrt p}\sum_{k=0}^{p-1}e^{-2\pi i k^2/p}
 =\frac{e^{-\pi i/4}}{\sqrt2}\left(1+e^{\pi i p/2}\right)=1
\]

when p=1 modulo 4. The hypotheses hold: 2 and p are positive and
coprime, and 2p is even. Counting roots of k^2=x shows that the sum
on the left equals the negative-exponent quadratic Gauss sum; the
trivial additive-character sum cancels. Complex conjugation gives the
positive-exponent evaluation as well. Hence both equal +sqrt(p).

For the unitary transform

\[
 (\mathcal Fv)(\xi)=p^{-1/2}\sum_xv(x)e^{-2\pi i\xi x/p},
\]

S has eigenvalue sqrt(p)chi(xi) at xi!=0 and zero at xi=0.
Thus P is indeed the projection onto the nonzero residue frequencies,
with exactly the plus sign displayed in the note.

This supplies the required proof-scope correction without changing its
mathematical statement. A suitable addition to the original source
scope would say that the Fourier-frequency identification uses the
classical positive quadratic Gauss evaluation, with a proof or source.
The reviewed code is consistent with this normalization but does not
independently certify it.

## 5. Mass normalization and the finite threshold

For u=1_C/sqrt(m), the exact projection identity gives

\[
 \|Pu\|^2=u^TPu
 =\frac12\left(1-\frac mp+\frac{R_C}{m\sqrt p}\right)
 =\frac{3p+5}{8p}+\frac{R_C}{2m\sqrt p}.
\]

The constant 3/8 refers to the fraction of **total** mass. The
zero-frequency mass is m/p; the nonzero mass is 1-m/p. Dividing the
residue mass by this latter quantity does give a limiting fraction 1/2.
The normalization and interpretation in the reviewed note are correct.

For p>=13, the lower quadratic-form bound implies
R_C/m>=-2-1/(p-5)>=-17/8. For p>=73 we have
sqrt(p)>=17/2, since 73>289/4. Therefore

\[
 q_C\ge\frac38-\frac{17}{16\sqrt p}\ge\frac14.
\]

The dropped term 5/(8p) is positive. The upper bound is strict because
R_C<0 and (3p+5)/(8p)<1/2 for p>5. Thus all constants in the stated
p>=73 conclusion hold with no hidden eventual-size assumption.
The complementary mass is consequently also bounded away from zero.

## 6. The one-vector limitation is real

The conclusion controls the Rayleigh quotient of P_(C,C) only at the
constant vector. It does not give a positive lower bound for its least
eigenvalue or an upper bound below one for its largest eigenvalue.
Basic contraction inequalities give some coarse second-moment bounds,
but do not supply the missing operator control or a vanishing residual.

For x in C, use the exact mask and Sv_i=p delta_i-1 to obtain

\[
 (S_C1)(x)=(L(x)-6)/4,
 \qquad L(x)=\sum_y\chi(y(y-1)(y-x)).
\]

Subtracting the mean row sum yields precisely

\[
 \|S_Cu-(R_C/m)u\|^2
 =\frac1{16m}\sum_{x\in C}(L(x)-6)^2-(R_C/m)^2.
\]

This verifies the diagnostic formula, not an upper estimate for its
right side. The stored example at p=37 has residual squared equal to
3, so even a genuine common-neighbor constant vector need not be an
eigenvector. Passing to three anchors or other subsets is not justified
by the one-vector theorem. The reviewed note states these limits
accurately and makes no spectral-edge inference from the small first
Rayleigh quotient.

## 7. Verifier assessment

The code uses only Python integers, Gaussian integer pairs and exact
Fractions. Its quartic character is generated from a checked primitive
root. The square agrees with the quadratic character, and the Gaussian
Jacobi norm is compared to p exactly. The polynomial double sum and
the common-neighbor quadratic form are evaluated independently. The
quadratic transformation is checked by both complete sums and direct
fiber counts; the residual identity has an independently evaluated
cubic character sum. The rational-plus-square-root sign comparison
handles all sign cases correctly. There is no fixed-width overflow or
floating-point acceptance issue.

The stored result reports 31,822 fiber counts, 1,418 complete
root-transform checks, 79 prime neighborhoods, 1,742 affine edges,
9,057 row identities and 72 p>=73 mass certificates. These are
historical finite checks, not new executions during this review and
not asymptotic proof evidence. The stored Fourier values are obtained
from the verified R_C and the stated surd formula; they do not test
the quadratic Gauss sign independently.

With the explicit Gauss normalization supplied above, the stated
constant-direction theorem survives this independent review. No
all-vector Fourier separation, longer complete-power estimate,
improved spectral edge, Paley proof or official-prize bridge follows.
