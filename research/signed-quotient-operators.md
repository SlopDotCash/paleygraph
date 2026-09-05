# Signed quotient operators and the loss from taking absolute values

**Status: the uniform subgroup estimate, classical Paley conjecture, and
prize reduction remain open.** We construct sparse symmetric integer
matrices for the nonprincipal periods and for their two-child products.
An unconditional estimate shows that discarding the matrix signs loses
almost all the desired cancellation in the working quartic window.
The construction does not itself bound the signed spectrum. All proofs
below are ordinary mathematics; no novelty or Lean-formalization claim
is made.

## A symmetric integer representation without the principal frequency

Let H≤F_p* have dyadic order n≥4, let g generate H, put
`m=(p−1)/n`, and let ψ:H→{±1} be given by `ψ(g)=−1`. Choose one
representative a_i from each nonzero H-coset. Consider a nonnegative
integer function w:F_p→Z with

\[
 w(0)=0,\qquad w(hx)=w(x)\quad(h\in H),\qquad
 d=\sum_x w(x)>0,\quad E=\sum_x w(x)^2.
\]

In particular w is even. Write

\[
 \lambda_w(b)=\sum_x w(x)e_p(bx).
\]

These numbers are real and constant on bH. Define functions

\[
 \phi_j(x)=
 \begin{cases}\psi(x/a_j),&x\in a_jH,\\0,&\text{otherwise},\end{cases}
\]

and an m×m integer matrix

\[
 A(w)_{ij}=\sum_z w(z)\phi_j(a_i-z).
 \tag{1}
\]

Then

\[
 \boxed{A(w)=A(w)^T,\qquad
       \operatorname{Spec}A(w)=\{\lambda_w(a_i):1\le i\le m\},}
 \tag{2}
\]

with multiplicities. In particular the principal eigenvalue d is absent.

**Proof.** Additive convolution by w is self-adjoint and preserves the
space

\[
 V_\psi=\{f:f(hx)=\psi(h)f(x)\text{ for all }h\in H\}.
\]

Every such f vanishes at zero. The φ_j form an orthogonal basis of
this space, all with squared norm n. Evaluation at a_i gives the
coefficient formula (1), so the matrix is symmetric and integral.
Additive Fourier transformation maps this space into the corresponding
covariant frequency space. On each nonzero frequency coset there is one
free coefficient, while the zero coefficient vanishes. Convolution
acts there by the constant λ_w(a_i), proving (2). Equivalently, the
nonzero functions `Σ_(h∈H) ψ(h)e_p(a_i h x)` are eigenvectors with
disjoint Fourier supports. This also proves completeness without an
assumption about distinct eigenvalues. ∎

For w=1_H, write S_H=A(1_H). Its eigenvalues are exactly the usual
nonprincipal periods η_H(a_i), once per H-coset. Thus

\[
 \|S_H\|_{\mathrm{op}}=\max_{b\ne0}|\eta_H(b)|.
\]

Unlike the earlier nonsymmetric integer representation and symmetric
matrix containing the principal eigenvalue in
[mixed-periods-and-shifted-energy.md](mixed-periods-and-shifted-energy.md),
S_H is symmetric, integral, and has the required spectrum directly.
Its dimension is m, and each row has at most n nonzero entries.
These structural facts supply no bound smaller than n by themselves.

## The two-child product uses the same construction

Write n=2k, K=⟨g²⟩ and L=gK, and take

\[
 w=1_K*1_L,\qquad R_H=A(w).
\]

Since −1∈K and K,L are disjoint, w(0)=0. Scaling by g exchanges K
and L, so w is H-invariant. It has mass d=k², and

\[
 E=\sum_x w(x)^2=B(K,L),
\]

where B is the mixed energy from
[dyadic-descent-and-mixed-energy.md](dyadic-descent-and-mixed-energy.md).
Consequently

\[
 \boxed{\operatorname{Spec}R_H
       =\{\eta_K(a_i)\eta_K(ga_i):1\le i\le m\}.}
 \tag{3}
\]

Here eigenvalues may coincide: for example the order-eight parent in
F_17 gives R_H=−I₂. Formula (3) does not assert simplicity.

Let c_r=(w^{*r})(0). Character orthogonality gives, for either weight
and every positive integer r,

\[
 \boxed{\operatorname{tr}A(w)^r=\frac{p c_r-d^r}{n}.}
 \tag{4}
\]

For R_H, c_r is exactly Z_(r,r), the count with r entries from each
child coset and total zero. For S_H, c_(2q)=E_q(H). Thus the product
criterion in [positive-product-moments.md](positive-product-moments.md)
has the precise operator form

\[
 \left(\frac1m\operatorname{tr}(R_H)_+^q\right)^{1/q}\le Cqk.
 \tag{5}
\]

The positive part is defined spectrally. For even q its stronger
absolute version replaces `(R_H)_+^q` by R_H^q. This recovers the
centered count exactly; it does not prove (5).

## Taking entrywise absolute values is asymptotically useless here

Define the nonnegative unsigned quotient

\[
 B(w)_{ij}=\sum_z w(z)1_{a_jH}(a_i-z).
\]

Then B(w) is symmetric and `|A(w)|≤B(w)` entrywise. Its row sums
are `d−w(a_i)`, so

\[
 \sum_{i,j}B(w)_{ij}=md-d/n.
 \tag{6}
\]

There is also an exact Frobenius identity:

\[
 \boxed{\sum_{i,j}\bigl(B(w)_{ij}^2-A(w)_{ij}^2\bigr)=d^2-2E.}
 \tag{7}
\]

**Proof.** On the H-invariant subspace, include the unit vector at
zero and normalize every nonzero coset indicator by √n. The convolution
matrix is

\[
 \begin{pmatrix}
  0&\sqrt n\,(w(a_i))_i^T\\
  \sqrt n\,(w(a_i))_i&B(w)
 \end{pmatrix}.
\]

Its spectrum consists of d and the same m nonprincipal eigenvalues as
A(w). Taking the trace of the square, and using
`nΣ_i w(a_i)²=E`, proves (7). In addition, (4) at r=2 gives
`Σ A(w)_ij²=(pE−d²)/n`. ∎

For each entry write `B_ij=u+v` and `A_ij=u−v`, where u,v are the
nonnegative integer sums over the two signs. Integrality gives

\[
 B_{ij}-|A_{ij}|=2\min(u,v)\le2uv
                =\tfrac12(B_{ij}^2-A_{ij}^2).
\]

Summing and testing |A(w)| on the constant vector yields

\[
 \boxed{
 d-\frac{d/n+d^2/2-E}{m}
 \le\rho(|A(w)|)\le d.
 }
 \tag{8}
\]

The upper bound uses the row sums; the lower bound is its Rayleigh
quotient. A negative lower bound outside the intended parameter range
is harmless. This proof uses integer weights, not arbitrary real ones.

Now impose the working quartic window `n⁴/4≤p≤n⁴`.
For S_H, d=E=n, so

\[
 \boxed{n-2/n\le\rho(|S_H|)\le n.}
 \tag{9}
\]

Indeed the deficit in (8) is
`n(1+n²/2−n)/(p−1)`, which is at most 2/n: replacing p by n⁴/4
reduces the required inequality to `n³−n²≥2`.

For R_H, d=k² and E=B(K,L)≥k², since w is a nonnegative integer
function of mass k². Since `p≥4k⁴`, (8) gives

\[
 \boxed{k^2-k/4\le\rho(|R_H|)\le k^2.}
 \tag{10}
\]

For the lower bound it suffices that

\[
 \frac{2k(k/2+k^4/2-k^2)}{4k^4-1}\le k/4,
\]

which follows from `8k²−4k≥1` for k≥2.

These are uniform statements about the actual matrices. They rule out
obtaining the desired cancellation by the bound
`||A(w)||≤ρ(|A(w)|)`: the latter is almost the trivial degree throughout
the quartic window. In particular, for the product matrix it is far
larger than the scale k log p needed by a pointwise sufficient version
of (5). This does not disprove (5), or exclude arguments that retain
signs in higher powers.

## Representative choices cannot make the signs independent

Replace a_i by h_i a_i, with h_i∈H, and set
D=diag(ψ(h_i)). Changing the basis in (1) gives exactly

\[
 A'(w)=D A(w)D,\qquad B'(w)=B(w).
 \tag{11}
\]

Hence all eigenvalues, the absolute matrix, and every product of
matrix entries along a closed walk are unchanged. Randomizing these
choices provides no independent random-edge-sign model. The two
factors of ψ(h_i) cancel at each visited vertex of a closed walk.
Other sources of cancellation are not excluded.

## Exact signed closed-walk counts

A quotient walk based at a_i lifts to an additive walk starting at a_i
and ending at u a_i for some u∈H. Call u its multiplicative endpoint
ratio. Count weights with their integer multiplicities. For length r,
let W_r(u) be the total number of these walks, over all representatives.
The additive walk count w^{*r} is H-invariant, so

\[
 W_r(u)=\sum_i(w^{*r})((1-u)a_i).
\]

For u≠1, multiplication by 1−u permutes the nonzero field elements.
Summing over H-cosets therefore proves

\[
 \boxed{
 W_r(1)=m c_r,\qquad
 W_r(u)=\frac{d^r-c_r}{n}\quad\text{for every }u\ne1.
 }
 \tag{12}
\]

The sign accumulated along the walk is ψ(u). Since Σ_(u∈H)ψ(u)=0,
(12) gives (4) again. Thus all nonidentity endpoint ratios have exactly
the same count, not independent fluctuations. The signed trace subtracts
that common count from the identity count. Controlling the remaining
imbalance at logarithmic depth is exactly the unproved centered-count
problem, rather than a new consequence of (12).

## Exact verification and scope

[The script](../experiments/signed_quotient_operators.py) constructs both
matrices directly from residues. It checks symmetry, the two Frobenius
formulas, the loss bound, and diagonal sign conjugation reconstructed
from different field representatives. Independent dense additive
convolution verifies 94 signed trace identities across 14 matrix cases.
The seven product cases also check 46 balanced counts independently by
convolving the child subgroup and correlating its scaled counts.
There are 32 checks of (12) in the two small fields.

The cases are `(p,n)=(17,8),(17,16),(97,8),(97,32),(1049,8),(2017,8)`
and `(17393,16)`. The last three are in their quartic windows. Matrix
dimensions reach 1087. Sparse int64 products are used only after the
explicit bound d^r≤2^63−1; all trace sums and principal subtractions
use Python integers. The order-32 product case stops at depth six to
respect that bound. No floating-point eigenvalue calculation is an
acceptance test, and the completed order-64 spectral maximum is not
recomputed.

For example, at p=17,n=8 the matrices are

\[
 S_H=\begin{pmatrix}-1&2\\2&0\end{pmatrix},\qquad R_H=-I_2.
\]

At p=17393,n=16, the subgroup matrix has absolute entry sum 17279,
whereas the unsigned entry sum is 17391. Its exact squared Frobenius
norm is 17377. The product matrix has degree 64 and absolute entry sum
67788 on 1087 vertices. These finite checks support the identities,
not a uniform bound on the signed spectrum.

Results and source hashes are saved in
[results/signed_quotient_operators.json](../results/signed_quotient_operators.json).
The next arithmetic issue is cancellation among signed closed walks
at logarithmic depth, or directly the positive spectrum of R_H. Neither
an unsigned spectral comparison nor random representative choices can
supply it. The full Paley and prize targets remain unproved.
