# An exact aggregate trace criterion for clique bounds

**Status: a quantitative reduction, not a proved aggregate bound or a
solution of Paley.** This note keeps the rank and anchor terms in the
localization argument exactly. It shows which growing-depth estimate
would imply a subpolynomial clique bound. That conclusion is weaker than
the full arbitrary-two-set character-sum conjecture tracked in this project.

[Kunisky, Sections 4 and 7](https://arxiv.org/html/2303.16475v1#S7),
already distinguish weak spectral convergence from control of extreme
eigenvalues, and identify the need for cancellation between necklaces.
The derivation below uses a direct two-projection calculation with all
normalizations displayed. No novelty claim is made.

## 1. The finite matrices

Let p be a prime congruent to 1 modulo 4 and set
S_xy=chi(x-y), chi(0)=0, J=11^T. The identities S1=0 and
S^2=pI-J imply that

\[
 P=\frac12\left(I-\frac Jp+\frac S{\sqrt p}\right)
 \tag{1}
\]

is an orthogonal projection of rank r=(p-1)/2.

Let I_0 be a clique of size a>=1. Let R be its common neighborhood,
excluding the anchors, let E be its diagonal indicator, and put m=|R|.
Then m<=r, since R is contained in the neighborhood of any one anchor.
For b=2^a, beta=1/b and c=sqrt(beta(1-beta)), define

\[
 D=bE-I,\qquad W=(\sqrt p S-J)D,
 \qquad X=P_{R,R}.
 \tag{2}
\]

Here the identity matrix I in (1)-(2) is distinct from the anchor set I_0.
The matrix W need not be symmetric; its traces, rather than a claimed
positive spectrum, are used below.

## 2. An exact two-projection identity

Let lambda_1,...,lambda_m be the eigenvalues of X, and let T_k denote
the usual Chebyshev polynomial, T_k(cos theta)=cos(k theta). For every
integer k>=1,

\[
 \boxed{
 \frac{\operatorname{tr}W^k}{[p\sqrt{b-1}]^k}
 =2\sum_{i=1}^m T_k\!\left(\frac{\lambda_i-1/2}{c}\right)
 +\frac{(p-r-m)+(-1)^k(r-m)}{(b-1)^{k/2}}.}
 \tag{3}
\]

**Proof.** This identity holds for any diagonal projection E of rank
m<=min(r,p-r); it does not require a common neighborhood. On each
two-dimensional canonical subspace of P,E associated to lambda, they
have matrices

\[
 E=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
 P=\begin{pmatrix}\lambda&\sqrt{\lambda(1-\lambda)}\\
 \sqrt{\lambda(1-\lambda)}&1-\lambda\end{pmatrix}.
\]

For 0<lambda<1, take a unit eigenvector v of EPE in the range of E;
the second unit vector is (Pv-lambda v)/sqrt(lambda(1-lambda)).
These pairs are orthogonal for an orthonormal eigenbasis. At endpoints
lambda=0 or 1, pair v with an available vector in the opposite P-eigenspace
inside the kernel of E. The rank assumptions guarantee enough such vectors.
After these m blocks there remain r-m directions with (P,E)=(1,0)
and p-r-m with (P,E)=(0,0).

On a two-dimensional block, U=(2P-I)(E-beta I) has trace 2lambda-1
and determinant c^2. Its power trace therefore equals
2c^k T_k((lambda-1/2)/c), by the characteristic-polynomial recurrence.
The remaining eigenvalues are -beta, repeated r-m times, and beta,
repeated p-r-m times. Finally U=W/(pb). Dividing its trace by c^k
proves (3). This argument also covers repeated eigenvalues and endpoints.

For even k the residual term is exactly
(p-2m)/(b-1)^(k/2), which is nonnegative. In particular the unwanted
rank terms do not need to be bounded by a constant before multiplication
by an exponentially growing normalization.

There is an equivalent identity with no Chebyshev evaluation and only
arithmetic in the real quadratic field Q(sqrt(p)). Put

\[
 Z=\sqrt p S_{R,R}-J_m,\qquad H_0=2I_m,\quad H_1=bZ,
 \quad H_k=bZ H_{k-1}-p^2(b-1)H_{k-2}.
\]

Then

\[
 \boxed{\operatorname{tr}W^k=\operatorname{tr}H_k
 +p^k[(p-r-m)+(-1)^k(r-m)].}
 \tag{4}
\]

This is the version checked using pairs of integer matrices in the verifier.

## 3. A nonasymptotic edge and clique implication

For even k>=2, suppose the explicit aggregate satisfies

\[
 \left|\frac{\operatorname{tr}W^k}{[p\sqrt{b-1}]^k}\right|\le B,
 \qquad B\ge0.
 \tag{AG}
\]

If m>0, (3) implies

\[
 \boxed{\max_i|\lambda_i-1/2|
 \le c\cosh\!\left(\frac{\log(B+2m)}k\right).}
 \tag{5}
\]

Indeed, write z_i=(lambda_i-1/2)/c. If |z_i|>1, then
T_k(z_i)=cosh(k arcosh |z_i|), since k is even. Every other T_k(z_j)
is at least -1. The nonnegative residual term in (3) gives
cosh(k arcosh |z_i|)<=B/2+m-1. The weaker bound using
arcosh(t)<=log(2t) proves (5). If every |z_i|<=1, (5) is automatic.

Let C be any clique containing I_0 and put ell=|C|-a. If ell>0, its
remaining vertices are in R. The indicator vector of these vertices has
Rayleigh quotient for X equal to

\[
 \frac12+\frac{\ell-1}{2\sqrt p}-\frac{\ell}{2p}.
\]

Using (5) therefore gives

\[
 \boxed{|C|\le a+
 \frac{1+2\sqrt p\,c\cosh(\log(B+2m)/k)}{1-1/\sqrt p}.}
 \tag{6}
\]

For ell=0, the trivial bound |C|=a suffices. This direct Rayleigh argument
uses neither a degree estimate for the localized graph nor a Weil estimate
for m. Uniform control of (AG) over all a-anchor cliques bounds every
larger clique by (6); cliques of size less than a already satisfy that bound.

For fixed 0<sigma<1, take a=floor(sigma log_2 p). If an even sequence
k/log p tends to infinity and (AG) holds uniformly with log(B+2p)=o(k),
then (6) proves

\[
 \omega(G_p)\le a+O\bigl(p^{(1-\sigma)/2}\bigr).
 \tag{7}
\]

Proving this aggregate hypothesis for every fixed sigma<1 would imply
omega(G_p)=p^{o(1)}. No such hypothesis is proved here. Even that clique
conclusion would not by itself establish the arbitrary-two-set Paley
conjecture, the thin subgroup target, or the prize connection.

At degree a=1, (AG) does hold unconditionally with B=p: D=2E-I and
(sqrt(p)S-J)/p are both orthogonal involutions, so W/p is orthogonal
and every power has trace of magnitude at most p. The resulting clique
bound remains at the square-root scale. This elementary case calibrates
the criterion but does not extend it to a>=2, where D is no longer an
involution and the missing exponential cancellation appears.

## 4. The exact weighted-necklace expression

For Z_0 a subset of the anchors, let D_(Z_0) have diagonal
product_(z in Z_0) chi(x-z). Since I_0 is a clique, its common-neighbor
indicator satisfies the exact identity

\[
 E=\frac1b\sum_{Z_0\subseteq I_0}D_{Z_0}-\frac12 E_{I_0},
 \qquad
 D=\sum_{\varnothing\ne Z_0\subseteq I_0}D_{Z_0}
       -\frac b2 E_{I_0}.
 \tag{8}
\]

At an anchor the product of all adjacency factors is 1/2; this is the
reason for the last term. It must not be omitted. Substituting (8) in W
expresses the aggregate (AG) using all nonempty necklaces, the anchor
projections, and the rank-one J terms. Equation (3) retains the complete
expression instead of discarding these terms separately.

For comparison, the term in tr(W^k) containing only sqrt(p)S and no
anchor projections is

\[
 p^{k/2}\sum_{Z_1,\ldots,Z_k\ne\varnothing}
 \Sigma(Z_1,\ldots,Z_k).
\]

Even hypothetical individual bounds
|Sigma|<=K^k p^((k+1)/2) would give, by the triangle inequality, only

\[
 \frac{p^{k/2}\sum |\Sigma|}{[p\sqrt{b-1}]^k}
 \le\sqrt p\,[K\sqrt{b-1}]^k.
 \tag{9}
\]

At fixed a>=2 and K>=1 this is exponential in k; it does not establish
log(B+2p)=o(k). This calculation concerns the strength of the available
upper estimate, not a lower bound on the actual aggregate. Cancellations
between the necklaces and the retained correction terms remain possible
and unproved. Thus the finite-depth results cannot be promoted to (AG)
by simply increasing the number of separately bounded words.

## Verification

[The exact verifier](../experiments/parallel2_spectral_transfer_2026_09_04.py)
checks (4), (8), the projection identities, and the clique Rayleigh
numerator using integer pairs representing u+v sqrt(p). It also includes
arbitrary coordinate projections to test the two-projection statement
independently of localization. [Results](../results/parallel2_spectral_transfer_2026_09_04.json)
record every case and source hashes. These tests check the algebraic
reduction; they do not supply the uniform bound (AG).

The final run checked 43 coordinate projections in five fields: 344
exact trace recurrences, 28 complete anchor expansions, and 122 clique
Rayleigh numerators. Empty projections and maximal allowed ranks are
included. The verification uses no floating-point acceptance condition.
