# Prime-field uncertainty: an exact gap and its quantitative limitation

Status: a prime-field nondegeneracy result with an explicit spectral gap,
and a proof that the determinant certificate developed here is too weak
at the required neighborhood sizes. No improved clique bound or full
Paley/prize proof is claimed. Standard Fourier uncertainty is credited
to its existing sources; no novelty claim is made for that theorem or
for concentration of interval Fourier projections.

## 1. The prime-field input

Let p≡1 mod4 be prime, r=(p−1)/2, S the Paley sign matrix, and

\[
P=\tfrac12(I-J/p+S/\sqrt p),\qquad
Q=\tfrac12(I-J/p-S/\sqrt p).
\]

These are the two rank-r Fourier projections, onto the nonzero
quadratic-residue and nonresidue frequencies. Let C⊂F_p have
1≤m=|C|≤r and write

\[
X=P_{C,C},\quad Y=Q_{C,C},\quad
D=X+Y=I_m-J_m/p,\quad \Delta=\det D=1-m/p>0.
\tag{1}
\]

Chebotarev's theorem says that every square minor of the prime-order
Fourier matrix is nonzero. We use the primary
[Tao paper, Lemma 1.3 and Corollary 1.4](../sources/tao-uncertainty-math0308286v6.pdf),
also available as [versioned HTML](https://arxiv.org/html/math/0308286v6).
Selecting any m rows from either frequency set gives an invertible
m×m matrix. Consequently its full Gram matrix is positive definite:

\[
X>0,\quad Y>0,\quad I-X=Y+J/p>0.                    \tag{2}
\]

In particular 0<λ_min(X)≤λ_max(X)<1. This excludes exact endpoint
eigenvalues using prime field order. The source theorem is qualitative;
the numerical estimates below require additional arguments.

## 2. A positive integer from the two conjugate determinants

There is an exact arithmetic invariant

\[
k(C):=p^{m+1}\det X\det Y\in\mathbb Z_{>0}.         \tag{3}
\]

To prove integrality, set B=2pX=pI-J+√p S_C. Since J has rank one,
multilinearity gives

\[
\det B=(\sqrt p)^{m-1}
\left(\sqrt p\det(\sqrt p I+S_C)
      -1^T\operatorname{adj}(\sqrt p I+S_C)1\right).
\tag{4}
\]

The term in parentheses lies in Z[√p]. Therefore the product of
det B and its √p↦−√p conjugate is an integer divisible by p^{m−1}.
Its conjugate is det(2pY).

The power of two can be removed. For ω=(1+√p)/2, the entries of
pX belong to Z[ω]: the off-diagonal entries are ω−1 or −ω and
the diagonal entries are (p−1)/2. This ring is closed under products
because ω²−ω=(p−1)/4 is an integer. Its conjugate norm is integral:

\[
(a+b\omega)(a+b\bar\omega)
=a^2+ab+(1-p)b^2/4\in\mathbb Z.
\]

Thus det(pX)det(pY) is an integer. Multiplying it by 4^m gives
det B det\bar B, which is divisible by p^{m−1}. Since p is odd,
det(pX)det(pY) is also divisible by p^{m−1}. This proves (3);
strict positivity follows from (2).

## 3. An explicit uniform spectral gap

Define the positive definite contraction

\[
R=D^{-1/2}XD^{-1/2},\qquad I-R=D^{-1/2}YD^{-1/2},
\]

and the exact rational number

\[
\theta(C)=4^m\det R\det(I-R)
=\frac{4^m k(C)}{p^{m+1}\Delta^2}\in(0,1].          \tag{5}
\]

For every eigenvalue μ of R, all the other factors
4μ_j(1−μ_j) are at most one. Hence 4μ(1−μ)≥θ(C), and

\[
\tfrac12(1-\sqrt{1-\theta})I\le R
\le\tfrac12(1+\sqrt{1-\theta})I.
\]

Using D≥ΔI and I−X≥Y gives

\[
g(C)I\le X\le(1-g(C))I,\quad
g(C)=\frac{\Delta}{2}(1-\sqrt{1-\theta(C)}).          \tag{6}
\]

In particular g(C)≥Δθ(C)/4 and k(C)≥1 give the fully explicit
uniform inequality

\[
\boxed{\quad\delta_{p,m}I\le X\le(1-\delta_{p,m})I,
\qquad \delta_{p,m}=\frac{4^{m-1}}{p^m(p-m)}.\quad}  \tag{7}
\]

The stronger rational bound retaining k(C) is k(C)δ_{p,m}.
There is no asymptotic claim that these lower bounds approximate
the actual endpoint gap. At m=1, (6) is exact; at larger m it can
discard substantial information about the remaining eigenvalues.

## 4. Why this determinant certificate cannot yield the needed edge

The limitation can be proved for actual Paley matrices. Set

\[
A=D^{-1/2}S_C D^{-1/2};\qquad
R=I/2+A/(2\sqrt p).
\]

If a_1,…,a_m are the eigenvalues of A, then |a_i|<√p by (2) and

\[
\theta(C)=\prod_{i=1}^m(1-a_i^2/p).                 \tag{8}
\]

The exact second moment is

\[
\operatorname{tr}A^2
=m(m-1)+\frac{2\|S_C1\|^2}{p-m}
       +\frac{(1^TS_C1)^2}{(p-m)^2}
\ge m(m-1).                                        \tag{9}
\]

Indeed D^{-1}=I+J/(p−m), so expanding tr(S_CD^{-1}S_CD^{-1})
gives the displayed expression. The first term uses the actual zero
diagonal and ±1 off-diagonal entries of S_C.

Let t=tr(A²)/p. Arithmetic-geometric mean and log(1−u)≤−u give
the exact rational and exponential bounds

\[
\theta(C)\le(1-t/m)^m\le e^{-t}
\le \exp[-m(m-1)/p].                               \tag{10}
\]

Consequently the *certificate in (6)*, even using the exact integer
k(C), satisfies

\[
g(C)\le\tfrac12\Delta\theta(C)
\le\tfrac12\Delta\exp[-m(m-1)/p].                   \tag{11}
\]

For m/p→c>0 it tends to zero exponentially in p. This is a bound
on the numerical certificate, **not an upper bound on the actual
spectral gap**. It therefore does not assert that the actual Paley
compression approaches an endpoint.

For our normalization b=2^a, N=b−1,
Z=(X−I/2)/(√N/b), an edge bound ‖Z‖≤1+o(1) would require
an endpoint gap approaching at least 1/2−√N/b, positive for a≥2.
The fixed-anchor common neighborhoods have m=p/b+O_a(√p), so
(11) shows that this determinant certificate cannot deliver that
positive constant. This conclusion applies to the certificate actually
derived here; it does not exclude using more information from minors,
determinant ratios, or arithmetic structure in another argument.

## 5. Prime-order Fourier nonvanishing also permits near concentration

For a complementary elementary example, write p=4h+1 and let
F={h+1,…,3h} be a set of r=2h Fourier frequencies. It is symmetric
under ξ↦−ξ and excludes zero. Its Fourier projection Π is real,
circulant, has rank r, and satisfies Π1=0. Every compression to
at most r coordinates and its complement are positive definite by
the same prime-order minor theorem.

Choose 1≤m≤r, n=m−1, C={0,…,n}, and the real vector

\[
v_j=(-1)^j\binom nj\quad(0\le j\le n),\qquad
\|v\|^2=\binom{2n}{n}.
\]

Under the unitary Fourier transform,

\[
\widehat v(\xi)=p^{-1/2}(1-e^{-2\pi i\xi/p})^n.
\]

For ξ outside F, |1−e^{-2πiξ/p}|²≤2, since its distance to
zero is at most h<p/4. Parseval therefore gives

\[
1-\frac{v^T\Pi v}{\|v\|^2}
\le\frac{2^n}{\binom{2n}{n}}
\le(2n+1)2^{-n}.                                   \tag{12}
\]

At m proportional to p this is exponentially small, even though
every relevant Fourier minor is nonzero. Moreover
S'=√p(2Π−I+J/p) has zero diagonal, S'1=0, and (S')²=pI−J.
Thus prime cyclic order, real translation invariance, the ambient
projection identities, and nonvanishing minors do not by themselves
give quantitative separation from the endpoints.

This interval-frequency family is **not a Paley family**. In fact
the geometric-series formula gives

\[
S'_{0,1}=\frac{1-2\sin(\pi r/p)/\sin(\pi/p)}{\sqrt p}
\sim-2\sqrt p/\pi,
\]

so its off-diagonal entries fail the Paley ±1 requirement for all
sufficiently large p. The small p=5 example happens to coincide
with a Paley projection; it does not affect this asymptotic distinction.
No Q(√p) entry condition is imposed, and C is not asserted to be
an anchor common neighborhood. This construction does
not preserve the elliptic character structure or sharp bounds from
the previous pass. No combined counterexample satisfying all Paley
constraints is claimed. Equation (12) is an illustration of the
qualitative theorem's limits, separate from the actual Paley result
(3)–(11).

## 6. Verification and the remaining target

The [exact verifier](../experiments/parallel17_prime_uncertainty_2026_09_05.py)
computes determinants over Z[√p], verifies positivity in both real
embeddings, the integer k(C), and an independent integer determinant
for θ(C). It checks (9) by matrix multiplication and (10) by rational
arithmetic. Selected matrices also receive exact certificates after
subtracting the rational gap. Boundary cases larger than rank are
checked to be singular. The interval example's binomial identity
and leakage certificates are checked separately. Finite computations
support, and do not replace, the uniform arguments above.

The prime-field input now supplies a proved explicit gap, but its
size is insufficient. The next argument must control how mass
distributes specifically across quadratic-residue frequencies, with
enough uniform strength to reach the required constant edge or
long-depth estimate. Neither support nonvanishing alone nor this
single product of all endpoint distances does that. Independent
mathematical review remains outstanding, and the full goal is open.
