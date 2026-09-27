# Exact lattice distances and subgroup periods

The full Paley conjecture and the Proximity Prize remain unproved. This
pass gives a constant-preserving reformulation of the subgroup period
target and a finite certificate that avoids high additive moments. It
does not supply the missing uniform estimate. The argument below is
ordinary mathematics with exact finite checks; separate-author review
and Lean formalization remain outstanding. No novelty claim is made.

## 1. Exact centering

Let p be an odd prime and A=-A be a nonempty subset of F_p^*, of size
n=2N. Choose one representative h_j from each pair {h,-h}. For a in F_p,
let r_j(a) be the unique integer in [-(p-1)/2,(p-1)/2] congruent to ah_j.
Set

    V(a) = sum_j r_j(a)^2/p^2,
    vbar = n(p+1)/(24p),
    eta(a) = sum_(h in A) exp(2 pi i ah/p).

The sum eta is real. Averaging over a!=0 gives exactly

    average V(a) = vbar,
    average eta(a) = -n/(p-1).

Indeed, multiplication by h_j permutes F_p^*, and
sum_(a!=0) r_j(a)^2 = p(p^2-1)/12. Additive orthogonality gives the
second identity. Define, on the nonzero frequencies,

    d(a)=V(a)-vbar,   f(a)=eta(a)+n/(p-1),
    D=max |d(a)|,     M_f=max |f(a)|,     M=max |eta(a)|.

For the operator formulas only, extend d(0)=f(0)=0. In particular,
d(0) is NOT V(0)-vbar, and f(0) is NOT eta(0)+n/(p-1).
The triangle inequality gives |M-M_f|<=n/(p-1)<=1.

## 2. An absolutely convergent inverse

The classical Bernoulli Fourier series, shifted by 1/2, gives

    ||x||_(R/Z)^2 = 1/12 + (1/pi^2) sum_(k>=1) (-1)^k cos(2 pi kx)/k^2.

This follows from [NIST DLMF 24.8.1](https://dlmf.nist.gov/24.8.E1)
with B_2(x)=x^2-x+1/6. The series is absolutely convergent, including
at half-integers. Pairing h with -h and subtracting the exact nonzero
frequency average therefore yields

    2 pi^2 d(a) = sum_(k>=1) c_k f(ka),    c_k=(-1)^k/k^2.       (1)

For multiples of p, the uncentered summand eta(ka)=n is constant and
disappears under centering. This is why the zero extensions above are
needed. Centering at N/12 instead of vbar would make (1) false.

Write k=2^t m with m odd, and let mu denote the Mobius function. The
Dirichlet-convolution inverse b of c is

    b_m = -mu(m)/m^2                         (t=0),
    b_(2^t m) = -mu(m)/(2^(t+1) m^2)         (t>=1).            (2)

To prove it, let u_k=1/k^2 and delta_j be the sequence supported at j
with value 1. Then

    c = -u * (delta_1 - delta_2/2),
    u^(-1)_k = mu(k)/k^2,
    (delta_1-delta_2/2)^(-1) = sum_(t>=0) 2^(-t) delta_(2^t).

Multiplying the last two sequences gives (2); the terms involving
mu(2m)=-mu(m) combine when t>=1. All series converge absolutely in
the l^1 Dirichlet-convolution algebra.

The elementary Euler product and the
[values of zeta(2) and zeta(4)](https://dlmf.nist.gov/25.6.E1) give

    sum_(m odd) mu(m)^2/m^2 = 12/pi^2,
    sum_k |b_k| = 18/pi^2,
    sum_(p does not divide k) |b_k| = 18 p^2/[pi^2(p^2+1)].      (3)

The last factor removes the odd prime's squarefree Euler factor
1+p^(-2). If U_k v(a)=v(ka), then U_k U_l=U_(kl). Absolute convergence
allows rearrangement of the double sum in (1), giving the exact inverse

    f(a) = 2 pi^2 sum_(k>=1) b_k d(ka).                        (4)

## 3. Comparison with no exponent loss

Taking absolute values in (1) and (4) proves

    12 p^2/(p^2-1) D <= M_f <= 36 p^2/(p^2+1) D.              (5)

For the forward inequality use
sum_(p does not divide k) 1/k^2 = (pi^2/6)(1-p^(-2)).

If 2A=A, in particular if A is a subgroup containing 2, then
f(2a)=f(a) and d(2a)=d(a). Grouping c_k by its odd part gives

    c_m + sum_(t>=1) c_(2^t m) = -2/(3m^2),   m odd.

Since sum_(m odd,p does not divide m) m^(-2)
=(pi^2/8)(1-p^(-2)), (5) sharpens to

    24 p^2/(p^2-1) D <= M_f <= 36 p^2/(p^2+1) D.              (6)

For every 1<=r<infinity the same two-sided comparisons hold between
the normalized l^r norms on F_p^*: use Minkowski and the fact that
U_k is a permutation when p does not divide k. Thus the reformulation
also preserves moment estimates at a selected growing degree.

In particular, an absolute bound

    D <= C sqrt(n log(p/n))

in the intended quartic subgroup family would imply the sought
square-root period scale, up to an absolute change of constant and
the additive shift at most 1. Conversely, that period bound implies
the same scale for D. The constant and actual sponsor parameters still
require the separate checks in [the target statement](subgroup-target.md).
No classical two-set quadratic-character reduction or official
extension-field prize transfer is supplied here.

## 4. The lattice and what its elementary norm does prove

For arbitrary representatives define

    L = {z in Z^N : sum_j z_j h_j = 0 mod p}.

This lattice has index p in Z^N, and

    L^* = Z^N + (1/p) Z(h_1,...,h_N).

The right-hand side lies in L^*, has p distinct cosets modulo Z^N,
and has the same covolume as L^*. The p cosets of Z^N in L^* are
therefore indexed by a in F_p. Their shortest squared norms are V(a):
each coordinate can be centered independently. The missing estimate
asks that every nonzero coset's shortest squared norm lie in a narrow
annulus around n(p+1)/(24p).

Now assume n=2N is a power of two, N>=2, and H=<g> has order n.
Set A=H and fix a in F_p^*. This domain was made explicit after the
[independent pass39 review](parallel39-independent-shell-review-2026-09-06.md).
Use representatives g^j, 0<=j<N, in that order, with their actual
signs after centering. Let

    F_a(X)=sum_(j=0)^(N-1) r_j(a) X^j,
    R=Z[X]/(X^N+1).

The polynomial X^N+1 is irreducible over Q: its translate by 1 is
Eisenstein at 2, since N is a power of two. In F_p it splits into
distinct roots g^k for odd k. At those roots,

    F_a(g^k) = a sum_(j=0)^(N-1) (g^(k+1))^j mod p.

This is zero at all N-1 roots except k=n-1, where it is aN!=0.
Multiplication by F_a in R modulo p consequently has rank one.
Its integer determinant, the cyclotomic norm of F_a, is nonzero and
divisible by p^(N-1). For the divisibility, integer row and column
operations lifting Gaussian elimination modulo p leave N-1 rows
divisible by p; equivalently use Smith normal form.

Parseval over the N primitive n-th roots of unity gives

    (1/N) sum_(k odd) |F_a(exp(2 pi i k/n))|^2 = sum_j r_j(a)^2.

The arithmetic-geometric mean inequality therefore implies

    p^(N-1) <= |Norm(F_a)| <= (sum_j r_j(a)^2)^(N/2),
    V(a) >= p^(-2/N).                                         (7)

There is also the exact congruence sum_j r_j(a)^2=0 mod p, by summing
the geometric progression g^(2j). Thus V(a) is a multiple of 1/p.

The limitation is quantitative: for growing N with p of quartic size
in n, the lower bound in (7) tends to 1, whereas vbar is about N/12.
It does not control deviations on the sqrt(n log n) scale, or even
provide the required dimension-scale lower center. This is a limit
of this norm argument, not a theorem excluding all lattice methods.

For the eligible example p=1153,n=8,g=75,a=1, the ordered centered
power coordinates are [1,75,-140,-123]. They have

    sum r_j^2=40355=35p,  V(1)=35/1153,
    Norm(F_1)=p^3,       vbar=1154/3459.

The mean is more simply kept as n(p+1)/(24p)=1154/3459.
The squared norm is only about 3.1% above the lower bound from (7).
This is one finite example, not a growing-family extremality theorem.
Permuting or changing signs preserves V, but generally changes the
polynomial norm; the ring calculation must use the ordered powers.

## 5. A rational finite certificate from truncated inversion

Choose L>=1, an integer Q>=1, and for 1<=k<=L, p not dividing k, set

    w_k = trunc_toward_zero(Q b_k),  W=sum_k |w_k|,
    Z=max_(a!=0) |sum_k w_k d(ka)|,
    R_L=(D W-Z)/Q >= 0.

Truncation toward zero keeps the sign of each residual coefficient.
Using (3) to bound all omitted coefficients, (4) gives

    M <= 36 p^2/(p^2+1) D - 2 pi^2 R_L + n/(p-1).              (8)

To see this without guessing a tail estimate, the retained operator
has norm at most Z/Q on this particular d. The total omitted l^1
mass is exactly 18 p^2/[pi^2(p^2+1)]-W/Q. Multiplying it by D and
2 pi^2 yields (8). Any proven rational P<=pi^2 may replace pi^2 in
the negative term, because R_L>=0. Every remaining quantity is rational.
The same estimate holds on a subset of frequencies if Z is maximized
only on that subset; the omitted term still uses the global D.

The implementation uses Q=2^40 and P=9869604401089/10^12. Its proof
P<pi^2 uses exact alternating bounds in Machin's identity, independently
of machine floating point. Cosine intervals use interval evaluation
of the degree-40 Taylor polynomial and remainder at most x^42/42!.

For a subgroup, let q=(p-1)/n and choose a primitive field generator u.
One loop through its powers accumulates

    C_j=sum_(x in u^j H) centered_residue(x)^2,
    d(u^j)=[12 C_j-np(p+1)]/(24p^2).

The exact mean is checked before inversion. Multiplication by a small
k is a cyclic shift of the q quotient indices. Therefore the algorithm
uses O(p+L p/n) bounded-integer arithmetic operations and O(p/n+L)
storage. This operation count is not a bit-complexity claim. Explicit
resource and integer-range checks restrict the C++ executable's inputs.
It does not compute any high additive moments or invoke Lean.

## 6. The known order-64 example, with attribution preserved

For p=6700417 and H=<2> of order 64, the L=4096 certificate proves

    M < 43.832,
    max_(a not in H) |eta(a)| < 39.838.

The L=256 version already gives M<44.151 and all other cosets below
40.150. A separate rational Taylor enclosure gives

    43.802482797626 <= eta(1) <= 43.802482797627.

Both truncations thus identify the unique maximizing coset as H.
The larger integer certificate was computed in about 0.63 seconds in
the recorded run; this is an observed runtime, not a universal estimate.

The project had ALREADY identified the maximum and enclosed it to
width 10^(-18) in the [period-polynomial certificate](period-polynomial-certificate.md).
The new result is an alternative derivation without those signed
moments, and a certified bound on the other cosets. It does not improve
the project's best enclosure of M. An initial comparison only with the
older sqrt(1970) moment bound missed that stronger prior result; that
comparison is corrected here and in the pass summary.

The [verifier](../experiments/parallel36_verify_2026_09_05.py) checks
4096 inverse coefficients by a separate divisor recurrence, general
symmetric sets with rational trigonometric intervals, all cosets in
12 small-field C++ comparisons, 28 dyadic norm cases, and direct
residue calculations at the large certificate's extremizers. It makes
41574 exact assertions. Most assertions concern the finite partition
and determinant arithmetic, not independent instances of the analytic
theorem. Its full run took about 25.4 seconds. This is same-author
verification by different computational paths, not independent review.

## 7. What remains for the actual goal

The positive moment upper bound in pass35 remains open. This pass
offers a geometric way to express the same missing cancellation:
prove an absolute thin-annulus estimate for these highly structured
dual-lattice cosets as n grows. The inverse is bounded uniformly, so
it does not itself lose an exponent. The elementary norm only supplies
(7), which is far too weak. Neither finite certificates, an average
distance identity, nor the lattice's index proves that every coset
has a typical distance. No period exponent, classical Paley bound,
or official Proximity Prize theorem is proved by this pass.
