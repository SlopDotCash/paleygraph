# Centered subgroup energies with the origin retained

**The Paley conjecture and the Proximity Prize remain unproved.**
This is an ordinary mathematical proof from a published incidence theorem.
It has an author audit and bounded computational checks, but no independent
review or Lean verification. No literature novelty is claimed.

Let p be an odd prime and H<=F_p^* a multiplicative subgroup of order n.
Write r_s=1_H^{*s}, E_s=sum_x r_s(x)^2, and

    T_s = E_s - n^(2s)/p = ||r_s-n^s/p||_2^2.

The main result of this note is an absolute constant D, independent of s,
p and H, such that for every integer s>=1,

    T_(2s) <= D n^(2s-1/2) T_s.                         (R0)

Unlike the [pass-28 full-energy recurrence](parallel28-moment-recurrence-2026-09-05.md),
this has no additive uniform term. It removes the principal-term cap in
that recurrence. The resulting deficit still grows too slowly to prove
the requested square-root bound. The mixed-moment exponent obtainable
from the recorded seeds remains 71/72 in the quartic window.

## 1. An exact-centered invariant-set estimate

Let Q be a nonempty H-invariant subset of F_p^*, with m=|Q|>=n, and let
E(Q)=#{a+b=c+d : a,b,c,d in Q}. We prove

    E(Q)-m^4/p <= K m^3/sqrt(n),                         (S0)

where K=3 C_R and C_R>=1 is an absolute incidence constant.

Use points P=Q x Q x H and planes

    x+h y-v z=u,        (h,u,v) in H x Q x Q.

Both sets have size N=n m^2. The planes are distinct, since the
coefficient of x is normalized to one. For each h,z in H, multiplication
by h and z permutes Q, so the incidence count is exactly

    I(P,Pi)=n^2 E(Q),       N^2/p=n^2 m^4/p.             (I)

Any affine line has an injective projection onto at least one coordinate.
Thus at most m of these product-set points lie on a line.

The sole non-elementary input is [Shkredov 2019, Theorem 8, equation (22)](https://msp.org/moscow/2019/8-1/moscow-v8-n1-p03-s.pdf),
printed page 20: for balanced point and plane sets in F_p^3, with maximum
collinearity k, I<=C_R(N^2/p+N^(3/2)+kN). The archived published PDF and
its exact hash are in the [pass-28 source ledger](../results/parallel28_source_scope_2026_09_05.json).
The present use needs this theorem, rather than subtracting the uniform
term from an upper bound whose coefficient on that term is unspecified.

For N<=p^2, N^2/p<=N^(3/2), and k<=m. Divide the incidence bound by n^2
and use n>=1 to obtain E(Q)<=3 C_R m^3/sqrt(n). Subtracting m^4/p
only improves this upper bound.

For N>p^2 we use the following elementary centered incidence bound.
For any M distinct affine planes in F_p^3, let m_d count the planes of
each parallel direction d and put F(x)=#{planes containing x}-M/p.
Each plane has p^2 points, two nonparallel planes intersect in p points,
and distinct parallel planes are disjoint. Expanding the square gives

    sum_(x in F_p^3) F(x)^2 = p^2 M-p sum_d m_d^2 <= p^2 M.

Cauchy--Schwarz therefore gives |I(P,Pi)-|P|M/p|<=p sqrt(|P|M).
In the balanced case this is at most pN<N^(3/2). Equation (I) then gives
E(Q)-m^4/p<=m^3/sqrt(n). This proves (S0) in both size ranges.
This variance calculation is supplied here in full; no theorem is
imported from the Vinh abstract located during the literature search.

## 2. Weighted functions and the exceptional value at zero

For a real function f on F_p, define

    E0(f) = sum_z (f*f)(z)^2 - (sum_x f(x))^4/p
          = p^(-1) sum_(a!=0) |f_hat(a)|^4.

Consequently E0^(1/4) is a seminorm. It is not monotone under taking
absolute values, so signed weights require a separate argument.

First let g>=0 be H-invariant and supported outside zero. Write
A=||g||_1, B=||g||_2^2, and P_t={x:g(x)>t}. Empty level sets contribute
zero. The layer decomposition and (S0) give

    E0(g)^(1/4) <= K^(1/4) n^(-1/8) integral |P_t|^(3/4) dt.

For A,B>0, |P_t|<=min(A/t,B/t^2). Splitting the integral at B/A yields

    integral |P_t|^(3/4) dt <=6 A^(1/2) B^(1/4),
    E0(g)<=1296 K A^2 B/sqrt(n).                        (W+)

The zero function is immediate. For real signed g supported outside
zero, write g=g_+-g_-. Both pieces are nonnegative and H-invariant.
The seminorm triangle inequality, (u+v)^4<=8(u^4+v^4), and
A_+^2 B_+ + A_-^2 B_- <= A^2 B give

    E0(g)<=10368 K ||g||_1^2 ||g||_2^2/sqrt(n).          (W+-)

Finally for arbitrary real H-invariant f, write f=g+f(0)1_{0}.
Here E0(f(0)1_0)=|f(0)|^4(1-1/p). Hence

    E0(f)<=82944 K ||f||_1^2 ||f||_2^2/sqrt(n)
             +8 |f(0)|^4.                              (W0)

The last term is essential in this general statement. For example,
f=1_0-1/p is H-invariant for H=F_p^*, with

    E0(f)=1-1/p,  ||f||_1=2(1-1/p), ||f||_2^2=1-1/p.

The ratio E0(f)/(||f||_1^2||f||_2^2/sqrt(n)) is
sqrt(p-1)/(4(1-1/p)^2), which tends to infinity. Thus (W0) does not
silently reintroduce the false origin-free assertion discussed in the
[earlier source qualification](analytic-bounds-and-amplification.md).

## 3. Applying the bound to actual centered convolutions

Put f_s=r_s-n^s/p. Its H-invariance follows by multiplying every
summand of r_s by the same element of H. Also

    ||f_s||_1<=2n^s,    ||f_s||_2^2=T_s,
    f_s*f_s=f_(2s),    E0(f_s)=T_(2s).

The convolution identity follows by expanding and using sum r_s=n^s.
For the origin, choosing the first s-1 summands determines the last,
so r_s(0)<=n^(s-1). Since n<=p, also n^s/p<=n^(s-1). Both numbers
are nonnegative, so

    |f_s(0)|<=n^(s-1),
    |f_s(0)|^2<=T_s,
    |f_s(0)|^4<=n^(2s-2) T_s.

Applying (W0) and using K>=1, n>=1 proves (R0) with

    D=331784 K,      331784=4*82944+8.

No constant here depends on s. Iteration gives, for j>=0,

    T_(s 2^j) <= D^j n^[2s(2^j-1)-j/2] T_s.            (Rj)

In particular, the cited E_3<=A n^4 log n bound gives, for r=3*2^j
and n in the quartic window,

    T_r <= A D^j n^(2r-2-j/2) log n.                    (R3)

The new proof establishes this centered consequence directly. It does
not validate every claim or intermediate step of the previously
qualified published centered-function argument.

## 4. What the improvement can and cannot give for the period

Let eta(a)=sum_(h in H) exp(2 pi i a h/p), M=max_(a!=0)|eta(a)| and
Delta=M/n. The already recorded centered mixed-moment inequality is

    Delta^(kl) <= 2/n + sqrt(p T_k T_l)/n^(k+l),          (G0)

for integers k>=2 and l>=1. For clarity, here is the origin bookkeeping.
Set A(x)=p^(-1) sum_t |eta(t)|^k exp(-2 pi i tx/p), a real, possibly
signed function. Then sum A=n^k, ||A-n^k/p||_2^2=T_k, and

    sum_y r_l(y)|eta(ay)|^k
      = sum_(x,y) (A(x)-n^k/p)(r_l(y)-n^l/p) exp(2 pi i axy/p)
          +n^l A(0)+n^k r_l(0)-n^(k+l)/p.

For k>=2, |eta|<=n and the exact second moment imply A(0)<=n^(k-1).
Also r_l(0)<=n^(l-1). The centered bilinear sum is at most
sqrt(p T_k T_l) in absolute value. Finally multiplicative invariance
gives sum_y r_l(y)eta(ay)=n eta(a)^l; Holder with weight r_l gives
n^k |eta(a)|^(kl)<=n^[l(k-1)] sum_y r_l(y)|eta(ay)|^k.
These facts prove (G0). The 2/n argument is not asserted for k=1.

Ignore logarithms only for the following power comparison. Define a
deficit d_s by T_s << n^(2s-d_s). The seeds are

    d_1=1, d_2=31/20, d_3=2.

Doubling now increases the deficit by 1/2 without a cap. Interpolation
of the nonzero Fourier moments supplies the concave linear envelope
through (1,1), (2,31/20), and (3*2^j,2+j/2) for every j>=0.
This envelope is closed under doubling and interpolation: doubling
is exact on [3,infinity); on [1,3], checking the linear breakpoints
1, 3/2, 2 and 3 gives an increment of at least 1/2. Thus it contains every bound
generated solely by those operations and seeds.

For k,l>=2, the second term of (G0) supplies saving
(d_k+d_l-4)/(2kl), with the first term also imposing 1/(kl).
On each pair of linear intervals, the former expression is bilinear
in 1/k and 1/l, so it is bounded by its endpoint values. All infinitely
many endpoint pairs reduce to the following three cases:

- k=3*2^i and l=3*2^j: the saving is
  (i+j)/(36*2^(i+j)), whose maximum is 1/72 at i+j=1 or 2.
- k=2 and l=3*2^j (or the transpose): the saving is
  (10j-9)/(240*2^j), whose maximum is 11/960 at j=2.
  It increases up to j=2 and decreases thereafter, by direct subtraction
  of consecutive terms.
- k=l=2: the saving is -9/80.

The first-term restriction in (G0) does not decrease the maximum at
(3,6), (3,12), (6,3), (6,6), or (12,3). The best power supplied by
this centered envelope and gate is therefore still 71/72.
This is a limitation of these estimates, not a theorem that the actual
period cannot satisfy a stronger bound or that all methods fail.

The direct coset estimate n M^(2r)<=p T_r gives from (R3) a saving
(j-2)/(12*2^j), maximized at 1/96 for j=3,4. At r of order log n,
the deficit in (R3) is only O(log r), while a square-root moment bound
would need T_r of order (C r n)^r. The gap in the power of n is still
of order r. Removing the uniform term has not removed this gap.

## 5. Feeding the amplitude bound back into the moments

There is another available operation, so the preceding logarithmic
deficit calculation should not be mistaken for a bound on everything
the existing amplitude estimate can supply. If M<<n^(1-delta), then
Fourier positivity gives

    T_(s+t)<=M^(2t) T_s.

Also exact coset repetition gives, for s,t>=1,

    T_(s+t)<=(p/n) T_s T_t.                              (P)

Indeed T_u=(n/p) sum_C |eta_C|^(2u); the diagonal terms in the product
of the two nonnegative coset sums contain the sum for s+t.
In the quartic window, the exponent rule in (P) is
d_(s+t)>=d_s+d_t-3. Feeding delta=1/72 back into the moments gives
d_(s+t)>=d_s+t/36. These operations improve some of the higher moment
bounds, but do not bootstrap the amplitude saving beyond 1/72.

Here is an explicit closed envelope proving this assertion. Let G(s)
equal the preceding deficit envelope on [1,24], and set

    G(s)=s/36+17/6                 for s>=24.              (F)

It is concave, with slopes everywhere at least 1/36. It is attained
as an exponent bound: beyond 24, use T_24 and the existing amplitude
estimate. Its closure under every operation just listed follows from:

- Interpolation: G is concave.
- Doubling (R0): G(2s)>=G(s)+1/2. Below s=12 this follows from the
  earlier envelope; on [12,24] the difference is s/72+1/3; for s>=24
  it is s/36.
- Amplitude feedback: G(s+t)>=G(s)+t/36, by its minimum slope.
- Coset product (P): G(t)-t/36<=17/6<3, so the same slope inequality
  gives G(s+t)>=G(s)+G(t)-3.

For the centered mixed gate, use the endpoint set {2,3,6,12,24,infinity}.
On every pair of pieces the saving expression is bilinear in 1/k,1/l.
At infinity with the other order k fixed its limit is 1/(72k); with
both orders tending to infinity it is zero. The finite endpoint maximum
is again 1/72, at the same five pairs identified above. The direct coset
bound has saving (G(r)-3)/(2r); on the tail it equals
1/72-1/(12r), approaching 1/72 from below.

Thus the seed estimates, interpolation, (R0), coset products, the direct
coset bound, the centered mixed gate, and feedback of its best amplitude
estimate reach a fixed point at saving 1/72 in this power ledger.
This is a closure statement about those upper-bound rules only. Constants
and logarithms are suppressed for fixed-order power comparisons; it is
not a proof of a sharp bound for actual energies, of a numerical constant,
or of an obstruction to other inequalities. With feedback, the deficit
does grow linearly at high orders, but its slope is only 1/36, far below
the slope near one needed for square-root cancellation.

## 6. Boundary with the official extension-field prize

The incidence input above is a prime-field statement. Replacing p
everywhere by the size q of an extension field would make (S0) false
with a uniform absolute constant. For example, embed Q=H=F_p^* in
F_(p^6), so m=n=p-1. Its additive energy is

    E(Q)=m^3-m^2+m,
    [E(Q)-m^4/p^6]/[m^3/sqrt(n)]
       =sqrt(n)(1-1/m+1/m^2-m/p^6),

which tends to infinity through primes p. This is the familiar subfield
obstruction, not a new counterexample to Paley or to the prize. It
prevents using this argument as an automatic transfer to the official
extension-field parameters. Arbitrary classical Paley input sets also
lack the large multiplicative invariance required by (S0).

## Verification and remaining input

The [verifier](../experiments/parallel29_verify_2026_09_05.py) independently
counts small point--plane incidences, checks the exact variance identity,
and checks centering, convolution, origin and exponent identities.
Floating-point Fourier checks are labeled separately. Its
[output](../results/parallel29_verification_2026_09_05.json) is bounded
evidence and does not formally verify the analytic theorem.

The next input must improve the estimates beyond the closed envelope
(F), or control the relevant signed spectral quantity by a different
argument. Deeper iteration and the listed feedback operations are
exhausted by the calculation above. The uniform square-root subgroup estimate,
arbitrary two-set Paley estimate, and scalar prize bounds remain open.
