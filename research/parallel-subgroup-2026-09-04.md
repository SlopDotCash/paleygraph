# Index-two projection bounds for mixed moments

**Status: neither Paley target is proved.** The new results in this pass are
an improved mixed-energy inequality, its extension at every convolution
depth, and two elementary restrictions on the balanced energy. They improve
how a hypothetical mixed-energy estimate propagates, but do not prove that
estimate, improve an exponent, or establish the positive-product condition
(PM+) or centered-moment condition (CM). These are ordinary proofs with exact
finite checks; no novelty or formal-verification claim is made.

The starting point is the notation and exact recurrence in
[dyadic-descent-and-mixed-energy.md](dyadic-descent-and-mixed-energy.md).
The high-moment target is stated in
[positive-product-moments.md](positive-product-moments.md).

## A sharper fourth-energy inequality

Let p be an odd prime, H a dyadic multiplicative subgroup of size 2k,
k≥2, K its subgroup of size k, and L=gK its other coset. Put

\[
 f=1_K*1_K,\qquad \widetilde f=1_L*1_L,\qquad w=1_K*1_L,
\]

where convolution is additive. Write

\[
 E=E_2(K)=\sum_x f(x)^2,\quad
 B=\sum_x w(x)^2=\sum_xf(x)\widetilde f(x),\quad
 T=\sum_x f(x)w(x).
\]

The last identity counts the three-in-K, one-in-L zero sums because both
sets are closed under negation. Multiplication by g exchanges f and
\(\widetilde f\), while leaving w invariant. Therefore

\[
 T=\sum_x u(x)w(x),\qquad
 u=(f+\widetilde f)/2,\qquad
 \sum_xu(x)^2=(E+B)/2.
\]

Since u(0)=k and w(0)=0, Cauchy–Schwarz on the nonzero elements gives

\[
 \boxed{2T^2\le B(E+B-2k^2).} \tag{1}
\]

The earlier bound was \(T^2\le B(E-k^2)\). The right side of (1),
after division by two, is smaller by exactly \(B(E-B)/2\). The saving
comes from projecting onto the H-invariant part before applying
Cauchy–Schwarz; it does not rely on a conjectural cancellation.

The exact recurrence remains

\[
 E_2(H)=2E+6B+8T.
\]

Consequently, suppose \(B\le Ck^2\Lambda\), with C,Λ≥1, at every
step of a fixed dyadic tower. Then at every level s,

\[
 \boxed{E_2(H_s)\le(7+4\sqrt3)C\Lambda s^2<14C\Lambda s^2.} \tag{2}
\]

To prove this, let D=7+4√3. The base group has energy 6. If
E≤DCΛk², drop the negative term in (1), obtaining
T≤CΛk²√((D+1)/2). The recurrence is then at most

\[
 C\Lambda k^2[2D+6+8\sqrt{(D+1)/2}]
 =4DC\Lambda k^2,
\]

because √((D+1)/2)=1+√3. This proves the induction. It replaces the
previous constant 22 with a constant below 14. **The hypothesis
B=O(k² log p) is still unproved; the change of constant does not discharge
it or move the asymptotic frontier.**

## Removing the residual mean as well

Let d=p−1. The nonzero sums of u and w are k²−k and k². Subtracting
their means on the nonzero field elements strengthens (1) to

\[
 \left(T-\frac{k^2(k^2-k)}d\right)^2
 \le
 \left(\frac{E+B}{2}-k^2-\frac{(k^2-k)^2}d\right)
 \left(B-\frac{k^4}d\right).
 \tag{3}
\]

This is an exact covariance bound. Both factors are nonnegative, since
they are sums of squares. In the quartic window, removing this mean does
not itself supply the missing order of magnitude.

## The same projection works at every convolution depth

For j≥1 write \(r_j=1_K^{*j}\), and set

\[
 E_j=\sum_xr_j(x)^2,\qquad \nu_j=r_j(0),
\]

\[
 Z_{a,b}=\#\{(x_1,\ldots,x_a,y_1,\ldots,y_b)\in K^a\times L^b:
                 \textstyle\sum x_i+\sum y_i=0\}.
\]

Thus \(Z_{j,j}=\sum_xr_j(x)r_j(gx)\). For j,s≥1 define

\[
 \begin{aligned}
 U_j&=\frac{E_j+Z_{j,j}}2-\nu_j^2
             -\frac{(k^j-\nu_j)^2}d,\\
 V_s&=Z_{2s,2s}-Z_{s,s}^2
             -\frac{(k^{2s}-Z_{s,s})^2}d,\\
 C_{j,s}&=Z_{j+s,s}-\nu_jZ_{s,s}
             -\frac{(k^j-\nu_j)(k^{2s}-Z_{s,s})}d.
 \end{aligned}
\]

Then

\[
 \boxed{U_j\ge0,\quad V_s\ge0,\quad C_{j,s}^2\le U_jV_s.} \tag{4}
\]

**Proof.** Set \(u_j(x)=(r_j(x)+r_j(gx))/2\) and
\(w_s=1_K^{*s}*1_L^{*s}\). Multiplication by g fixes w_s. Hence
\(\langle r_j,w_s\rangle=\langle u_j,w_s\rangle=Z_{j+s,s}\).
Their squared norms, values at zero, and total masses are respectively

\[
 \begin{array}{c|ccc}
       &\text{squared norm}&\text{value at zero}&\text{mass}\\\hline
 u_j&(E_j+Z_{j,j})/2&\nu_j&k^j\\
 w_s&Z_{2s,2s}&Z_{s,s}&k^{2s}.
 \end{array}
\]

These identities follow from convolution associativity, negation
invariance, and the fact that multiplication by g acts as an involution
on K-invariant functions. Restrict to the nonzero elements, subtract each
mean, and apply Cauchy–Schwarz. ∎

A useful specialization controls a three-to-one mixed moment without
asking for moments of higher total degree. Put j=2r and s=r. Since
\(\nu_{2r}=E_r\), dropping only the residual mean gives

\[
 \boxed{(Z_{3r,r}-E_rZ_{r,r})^2
 \le\left(\frac{E_{2r}+Z_{2r,2r}}2-E_r^2\right)
          (Z_{2r,2r}-Z_{r,r}^2).} \tag{5}
\]

All moments on the right have total degree at most 4r, the degree of
Z_(3r,r). Formula (4) supplies its stronger centered version. At r=1,
(5) is exactly (1), since E₁=k and Z_(1,1)=0.

This controls an unbalanced mixed count in terms of a pure moment and
balanced counts. It does **not** bound the balanced quantity
\(pZ_{q,q}-k^{2q}\) appearing in (CM). The all-depth formula therefore
sharpens propagation once such input is available, without creating that
input.

## Two additional unconditional restrictions

First,

\[
 \boxed{B\le E-k.} \tag{6}
\]

Indeed f(x) is odd exactly when x∈2K: off-diagonal ordered pairs occur
in swapped pairs, and the diagonal pair is unique. The analogous
statement for \(\widetilde f\) uses 2L. Hence f−\(\widetilde f\)
is odd at each of the 2k elements of 2H, so

\[
 2(E-B)=\sum_x(f(x)-\widetilde f(x))^2\ge2k.
\]

Second, w is nonnegative, vanishes at zero, has mass k², and is constant
on each H-coset of size 2k. If its coset values are c_i≥0, then
\(\sum_i c_i=k/2\). Therefore

\[
 \boxed{B=2k\sum_i c_i^2\le2k(\sum_i c_i)^2=k^3/2.} \tag{7}
\]

This is cubic and does not give the needed near-quadratic bound.

There is a parity extension of (6). If q is a power of two, repeated
squaring in the group algebra over F₂ gives
\(r_q(x)\equiv1_{qK}(x)\pmod2\). Thus

\[
 \boxed{E_q-Z_{q,q}\ge k\qquad(q\text{ a power of two}).} \tag{8}
\]

The saving k is negligible compared with the high moments required by
(CM). It records a forced antisymmetric component, rather than asserting
that this component is large enough.

## Exact checks and their limits

Run `python3 experiments/parallel_subgroup_2026_09_04.py`. The script uses
Python integer counters and exact rational arithmetic. It checks:

- 52 fourth-energy cases, including six with nonminimal child energy;
- 608 all-depth inequalities and independent convolution identities, with
  1≤j≤8 and 1≤s≤4 across 19 field/subgroup pairs;
- 76 parity identities at depths 1, 2, 4, and 8;
- the two quartic resonances (p,n)=(6700417,64) and (67403009,128).

For the first resonance, (E,B,T)=(2976,1024,96): the old upper bound
for T² is 1,998,848, and (1) gives 999,424. For the second,
(E,B,T)=(12096,4608,0): the bounds are 36,864,000 and 19,611,648.
The inequalities remain far from equality in these examples. Dense
cases such as (p,n)=(17,16), where (E,B,T)=(264,256,224), ensure the
checks are not restricted to minimum-energy children or vanishing T.

The result file records the script hash. The finite checks support the
implementations and arithmetic identities, not any uniform Paley bound.
