# General level concentration and additional matrix identities

The [supergroup-coset mass bound](parallel44-structured-levels-2026-09-06.md)
excludes the exact subgroup model that made the earlier functional
inequality sharp. It does not improve that functional power for arbitrary
level sets: the explicit interval family below obeys every new coset-mass
bound and still has the same correlation power. This is an abstract
function family, not actual subgroup difference data or a Paley
counterexample. Further arithmetic constraints remain available.

## 1. An interval family retaining the old functional power

Let t>=2 be a power of two and set

    n=t^5, A=t^3, d=t, I={1,...,d},
    l=(n-2-Ad)/2, B=2l^2+1,
    S={Bj+j^2+10d+1:1<=j<=l}.

Choose a prime m>=n^3 large enough that m>2max(S), and work additively
in G=Z/mZ. In fact 2max(S)<n^3 for these parameters, so one can choose
n^3<m<2n^3 by Bertrand's postulate. This does not assert that nm+1 is
prime. The sets I,S,{d+1} are disjoint and avoid the identity 0. Define

    a=A*1_I+2*1_S+1_{d+1}.

Then sum a=n-1, and exactly one positive value is odd. Moreover

    A2=sum a^2= A^2d+4l+1 <=3n^(7/5),
    max a=n^(3/5), Q=max_(dyadic levels) T^3|D_T|=n^2.      (1)

In particular the second-moment size is compatible with the current
upper input A2<<n^(29/20), with room in the exponent.

The filler S is Sidon. Equality of two pair sums first forces equality
of the sums of their indices, since B>2l^2, and then equality of the
unordered index pairs. There is no wrap modulo m. Therefore
E(S)=2l^2-l. With

    K=sum_q (sum_x a(x)a(x+q))^2,

the interval component alone gives K>=A^4(2d^3+d)/3>=(2/3)n^3.
For an explicit upper bound, expand a*a into its nine component
convolutions and use Cauchy on their L2 norms. The mixed energy satisfies
E(I,S)<=d^2 l. Consequently

    K<=9[n^3+8n^2+1+4n^(13/5)+2n^(7/5)+4n]
      <=180n^3.                                           (2)

Because G has prime order, its only subgroup cosets are singletons and
the whole group. Both satisfy the new mass bound with constant 1:

    (max a)^3<=n^2,
    (sum a)^3<=n^2 |G|.                                    (3)

Thus (1)-(3) retain the actual mass, parity pattern, second-moment scale,
weak cubic tail, quotient-energy scale, and all subgroup-coset mass caps
being compared. They do not impose all identities of actual difference
multiplicities.

Put T(u,v)=sum_x a(x)^2 a(x+u)^2 a(x+v)^2. If
|u|,|v|<=floor(d/4), the three intervals intersect in at least d/2
points. There are at least d^2/4 such pairs, and each has T>=A^6d/2.
It follows exactly that

    ||T||_(3/2)^3 >= A^18 d^7/128 = n^(61/5)/128,
    ||T||_(3/2) >=128^(-1/3)n^(61/15).                     (4)

The previous functional upper expression
K^(1/5)Q^(26/15) has the same power n^(61/15). Therefore these extra
subgroup-coset caps alone cannot improve that expression by a factor
n^epsilon for any fixed epsilon>0 on all abstract functions satisfying
the listed constraints. No incidence matrix rho is supplied for this
family, so (4) is not a lower bound for an actual weighted triangle W.
It does not rule out an improved estimate using more subgroup arithmetic.

## 2. A nonlinear identity that actual incidence matrices must satisfy

There is additional structure beyond those scalar constraints. Let H=-H,
write its quotient additively, and let

    C_ij=#{x in H_i:x+1 in H_j}, a_i=C_0i.

Define the simultaneous shift (S_k C)_ij=C_(i-k,j-k), and for a real
weight vector w let L(w)=sum_k w_k S_k C. Then

    C^2=n I-n e_0 e_0^T+L(a).                              (5)

One direct proof uses the exact subgroup-coset convolution algebra
already established in
[the mixed-period note](mixed-periods-and-shifted-energy.md). If A_i is
the characteristic function of H_i and delta is the additive identity,
then

    A_i*A_j=n 1_(i=j) delta+sum_k C_(j-i,k-i) A_k.

Compare the coefficient of A_j in (A_0*A_0)*A_i and
A_0*(A_0*A_i). Their respective values are

    n 1_(i=j)+sum_k a_k C_(i-k,j-k),
    n 1_(i=0,j=0)+sum_k C_ik C_kj.

Associativity proves (5). All terms, including the special identity
coefficient, are retained. In particular its diagonal gives the exact
row second moment

    sum_j C_ij^2=n-n 1_(i=0)+sum_k a_k a_(k-i).              (6)

## 3. Inner products of simultaneous shifts

The same actual matrices satisfy

    <C,S_k C>_F=n(n-1)+(p-2n)1_(k=0).                      (7)

The left side counts x,y!=0,-1 whose two ratios x/y and
(x+1)/(y+1) lie in the same specified coset H_k. If the ratios r,s
are distinct, they uniquely determine y=(s-1)/(r-s), x=ry.
For k!=0 every ordered distinct r,s in H_k works, giving n(n-1).
For k=0 omit r=1 or s=1, giving (n-1)(n-2), and add p-2 pairs x=y.
This is exactly (7). Hence for every real weight w,

    ||L(w)||_F^2
      =n(n-1)(sum_k w_k)^2+(p-2n)sum_k w_k^2.               (8)

These are exact identities, not new bounds on the spectral edge. For
b=a^2, the triangle still equals W=n*b^T L(b)b. The identities provide
constraints for further work on that quadratic form, but the present
argument does not turn them into a stronger uniform estimate.

## 4. Validation and limits

The [checker](../experiments/parallel44_structured_levels.py) verifies the
interval construction at n32 and n1024, including exact weighted quotient
energies, the filler Sidon property, every subgroup-coset cap, and the
lower bound from the displayed target pairs. The uniform family and its
energy bounds follow from the proof above, not those two examples.

It independently checks (5) at every matrix entry, (7) at every shift,
and (8) for three weight vectors in each of (p,n)=(97,8),(353,16),
(1153,8). The weights include a signed vector. Results are saved in
[the certificate](../results/parallel44_structured_levels_2026_09_06.json).

Root completed ordinary proof review and distinct exact calculations.
No separate-agent, Lean, or external peer review is claimed. The generic
interval example is not asserted to satisfy (5) or to have any actual
prime-field realization. The high-level nested-coset partition remains
an unproved sufficient input, and the full Paley and prize goals remain
unchanged.
