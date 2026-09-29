# Anchored rows and exact arithmetic recognition of the difference profile

The full Paley and Proximity Prize goals remain unproved. Known entries
in the incidence matrix give a simpler obstruction to both earlier
uncentered interval models. Centering repairs that obstruction, and a
product congruence can also be satisfied formally. Power sums in the
actual prime field impose further constraints: the first n-1 determine
the entire difference profile. This is an exact recognition theorem,
not a bound on those coefficients or a solution of the conjecture.

## 1. Known entries constrain individual rows

Use the established notation C_ij, a_i=C_0i, and
h_i=sum_j C_ij(C_ij-1). The incidence symmetries give

    C_i0=a_i,   C_ii=a_(-i).

For i!=0 these are distinct entries in row i. Therefore

    a_i+a_(-i)<=n,
    h_i>=a_i(a_i-1)+a_(-i)(a_(-i)-1).                       (1)

At i=0 one has h_0>=a_0(a_0-1). As in the preceding pass,

    h_i=sum_j a_j a_(j+i)-(n-1)1_(i=0).                    (2)

These conditions require no estimate on X. In the previous long-filler
model a=A1_{1,...,d}+2 1_J+1_{d+1}, with n=t^5,A=t^3,d=t,
l=|J|=(n-2-Ad)/2 and J far from the short interval, direct counting gives

    h_d=4(l-d)+A < A(A-1)=a_d(a_d-1)                       (3)

for all dyadic t>=2. For t>=4, use h_d<=2t^5+t^3 and
t^6-t^3>2t^5+t^3; t2 is a direct check. With the Sidon filler one has
instead h_d=A, so the same violation is immediate. Thus both uncentered
families fail actual anchored-row conditions, even though the long filler
passed the weaker total-X subset conditions. The previous necessary
inequality remains valid; (1) is a simpler exclusion for these examples.

## 2. The exact product and power sums

Let p be an odd prime, let H have even order n dividing p-1, let
m=(p-1)/n, and choose q of order m in F_p*. Index the cosets by

    H_i={x!=0:x^n=q^i},
    a_i=#{h in H\{1}:(h-1)^n=q^i}.

Suppose n=2^nu. The single odd coefficient of a is at the class gamma
of2, so q^gamma=2^n. This parity assertion follows from the involution
h->h^(-1) on H\{1}; its only fixed point is -1, whose difference is -2.
The derivative of X^n-1 at1 gives product_(h!=1)(1-h)=n. Since -1 is
in H, taking quotient classes proves

    sum_i i a_i == nu gamma  (mod m).                       (4)

There are also exact identities in F_p, for every integer k>=1:

    sum_i a_i q^(ki)
      = n sum_(j=0,...,k) binom(nk,nj)  (mod p).             (5)

Expand sum_(h in H)(h-1)^(nk). The subgroup sum of h^r is n when
n divides r and0 otherwise. All surviving signs are positive because n
is even. The h=1 term vanishes, proving (5). In particular the first
two residues are2n and n[2+binom(2n,n)]. Equation (4) is only one
necessary congruence; satisfying it does not imply (5).

## 3. Recognition from n-1 power sums

Let u_i be nonnegative integers with sum u_i=n-1. Then u is the actual
profile a if and only if (5), with u in place of a, holds for
1<=k<=n-1.

Indeed, consider the degree n-1 monic polynomial

    P_u(Y)=product_i (Y-q^i)^(u_i).

The first n-1 power sums of its roots determine its elementary symmetric
coefficients by Newton's recurrence

    r e_r=sum_(j=1,...,r)(-1)^(j-1)e_(r-j)s_j,  e_0=1.

All r<=n-1 are invertible in F_p. Equal power sums therefore give
P_u=P_a. Unique factorization over F_p, with the q^i distinct, makes
their root multiplicities equal. The converse is (5).

This provides a complete arithmetic check for proposed profiles without
constructing the full m by m incidence matrix. It does not give an upper
estimate on sum a_i^r, shifted excess, W, or the spectral edge. Using
only a few of the congruences is a necessary filter, not this complete
recognition criterion.

## 4. A repaired family and the limits of the coarse constraints

Here is an abstract repair that shows (1) and (4) do not themselves yield
a power saving. Let t>=4 be dyadic, n=t^5,A=t^3, nu=log_2 n, and

    I={1-t,...,t},   l=(n-2-2At)/2.

Work in a prime cyclic group with n^3<m<2n^3. For an integer R in
{10n+1,...,12n}, choose gamma modulo m by

    (nu-1)gamma == At+2lR+l(l-1) (mod m),                  (6)

and require gamma outside I and J_R={R,...,R+l-1}. Set

    a=A1_I+2 1_(J_R)+1_gamma.

Such R exists. The map R->gamma is injective modulo m, so at most2t
positions put gamma in I. For each j in {0,...,l-1}, the equation
gamma=R+j has at most one solution R: its coefficient2l-(nu-1) is
nonzero modulo m. Thus at most l+2t<2n choices are forbidden. Both
intervals are disjoint and have no wrap; gamma may lie elsewhere in the
cyclic group. The unique odd label is gamma. Since sum_(i in I)i=t,
(6) proves the formal product congruence (4).

Now a_0=A=max a. For i!=0, the two correlation terms using index0 give
r_i>=A(a_i+a_(-i)), which implies (1). The zero-row condition also holds:
A2>=2t A^2 and n=t^5 ensure A2-n+1>=A(A-1). The two known row entries
have total at most2A<=n.

The mass is n-1 and A2<=4n^(7/5). The exact dyadic value convention
gives Q=2n^2. The singleton and whole-group coset mass caps hold with
constant1. Expanding the nine component convolutions gives

    K<=9[10n^3+16n^(13/5)+4n^(7/5)+4n+1]<=315n^3.          (7)

For |u|,|v|<=t/2, at least t indices remain in all three translates of
I. There are at least t^2 such pairs, so the triple-correlation norm
still satisfies ||T||_(3/2)>=n^(61/15).

For all t>=64, the total-X subset inequalities from the preceding pass
hold with a formal X=2n^2. The same proof applies with h_max<=4n^(7/5):
a positive numerator requires R_D<16n^(4/5)<n/2, forcing at least n/4
central shifts of the long filler into D. Explicitly, for
0<|i|<=n/8, h_i>=4(l-|i|)>=3n/2-4t^4-4>=n/2; also
h_0>=4l-n+1=n-4t^4-3>=n/2. These are at least n/4 shifts.
The total h mass is (n-1)(n-2)<n^2, while
6^(1/3)(n/4)^(2/3)(2n^2)^(2/3)=(3/2)^(1/3)n^2>n^2.
Thus the cubic term alone suffices whenever the numerator is positive.
A budget X_model=384n^2 also satisfies divisibility by6, the compared
diagonal-excess constraint, and K<=n[2(n-1)^2-(n-1)+X_model].

These are formal scalar and anchored-row constraints. No local X_D
allocation, full multiplication law, field embedding of the labels, or
incidence matrix realizing this family is supplied. In particular gamma
is only a proposed class of2 until an actual field embedding is checked.

## 5. One certified quartic field rejects the tested repaired profiles

The exact example is

    n=1024, m=1094909953, p=1121187791873=nm+1,
    n^4<p<2n^4.

Both primalities have short certificates in the checker. Write
m=16707*2^16+1. The residue5^((m-1)/2)=-1 mod m forces every prime
divisor of m to be1 mod2^16, hence larger than sqrt(m). Thus m is prime.
For p, 2^(p-1)=1 mod p and gcd(2^n-1,p)=1 force every prime divisor
of p to have multiplicative order divisible by m, hence to exceed
sqrt(p). This proves p prime without a probable-prime assumption.

At t4 all2,048 positions R in the stated interval give disjoint odd
labels under (6). Since 2^n!=1 and m is prime, the proposed class gamma
uniquely fixes q=(2^n)^(gamma^(-1) mod m). Direct geometric sums show that
none of these profiles has sum_i a_i q^i=2n mod p. This is exhaustive
only over this specified range at this one n,p. No uniform rejection of
the repaired family follows.

The actual order1024 subgroup has a0=0 and max a=2, so it also directly
differs from these centered high-interval profiles. Its product identity,
first three power sums, and1,024 anchored-index checks pass. The finite
sieve is an illustration of the arithmetic filter, not a new uniform
bound or evidence of a proportion of exceptional primes.

## 6. Validation and next step

The [checker](../experiments/parallel47_anchored_arithmetic.py) records
the exact prime certificates, direct subgroup counts, model tests, and
all2,048 modular rejection residues through a reproducible hash. The
first three geometric-sum evaluations are independently checked term by
term. For (p,n)=(97,8),(353,16),(278177,32), Newton reconstruction from
all n-1 target power sums agrees with the directly multiplied root
polynomial, covering53 moments and three full reconstructions.
Results are in the
[certificate](../results/parallel47_anchored_arithmetic_2026_09_06.json).

Root completed the ordinary derivations and exact checks. No separate
agent, Lean, external review, or novelty claim is made. The next task is
to extract a uniform bound from the arithmetic moment conditions or the
full off-diagonal multiplication law. Merely satisfying the coarser
constraints or recognizing the true profile does not bound it. All
existing uniform exponents and full goals remain unchanged.
