# Exact carry and norm-compression identities

All statements concern R=Z[X]/(X^N+1) with N>=4 a power of two, n=2N, p an odd prime with n|(p-1), and g of exact order n. Conjugation sends X to X^-1. Centered residues lie in[-(p-1)/2,(p-1)/2]. No principal-ideal assumption is needed.

## The distinguished split-root support

Modulo p,

```
F_a(X) = a*sum_(j=0)^(N-1) g^j X^j.
```

At any root t of X^N+1, (g*t)^N=1. The geometric sum is zero except when t=g^-1, where it equals a*N. Since p does not divide N, F_a has exactly one nonzero split-root value for a!=0. The multiplication matrix modulo p therefore has rank one. Its determinant is divisible by p^(N-1), giving the integer norm defect tau. The cyclotomic field is totally imaginary with conjugate pairs of embeddings, so the norm is positive for nonzero F_a.

## A ternary carry cocycle

The identity

```
1+X+X^-1 = (X^3-1)/(X*(X-1))
```

shows that epsilon has norm one: X->X^3 is a field automorphism, the numerator and X-1 have the same norm, and Norm(X)=1. Multiplication by epsilon is an integral invertible map. Its coefficient action is the sum of the original vector and its two signed neighbouring shifts. Each output coefficient has magnitude at most3(p-1)/2.

Modulo p, epsilon acts on the distinguished root by m=1+g+g^-1. Thus recentering epsilon*F_a gives F_(ma). Their difference is p*C_a. The centered reduction of a number in this range subtracts only -p,0 or p; hence every carry coefficient lies in{-1,0,1}. For two successive steps,

```
epsilon^2*F_a = F_(m^2*a) + p*(epsilon*C_a + C_(ma)).
```

If raw=epsilon*F_a and next=raw-p*C_a, then

```
||next||^2 = ||raw||^2 - 2p*<raw,C_a> + p^2*||C_a||^2.
```

This exactly locates the nonlinear section change. Norm(raw)=Norm(F_a), but subtraction of p*C_a has no general norm-monotonicity property. The finite witnesses test that failure directly on actual cosets. The carry mechanism is a standard choice-of-section cocycle for an integral toral action; its name is not a novelty claim.

## Height reduction with a kernel relation

Take any nonzero f in R with f(g^-1)=0 modulo p. All split-root values of f*F_a vanish, so f*F_a is zero in R/pR and D=f*F_a/p is integral. Also p divides Norm(f), so write Norm(f)=k*p with k a positive integer. Multiplicativity gives

```
Norm(D) = Norm(f)*Norm(F_a)/p^N = k*tau(a).
```

No assumption that p does not divide k is needed. Each coefficient in a negacyclic convolution is a signed sum. Therefore

```
|D_i| <= ((p-1)/(2p))*sum_j |f_j|.
```

Taking the floor gives the integer digit-height bound. In particular, coefficient absolute sum3 forces ternary digits. The reduction saves coefficient height when f is short; the cost of finding f, its norm cofactor and its possible dependence on n remain explicit.

For a fixed f this encoding is injective: f*F=p*D has a unique solution in the cyclotomic field. This suggests an exact realizability checker. If adj(f) satisfies f*adj(f)=Norm(f)=k*p, then

```
F = adj(f)*D/k.
```

Integral centered coefficients and the congruences F_j=a*g^j modulo p characterize whether D belongs to the actual centered section. A norm constraint alone does not impose those conditions. The complete N8 comparison is now implemented and separately checked:257 of6561 ternary words are actual encodings, including zero. Even knowing the exact actual norm-energy pairs admits945 words.

## The exact lattice and section

Let I be the inverse image in R of the one-dimensional subspace spanned by F_1 modulo p. Then [R:I]=p^(N-1), and pR is contained in I. Define J=(f/p)I. The split-root argument above shows J is an integral lattice, and determinants give

```
[R:J] = Norm(f)*[R:I]/p^N = k,
fR is contained in J,
[J:fR] = p.
```

Use the half-open coefficient cube C=[-p/2,p/2)^N. Because p is odd, its integer points have precisely the centered coefficients used here. The set I intersect C is one representative of each of the p classes I/pR. Multiplication by f/p sends C to the parallelotope P=f[-1/2,1/2)^N, which is a fundamental domain for the lattice fR. Consequently the actual digit words, including zero, are exactly

```
J intersect P,
with one representative for each of the p classes J/fR.
```

This is classical lattice geometry specialized to the checked relation. It is not an additional small-index estimate: J has index k, which can be large, and P can be skewed. Replacing P by its coordinate bounding box is precisely the relaxation that introduces6304 extra small-cube points in the N8 example.

When p does not divide k, f has exactly one vanishing split-root evaluation modulo p: its multiplication determinant has p-adic valuation one, so its reduction has nullity one. Thus an integral inverse F satisfying fF=pD automatically lies in I. When p divides k this argument does not apply; the separate scalar-congruence checks are necessary. The f=p control makes that distinction concrete.
