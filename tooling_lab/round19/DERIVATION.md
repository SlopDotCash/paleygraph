# Exact limits of relation bounds and a projection codec

Throughout the field arguments, p is odd prime, N>=4 is a power of two, g has order2N, u=g^-1 mod p, and R=Z[X]/(X^N+1). Write m=(p-1)/2 and F_a for the centered coefficient vector congruent to a*(1,g,...,g^(N-1)). Every nonzero f in R with f(u)=0 mod p gives the integral encoding E_f(a)=fF_a/p. Its coefficient bound is floor(m*||f||_1/p). These facts were established in round17.

## 1. Fixed windows lose the sign memory

First use Q=2^N+1 as the modulus, without assuming it prime. For a ternary digit vector D, extend D_(j+N)=-D_j and put

```
W_(s,j) = sum_(t=0)^(s-1) 2^(s-1-t)*D_(j+t).
```

The relation2^s-X^-s gives |W_(s,j)|<=2^(s-1). Inverting the generator-2 encoding gives integral coordinates F_j with F_(j+1)=2F_j-QD_j and F_N=-F_0. Telescoping yields W_(s,j)=(2^s F_j-F_(j+s))/Q. At s=N this is exactly F_j. Thus all depth-N windows are equivalent to centered inverse membership.

A collection of windows up to depth d forbids a pair of consecutive nonzero digits of the same sign at distance less than d, including the antiperiodic boundary. To see necessity, a window beginning at that pair's first digit and ending at the second has value of magnitude2^distance+1. Conversely, a window with alternating nonzero signs has magnitude at most its leading power, hence at most2^(s-1).

Without a parity condition, D=(1,-1,0,...,0) passes every proper window and first fails at depth N. Its inverse has F_1=-2^(N-1)-1, one unit outside the centered range. There is no uniform fixed proper depth that recognizes this family.

With odd nonzero support, depth N/2 suffices. Around the antiperiodic circle, the number of sign-changing gaps is odd. With odd support, the number of same-sign gaps is therefore even. If there are any, there are at least two. If depth N/2 windows pass, each such gap has length at least N/2, leaving no positive length for the at least one other gap. This is impossible. The signs must alternate globally, which together with odd support is exactly the actual Q-section from round18.

The threshold is sharp even after requiring digit sum1: put+1 at positions0 and N/2-1, and-1 at N/2, with zeros elsewhere. Its first failing window has depth N/2. This gives an explicit separation between a growing raw digit window and constant sign memory. It is consistent with the classical distinction between finite-type and sofic descriptions. The infinite witness family concerns the moduli Q; it does not assert infinitely many prime Fermat numbers.

## 2. A certificate against every small relation portfolio

Set a=m+1, so F_a[0]=-m. Form the noncentered vector V=F_a+p*e_0. Then V_0=m+1, all other coordinates remain centered, and V has the correct scalar congruences. Define

```
delta = m-max_(j>0)|V_j|,
B = delta+1.
```

The distinct powers g^j, 0<j<N, are never+1 or-1. Hence none of these coordinates equals+m or-m, so delta>=1 and B<p.

**Claim.** For every nonzero relation f with f(u)=0 mod p and ||f||_1<=B, every coefficient of fV/p is integral and satisfies the same digit bound floor(m*||f||_1/p) as actual centered encodings. This includes every number of such relations and all their coefficient inequalities.

Each coefficient of fV is a linear functional w·V. Its normal w is a signed permutation of f, satisfies sum_j w_j*g^j=0 mod p, and has the same L1 norm. If w is supported only at coordinate0, its nonzero coefficient must be divisible by p, contradicting ||w||_1<p. Otherwise write t=|w_0| and s=sum_(j>0)|w_j|. Then s>=1 and t+s<=B, and

```
|w·V| <= (m+1)t+(m-delta)s
        = m(t+s)+t-delta*s
        <= m(t+s).
```

The last step has the explicit nonnegative slack certificate

```
delta*s-t = (delta+1)*(s-1) + (B-t-s) >= 0.
```

Finally w·V is divisible by p because V has the correct scalar congruences. Dividing and rounding down gives precisely the coefficient bound. The verifier checks the arithmetic hypotheses, the polynomial identity in t,s, and the excluded axis case. This certifies a universal inequality-family obstruction, not just failure of the16 relations sampled by the experiment.

This obstruction concerns the stated L1-derived coefficient bounds. It does not rule out nonlinear membership algorithms, arbitrary arithmetic constraints, or every possible norm method. The squared coefficient energy changes by exactly p, since (m+1)^2-m^2=p.

## 3. A factor-two bracket for the separator cost

Let r be the least odd residue in the multiplicative subgroup H=<g>, excluding1. In each pair{h,-h}, exactly one representative in{1,...,p-1} is odd. If that representative is t, the centered magnitude of h/2 is(p-t)/2. It follows that

```
r = 2*delta+1.
```

Choose e with g^e=r. The relation q=r-X^-e has L1 norm r+1=2B. Its digit bound is(r-1)/2=delta. At the escaped V, its coefficient0 divided by p is(r+1)/2=delta+1, so it rejects. Thus the least L1 norm of a relation whose coefficient bound rejects this V lies in

```
[B+1, 2B].
```

This is a bracket, not an exact optimum. At the largest tested prime it is[11607942,23215882]. The budget depends on the subgroup's least nonidentity odd residue, so changing the order of the same subgroup generator does not improve this certificate.

## 4. Extract one scalar, then check the section

The previous obstruction suggests retaining the section operation explicitly. Fix f and calibrate

```
beta = E_f(1)(u) mod p.
```

For every scalar a, F_a-aF_1=pZ for an integral Z. Therefore

```
E_f(a)-aE_f(1)=fZ,
E_f(a)(u)=a*beta mod p.
```

If beta is nonzero, a candidate digit word D determines the scalar

```
a(D)=D(u)*beta^-1 mod p.
```

The complete criterion is D=E_f(a(D)). Necessity follows from the calibration identity. Sufficiency follows from the equality itself. The final comparison reconstructs the centered F_a and multiplies by the sparse f; it retains the nonlinear section that the relation bounds miss. The projection alone is not a membership certificate: the near-boundary counterexamples have exactly the same projected scalar as their centered partners.

This branch needs no field norm, adjugate or cofactor. It works even when p divides the norm cofactor: f=p has beta=N, which is nonzero. At the general N64 input, beta=272772816 mod2013265921. A word uses64 projection terms,64 reconstructed coordinates, and576 coefficient products because f has9 nonzero coefficients. These are arithmetic counts, not a measured runtime speedup.

## 5. A complete fallback when calibration vanishes

If beta=0, the fast projection cannot be inverted. Let Norm(f)=kp and let A=adj(f), so fA=kp. Define the single coefficient row

```
w=(A_0,-A_(N-1),...,-A_1).
```

For an actual digit word, w·D=k*F_a[0]. Compute r=w·D modulo kp in[0,kp). If r is not divisible by k, reject. Otherwise a=r/k belongs to{0,...,p-1}, and the same exact re-encoding comparison decides membership. This is complete because for actual D, r/k equals F_a[0] modulo p. False positives are excluded by re-encoding regardless of how arbitrary nonactual words project.

Only one adjugate row is used per word. Its construction and its potentially large integer weights remain explicit costs; no cofactor-residue state table is constructed. Squared relations and f=p^2 exercise this fallback. It is different from assuming that beta=0 means the encoding itself is invalid.

## 6. What remains missing

An exact membership test does not count or search the whole section. The general scalar space still has p elements, and no uniform norm envelope follows. The next tool should complete partially specified digit words, or aggregate over their completions, while retaining this exact projection/re-encoding criterion. Its completeness and cost need independent tests; a fast check of one supplied word does not establish fast enumeration or solve either prize.
