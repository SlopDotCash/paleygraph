# Fixed deletion pair, exact two-insertion moments

Let S be symmetric, have zero diagonal and entries ±1 off the diagonal, and satisfy S1=0 and S²=qI−J. The matrix adapter checks these identities. The fast backend constructs S from the quadratic character of a prime q≡1 mod4. Fix C, |C|=n, and distinct a,c∈C. Put A=C\{a,c}, O=Fq\C, m=q−n≥2 and N=m(m−1). All sums over insertion pairs below are ordered and exclude b=d. Each final set occurs twice.

For degree D, define elementary symmetric functions of a row restricted to A:

```
t0 = Σx eD(S[x,A])
g(x) = e(D−1)(S[x,A]); h(x) = e(D−2)(S[x,A])
f = Sg; v = Sh; w = S(g h); K = S diag(h) S.
```

Negative elementary degrees mean zero. Multiplication of the two insertion factors gives the exact identity

```
TD(A∪{b,d}) = t0 + f(b) + f(d) + K(b,d).
```

The computational input consists of t0, G0=Σg, G2=Σg², H0=Σh, H2=Σh², Q=gᵀSh, the selected entries of g,h,f,v,w on C, and K restricted to C×C. Set

```
F = −ΣC f
F2 = qG2−G0²−ΣC f²
H = H0−ΣC h
H2O = H2−ΣC h²
FH = Q−ΣC f h
KsumO = Σi,j∈C Kij
K2O = (q²−2q)H2+H0²
       −2Σi∈C [q(H2−hi²)−vi²] + Σi,j∈C Kij²
L = −qΣC w + G0ΣC v + Σi,j∈C fi Kij.
```

These are the outside totals Σf, Σf², Σh, Σh², Σfh, ΣK, ΣK², and Σb,d∈O fb Kbd, respectively. The identities follow from S²=qI−J:

```
K1=0
||K||F²=(q²−2q)H2+H0²
Σj Kij²=q(H2−hi²)−vi²
Kbb=H0−hb
(fᵀK)j=q wj−G0 vj.
```

Writing Zbd=fb+fd+Kbd, subtract its diagonal before forming the moments:

```
Z1 = Σb≠d∈O Zbd
   = 2(m−1)F+KsumO−mH0+H

Z2 = Σb≠d∈O Zbd²
   = 2(m−2)F2+2F²+K2O+4L
     −mH0²−H2O−4H0F+4FH+2H0H

ΣY = N t0+Z1
ΣY² = N t0²+2t0 Z1+Z2.
```

Divide by N and subtract the square of the mean for the variance. Every final operation in Python uses unbounded integers or exact fractions. The literal pair census and the older marked-moment compiler with inclusion-exclusion independently check the exclusions.

Only the term 4FH contains Q; its variance coefficient is exactly 4/N. Centering g and using the conference identity gives

```
Q² ≤ (qG2−G0²)(qH2−H0²)/q.
```

Consequently the variance is `variance_without_Q + 4Q/N`, with a certified squared radius stored as `Q_free_variance_radius_squared`. The value at Q=0 is an algebraic baseline, not an assertion that another graph realizes Q=0. At D=6,n=6, g is identically zero, so this contraction vanishes exactly. That case cannot demonstrate its necessity for larger inputs.

The first moment has a useful further reduction. For R(x)=Σy∈C Sxy and iC(x) the indicator of C, its ordered numerator is

```
Σx [N eD(S[x,A]) −2(m−1)R(x)e(D−1)(S[x,A])
     +(R(x)²−(m−1)−iC(x))e(D−2)(S[x,A])].
```

This expression needs no global convolution. Its efficient compilation for every deletion pair is a next experiment, not an implemented interface in this round.

## Arithmetic and implementation bounds

The accepted prime backend domain is q≤10,000,000, n≤64, 0≤D≤min(6,n), 2≤n≤q−2. A has at most62 elements. Thus |g|≤C(62,5)=6,471,002 and |h|≤C(62,4)=557,845 for the largest degrees; lower degrees fit these bounds. Selected transforms, K entries and t0 fit signed64 bits. Squares, products and bilinear totals use signed128 bits; final moment products use Python integers. The Python dense adapter is limited to q≤257, where its matrix intermediates also fit signed64 bits. Independent NumPy review reduces bounded4096-row blocks to Python integers before accumulating across the field.

The producer uses two exact NTT primes, 2013265921 and1811939329, with roots31 and13. Both support transform length2²⁴, required at q=6,700,417. Linear convolution is folded modulo X^q−1. Each signed Sh coordinate is reconstructed before taking its dot product with g: reconstructing Q itself from only these two primes would have insufficient general range. The product exceeds twice the certified coordinate bound (q−1)max|h|. The producer checks every selected Sh coordinate, its full zero sum and its full squared norm.

The separate reviewer uses coefficient counts instead of row DP, decimation in frequency instead of the producer's decimation in time, computes Sg instead of Sh, and reconstructs the final scalar using a third prime469762049. Its CRT product is1,713,652,354,748,588,808,931,901,441, larger than twice q(q−1)max|g|max|h| throughout the stated domain. This makes its scalar reconstruction exact, not a modular spot check. The row pass costs O(qn²), the transforms O(q log q), and storage O(q+n²). No q×q matrix is allocated by either fast convolution backend.
