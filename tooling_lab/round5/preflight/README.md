# Global third-moment compiler preflight

2026-09-05. The proposed bridge is sound. The tiny implementation already recovers the two required values without enumerating six-column sets:

| Matrix | Exact global E[T6(C)^3], |C|=6 |
|---|---:|
| Paley49 | 44555/35673 |
| Peisert49 | 6145/3243 |

[third_moment_preflight.py](third_moment_preflight.py) independently constructs both GF49 matrices, verifies the normalization maps, and compares normalized results against literal unnormalized row-triple calculations. [third_moment_preflight.json](third_moment_preflight.json) preserves all normalized trace/correlation histograms and the exact repeated-row contributions. This is a preflight for a new local instrument, not an invention of elliptic trace moments or a prize proof. Frozen rounds1–4 remain unchanged.

## Exact n=d=6 identity

Write `E(C)=sum_a product_(x in C) S_ax` for |C|=6. Interchanging finite sums gives

```
sum_(|C|=6) E(C)^3
 = sum_(ordered a,b,c) e6((S_ax*S_bx*S_cx)_x).
```

Let `K6(m,tau)` be e6 of a sign vector with m nonzero entries and sum tau. It is a standard elementary-symmetric/Krawtchouk coefficient:

```
K6(m,tau) = [z^6](1+z)^((m+tau)/2)(1-z)^((m-tau)/2)
 = [tau^6 -(15m-40)tau^4 +(45m^2-210m+184)tau^2
    -15m(m-2)(m-4)] / 720.
```

The code verifies the polynomial against the literal binomial coefficient formula. In particular K6 is even in tau.

For all-equal rows, the contribution is `-q*binom((q-1)/2,3)`. For exactly two equal rows there are `3q(q-1)` ordered triples. Their product vector is the remaining row with one additional position zeroed. Its positive/negative counts are `(q-3)/2,(q-1)/2` in either order, so their total contribution is `-3q(q-1)*binom((q-3)/2,3)`.

For distinct Paley rows normalize `(a,b,c)` to `(0,1,t)`, t≠0,1. The affine change may flip all three row signs, which changes their product sum by a sign but leaves K6 unchanged. Define

```
tau(t)=sum_x chi(x(x-1)(x-t)).
```

There are q(q-1) choices of the first ordered pair for each normalized parameter t. Therefore

```
E[T6(C)^3] = {
 -q*binom((q-1)/2,3)
 -3q(q-1)*binom((q-3)/2,3)
 +q(q-1)*sum_(t!=0,1) K6(q-3,tau(t))
} / binom(q,6).
```

No extra factor6 belongs in the normalized t-sum: it already parametrizes ordered triples. The separate unnormalized check enumerates unordered row triples and multiplies their contributions by6.

At q49 the shared repeated contributions are −99,176 and −12,496,176. The distinct contributions are30,060,912 for Paley and39,092,592 for Peisert, yielding the exact acceptance values above.

## Peisert qualification

The Peisert sign function is not a multiplicative quadratic character. Arbitrary affine scaling is not an allowed sign symmetry. The q49 normalization nevertheless works by an explicit semilinear map.

Use F49=F7[X]/(X²+1), generator g=2+X, and Peisert positive exponents0,1 modulo4. For a nonzero d=g^j choose exponent e=1 if j is even and e=7 if j is odd. Then

```
x -> d^(-e)*x^e
```

sends0,d to0,1 and changes the Peisert sign function by the constant sign psi(d). Translation supplies the general first point. The script checks every nonzero d and every field element, hence all48 normalization maps. This proves the required reduction for this actual matrix; it is not an assumption that all conference matrices have affine symmetry.

Its normalized tau is `sum_x psi(x)psi(x-1)psi(x-t)`, a signed triple correlation. It is **not** the trace of the Legendre elliptic family. Indeed the Peisert49 histogram includes tau=22, outside the Legendre Hasse interval[-14,14]. No elliptic bound should be imposed on that control.

| Even power sum over normalized parameters | Paley49 | Peisert49 |
|---|---:|---:|
| sum1 | 47 | 47 |
| sum tau² | 2300 | 2300 |
| sum tau⁴ | 251120 | 552176 |
| sum tau⁶ | 35149760 | 233600960 |

Thus retaining row triples demonstrably distinguishes a pair whose global Johnson L2 spectrum agrees. This does not imply that the third moment determines the full distribution, extrema, or useful uniform tail bounds.

## Generic n>6: sufficient inventory

For distinct rows a,b,c let e01,e02,e12 denote their mutual signs. The three zero-containing column patterns are

```
(0,e01,e02), (e01,0,e12), (e02,e12,0).
```

On the remaining q−3 columns, write the sign variables as A,B,C. Their Walsh sums are

```
sum1=q-3;
sumA=-e01-e02, sumB=-e01-e12, sumC=-e02-e12;
sumAB=-1-e02*e12, sumAC=-1-e01*e12, sumBC=-1-e01*e02;
sumABC=tau.
```

Walsh inversion determines all eight bulk type multiplicities: for a sign pattern u, its count is one eighth of the sum of these eight moments times the matching parity monomial in u. These formulas are verified for every distinct row triple in the tiny matrices, in addition to the normalized cases. Therefore the generic global third moment requires the joint normalized inventory `(tau,chi(t),chi(t-1))` in the Paley setting. A tau-only histogram is not justified for n>6 until a further elimination is proved.

For three row entries a,b,c at a column, use

```
H=(1+a*t1)(1+b*t2)(1+c*t3)-1,
product_columns(1+z*H).
```

Extract degree6 in each ti and apply the inclusion functional `z^k -> (n)_k/(q)_k`. It counts k distinct columns, even when one column contributes to all three factors. Use separate all-equal and exactly-two-equal row templates; the n=6 collapsed product shortcut cannot be substituted for those templates when n>6.

### Why only seven tau evaluations per sign class suffice

The parent's proposed stronger compression has a valid structural proof. At fixed q and mutual signs, take the formal logarithm

```
log product_y(1+z*H_y) = sum_y sum_(r>=1) (-1)^(r-1)*z^r*H_y^r/r.
```

The coefficient of `t1^i t2^j t3^k` in `H_y^r` is a universal integer times `a_y^i b_y^j c_y^k`. For bulk signs, only the case in which i,j,k are all odd depends on tau. Every tau-dependent logarithm term consequently has degree at least1 in **each** ti. In the exponential at most6 such factors can contribute to the target `(6,6,6)` coefficient. The fixed zero-boundary factors and final inclusion functional do not increase tau degree. The result is a polynomial of degree at most6 in tau, not18.

Seven exact values determine that polynomial. After normalizing e01=+1, four mutual-sign classes remain, so28 paired-type coefficient evaluations are sufficient in principle. This is interpolation with an independently justified degree bound, not a fit inferred from numerical agreement.

The interpolation nodes need not be actual elliptic traces. Seven admissible nonnegative integral synthetic type inventories may be used, with the lack of graph realizability stated explicitly. For a fixed mutual-sign class, integrality forces a congruence class modulo8 and nonnegativity bounds the available interval. Small q can provide fewer than seven admissible values. The safe fallback is direct evaluation at its actual finite set of tau values. Using generalized binomial powers outside the nonnegative domain would require a separate formal-series implementation and justification.

## Implementation constraints and practical architecture

1. **Reuse the convolution already available.** In the prime Paley case, `tau(t)=sum_x chi(x(x-1))*chi(x-t)` is one character convolution. Its full array is computable with the existing exact NTT infrastructure in O(q log q) arithmetic and O(q) memory. O(q²) is the elementary preflight cost, not an inherent cost of the proposed instrument. Keep modulus, transform-length, cyclic-wrap and integer-recovery guards. A prime-field cyclic backend does not automatically implement the additive group of an extension field.
2. **Keep signs and singular parameters explicit.** Omit t=0,1 from the distinct-row family; restore all repeated rows separately. The usual Frobenius trace is `a_t(q)=-tau(t)`, because the projective point count is `q+1+tau(t)`. Even moments ignore this sign, but generic weighted odd moments may not.
3. **Accumulate the sufficient statistics.** At n=6 only tau powers0,2,4,6 are needed. For generic n keep the four mutual-sign classes and powers tau^j through j=6. Retain the full tiny histograms as audit data. The Hasse support can reduce the number of Paley histogram bins, but it does not apply to Peisert.
4. **Use exact rational coefficients.** The union polynomial uses at most6 degrees in each of three variables and union degree at most18. A naive grouped convolution can still be expensive despite this bounded state space. Cache type powers and sign symmetries; evaluate the seven-node polynomial by exact arithmetic and compare an extra admissible node as a regression check. The proof, not that extra check, establishes its degree bound.
5. **Separate curve-family ensembles.** A formula summing elliptic isomorphism classes with automorphism weights cannot be inserted as an unweighted sum over Legendre parameters without multiplicity and exceptional-j corrections. A theorem for prime p does not automatically handle q=p^e. Preserve the parameter ensemble and moment order in source adapters.
6. **Acceptance ladder.** Keep the current n6 twins as exact acceptance fixtures. Independently enumerate small n7/n8 slices for the generic compiler. Verify full triple-type histograms before any moment compression. Only then scale prime-field tau arrays and record coefficient costs, moment size and whether the third-moment distinction remains quantitatively informative.

## Prior art and local overlap

Classical trace moments are not new. [Kaplan–Petrow, *Elliptic curves over a finite field and the trace formula*](https://arxiv.org/abs/1510.03980),2015 preprint/2017 publication, gives finite-field point-count power moments with prescribed subgroups through Hecke traces. Its abstract explicitly identifies the full2-torsion and Γ0(4) predecessors. [Gallagher–Li–Sweeting–Vassilev–Woo, *Generating functions for power moments of elliptic curves over Fp*](https://arxiv.org/abs/1807.00749),2018, packages known Birch/Ihara/Kaplan–Petrow moment formulas into rational generating functions. These sources provide established context; their full ensemble normalizations must be audited before becoming an exact computational adapter.

[Grove, *Hypergeometric moments and Hecke trace formulas*](https://link.springer.com/article/10.1007/s11139-026-01372-y), publishedApril6,2026, gives the Legendre point-count normalization in §1.1 and records the Ahlgren and Ahlgren–Ono Γ0(4)/Γ0(8) even-weight formulas in Theorem2.3. They express sums of Chebyshev-type trace polynomials through Hecke traces and can recursively supply exact low power moments in their stated prime-field ensembles. This is a possible independent arithmetic oracle, not a new theorem supplied by the proposed compiler. [The Stacks Project's Legendre-family section](https://stacks.math.columbia.edu/tag/03VA) provides the standard geometric family. [Brouwer–Martin's triple-intersection paper](https://arxiv.org/abs/2109.03654),2021, connects Paley triple intersection counts to the same cubic character sums.

The workspace already contains close components:

- [research/parallel2-classical-2026-09-04.md](../../../research/parallel2-classical-2026-09-04.md), §§1–2, derives an exact elliptic reduction of M4 and the full tau convolution and second-moment identities. [Its executable](../../../experiments/parallel2_classical_2026_09_04.py) has an `elliptic` function. The cubic trace vector and exact graph/elliptic bridges are therefore already local tooling.
- [research/parallel23-spectral-coupling-2026-09-05.md](../../../research/parallel23-spectral-coupling-2026-09-05.md) already uses Grove's theorem, square-restricted Legendre moments, and a source-audited Hecke identity. It explicitly distinguishes algebraically assigning a Hecke trace from computing it independently.
- The [round3 twin census](../../round3/pair_type_twins/README.md) already records the unequal third moments and proposes retaining row triples. It obtains those moments from full six-set distributions. In the bounded inspected Python/Markdown scope, no prior executable implementing the present normalized triple-row-to-global-third-moment bridge or the degree-six tau interpolation architecture was identified. That is not proof of absence elsewhere.
- Row-product expansion, Walsh type inversion, finite-population inclusion and graph contractions were already reviewed in rounds2–4. They should retain that attribution.

The specific potential contribution is a **certified global third-moment compiler** with trace/correlation backend adapters, a proved finite input budget, explicit repeated-row corrections, and actual twin acceptance tests. Its insight is that the missing third-order information is computable from a small arithmetic inventory and can be tested separately from the global L2 data. Neither a new classical trace-moment theory nor a uniform exceptional-set bound follows.

## Search log

All searches and primary-page reads occurred2026-09-05:

- `Legendre elliptic curve family even moments sixth moment Hecke trace formula Kaplan Petrow`
- `Legendre family trace distribution finite fields class numbers moments Birch sixth`
- `"Hypergeometric moments and Hecke trace formulas"`
- `"Paley" "third moment" character graph`
- Targeted primary reads: Kaplan–Petrow arXiv abstract and bibliographic record; Gallagher et al. arXiv abstract; Grove full publisher HTML, §§1.1,2.2/Theorems2.3–2.4; Stacks Legendre section; Brouwer–Martin abstract, with its paper already inspected in round4.
- Local `rg` over research notes, selected experiment implementations, the twin census, and the Proximity knowledge base for third moments, Legendre traces, Hecke, and elliptic moment terminology; followed by the focused reads listed above.

Replay the bounded acceptance prototype with `python3 tooling_lab/round5/preflight/third_moment_preflight.py`. It performs zero six-column subset enumerations. The separate full row-triple checks deliberately cost more than the proposed normalized backend, to audit the reduction independently.

An additional [root acceptance check](root_acceptance_check.json) uses fresh complete six-set enumeration at q13 and q17 and the independently verified full round3 histograms for both q49 twins. All four moments agree exactly. Its [separate implementation](root_acceptance_check.py) deliberately uses the older, more expensive viewpoint to verify the proposed reduction.
