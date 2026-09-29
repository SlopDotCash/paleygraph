# Independent checks of the exact critical-size moment tool

The root review passes the exact q-polynomial, critical-exponent and inventory-free enclosure checks. The separately assigned subagent proof review stopped on a usage limit before creating any files; this report records the review actually completed locally. It is not an additional completed subagent review or a Lean formalization.

[review_polynomials.py](review_polynomials.py) imports no production compiler or interpolation helper. It uses the previously independent literal-column expansion from round5 and independently builds Paley matrices at q13,17,29. Every actual family/theta type at those orders, plus repeated-row templates, agrees with the new q-polynomials:19 triple cases and361 raw union coefficients. A separately implemented two-row column recurrence checks39 second-moment union coefficients. [Exact results](results.json).

The script also scans **every nonzero monomial**, rather than only the leading coefficient of each union polynomial. With n=alpha*q^(1/4)+O(1), a term q^d at union size k has weight d−3k/4. A measured trace power of order j has Hasse-supported weight at most1+j/2. The unique baseline maximum is q^10/216 at union9, weight13/4; the next weight is5/2. The second moment has unique maximum q^7/720 at union6, weight5/2; the next is7/4. All arithmetic terms have weight at most3/2. After fixing ordinary powers0–6 the remaining three directions have weights at most−1/4.

The finite interpolation budget has the required algebraic justification. Triple-column multiplicities are affine in q and theta. In the logarithm of the generating product, a q-dependent sign-power sum consumes at least two total target degrees, because every row is balanced and the bulk q contribution vanishes unless all appearing sign exponents are even. A theta-dependent contribution consumes at least one degree in each of three variables. Target total degree18 therefore permits q degree at most9, and q degree at most floor((18−3j)/2) at theta power j. Raw union degree k gives the additional multiplicity-degree bound k. The repeated-row templates obey the same two-degree lower bound. The pair target has total degree12 and degree bound6 before global row weights. These are formal rational polynomial statements, valid across the two residue branches; graph realizability of synthetic interpolation nodes is unnecessary.

The known universal theta6 coefficient and vanishing theta5 term reduce the two family polynomials to quartics. Five theta values and ten q values determine the largest relevant polynomials, with additional values used as checks. This makes the exact interpolation an identity certificate under the displayed degree argument. Numerical agreement at a few primes alone would not suffice.

The leading constants have an independent elementary count. Three labelled six-sets with nine distinct columns appearing twice split into three labelled disjoint blocks of size3, giving (n)_9/(3!)³. Identical six-sets give the leading second-moment count (n)_6/6!. Thus the constants are1/216 and1/720, and their standardized ratio is40*sqrt(5). The full polynomial scan excludes other row partitions; no independence assumption on actual character rows enters this conclusion.

The exact mean is `−(n)_6/[48(q−2)(q−4)]`. Its centering correction has weight2, below the raw third-moment remainder weight5/2. Squaring the mean does not affect the variance leading term. Therefore, uniformly for fixed alpha in compact subsets of positive real numbers,

```
variance = q*(n)_6/720 + O_alpha(q^(7/4)),
central third moment = q*(n)_9/216 + O_alpha(q^(5/2)),
skewness ~ 40*sqrt(5)/sqrt(q).
```

The retained finite-n ratio `(n)_9/(n)_6^(3/2)` has expansion1−27/(2n)+O(n⁻²). It explains why fitting an additional factor n to the early finite samples would be wrong. Neither these statements nor the first three moments imply a central limit theorem or a useful extreme-set bound.

## Sensitivity constants and finite enclosures

For arbitrary permitted trace arithmetic, the leading union6 term is proportional to `M6−(15q−85)M4`. Set x=theta², N=q−2 and m2=(q−3)(q+1). The function f(x)=x³−(15q−85)x² is concave on[0,floor(2sqrt(q))²] for q≥29. Its chord and Jensen bounds enclose the sum with exact second moment m2. Their asymptotic width is30q⁴. Multiplication by the inclusion and ordered-row factors gives raw sensitivity at most `(alpha^6/24+o(1))*q^(3/2)`. All other arithmetic terms have weight at most1. Dividing by the variance to power3/2 gives standardized sensitivity at most `(360sqrt(5)*alpha^−3+o(1))*q^−9/4`.

If ordinary theta moments through6 agree, only changes in A2,A3,A4 remain. The leading term is `(n)_7/q^5` times the change in `2q*A2−A4/3`. On0≤u=theta²/q≤4, the function2u−u²/3 has range[0,3]. The monochromatic mass is(q−5)/4, proving raw difference at most `(3alpha^7/4+o(1))*q^−1/4`, and standardized difference at most `(6480sqrt(5)*alpha^−2+o(1))*q^−4`. The fixed-ordinary-moment hypothesis cannot be removed from that stronger statement. The remaining raw terms have weight at most−1/2.

[review_enclosures.py](review_enclosures.py) reconstructs the finite query with binomial inclusion probabilities, independent polynomial evaluation, and the exact chord/Jensen endpoints. It imports no production helper. It independently proves primality of both large orders, checks every exact interval/range/coefficient field at q6,700,417,n50 and q2,013,265,921,n211, and verifies containment of six known exact moment controls. [Evidence](enclosure_results.json).

The large enclosures require no measured trace moments or character arrays. They rely on the proved coefficient identities, prime-Paley normalization, Hasse, conference identities and positive variance. Rounded standardized endpoints are display values; exact rational central endpoints and variance define the enclosure, and width squared is checked directly. Degenerate n=q or q−1 does not define skewness and is rejected by the producer's positive-variance assertion.

All sources and input hashes are recorded. The conceptual foundations—finite-population moments, parity diagrams, interpolation, Hasse and Jensen—are classical, as the producer's attribution states. The exact local coefficient identity and sensitivity tool are validated here; historical originality, graph-realizability of relaxed witnesses, a distributional limit theorem and either prize remain unestablished.
