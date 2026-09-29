# Independent audit of the trace-histogram envelopes

All eight saved cases pass an independent exact rational audit: four prime orders, each with the conference-only and fixed-ordinary-moment models. The reviewer reconstructs the finite support, equality constraints, objective, and actual inventory weights directly from frozen round5 data. It imports no production model or verifier, uses only the Python standard library, calls no optimizer, and does not rerun the old character convolution or third-moment backend.

The reviewed production snapshot is fixed:

- `trace_envelopes.py`: `0315f4367914d1cad6ceef2e08647cc706357762df0c00c5b68175aab41ee6ed`
- `results.json`: `10f8f77d25261422c9c0077a9cff4fe23843a8d08fe2b7688677ad3b3534689d`

[audit_certificates.py](audit_certificates.py) checks these hashes, all production input hashes, and each folded inventory hash. [review_results.json](review_results.json) records the exact output values and the review source hash. Earlier rounds remain unchanged.

## Model and normalization

The support contains only integers theta with theta²≤4q, the correct family residue modulo8, and nonnegative integral bulk sign-type counts. The reviewer independently reconstructs the eight bulk counts by parity inversion from the three boundary columns, row sums, pair inner products, and canonical triple sum. It confirms equality with every exported support point, including exclusions caused by type nonnegativity.

Let N_m(theta) be the monochromatic count and N_u(theta) the combined count for all three mixed normalized edge classes. The production variables are

```
w_m(theta)=N_m(theta)/q,
w_u(theta)=N_u(theta)/(3q).
```

Thus the mixed variable represents one third of the combined mixed count. The exact count identities are

```
sum N_m=(q−5)/4,       sum N_u=3(q−1)/4,
sum theta*N_m − sum theta*N_u/3 = 2,
sum theta²*(N_m+N_u)=(q−3)*(q+1).
```

The additional model fixes the actual ordinary sums of theta powers1,3,4,5,6; powers0 and2 are already fixed. The reviewer checks the scaled rows and right sides after setting R=floor(sqrt(q))+1. It then reconstructs the objective from the two frozen degree-six family polynomials:

```
constant repeated-row contribution
+ q²(q−1)*sum_theta(P_m(theta)*w_m(theta)+3*P_u(theta)*w_u(theta)).
```

The q(q−1) ordered-triple factor is applied once. The actual folded prime inventory is nonnegative, satisfies every equality exactly, and reproduces the saved exact global third moment in each case.

These variables are nonnegative **real/rational histogram weights**. Integrality, curve realizability, orbit stabilizer multiplicities, and existence of an underlying graph are not asserted. The support and constraints provide a relaxation containing the actual input, not an inverse construction theorem.

## Certificate audit

For each lower or upper certificate the reviewer checks every rational dual inequality on the entire finite support, the exact dual bound, every primal weight and equality, primal objective, primal/dual gap, minimum slack, and all point indices. The signed lower-bound convention is checked explicitly.

It also reverses the two exported class-constant repairs and independently recomputes their maximum deficits. Each support point belongs to exactly one constant row, so these nonnegative repairs suffice to make the proposed dual feasible. The checked inequalities establish the bounds irrespective of the optimizer's numerical accuracy or messages.

The exported known-constraint component is subtracted with exact arithmetic. The reviewer recomputes the maximum absolute residual and verifies that adding its fixed constraint contribution back leaves the exact objective unchanged. Floating-point HiGHS proposals are used only to discover candidate supports and dual multipliers. The reconstructed rational primal equations and repaired dual inequalities are the certificates.

All5,208 pointwise dual inequalities and104 primal equalities pass. Fifty-six corrupted certificates are rejected: negative dual directions, negative primal weights, positive weight perturbations, duplicate indices, incorrect bounds, incorrect minimum slacks, and incorrect quotient magnitudes. These are checks of the independent verifier, not evidence that every malformed input has been enumerated.

Means, second moments, variance, central third moments, and standardized third-moment squares are checked from exact fractions. For an interval of raw width W, the reviewer verifies

```
standardized_width_squared = W² / variance³
```

directly. It never obtains the width by subtracting rounded standardized endpoints; those endpoints can coincide numerically while their exact difference is nonzero.

## Confirmed ambiguity and a stronger small-case result

For the ordinary-moment model, both rational primal histograms have exactly the same theta power sums through degree6 as each other and the actual inventory. Their exact third moments differ at q1297,65537,1,000,033:

| q | Difference between two feasible raw third moments | Certified outer width | Standardized outer width |
|---:|---:|---:|---:|
|101|0|7.52386e−14|6.60396e−19|
|1297|0.0165967497575|0.0165967497583|2.44875e−9|
|65537|0.00986133938049|0.00986133938049|8.21440e−16|
|1,000,033|0.00981083030144|0.00981083376020|1.55318e−20|

These decimals summarize exact fractions stored in the review output. The certified bounds need not be exactly optimal. At the three larger orders, the unequal **primal** values establish ambiguity inside the relaxation; a positive dual interval alone would not establish it.

The q1297 witness is compact and fully explicit in [rational_ambiguity_witness_1297.json](rational_ambiguity_witness_1297.json). Its lower histogram has eight nonzero entries and its upper histogram has nine. Entries are N_m or N_u/3, with this convention included in the file. Both have ordinary powers0…6

```
1295, −2590, 1679612, −3814648,
4393553648, −14040602080, 14383212251072,
```

but their raw third moments differ by exactly

```
579595511175839744 / 34922229933188280115.
```

These are fractional trace histograms, not two constructed graphs.

At q101 the independent audit proves more than coincident optimizer values. The nine equality rows have exact rank9 on ten support points. Consequently every solution is the actual histogram plus a scalar multiple of one exported null vector. The actual monochromatic weights at theta=−14 and theta=2 are both zero, while the corresponding null-vector entries are−3 and42/5. Nonnegativity forces the scalar to be both nonpositive and nonnegative. Thus the feasible histogram is unique. [q101_uniqueness.json](q101_uniqueness.json) records the full null direction, support, actual weights, and the exact feasible parameter interval[0,0]. The tiny positive outer width at this order comes from the repaired certificates and is not evidence of ambiguity.

The result isolates a specific information gap: ordinary moments through degree6 can leave the two residue-class polynomial contributions undetermined in the stated relaxation. Their additional ambiguity becomes extremely small after standardization in these scaled examples. It is not a uniform theorem in q, a counterexample pair of actual graphs, a tail estimate for column subsets, or progress toward a prize proof by itself.

Replay only the saved-certificate audit with

```
python3 tooling_lab/round6/trace_envelopes_review/audit_certificates.py
```

The arithmetic and model here are independently coded. They still rely on the already reviewed frozen round5 polynomial and actual-inventory data; that dependency is explicit rather than counted as a fresh full backend verification.
