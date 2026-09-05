# Exact critical-size third-moment laws and inventory-free bounds

The exact coefficient calculation supports a different asymptotic from the apparent finite-size trend. For the prime Paley statistic

```
T6(C)=sum_x e6((S_xy)_(y in C)),     |C|=n,
```

with q≡1 modulo4 prime and n=alpha*q^(1/4)+O(1), fixed alpha>0, its standardized third central moment satisfies

```
skewness(T6) ~ 40*sqrt(5)/sqrt(q).
```

It does not acquire a further factor n. The exact finite-n falling-factorial correction explains the previously observed values. This lane also bounds how much trace arithmetic can change this particular global average and supplies a finite enclosure at q=6,700,417 and q=2,013,265,921 without constructing a character array. These are global-moment statements, not bounds on exceptional column subsets, spectral extrema, or either prize. The classical foundations and historical novelty limits are stated below.

## Exact algebra and proof budget

The frozen [round5 symmetry reduction](../../round5/trace_invariants/README.md) expresses the global third moment using seven measured statistics

```
A1,A2,A3,B3,A4,B4,M6,
```

where A_j and B_j are theta power sums over monochromatic and mixed normalized parameters and M_j=A_j+B_j. The present compiler exports the identity

```
E[T6^3] = sum_(k=0)^18 (n)_k/(q)_k *
          (B_k(q)+sum_(i=1)^7 C_(k,i)(q)*statistic_i).
```

The symbol B_k(q) here is a baseline polynomial, distinct from the mixed moment B_j. Every polynomial coefficient is rational and is saved in [q_polynomials.json](q_polynomials.json). The baseline includes all three repeated-row templates and the conference constraints; their normalization is inherited unchanged from round5. The global second moment is represented by thirteen further polynomials, for union sizes0…12.

The finite q evaluation budget has a structural proof. In the formal logarithm of the column product, each multiplicity is affine in q and theta. For distinct rows, a q-dependent sign-power sum has even exponent on every coordinate that occurs, so each q-dependent logarithm factor consumes at least two target degrees. The target is(6,6,6), with total degree18. Thus every raw coefficient has q degree at most9. A theta-dependent factor consumes at least one degree on all three variables; the coefficient of theta^j has q degree at most floor((18−3j)/2). At fixed union size k, the elementary coefficient is also polynomial of degree at most k in the multiplicities, giving the additional bound k. The repeated-row templates have balanced signs as well; their q-dependent logarithm terms again consume at least two degrees. For the two-row target(6,6), the corresponding raw q degree is at most6, and the global row weights increase the bound to8.

These arguments concern formal polynomial identities. They do not require every sampled q,theta pair to come from a graph. The affine histogram formula defines the same polynomial over the rationals on both q residue branches. Consequently interpolation on one admissible q congruence class determines the polynomial; it is not a separate empirical fit for each residue branch.

[compile_q_polynomials.py](compile_q_polynomials.py) uses q=81,89,…,161. At each q it evaluates six admissible theta nodes in each of the two families, plus three repeated-row templates. The already proved universal theta^6 coefficient and vanishing theta^5 coefficient leave two quartic remainders. Five nodes determine each remainder; the sixth checks it. Ten q nodes suffice for the largest degree9 coefficients; the eleventh checks the result. Smaller degree bounds use fewer determining nodes and check all remaining ones.

This costs165 exact third-moment coefficient evaluations. A further75 evaluations verify every raw polynomial at q49,101,1297,65537,1,000,033, including the other q residue branch and values far outside the determining interval. The second-moment compiler is queried for three pair templates at eleven q values. Compilation took about23 seconds in this run; query costs are much smaller. Timing is incidental to the proved degree budget. The source and all imported backend hashes are stored with the polynomial output.

## Leading terms and excluded competitors

If a polynomial term has q degree d and union size k, the inclusion factor contributes q^(−3k/4) in the critical regime. A statistic of order j has Hasse size at most O(q^(1+j/2)). [scaling_results.json](scaling_results.json) records the degree/exponent scan for every nonzero coefficient polynomial and all seven arithmetic directions; the polynomial arrays retain the zero coefficients explicitly.

| Component | Unique leading union size | Leading q coefficient | Critical exponent | Largest competing exponent |
|---|---:|---:|---:|---:|
| Baseline third moment |9|q^10/216|13/4|5/2|
| Second moment |6|q^7/720|5/2|7/4|
| Arithmetic part of third moment |6|depends on M4,M6|at most3/2|at most1 after its leading union6 part|

The leading first row is also transparent combinatorially. For one balanced row, the leading contribution to the cube of e6 uses three labelled six-sets in which each selected coordinate occurs twice. Their three pairwise-only blocks each have size3. Choosing the three labelled disjoint blocks gives (n)_9/(3!^3)=(n)_9/216. There are q choices of the row. In the leading variance term, the two six-sets are identical, giving q*binom(n,6). This argument identifies the constants; the exact polynomial scan is what excludes contributions from other row partitions and trace terms. No independence of the actual character rows is assumed.

The mean simplifies exactly:

```
mu = -(n)_6/[48*(q−2)*(q−4)].
```

The resulting laws, uniformly for alpha in a fixed compact subinterval of(0,infinity), are

```
variance = q*(n)_6/720 + O_alpha(q^(7/4)),
raw third moment = q*(n)_9/216 + O_alpha(q^(5/2)),
central third moment = q*(n)_9/216 + O_alpha(q^(5/2)).
```

The centering correction has size O_alpha(q²), so it lies below the stated third-moment error. The variance correction mu² is smaller still. Equivalently,

| Quantity | Leading critical law |
|---|---|
|Mean|−alpha^6*q^(−1/2)/48|
|Variance|alpha^6*q^(5/2)/720|
|Raw and central third moment|alpha^9*q^(13/4)/216|
|Standardized third moment|40*sqrt(5)*q^(−1/2)|

Retaining the finite-n correction gives

```
skewness = [40*sqrt(5)/sqrt(q)] *
           [(n)_9/(n)_6^(3/2)] * (1+O_alpha(q^(−3/4))).
```

The falling-factorial ratio tends to1, with expansion1−27/(2n)+O(n^(−2)). These are asymptotic statements along n growing in the critical regime; they do not predict the sign at every small n. For n<9 the displayed leading proxy is zero, and the lower terms can determine the entire result.

| q,n | Exact standardized third moment | Finite-n leading proxy | Limiting proxy40sqrt(5/q) |
|---|---:|---:|---:|
|65537,16|0.1047220069|0.1047625186|0.3493829559|
|1,000,033,31|0.0536038970|0.0536079092|0.0894412433|

Thus the finite values are consistent with a substantial falling-factorial correction; extrapolating their near proportionality to n would have given the wrong limiting law. All14 saved finite mean, variance, raw third, and central third values, including the q49 twins, match the new polynomial identity exactly. The Hasse-based asymptotic and bounds below apply to prime Paley, not to the Peisert control.

## Uniform arithmetic sensitivity

These statements concern two nonnegative Hasse-supported normalized histograms with the displayed conference constraints, at the same q,n. They include actual prime Paley inputs and the larger fractional relaxation. The mean and variance are fixed by q,n, so differences of raw and central third moments agree.

Without fixing the ordinary theta moments, the complete arithmetic-dependent part has size O_alpha(q^(3/2)). A more informative uniform difference bound follows from its leading union6 term. Let N=q−2, x=theta², and m2=(q−3)(q+1). The exact union6 arithmetic contribution is

```
(n)_6/(q)_6 * q(q−1)/720 *
  [M6−(15q−85)*M4].
```

With H=floor(2sqrt(q)), the polynomial f(x)=x³−(15q−85)x² is concave on[0,H²] for q≥29. Therefore every admissible nonnegative histogram satisfies the chord and Jensen bounds

```
m2*[H^4−(15q−85)H²]
 <= sum f(theta²)
 <= N*f(m2/N).
```

The difference between these endpoints is(30+o(1))*q^4. All remaining arithmetic terms have critical exponent at most1. It follows that

```
limsup q^(−3/2)*|Delta E[T6^3]| <= alpha^6/24,
limsup q^(9/4)*|Delta skewness| <= 360*sqrt(5)*alpha^(−3).
```

These bounds quantify worst sensitivity under the stated Hasse and conference constraints. They do not describe the actual variation of an arithmetic family, which may be smaller.

If ordinary theta powers0…6 are fixed as well, A1 is determined, and changes in B3 and B4 are the negatives of changes in A3 and A4. Only three directions remain. The exact degree scan gives maximal critical exponents−1/4,−1/2,−1/4 for A2,A3,A4 respectively. The leading union7 expression is

```
(n)_7/q^5 * Delta[2q*A2−A4/3].
```

The exact coefficient before taking leading terms is

```
(n)_7/(q)_7 * q(q−1) *
 [(2q−134/3)*DeltaA2 +(32/3)*DeltaA3 −DeltaA4/3].
```

For u=theta²/q in[0,4], the function2u−u²/3 lies in[0,3]. The monochromatic mass is(q−5)/4. Consequently

```
limsup q^(1/4)*|Delta E[T6^3]| <= 3*alpha^7/4,
limsup q^4*|Delta skewness| <= 6480*sqrt(5)*alpha^(−2).
```

The lower-order raw remainder here is O_alpha(q^(−1/2)). This is a bound conditional on agreeing ordinary trace moments; it cannot be substituted for the unconditional-in-prime-Paley bound without those measurements.

There is also a simple finite, exact, no-optimizer bound. Compute the three actual coefficients D2,D3,D4 of the remaining directions and evaluate D2*theta²+D3*theta³+D4*theta^4 on the finite monochromatic Hasse/residue support. Its range times(q−5)/4 bounds every permitted difference. This independent range bound contains all eight previously exhibited envelope witnesses. For fixed ordinary moments it gives0.0299335 at q1297,0.0116281 at q65537, and0.0106739 at q1,000,033, compared with feasible LP witness widths0.0165967,0.00986134,0.00981083. It is an inexpensive conservative control, not an optimality claim. At q101 it remains positive because it omits further feasibility constraints; the separate exact envelope review proved uniqueness there.

## Applications without a trace inventory

[inventory_free.py](inventory_free.py) uses only q,n and the polynomial file. It proves primality by trial division, computes the exact mean and variance, applies the preceding chord/Jensen bound to the union6 contribution, and bounds every remaining arithmetic term using its exact coefficient and the class-mass Hasse ranges. It computes no character entries and assumes no ordinary trace moments have been measured.

| Prime q | n | Certified standardized interval, rounded | Width computed from exact fractions |
|---:|---:|---|---:|
|6,700,417|50|[0.025669519976984242,0.025669519977420386]|4.36146425765e−13|
|2,013,265,921|211|[0.0018672388434094857,0.0018672388434094866]|9.83498625501e−19|

[inventory_free_results.json](inventory_free_results.json) preserves the exact rational raw and central interval endpoints, exact variance, and exact standardized width squared. Rounded endpoint subtraction is never used to obtain the widths. Each query took about0.002–0.003 seconds in this run, including trial division. The first order is the original small-prize-scale parameter, but this calculation certifies only a global average at that order; it gives no required maximum or exceptional-set bound.

Once the polynomial file is compiled, the rational arithmetic has a fixed coefficient budget and no q-sized allocation. Trial-division primality checking uses O(sqrt(q)) small divisions; the moment query itself uses a fixed number of polynomial evaluations whose integer bit lengths grow with log(q) and log(n). Falling factorials use at most18 factors. The `enclose` API requires positive variance: it rejects degenerate choices such as n=q or n=q−1, where skewness is undefined. The two critical applications satisfy this requirement. Six existing prime cases supply scalar readback controls for the inventory-free enclosure, without their arrays being read.

## Attribution, limitations, and replay

Finite-population inclusion, Newton identities, parity/diagram counting, polynomial interpolation with proved degree bounds, Hasse's bound, and Jensen's inequality are established mathematics. The local [parallel23 classical-upper note](../../../research/parallel23-classical-upper-2026-09-05.md) already uses balanced sampling and Rademacher pairing counts; that note's power statistic is different from the present cube of a degree-six elementary symmetric statistic. The [round4 moment compiler](../../round4/marked_moments/README.md), [round5 preflight](../../round5/preflight/README.md), and [trace-invariant reduction](../../round5/trace_invariants/README.md) supply the exact interfaces reused here.

[Peccati–Taqqu's survey on moments, cumulants, and diagram formulae](https://arxiv.org/abs/0811.1726), checked2026-09-05, documents the established combinatorial moment/chaos framework. The constants here have a direct finite combinatorial proof, rather than relying on a new Gaussian-limit theorem. No exhaustive literature search establishes historical originality of these laws or this computational reduction. The specific local output is an executable q-polynomial certificate, a failed-fit control, explicit sensitivity bounds, and a scalable global-average enclosure.

The saved identities are exact rational computations with a mathematical degree argument; they have not been formalized in Lean. An independent root review is recorded separately under `round6/critical_third_review/`. No probabilistic central limit theorem, distributional convergence, pointwise control, or prize proof follows merely from the first three moment asymptotics.

```
python3 tooling_lab/round6/critical_third/compile_q_polynomials.py
python3 tooling_lab/round6/critical_third/derive_scaling.py
python3 tooling_lab/round6/critical_third/inventory_free.py
```

The compiler queries frozen prior sources read only. Source/input hashes are recorded in each output, and `manifest.json` binds this completed lane.
