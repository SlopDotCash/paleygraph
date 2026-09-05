# Symmetry-folded trace inventory

This lane compresses the existing exact third-moment compiler from four normalized edge classes to two polynomial families, then to seven measured statistics and two quartic remainders. Both iterations match every one of the 14 saved small, twin, and scaled results, including prime order 1,000,033 and subset size31. The generic coefficient budget falls from35 to19, then15. It uses existing row permutation, global complementation, Walsh inversion, Legendre-family, and exact polynomial machinery; historical originality is unestablished.

The useful refinement is that the complete scalar histogram of **theta**, defined below, already determines the triangle label through a congruence modulo8. Two polynomial families remain necessary for the tested bounded-degree interface. This is different from saying that the complete scalar histogram loses information.

## The invariant and its precise symmetry

Let S be a symmetric zero-diagonal sign matrix satisfying S1=0 and S²=qI−J, and let a,b,c be distinct rows. Write

```
edges = (e01,e02,e12) = (S_ab,S_ac,S_bc),
tau = sum_x S_ax S_bx S_cx,
theta = e01 e02 e12 tau.
```

Walsh inversion reconstructs the entire histogram of the three row entries from these four numbers; the exact formula and zero-containing columns are in the [preflight](../preflight/README.md). A row permutation permutes the three edges and preserves tau. Global sign reversal negates all three edges and tau, preserving theta. The eight sign triangles have two orbits under these operations:

| Family | Canonical edges | Canonical tau |
|---|---|---|
| Monochromatic, all three edges equal | (1,1,1) | theta |
| Mixed | (1,1,−1) | −theta |

These are literal equivalences of column-type histograms up to the stated operations. If E_d(row,C) is the row's degree-d elementary symmetric coefficient on C, the triple product of the three E_d values is invariant under row permutation. Complementing all entries multiplies the product by (−1)^(3d), which is1 for even d. Thus both the union coefficients and their finite-population inclusion averages depend only on the family and theta for even d. The adapter rejects odd degree; no odd-degree quotient is implied.

For the canonical monochromatic triple, the count of the (+++) bulk pattern is (q−15+tau)/8, forcing theta≡15−q modulo8. For canonical mixed edges it is (q−7+tau)/8, forcing theta≡q−7 modulo8. These residues differ by4 when q≡1 modulo4. Therefore the **entire theta histogram alone identifies the family of each bin**. The separate Boolean is redundant audit metadata. This does not make the two degree-six coefficient polynomials equal.

The global normalized-pair contract remains essential. For each actual normalized parameter t there are q(q−1) ordered distinct row triples. The folded record counts have not yet been multiplied by that factor; the moment consumer applies it exactly once. There is no extra factor6. Arbitrary conference matrices without a verified normalized-pair symmetry require a separate inventory adapter.

Under that transitivity contract, each mixed equivalence class splits equally among the three edge patterns with e01=+1. Its folded count is divisible by3. `unfold_inventory` reconstructs those three entries, with tau=(edge product)*theta and one third of the folded count; the monochromatic entry is unchanged. This recovers all nine tested original inventories byte-for-byte at the record-value level. Folding is consequently lossless for these normalized joint histograms, although taking finitely many power sums is a further compression.

## Prime parameter action and exceptional orbits

For a prime p≡1 modulo4 and t≠0,1,

```
tau(t) = sum_x chi(x(x−1)(x−t)),
theta(t) = chi(t(t−1))*tau(t).
```

The generators R(t)=1−t and I(t)=1/t satisfy

```
tau(1−t)=tau(t),
tau(1/t)=chi(t)*tau(t).
```

For reflection substitute x=1−y; the cubic changes by a factor−1, whose character is1. For inversion substitute x=y/t; the cubic changes by t^(−3), whose character is chi(t). Reflection preserves chi(t(t−1)); inversion changes its value to chi(t−1). Multiplying by the corresponding tau identities proves theta invariance. The triangle family is preserved as well: reflection exchanges the two variable edge signs, while inversion sends the sign pair (u,v) to (u,uv). In particular the pair(+,+) is preserved.

The generated S3 orbit consists of the distinct values among

```
t, 1−t, 1/t, 1/(1−t), (t−1)/t, t/(t−1).
```

For p>3 its special orbits are:

- The harmonic orbit {−1,2,1/2}, of size3 and stabilizer size2.
- The two roots of t²−t+1, when they exist, forming one size2 orbit with stabilizer size3. They exist exactly when p≡1 modulo3.
- Every other orbit has size6 and trivial stabilizer.

One obtains these exceptions by solving the fixed-point equations of the three transpositions and two3-cycles. The harmonic and equianharmonic solutions intersect only in characteristic3, excluded here. Hence the weighted orbit census is exact; no assumption that every orbit has six elements is made.

At p65537 the census is one size3 orbit and10,922 size6 orbits. At p1,000,033 it is one size2 orbit, one size3 orbit, and166,671 size6 orbits. The source verifies both generator identities, theta and family preservation, every orbit, all exceptional root equations, and total multiplicity p−2 against the complete saved arrays at all seven primes. This is a symmetry check on previously certified arrays, not a second independent full convolution at the largest primes.

The prime action above is not asserted for the Peisert character. The generic row-permutation/complement proof still applies to the actual Peisert49 matrix and its separately verified semilinear normalized-pair adapter. Its tau is a triple correlation, not a Legendre elliptic trace.

## Input and coefficient budget

[folded_moment.py](folded_moment.py) provides `fold_inventory`, `unfold_inventory`, and `folded_global_third_moment`. The last function accepts the original joint record or a folded record. Its arithmetic consumer needs only q,n,d, the normalization contract, and the seven power sums of theta for each of the two families. This is14 integer statistics, two of which (the zeroth powers) are fixed by q. It does not enumerate or consume individual parameter records; the replay deliberately removes them before evaluation. Fourteen is a proved sufficient budget, not a claim of minimality.

For fixed canonical edges, the degree-d triple-product coefficient is polynomial in tau of degree at most d. The [preflight logarithmic argument](../preflight/README.md#why-only-seven-tau-evaluations-per-sign-class-suffice) proves this independently of interpolation: a tau-dependent logarithm term consumes at least one degree in each of three variables. Replacing tau by ±theta preserves the bound.

At d6, seven exact admissible nodes determine each family polynomial; an eighth is an implementation check. Add the three repeated-row templates: 2×8+3=19 coefficient evaluations. A minimal interpolation computation would use17, omitting the two extra checks. Repeated-row complement symmetry could remove another redundant check, but this implementation retains it. If fewer than seven admissible values exist at small q, the adapter evaluates all possible values and uses a polynomial agreeing on that finite domain. Such a polynomial is not promoted to a degree-six identity outside the admissible domain. The present API explicitly requires q≥13 and d∈{0,2,4,6}.

The work before folding is the frozen trace backend's O(p log p) exact modular arithmetic and O(p) memory. Folding costs O(h), with h the number of input histogram bins. The consumer uses at most19 bounded-degree coefficient problems for d6; arithmetic cost also depends on integer bit lengths, which grow with q. Its union degree is at most18 and it allocates no p² matrix. The independent full-array orbit verifier is linear after character and inverse tables (its Euler-criterion character oracle deliberately costs O(p log p)). Literal tiny-matrix verification costs O(q^4) and is not used for scaling.

### Proved revision: seven statistics and fifteen evaluations

[residue_moment.py](residue_moment.py) strengthens the first adapter by identifying its two highest polynomial coefficients before interpolation. Individual row sign changes preserve the degree-six row-product coefficient, even though such changes need not preserve the whole conference matrix. For each canonical family choose these signs so that the three-row product sum is theta. All three rows remain balanced. If F is the column product used by the compiler, its formal logarithm has the form

```
log F = G + theta*L.
```

G is independent of theta. Every monomial of L has positive odd degree in all three target variables. Its lowest term is

```
t1*t2*t3*K(z),       K(z)=z−3z²+2z³.
```

Indeed the coefficient of t1*t2*t3 in H,H²,H³ is respectively1,6,6 times the column product, and higher powers cannot contribute. The logarithm coefficients therefore give K. To extract theta^6 at target degrees(6,6,6), all six factors of L must use their lowest term, and exp(G) must contribute its constant term. Consequently the **raw union polynomial** coefficient of theta^6 is exactly

```
K(z)^6 / 720.
```

For theta^5, the five factors of L must again use their lowest terms: the next positive odd exponent would exceed6. The remaining degree(1,1,1) would have to come from exp(G). The singleton coefficients of G vanish because each row sums to zero. Its (1,1,1) coefficient is zero because that is precisely the theta-dependent term already removed. Every decomposition of(1,1,1) into two or three nonzero exponent vectors uses a singleton. Hence [t1*t2*t3]exp(G)=0, proving that the theta^5 coefficient vanishes. This proof holds before the inclusion functional and for every valid q; it is not an inference from the quartics observed at n7 and n8.

After applying inclusion, write

```
P_family(theta) = c6(q,n)*theta^6 + R_family(theta),
degree(R_family) <= 4,
c6(q,n) = inclusion_(q,n)(K(z)^6/720).
```

Five nodes determine each remainder, a sixth checks the implementation, and three repeated-row templates remain:2×6+3=15. Small-q domains use all available admissible values as before. The [raw coefficient audit](residue_results.json) checks both highest-degree identities at all19 union degrees in both families for q49,101,1297, before choosing n. All14 end-to-end moment values again agree with both earlier compilers.

Let A_j denote monochromatic theta power sums, B_j mixed sums, and M_j=A_j+B_j. The conference identities imply

```
A0=(q−5)/4,  B0=3(q−1)/4,
B1=3*A1−6,
B2=(q−3)*(q+1)−A2.
```

For the first-power relation, on normalized rows set u=S_0t,v=S_1t. The monochromatic indicator is(1+u)(1+v)/4 and theta=uv*tau. Thus4A1=M1+sum(tau)+sum(u*tau)+sum(v*tau)=M1+6. The three last sums equal2 by row balance and S²=qI−J, with the two repeated parameters removed. The second-power identity is the corresponding squared norm identity and theta²=tau².

It follows that the following **seven measured integers** suffice, together with q,n and the normalization contract:

```
A1, A2, A3, B3, A4, B4, M6.
```

`seven_statistics(record)` exports them and `seven_stat_global_third_moment(summary,n)` consumes only that summary. The summary is an exact sufficient-statistic contract from a separately certified inventory; it is not a graph-realizability certificate. The full-record path checks admissible support through the original compiler; a supplied powers-only summary is trusted provenance, with inexpensive identities checked but no inverse reconstruction of an actual graph. Neither adapter silently claims such an inverse theorem.

Equivalently, with kappa=(15−q) modulo8, define the real-valued mod8 character

```
s(theta)=exp(2*pi*i*(theta−kappa)/8),
```

which is+1 on the monochromatic residue and−1 on the mixed residue. A_j=(M_j+sum(s(theta)*theta^j))/2. The seven-stat interface can therefore be expressed as ordinary moments M1,M3,M4,M6 and the three character-weighted moments of orders2,3,4, with the lower quantities eliminated by the displayed identities. This is a residue-weighted trace-moment calculation based on a known finite character. It does not invent character twists or Legendre trace moments.

The raw-coefficient audit additionally computes the exact rank of the seven statistic coefficient columns after the stated conference identities: it is7 at each of q49,101,1297. Thus none of the seven formal statistic directions can be discarded while retaining all union coefficients at those q. Since the inclusion polynomials are independent on the valid subset sizes, the same statement holds for simultaneous queries over all n there. These ranks do not prove independence over actual graph families or a universally minimal data representation.

## Results and rejected shortcut

[results.json](results.json) records exact equality with all14 existing general-compiler cases, including the three subset sizes6,7,8 for both q49 twins. Degree0,2,4 checks also match at q13; degree3 is rejected. At p1,000,033,n31 the exact third moment is unchanged, approximately3.3859584991473612e16. The folded backend used19 evaluations and took about1.83 seconds in the first measured run. Timings vary with concurrent load; the operation-count reduction is the dependable comparison.

| Prime order | Original joint bins | Folded theta bins |
|---:|---:|---:|
|13|7|3|
|17|7|3|
|29|10|4|
|101|17|7|
|1297|63|27|
|65537|448|192|
|1,000,033|1750|750|

The first proposed further shortcut was one degree-six polynomial across both families. At q101,n6 it works: both reduce to the same even Krawtchouk polynomial. At n7 it fails exactly. The monochromatic polynomial minus the mixed polynomial is

```
−4/11315535 −(542/3224927475)*theta
+(59/6449854950)*theta² +(2/3224927475)*theta³
−theta⁴/51598839600.
```

At n8 the difference is also a nonzero quartic, exported exactly. Both canonical families have at least seven admissible values at q101, so their degree-six polynomial extensions are uniquely determined. Their non-equality rules out a common degree-six polynomial on the full admissible synthetic-type domain. The actual theta residues are disjoint: this is not an example of two actual triples with the same theta but different labels. We have not produced two actual graphs with identical first seven unlabelled theta moments but different global third moments. No such stronger insufficiency claim is made.

[symmetry_validation.json](symmetry_validation.json) also reports literal canonical type checks for every unordered distinct triple in Paley13,17,29,49 and Peisert49. All canonicalizing row permutations and complements are tested, not just one representative; the normalized q(q−1) weights agree with the full ordered-triple inventory.

The earlier finite controls still apply unchanged. At n6 the q49 twins have different third moments but the same maximum27. At n7 the Peisert third moment is larger (−1672321/1533939 versus −99128449/1533939), while its maximum is smaller (61 versus77). At n8 its third moment and maximum are both larger (maximum132 versus116). These [independent full-distribution controls](../review/README.md) refute a monotone ordering of maxima by the raw third moment in these examples. They do not prove that every use of the full invariant inventory must fail. No exceptional-set bound, signed spectral-edge bound, or prize proof follows from this compression.

## Attribution and replay

The [Stacks Project's Legendre-family section](https://stacks.math.columbia.edu/tag/03VA) gives the standard family y²=x(x−1)(x−lambda). [Brouwer–Martin](https://arxiv.org/abs/2109.03654) studies Paley triple intersections; both sources were freshly checked2026-09-05. The affine parameter permutations, quadratic twists, row-sign symmetries, and Walsh inversion are established foundations. The local [preflight](../preflight/README.md#prior-art-and-local-overlap) documents earlier workspace trace tooling and the broader trace-moment literature. The specific local addition here is the exact symmetry adapter, its finite input budget, residue-aware audit, exceptional-orbit checks, and replay against the existing global-moment compiler. No exhaustive historical search or theorem of originality is claimed.

```
python3 tooling_lab/round5/trace_invariants/run_experiments.py
python3 tooling_lab/round5/trace_invariants/run_residue_experiments.py
python3 tooling_lab/round5/trace_invariants/verify_symmetry.py
python3 tooling_lab/round5/trace_invariants/write_manifest.py
```

The source imports the existing exact coefficient compiler; no source in rounds1–4, preflight, or trace_backend is modified. `manifest.json` binds this lane's final files. External source and input hashes are included in result records. The scope of the independent reviewer is recorded separately under `round5/review/fold_review*`.
