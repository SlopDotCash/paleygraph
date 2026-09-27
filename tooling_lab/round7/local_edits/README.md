# Exact local variance at an actual arithmetic set

The prototype computes the first and second moments of `T_d(C')` over every one-point swap of a fixed set C, without evaluating those neighbouring sets. The new retained object is an n-by-n internal contraction of the deletion derivatives. A direct toy census shows that the complete signed row-count histogram and the swap mean do not determine this local variance.

Both implementations and the scale comparison are completed. Every retained contraction agrees across eight critical-size inputs through q=6,700,417,n=50. A further ablation at n=7 and8 shows why the n=6 witness alone would have overstated the value of the new representation.

Let S be a real symmetric matrix with zero diagonal, signs ±1 off the diagonal, `S1=0` and `S²=qI-J`. For C of size n and `a in C`, define

```
r_a(x) = e_(d-1)(S[x,C without a])
v_a(b) = sum_x S[x,b] r_a(x)
u_a = v_a(a)
V[a,b] = v_a(b), for a,b in C.
```

Elementary symmetric polynomial insertion gives the exact identity

```
T_d(C without a plus b) - T_d(C) = v_a(b) - u_a,  b outside C.
```

The balanced column sum and conference identity imply

```
sum_b v_a(b) = 0
sum_b v_a(b)^2 = q sum_x r_a(x)^2 - (sum_x r_a(x))^2.
```

Thus we can remove the n forbidden insertion points using just V. Write `m=q-n`, `t_a=sum_(b in C) V[a,b]`, `w_a=sum_(b in C) V[a,b]^2`, and `g_a=q||r_a||²-(sum r_a)²`. If delta is the change under a uniform swap, then

```
n*m E[delta]   = sum_a (-t_a - m*u_a)
n*m E[delta²]  = sum_a (g_a - w_a + 2*u_a*t_a + m*u_a²).
```

Subtracting the square of the mean gives the exact local variance. Degree0 uses r=0. No independence of character values, random-seed assumption or neglected boundary terms enters these identities. They hold at each admissible actual C.

The full signed/zero row-count histogram determines every sum `sum_x f(sum_(a in C) S[x,a], 1_(x in C))`, all elementary targets T_d, and all their one-swap means: for a fixed row, the inside and outside sign populations follow from those two counts and row balance. It does not keep the correlations between different rows needed for squared insertion sums.

The strongest saved Paley17 witness is:

| Quantity | C={0,1,2,3,4,9} | D={0,1,2,3,4,8} |
|---|---:|---:|
| Complete signed/zero row-count histogram | Identical | Identical |
| T6 | -1 | -1 |
| Mean T6 after one swap | 1/33 | 1/33 |
| Variance after one swap | 9140/1089 | 12308/1089 |
| Minimum / maximum neighbouring T6 | -5 / 3 | -5 / 11 |
| Complete insertion-transform squared energy, sum g_a | 1144 | 1144 |
| Internal contraction squared energy, sum w_a | 590 | 398 |

Here every u_a=-1, the sum of the t_a is -2 in both cases, all six derivative norms are12, and the multisets of derivative sums also agree. The entire variance difference is therefore `(590-398)/(6*11)=32/11`. It comes from correlations in the internal contraction, not the available complete-transform energy. All66 neighbours of each set were evaluated directly. This toy is outside the critical-size constraint n^4<=q; it is an exact information-loss witness, not a critical-scale counterexample.

The normalized census covers all1365 six-sets containing0 and1, checks90090 neighbour values against a cache of all12376 six-set targets, and finds12 ambiguous variance fibres among70 complete signed-histogram fibres. These are normalized-set fibres, not a count of affine equivalence classes. A separate degree/boundary suite compares120 cases and6432 direct neighbour evaluations across Paley5,13,17,29, including degree0,1,3,6, n=1 and n=q-1.

The next ablation retained the **unordered list of every quartic correlation** as well as the full histogram. At n=6 it resolved every observed variance ambiguity, revealing the shortcut: for b distinct from a, V[a,b] reduces to the quartic correlation on C without a,b, minus its value at the excluded zero row. At n=7 and8 the quartic deck no longer suffices.

| Size in Paley17 | Normalized inputs | Histogram ambiguity fibres | Histogram + quartic-deck ambiguity fibres | Histogram + deletion-summary ambiguity fibres |
|---:|---:|---:|---:|---:|
| 6 | 1,365 | 12 | 0 | 0 |
| 7 | 3,003 | 30 | 6 | 0 |
| 8 | 5,005 | 47 | 6 | 0 |

The seven-point sets `{0,1,2,3,5,6,7}` and `{0,1,2,4,5,7,9}` have the same complete row histogram, the same unordered quartic deck, T6=-3 and mean swap value -17/35. Their local variances are57216/1225 and75136/1225, an exact difference512/35. The eight-point witness has variance difference64/9. [Ablation results](ablation_results.json) contain both full witnesses. The census compares660,660 direct neighbour values over9,373 normalized sets; it is outside the critical-size regime.

The deletion summary here is the unordered list of `(sum r_a, ||r_a||², u_a, t_a)` alongside the full histogram. It omits each row's squared internal energy w_a. Its zero ambiguity on this particular finite domain is an observation, not a proof that these summaries always determine the variance. A further orbit audit below tests whether this apparent sufficiency merely reconstructs the full input up to symmetry.

The Python backend groups each row by its numbers of positive, negative and zero signs. A short binomial table evaluates every r_a. It streams row blocks, calculates `R^T S_C`, and accumulates derivative sums and squared norms. Its work is O(q*n²) plus fixed-degree row work; storage is O(q+block*n+n²), including the prime character vector. This prototype restricts q<=10,000,000, n<=64 and d<=6. NumPy integer operations have explicit absolute-sum bounds; norms are accumulated as Python integers across blocks and final moments are rational.

The separate [Gram backend](gram_oracle.cpp) computes full elementary coefficients by multiplication, divides synthetically by `1+s*t`, and writes the deletion field as `A(x)+B(x)S[x,a]` with an explicit correction at x=a. It then uses the symmetric weighted Gram matrix `S_C^T diag(B) S_C`. This is a separately implemented arithmetic check of every retained contraction, not an independent mathematical theorem or a direct large-neighbour census.

The relation with established tools is explicit. With the uniform swap operator P, `2*Gamma(T,T)=E[delta²|C]` is the usual discrete carré du champ. Local variance is this quantity minus `(PT-T)²`. No new notion of conditional variance, exchangeable pairs or discrete derivative is claimed. The [prior-art audit](prior_art.md) also identifies the previous workspace drift and variance calculations.

An exact local variance still does not bound every exceptional set. The prior thin-shell persistence obstruction remains relevant. The immediate research use is to evaluate and falsify proposed local estimates at actual inputs, while exposing exactly which correlations a coarse feature map discards. A uniform estimate for these retained contractions, or a more effective use of them, remains missing.


## Scale, exact finite bounds and use

The independently compared critical-size cases use arithmetic progressions and one seeded uniform set at each q. They preserve every input label and satisfy n^4<=q.

| q | n | Neighbours covered by each exact moment query | First backend, AP / random (seconds) | Gram backend, AP / random (seconds) |
|---:|---:|---:|---:|---:|
| 1,297 | 6 | 7,746 | <0.002 / <0.002 | <0.01 / <0.01 |
| 65,537 | 16 | 1,048,336 | 1.85 / 1.02 | 0.019 / 0.019 |
| 1,000,033 | 31 | 31,000,062 | 50.8 / 58.9 | 2.47 / 2.43 |
| 6,700,417 | 50 | 335,018,350 | 224 / 243 | 43.1 / 36.8 |

The first backend reports contraction time after character construction; the Gram time includes process and field setup. These are observed timings on a shared, busy machine, not controlled performance benchmarks. The Gram method's contraction work remains O(q*n²). It is now exposed through [query.py](query.py), which returns exact rational moments and optionally the complete contraction record.

```
c++ -O3 -std=c++17 tooling_lab/round7/local_edits/gram_oracle.cpp -o tooling_lab/round7/local_edits/gram_oracle
python3 tooling_lab/round7/local_edits/query.py --q 65537 --selected 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 --upper-threshold 50000
```

An optional threshold uses classical Cantelli with exact rational arithmetic. If the neighbourhood mean is mu, variance v and threshold t>mu, the proportion with value at least t is at most `v/(v+(t-mu)^2)`. Multiplying by the exact neighbour count and taking its integer floor yields a valid finite count bound. The lower-tail version applies the same calculation to the negative target. For the wrong side of the mean, the API returns the trivial bound.

For the progression `{0,...,49}` at q=6,700,417, T6=-98,772,836. The exact moments certify that at least330,335,111 of its335,018,350 neighbours retain the negative sign and more than half this magnitude: over98.602%. Those neighbours were not individually evaluated. The uniformly seeded input at the same q has only a29.082% lower certificate for retaining its own sign and half-magnitude. These are input-specific one-step lower bounds, not exact frequencies, multi-step persistence, Sidon-restricted neighbourhoods or an exception-removal theorem. They do not overcome the prior thin-shell obstruction.

The known bound is included to make the local moments operational, not as a claim to invent a concentration inequality. [Threshold validation](threshold_results.json) checks8,298 exact finite-distribution event bounds, rejects11 malformed queries, checks all four supported test degrees at the maximum n=64, and compares the fast consumer on all eight saved scale cases.

[Independent readback](verification.json) verifies7,506 internal entries and206 derivative sums/norms against the separate Gram implementation, literally rebuilds ten witness sets and700 neighbour targets, checks both equal-quartic-deck pairs, and rejects seven corrupted records. Its Python code imports no production module. The producer's exhaustive direct oracles remain separate from this narrower saved-witness review. All work was completed locally after the parallel agents became unavailable; these are separate implementations and root-agent checks, not human or Lean formalization.


## A saturation control for apparent sufficiency

The [orbit audit](orbit_results.json) and [independent Burnside/feature review](orbit_verification.json) explain much of the deletion-summary result. At n=6 its98 summary fibres each identify one of the98 square-affine orbits. At n=8 the same is true for190 fibres and190 orbits. A feature that identifies the entire set up to symmetry automatically determines every invariant on that finite domain. Those two zero-ambiguity results therefore do not supply evidence of a useful general moment closure.

At n=7 there are150 square-affine orbits but148 summary fibres. Two fibres each join two genuinely different full-affine orbits, while retaining the same local variance. Every orbit count is independently confirmed by Burnside's lemma, and the reviewer rebuilds derivative profiles at every canonical representative using literal products. The finite normalization coverage is also checked on every six-, seven- and eight-set of F17. Direct enumeration shows that both residual pairs also have identical complete one-swap target histograms: values -19,-11,-3,5 with multiplicities2,10,40,18. Consequently no higher unlabelled moment of this one-swap distribution can separate either pair. The next distinction must retain how edits share their deleted or inserted points, or inspect a different neighbourhood.
