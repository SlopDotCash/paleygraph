# Covariance between edits that share an inserted point

This completed prototype retains a joint distribution that the previous tool flattened. It distinguishes actual Paley sets with the same complete one-swap target histogram, the same unlabelled collection of complete per-deletion distributions, and the same per-deletion variance multiset. The distinguishing quantity is the covariance between different deletion choices when the inserted point is held in common.

The computation extends through q=6,700,417,n=50 using a50-by-50 derivative Gram matrix. It does not evaluate all335,018,350 neighbouring targets. The extra matrix is calculated by streaming character rows; the API remains limited to prime q congruent to1 modulo4, q<=10,000,000, n<=64 and degree<=6. Its numerical experiments below use degree6. These are exact arithmetic tools at specified inputs, not a uniform bound for either prize.

## The finite separation

Consider C={0,1,2,3,4,6,10} and D={0,1,2,3,4,6,15} in F17. The [previous independent orbit audit](../../round7/local_edits/orbit_verification.json) establishes that they are not affine-equivalent and that their complete one-swap T6 distributions agree: values -19,-11,-3,5 with multiplicities2,10,40,18. Every higher moment of that flattened distribution must therefore agree as well.

The [preflight](results.json) additionally checks the unlabelled per-deletion distributions and simpler possible separators. The unlabelled per-insertion distributions already differ; covariance is not the only possible distinguishing statistic. Its benefit is an exact streaming contraction that avoids constructing the complete edit matrix.

Write Sigma for the covariance matrix of the deletion responses under a shared uniform outside insertion. Let P=I-11^T/n remove the common deletion mode.

| Quantity | C | D |
|---|---:|---:|
| Complete one-swap target distribution | Identical | Identical |
| Unlabelled collection of per-deletion distributions | Identical | Identical |
| Multiset of diagonal entries of Sigma | Identical | Identical |
| Sum of all covariance entries | 3584/25 | 3584/25 |
| Trace of Sigma | 4928/25 | 4928/25 |
| Trace of Sigma² | 10317824/625 | 9211904/625 |
| Trace of (P Sigma P)² | 68308992/4375 | 60321792/4375 |

A second non-affine pair has the same separation. Removing the response common to all deletions preserves the distinction. Translation and nonsquare-affine changes preserve each covariance matrix after the corresponding row permutation. No ordering convention is responsible for the separation.

This is a toy input outside n^4<=q. It establishes a precise information requirement on actual arithmetic inputs; it does not establish a critical-scale collision or an inverse theorem.

## Exact shared-insertion identity

For a fixed actual set C, let

```
r_a(x) = e_(d-1)(S[x,C without a])
V[a,j] = sum_x S[x,j] r_a(x), for a,j in C
R[a,c] = sum_x r_a(x) r_c(x)
s_a = sum_x r_a(x),  t_a = sum_(j in C) V[a,j],  u_a = V[a,a]
m = q-n.
```

For b outside C, the edit response is `delta_a(b)=sum_x S[x,b]r_a(x)-u_a`. The balanced symmetric conference identity gives

```
G[a,c] = sum_(b outside C) delta_a(b) delta_c(b)
       = q R[a,c] - s_a s_c - sum_(j in C) V[a,j] V[c,j]
         + u_c t_a + u_a t_c + m u_a u_c.
```

The insertion sums are `h_a=-t_a-m*u_a`. Consequently

```
Sigma = G/m - h h^T/m²
      = (q R - s s^T - V V^T)/m - t t^T/m².
```

The second line makes the cancellation of constant deletion offsets explicit. The only additional contraction beyond round7 is the full derivative Gram matrix R. The diagonal entries were already used by the marginal variance tool; the off-diagonal entries retain dependencies between deletion choices.

To remove the common response, write c_a for a covariance row sum and c for the sum of all entries. Then

```
(P Sigma P)[a,b] = Sigma[a,b] - c_a/n - c_b/n + c/n².
```

The implementation uses a single exact denominator for each matrix. The [projection review](coupling_verification.json) independently centers the literal edit matrix in both directions and forms its covariance, checking every entry on all four twin inputs.

## Streaming the derivative Gram matrix

The [compiled backend](covariance_backend.cpp) extends the frozen round7 synthetic-division method. For d=6, put E_j(x)=e_j(S[x,C]). Since S[x,a] is±1 except at x=a,

```
r_a(x) = A(x) + B(x) S[x,a] + delta_a 1_(x=a)
A = E5+E3+E1
B = -(E4+E2+E0)
delta_a = E5(a)-A(a).
```

Expanding the products r_a r_c reduces R to sums of A², transforms of AB at selected points, a symmetric weighted Gram matrix with weight B², and explicit zero-row corrections. V uses the corresponding weight B. The backend computes both weighted matrices in the same field-row pass. It uses O(q*n²) arithmetic and O(q+n²) storage, including the character vector. It is not an asymptotically optimality claim; other representations, including exact character convolutions, remain possible alternatives.

The prototype checks primality by trial division. At n<=64,d<=6, the coefficient bound `|B|<=637393` gives `q*B²<=4.063e18`, within signed64-bit range. A² and AB sums use signed128-bit accumulation. Final rational covariance and squared norms use Python integers. The boundary review includes an actual64-point quadratic-residue subset of F257 whose row at0 attains `|r_a(0)|=binom(63,5)=7028847`, testing the large-coefficient path.

[Scale results](scale_results.json) contain eight inputs reused from round7, with full inherited marginal agreement. The new backend took about25 and32 seconds on the two largest inputs in this run. Timings were collected on a shared machine and are not controlled comparisons with earlier runs.

## A necessary normalization control

Most covariance energy on the larger examples comes from a response shared by every deletion choice. The contrast squared norm is only about0.00031 or0.00025 of the full squared norm for the progression or seeded input at q=6,700,417,n=50. Interpreting the large raw off-diagonal share as special arithmetic structure would therefore be misleading.

The [reference model](reference_results.json) derives an exact benchmark for independent Rademacher signs, not a conference matrix. Orthogonality of distinct squarefree monomials gives

```
E[r_a r_c] = binom(n-2,5) + 1_(a=c)*binom(n-2,4).
```

Its projected-to-full squared Frobenius norm ratio is

```
25 / [25 + (n-1)(n-5)²].
```

Thus n³ times this independent-sign ratio tends to25. This is a proved identity for the stated reference model; no Paley asymptotic follows. Complete enumeration over sign vectors verifies the formula at n=6,7,8. At the largest saved field, the measured ratio is1.2294 times the reference for the progression and1.0018 times the reference for the single seeded set. These finite comparisons neither characterize arithmetic progressions nor quantify all random inputs.

## Reproduction and review scope

```
c++ -O3 -std=c++17 tooling_lab/round8/shared_insertions/covariance_backend.cpp -o tooling_lab/round8/shared_insertions/covariance_backend
python3 tooling_lab/round8/shared_insertions/query_covariance.py --q 17 --selected 0,1,2,3,4,6,10 --full-certificate
```

The [query command](query_covariance.py) returns exact marginal moments, covariance squared norm, and the projected covariance squared norm. `--full-certificate` includes both matrices and the contraction records. The lower-level algebra consumer requires genuine producer data; consistency checks alone do not prove that arbitrary supplied Gram data come from a graph.

Verification is separated by coverage:

- [Literal preflight](results.json): every entry for four actual twin inputs agrees with direct neighbour-target calculations.
- [Degree and boundary review](boundary_verification.json):25 cases,33,480 exact derivative-Gram and covariance entries, all degrees0–6, full outside insertion transforms, and six rejected invalid queries.
- [Streamed second implementation](verification.json): all584 new Gram entries at q<=65537;12 selected off-diagonal entries on each larger input,48 total; every one of206 Gram row sums; the pointwise derivative Euler identity on15,534,568 rows; and reconstruction of all7,506 large cross-moment formulas. The largest new Gram matrices were **not** completely recomputed by a second implementation.
- [Coupling controls](coupling_verification.json):392 covariance entries across eight affine controls and196 literal projected-covariance entries.

All original round7 artifacts remain frozen. These are separate arithmetic implementations and root-agent reviews, not human refereeing or Lean formalization. The exact identities are derived above; the [prior-art record](prior_art.md) distinguishes them from established covariance and contraction methods.
