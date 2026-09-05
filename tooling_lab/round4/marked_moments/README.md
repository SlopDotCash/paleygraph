# Exact marked moments and the missing contraction

This instrument computes exact mean and second moment of `T_d(C)=sum_x e_d(S[x,C])` over all n-element sets C containing specified marked columns M. Here S is a symmetric balanced conference sign matrix: zero diagonal, off-diagonal entries±1, `S1=0`, and `S²=qI−J`. Prime Paley matrices and the checked Paley49/Peisert49 matrices meet this contract.

The useful new local finding is more precise than “keep arithmetic labels.” For three marks, all additional relation information needed by this second moment reduces to **one integer contraction**. Actual Paley triples have identical marked cell sizes and different values of that contraction, with different conditional moments. The contraction can be computed with one exact character convolution. A separate norm calculation shows that its contribution to this particular conditional average becomes extremely small in the selected large critical-size examples. These are finite identities, counterexamples and diagnostic bounds, not a prize proof or a certification that the mathematics has never appeared before.

## General compiler

Rows are partitioned by their complete sign vector at M. The general input contains the size of each cell and, for each ordered pair of cells, the number of row pairs having mutual sign0,+1,−1. These are produced by the [exact counting backend](../marked_counts/README.md), or directly from a small matrix.

For two rows with entries a_y,b_y, define

```
G_y(t,w) = (1+a_y*t)(1+b_y*w).
```

The product of their degree-d elementary coefficients is `[t^d w^d] product_(y in C) G_y`. Marked columns contribute their G factors deterministically. Outside M, form

```
product_(y outside M) [1+z*(a_y*t+b_y*w+a_y*b_y*t*w)].
```

Replace each z^k coefficient by `(n−|M|)_k/(q−|M|)_k`. The k counts **distinct outside columns**; an `a_y*b_y*t*w` term consumes one column. No factor choosing additional unused columns is missing: the falling-factorial ratio is exactly the inclusion probability of those k specified columns in the uniform completion of M.

The paired-row sign-type histogram follows from the conference identities and the mutual row sign. Subtract the known pairs at marked columns. Group equal remaining types, use multinomial powers, and truncate at t,w degree d and z degree at most2d. Every coefficient is an integer; only the final inclusion probabilities require rational arithmetic. Cache classes under simultaneous mark permutations, row exchange and global sign reversal.

The single-row version computes the conditional mean. The compiler validates compressed row/column marginals, symmetry, marked-column Gram identities, diagonal counts and signs to the marked singleton rows. These checks do not replace provenance of the underlying conference matrix. `matrix_counts` checks the full matrix identities; the prime backend independently validates its Paley character construction.

### A necessary special-case control

When n=d, the target is a single d-column kernel. For two distinct rows, the products `a_y*b_y` have `(q−3)/2` positive entries, `(q−1)/2` negative entries and two zeros, regardless of the rows' mutual sign. After subtracting marks, their conditional product expectation depends only on the marked patterns. The diagonal case is also determined by cell sizes. Thus **cell sizes alone determine the second moment when n=d**.

The implementation includes this separate shortcut and asserts its exact equality with the general compiler. Testing only n=d would have falsely suggested that the relation backend was adding useful information. The n7/n8 tests were essential.

## Three marks leave one unknown contraction

Write the three marks as m0,m1,m2. Define two vectors, both set to zero on the marked rows:

```
h2(x) = S[x,m0]*S[x,m1] + S[x,m0]*S[x,m2] + S[x,m1]*S[x,m2]
h3(x) = S[x,m0]*S[x,m1]*S[x,m2]
Q(M)  = h2^T S h3.
```

For every allowed q,n,d in this implementation with q≥17 and three marks,

```
E[T_d(C)^2 | M subset C] = B(cell sizes; q,n,d) + c(q,n,d)*Q(M).
```

B depends on the full signed cell-size profile, including the three marked singleton rows. It does not require any cell-pair relation counts. Q and c can both be zero. The `q≥17` implementation guard ensures that the hypothetical positive and negative paired-row histograms used by the Boolean expansion are nonnegative for every bulk pattern; it is a domain guard, not a claim that no smaller-field identity exists.

To derive the reduction, let K_s(u,v) be the conditional paired-row coefficient at bulk sign patterns u,v and mutual sign s. Expand `K_+(u,v)−K_−(u,v)` in Walsh characters of the two three-bit patterns. Global sign reversal forces odd total Walsh degree. Row exchange gives symmetry, and permutations of the three marks identify coefficients.

For a subset I of the marks, let f_I be its sign-product character outside M and zero on M. If one side of `f_I^T S f_J` has degree0 or1, the conference identities determine it from the cell sizes. Specifically,

```
S f_empty = −sum_i S e_(m_i)
S f_{i}   = q e_(m_i) − 1 − sum_j S[m_j,m_i] S e_(m_j).
```

The boundary term `q e_(m_i)` vanishes when paired with f_J because f_J is zero on M. The remaining sums are cell-size Walsh moments. Among terms with both degrees at least2, odd parity leaves only degrees(2,3) and(3,2). All six ordered coefficients are equal. Symmetry pairs the two orientations, giving exactly cQ; the factor1/2 in the sign split cancels the two orientations. Diagonal and marked-row contributions are already known from cell data.

[contraction_reduction.py](contraction_reduction.py) reconstructs the entire Walsh coefficient table exactly, checks these symmetries, eliminates the known terms, and compares the final formula with the general compiler. Its `reduce_from_cells` API accepts only q, marks, cell sizes, Q, n and d. The [optimized backend](../marked_counts_ablation/README.md) obtains Q with three NTT transforms, including one inverse, instead of the full eleven-cell counting pipeline's seventeen transforms.

## Actual information-loss witnesses

In Paley29, the ordered triples `(1,4,0)` and `(1,9,0)` have **identical literal signed cell-size profiles**, including every zero pattern. At n7,d6:

| Quantity | Marks(1,4,0) | Marks(1,9,0) |
|---|---:|---:|
| Cell-only term B | 941207/7475 | 941207/7475 |
| Coefficient c | −2/7475 | −2/7475 |
| Contraction Q | −42 | 86 |
| Conditional second moment | 72407/575 | 188207/1495 |

Their moments differ by `256/7475`, explained exactly by the contraction. At n6 the moments agree, as the shortcut requires. Every completion at n6,n7,n8 was independently enumerated with direct integer target evaluation.

Two further ablations use Paley49 marks(0,1,3) versus Paley49 marks(0,1,7), and Paley49 marks(0,1,3) versus Peisert49 marks(0,1,2). In each pair the signed cell-size profiles agree, n6 second moments agree, and n7/n8 moments differ. The ablation lane verified **6,363,136 completions** across all three pairs. This is an actual-matrix sufficiency failure, not an arbitrary alteration of counts that might fail to describe a graph.

The first attempted ablations at Paley13 and Paley17 were negative: identical cell profiles also had identical moments in those scans. Larger small cases were needed. See the [complete witnesses and independent validation](../marked_counts_ablation/README.md).

## Scale and a limitation the tool exposed

The general exact compiler evaluates the million-prime case after its compressed counts have been supplied; it does not enumerate field elements or n-sets. That example uses23 canonical relation classes and integers up to149 bits in its aggregated union coefficients. The arithmetic backend still scales with the field. The optimized one-contraction backend retains arrays of length proportional to q; this is not a constant-memory end-to-end algorithm.

| Prime q | Set size n | Q for marks(0,1,2) | Relative change of conditional second moment from the unmarked average |
|---:|---:|---:|---:|
| 1297 | 6 | 3330 | −1.78196×10⁻⁶ |
| 65537 | 16 | 196626 | −8.16764×10⁻⁶ |
| 1,000,033 | 31 | −571662 | −1.66823×10⁻⁷ |

The last column includes changes in B as well as the relation correction cQ. They must not be confused.

On the orthogonal complement of the constant vector, S has norm √q. Ordinary centered Cauchy–Schwarz therefore gives the exact bound

```
Q² ≤ [(q*||h2||²−sum(h2)²)*(q*||h3||²−sum(h3)²)] / q.
```

Multiplying by c² gives an exact squared radius around B for every conference matrix realizing the same marked cell data. [certified_envelope.py](certified_envelope.py) records the rational squared radius and rounds its square root outward with integer arithmetic.

For q65537,n16 the radius divided by |B| is below `3.4×10⁻¹²`; for q1,000,033,n31 it is below `1.2×10⁻¹⁵`. At n=d it is exactly zero. The additional relation statistic is real and necessary for exact answers, but in these cases even its conference-norm outer bound contributes very little to this conditional average. This does **not** bound individual T6 values, rule out another use of Q, or show that all arithmetic conditioning is weak. It identifies the limited effect of one precisely defined diagnostic.

The subsequent [exact coefficient study](../coefficient_structure/README.md) derives a closed formula for c(q,n,6) by a general-q polynomial recurrence. It explains the sign reversal between n7 and n8 in Paley29 and separates fixed-size, critical-size and positive-density regimes. In the critical regime n/q^(1/4) tending to a fixed alpha>0, it proves `c ~ −alpha^4/q²`; the norm bound then gives `|cQ| <= (sqrt(3)*alpha^4+o(1))/sqrt(q)` uniformly over genuine three-mark conference inputs. This uniform statement concerns the absolute correction, with no inferred lower bound on B or pointwise bound on T6.

## Verification and reproduction

The [independent compiler review](../novelty/conditional_compiler_review.json) checked2,454 saved complete conditional records, four full-mark boundaries,316 degree-general checks over Paley5, and eight mark-order controls. It includes odd degrees, degree0, n=q and fully fixed sets. An independent contraction review also checks direct bilinear sums, the boundary eliminations and direct completion moments. These are independent programs and agent review, not human refereeing or a Lean formalization.

```sh
python3 tooling_lab/round4/marked_moments/run_experiments.py
python3 tooling_lab/round4/marked_moments/certified_envelope.py
python3 tooling_lab/round4/novelty/review_conditional_compiler.py
```

- [conditional_moments.py](conditional_moments.py): general exact compiler, count checks, small matrix adapter and n=d shortcut.
- [results.json](results.json): scale cases, affine/null controls, degree and boundary checks.
- [contraction_reduction.py](contraction_reduction.py): the one-contraction identity and relation-free consumer API.
- [contraction_envelopes.json](contraction_envelopes.json): exact correction envelopes and source hashes.

The foundations are established: finite-population inclusion probabilities and conditional expectations; Walsh expansion; subconstituent/marked graph algebras; conference identities; exact character convolution. Earlier local work already calculated conditional fourth character moments with prescribed points. The specific executable second-moment compiler, the three-mark reduction and its finite witnesses are the local research outputs. Historical novelty remains unestablished; the [round4 source audit](../novelty/README.md) records the overlap and searches.
