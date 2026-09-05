# Independent round4 mathematical and implementation review

Reviewed 2026-09-05. These are exact finite tooling results, not prize proofs or a certification of historical originality. The candidate sources belong to the other lanes; this directory contains independent oracles and reviews only.

The conditional compiler passed all **2,454** previously generated direct conditional records, four full-mark boundaries, 316 additional generic-degree cases over F5, and eight mark-permutation checks. The records include complete Paley13 and Paley17 subset censuses and actual GF49 Paley/Peisert constructions. The generic cases include odd degrees, degree zero, all columns selected, and all selected columns marked. See [conditional_compiler_review.json](conditional_compiler_review.json) and [review_conditional_compiler.py](review_conditional_compiler.py).

The direct data generation is independently implemented in [direct_conditional_oracle.py](direct_conditional_oracle.py), [origin_marked_oracle.cpp](origin_marked_oracle.cpp), and [selected_marked_oracle.cpp](selected_marked_oracle.cpp). It does not import a conditional-moment or coloured-cell compiler. Prime cases enumerate every n-subset and aggregate every contained mark of size at most three. The GF49 sixth-kernel oracle uses ordinary row products on origin-containing sets, with the exact 49/6 translation correction only for the unmarked baseline. The n=7 GF49 cases directly enumerate all 163,185 completions of each selected three-mark set. Complete input/output values are in [direct_conditional_oracle.json](direct_conditional_oracle.json).

## Conditional coefficient audit

For a row pair, write `G_y=(1+a_y t)(1+b_y u)`. Once M is fixed, its factors are deterministic. Expand the remaining-column product as

```
product_(y outside M) [1+v*(G_y-1)].
```

A term of v-degree k uses exactly k distinct outside columns. Its probability of inclusion is `(n-m)_k/(q-m)_k`, where m=|M|. In particular, `a_y*b_y*t*u` uses one column and has v-degree one. There is no extra binomial or factorial multiplier. The compiler handles the cases n=m and n=q consistently; truncation avoids division beyond the outside population size.

The ordered cell-pair convention is also correct. If U,V are cells and `B(U,V)=sum_(x in U,y in V) S_xy`, then

```
N_plus(U,V) = (|U||V| - diagonal(U,V) + B(U,V))/2,
diagonal(U,V) = |U| if U=V, else 0.
```

The moment expansion sums ordered row pairs. Applying another factor 1/2 to that sum would be incorrect. Marked rows are singletons, and their mutual signs to other rows are already determined by those rows' mark patterns.

One- and two-mark nulls hold for even-degree Paley targets because the full affine group preserves the target, including nonsquare dilations. They also hold for the particular Peisert49 graph: graph automorphisms together with a complementing isomorphism act transitively on ordered distinct pairs, and an even-degree target is unchanged by complementation. These symmetry explanations should not be asserted for every conference matrix without checking its automorphisms. The direct oracles verify the relevant nulls rather than assuming them.

When n=d, `T_d(C)` is a single row-product sum. For distinct rows the products `S_xy*S_zy` have `(q-3)/2` positive, `(q-1)/2` negative, and two zero positions before removing M. Thus its conditional second moment depends on the mark cell sizes alone; the mutual row sign disappears. Consequently the n=d Paley/Peisert differences do not demonstrate the value of coloured edge counts. The n=7 and n=8 tests are essential controls.

## Three-mark contraction reduction

The final [contraction_reduction.py](../marked_moments/contraction_reduction.py) specialization is mathematically sound for its stated input contract. Let

```
f_I(x) = product_(i in I) S_(x,m_i), for x outside M; 0 on M.
h2 = f_{01}+f_{02}+f_{12}; h3=f_{012}; Q=h2^T S h3.
```

For a bulk pattern pair u,v let `D(u,v)=K_+(u,v)-K_-(u,v)`, the difference of its two hypothetical off-diagonal moment coefficients. Global sign reversal makes D odd under `(u,v)->(-u,-v)` for every degree d, since two degree-d factors contribute the even sign `(-1)^(2d)`. Therefore its Walsh expansion on the six signs has odd total degree. Row exchange and simultaneous permutation of the three marks preserve its coefficients.

An irreducible term with degree at least two on both sides must have bidegree (2,3) or (3,2); all six such coefficients are equal, say c. Every other contraction is determined by the cell data using

```
S f_empty = -sum_j S e_(m_j),
S f_{i} = q e_(m_i) - 1 - sum_j S_(m_j,m_i) S e_(m_j).
```

Pairing with a vector that vanishes on M removes the `q e_(m_i)` term. The remaining sums are Walsh cell sums with the displayed boundary corrections. Bulk diagonal pairs and pairs touching a marked row are handled separately. The six high-degree terms, each initially multiplied by 1/2 from the +/- decomposition, pair under symmetry of S, leaving **c*Q**, not c*Q/2. Hence

```
E[T_d(C)^2 | M subset C] = cell_only(q,n,d,M) + c(q,n,d)*Q(M).
```

The guard q>=17 is sufficient for both hypothetical bulk signs: the smallest nonzero paired type multiplicity before mark subtraction is `(q-5)/4`, at least three. This is a domain guard for the current Walsh calculation; it is not a statement that no identity can be obtained for smaller q.

The cell-only API still requires genuine conference-matrix provenance and a certified Q. Its consistency checks do not establish realizability of an arbitrary supplied histogram or contraction. The complete relation-table wrapper checks all three returned rational moments against the previously audited compiler.

[review_contraction_reduction.py](review_contraction_reduction.py) checks 384 low-degree contractions by literal matrix sums and 18 complete conditional moment cases over 31,372 enumerated subsets. It includes both GF49 graphs, odd and even degrees, and n=d degeneracies. It also compares the relation-free API's mean, second moment and variance with direct enumeration. The two F29 marks `(1,4,0)` and `(1,9,0)` give

```
cell_only = 941207/7475,  c = -2/7475,
Q = -42 and 86,
second_moment = 72407/575 and 188207/1495,
difference = 256/7475.
```

This is a genuine finite information separation: the same literal cell histogram can require one further signed contraction. It does not bound worst-case subsets, prove Q is independent of every simpler possible invariant, or establish historical novelty. See [contraction_review.json](contraction_review.json).

There is also a useful negative scale result. Center the two vectors by subtracting their means. Since S annihilates constants and has norm sqrt(q) on their orthogonal complement,

```
Q^2 <= q * (||h2||^2 - sum(h2)^2/q) * (||h3||^2 - sum(h3)^2/q).
```

Multiplication by c² gives a squared radius around the cell-only term. This is ordinary centered Cauchy–Schwarz applied after the exact reduction. The outward rational square-root calculation is correct. [review_envelope.py](review_envelope.py) passes 424 exact rounding checks, six saved envelope checks, and three direct matrix checks. The bound's relative radius is at most `16845159/(5*10^18)=3.3690318e-12` for q65537,n16 and `11553/10^19=1.1553e-15` for q1000033,n31. These refer to the fixed mark `(0,1,2)` and the reported cell data. They bound the relation correction to a conditional second moment, not T6 on individual subsets. Thus the genuine finite separation should not be described as a large correction at these critical-scale samples. See [envelope_review.json](envelope_review.json).

## Symbolic interpolation batch certificate

The support-batch proof is sound under distinct canonical field coordinates, `1<=k<=s<=n`, and degree less than k. If two polynomial pairs `(A,B)` and `(A',B')` had k coordinates in the intersection of their common input-agreement sets, then both differences vanish at those k points. Polynomial uniqueness forces the two pairs to coincide. Thus distinct tracks have disjoint families of selected k-bases.

Each track is backed by one selected input basis; interpolation readback proves its degree and its agreement with the input there. Its complete selected-base multiplicity is the sum over disjoint cover blocks of `binom(|block intersect common|,k)`. Distinctness and equality of the total with the entire selected-cover cardinality prove that no selected basis is missing. This argument does not trust the ZDD or enumerate skipped bases.

Any s-coordinate agreement support contains a selected k-base by the cover certificate. At a qualifying scalar, the reconstructed track and the qualifying codeword agree on that basis, so they are equal. Every coordinate equation is either automatic, impossible, or satisfied at exactly one field scalar. The bucket calculation is therefore complete, including whole-field tracks. A symbolic whole-field export describes the whole family; an unmaterialized `nodes` list by itself does not.

`verify_batches` proves track coverage only. The final `verify_export` additionally binds field/problem metadata, checks canonical input encodings, recomputes every scalar and support branch, and checks node and track exports. These are distinct responsibilities; a call to the narrower routine is not a full export audit.

[review_support_batches.py](review_support_batches.py) uses independent F3, F5 and GF9 arithmetic and exhaustive coefficient/scalar enumeration. It passed **261** full node comparisons, **1,335** explicit selected-base interpolations, **162** whole-field cases including **87** symbolic ledgers, and **2,200** sequential ZDD-subtraction checks. It exercises partition, anchor and full-k covers, k=1, k=s, and s=n. See [support_batches_review.json](support_batches_review.json).

The speedup parameter is the number of distinct interpolation tracks, not the number of final codewords. Tracks may meet at one scalar while differing elsewhere. The recorded interleaved failure is therefore substantive. The large planted cases are certified finite outputs; they do not show a uniform small-track bound for arbitrary or prize-obstructing inputs.

## Replay

```
python3 tooling_lab/round4/novelty/review_conditional_compiler.py
python3 tooling_lab/round4/novelty/review_contraction_reduction.py
python3 tooling_lab/round4/novelty/review_envelope.py
python3 tooling_lab/round4/novelty/review_support_batches.py
```

The first replay reads the saved independent oracle data and verifies its source hashes. To regenerate that more expensive data, run `direct_conditional_oracle.py`. Prior rounds remain unchanged. The primary-source and local-prior-art audit is in [prior_art.md](prior_art.md).
