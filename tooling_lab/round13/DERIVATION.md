# Coverage is a different question from minimum pinning density

Let H=f0+span(f1,f2,f3) and v_i=(f1(x_i),f2(x_i),f3(x_i)). For any received word r, a member f in H has agreement set A_f={i:f(x_i)=r_i}. If an independent triple Q is contained in A_f, solving its three coordinate equations recovers f uniquely. Verifying the resulting polynomial against all coordinates filters out all spurious candidates.

Thus it suffices to find a family Q of independent triples such that every set of size s contains one member of Q. This is a Turán covering problem with a restriction to triples allowed by the evaluation matroid. It is not necessary to compute the minimum total number of independent triples inside an s-set.

Partition [n] into disjoint blocks B_j. Choose independent query triples Q_j inside each block. Define

```
alpha_j = max {|A|: A is a subset of B_j and contains no query in Q_j}.
```

For any set A avoiding every query, its intersections with the blocks avoid each local query family, so

```
|A| = sum_j |A intersect B_j| <= sum_j alpha_j.
```

Conversely, choose a maximizing avoiding set separately in each block. Their union avoids all queries and has exactly sum_j alpha_j coordinates. Consequently the query family's independence number is exactly that sum. If it is less than s, every s-set contains a query. This is a direct disjoint-union argument from ordinary hypergraph theory.

For a block of size b>=3 in which every triple is independent, querying all triples gives alpha=2 and costs C(b,3) queries. Singleton and two-coordinate blocks have no triples and alpha=b. The actual large partition uses

```
n = 201*5 + 3*6 + 1 = 1024,
queries = 201*C(5,3) + 3*C(6,3) = 2070,
sum alpha = 201*2 + 3*2 + 1 = 409 < 410.
```

The prototype does not assume blocks are of this form: it enumerates local query-free subsets, so it also handles zero columns, parallel columns and dependent triples. Blocks are currently capped at 12 coordinates to bound this explicit verification. The constructor uses integer bitmasks; the verifier uses set containment over literal combinations and Gaussian rank rather than determinant formulas.

Search uses a recorded seed and bounded shuffled partitions. Only the returned certificate and its verifier are needed to establish coverage. Search randomness has no role in the decoded-list guarantee. Failure to find a good partition establishes no lower bound on the optimal cover size.

The family census checks the partition theorem against all 4,368 eleven-subsets of each actual small space. The complete received-word census checks the interpolation consequence on every F5-valued length-five word. These are distinct checks: set coverage is uniform over received words, while the list census tests the implementation's reconstruction and filtering.

For dimension d, replace triples by d-element coordinate bases and require sum alpha_j<s. On a block whose every d-set is independent, alpha=min(d-1,|B_j|). The certificate logic extends immediately, but the existence and cost of a suitable partition do not. Enumerating all local subsets and all d-queries can grow exponentially in d. Existing general pruning theorems already address the underlying derandomization problem; this interface needs a practical and precisely delimited comparison.
