# Coding transfer audit: preserve which roots a block reuses

The existing local work already retains codewords, maximal agreement supports, scalar specializations and affine tracks. The [overlap certificates](../round6/overlap_preflight/README.md) also allow supports larger than a basis. A generic syndrome, matroid or overlap construction would repeat that work. The read-only [syndrome-chord note](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-rate-quarter-syndrome-chord-2026-07-10.md) further warns that abstract signatures can violate the geometry of actual quotient columns; its referenced Lean results were not rebuilt here.

A current source adds a useful hypothesis to audit. Jeronimo, Liu and Rajpal's [January2026 v1 paper](https://arxiv.org/html/2601.10047v1) develops proximity-gap machinery for folded codes using evaluation-kernel dimensions and pruning. Its Lemma4.1 on PDF page10 displays a dimension bound summed over all nonzero basepoints. The original [Guruswami–Kopparty construction](https://www.math.toronto.edu/swastik/subspace-designs.pdf), Section4.2, chooses basepoints whose associated evaluation sets are disjoint. These are different contracts. This audit addresses the displayed v1 lemma and transfer assumptions; it does not determine the status of the newer paper's main result or proposed repairs.

An exact counterexample to that displayed all-basepoints inequality is tiny. Work over F17, take γ=3 of order16, block length m=2, degree limit k=4, and U spanned by

```
f(X)=(X−1)(X−3)(X−9)=X³+4X²+5X+7 mod17.
```

The evaluation kernel H_a requires vanishing at both a and3a. For all a≠0, U∩H_a has dimension1 at a=1 and a=3 and dimension0 elsewhere. The sum is2. The displayed right-hand side at dimension d=1 is (k−1)/m=3/2. The two successful blocks share the root3. Restricting to the disjoint blocks starting at γ^(2j) leaves only one successful block. The polynomial from roots{1,2,4} has the same degree and number of roots but no successful overlapping block.

[coding_block_audit.py](coding_block_audit.py) verifies this through finite-field ranks and Horner evaluation. A separate power-evaluation census checks every5,220 monic polynomial of degree below4: the maximum is2 over all basepoints and1 over the chosen disjoint family. These statements include every one-dimensional polynomial subspace, since each nonzero polynomial has a unique monic representative. The source PDF and HTML text agree on the all-basepoints quantifier; no screenshot inspection was available.

The overlap correction is elementary and already within root-counting mathematics. If each coordinate appears in at most μ blocks of size m, any nonzero degree-<k polynomial has at most μ(k−1)/m vanishing blocks. For consecutive windows on a multiplicative cycle with a proper root set of size t, the sharper maximum is max(0,t−m+1): decompose the root set into runs. At most k−m windows can vanish when k−1 is below the cycle length. This explains the counterexample without inventing a new decoding bound.

There is also an immediate scalar-alphabet obstruction to a direct transfer. For a d-dimensional polynomial subspace, evaluation at one coordinate has rank at most1, so its kernel has dimension at least d−1 at every coordinate. The normalized average loss is at least(d−1)/d. The exact U=span{1,X} example at length16,k4 has loss1/2, exceeding rate1/4. Folding changes the observation map and the distance notion; it cannot silently be supplied to an ordinary-code problem.

The next candidate tool is an incidence-preserving block verifier on **actual candidate polynomial clusters**, with an explicit cost for reused coordinates and for lost agreement under blocking. It would have to improve a measured instance using more than the generic root-count bound. The current artifact diagnoses failed imports; it is not that new general tool. Historical novelty is unestablished, and no external message or submission was made.
