# Pinning a declared polynomial space

Let U be a d-dimensional space of polynomials of degree less than k over Fp, evaluated at n distinct domain points. Fix a basis f1,...,fd and write vi=(f1(xi),...,fd(xi)). A d-coordinate set B identifies a unique member of any fixed affine translate of U precisely when its evaluation columns are independent. If a received word agrees with a particular member on a set A of size at least s, every independent d-subset of A is a successful pinning choice for that member. This statement supplies neither the affine clusters nor their coverage.

For a fixed s, define

```
b_U(s) = min over |A|=s of #{B ⊆ A : |B|=d and rank(vB)=d}.
```

Uniform sampling of a d-subset of the whole domain succeeds with probability at least b_U(s)/C(n,d). Sampling only inside A would use a different denominator and is not the interface here. The origin of the affine space does not enter the evaluation ranks. The proposed object is a restricted-basis minimum of a representable matroid; its general ingredients are established mathematics.

## Rank two: exact optimization

Separate zero columns from nonzero projective classes of capacities c1≥...≥cm. In a chosen set, let z be the number of zero columns, aj the number chosen from class j, and t=s−z. Its number of independent pairs is

```
Σi<j ai aj = (t²−Σj aj²)/2.
```

Replacing a selected nonzero column by an available zero cannot increase this count. Thus some minimizer takes min(s,z_total) zeros. For fixed t, the allocation that fills capacities in descending order majorizes every feasible occupancy vector: the sum of its largest h entries is min(t,Σj≤h cj), which bounds the largest-h sum of every other allocation. Convexity of x² now proves optimality. Class construction, sorting and filling cost O(n log n) ordinary operations; finite-field arithmetic costs are additional.

The separate verifier uses dynamic programming. If the already processed classes contain t selected points, adding a points from a new class creates exactly a·t pairs. The recurrence minimizes that count at each total occupancy. It also minimizes over every feasible number of selected zeros rather than assuming the constructor's zero-first rule. Direct determinants reconstruct the classes and check the attaining set.

## A universal degree-bound baseline

An h-dimensional subspace of degree-<k polynomials can vanish simultaneously at at most k−h distinct points: otherwise every polynomial would share a factor of degree greater than k−h, leaving fewer than h degrees of freedom. After j independent evaluations are selected, the subspace vanishing on them has dimension d−j. At least s−k+d−j coordinates in A extend independence. Counting ordered sequences and dividing by d! yields

```
b_U(s) ≥ C(s−k+d,d),  for s≥k.
```

The space P(X)·{polynomials of degree<d}, where P has k−d distinct roots in A, attains it. This is an elementary MDS/common-root bound, not a new theorem claim. The rank-two controls check 38 such extremal parameter choices.

## Rank three: retain shared coordinates on lines

This version requires a simple rank-three evaluation matroid: all columns are nonzero and every pair is independent. Each dependent triple lies on exactly one maximal projective line. If L ranges over nontrivial lines (at least three domain points), then

```
D(A) = ΣL C(|A∩L|,3),
b_U(s) = C(s,3) − max over |A|=s D(A).
```

The line-size list alone loses which coordinates the lines share. The explicit F17 pair below demonstrates that this loss affects b_U(s). Full line incidence determines the rank-three matroid; it is not a newly invented invariant or a compressed substitute for the entire matroid.

The optimizer uses exact integer branch-and-bound on chosen coordinates S and undecided coordinates V, with r further choices required. For each line put a=|S∩L|, b=|V∩L|. One upper bound on the final dependent count is ΣL C(a+min(b,r),3).

A second upper bound charges each newly completed triple equally to its h newly selected coordinates. In units of six, each such coordinate receives 6/h, for h=1,2,3. For a fixed undecided coordinate v on a line, at most hL=min(b−1,r−1) other new coordinates on that line can be selected. Therefore its charge from that line is at most

```
6 C(a,2) + 3 a hL + 2 C(hL,2).
```

Sum these charges over the lines through v to obtain wv. Every completion has at most

```
D(S) + floor((sum of the r largest wv)/6)
```

dependent triples. Taking the minimum of both bounds and C(s,3) is valid. This is a fractional incidence upper bound used in an ordinary finite search; no algorithmic originality is asserted.

Each tree node records its split, a proved upper-bound leaf, a unique completion, or an unresolved frontier. The verifier reconstructs the field columns independently, checks every triple by row reduction, verifies every line's maximality, and replays the full tree with integer arithmetic. A node budget can leave a certified interval. The budget is a stopping threshold; pending siblings needed to cover the search space can make the final node count slightly exceed it. No interval is labeled exact unless its endpoints coincide. There is no polynomial runtime guarantee.

## What the finite collision proves

On the nonzero elements of F17, consider the degree-<8 subspaces

```
U6 = span(1, X, X⁶+X²),
U4 = span(1, X, X⁴+X²).
```

Both have 54 two-point lines, eight three-point lines and seven four-point lines. For any j≥2, their number of rank-two j-subsets is ΣL C(|L|,j). Singletons have rank one, and every remaining nonempty subset has rank three. Consequently they have identical full subset-rank distributions and identical Tutte polynomials. This implication is known; the current primary-source audit records it explicitly.

Their codeword weight distributions are also equal:

```
weight       0   12   13   14    15    16
codewords    1  112  128  864  2048  1760
```

Their generalized Hamming weights are (12,15,16). A separate enumeration checks all 65,536 subsets and all 4,913 codewords of each space. Nevertheless b_U6(11)=147 and b_U4(11)=148. Even the means and variances of the two complete eleven-set pinning distributions agree; the third moments differ. Thus these global data do not determine this worst-case pinning query, already on these actual polynomial inputs.

This is a verified finite obstruction to a specified information reduction. It does not prove historical novelty, supply a decoder, establish a uniform asymptotic estimate, or transfer automatically to a Paley spectral problem. The connection to earlier Paley tooling is methodological: both studies must test what joint incidence their summaries discard.
