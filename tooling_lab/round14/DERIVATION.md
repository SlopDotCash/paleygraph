# What a complete query transcript tells us

Fix H=f0+span(f1,...,fd), an affine polynomial space over Fp, and distinct evaluation coordinates. Write v_i=(f1(x_i),...,fd(x_i)). Partition the coordinates into blocks B_j such that every d-element subset of a block is independent. Blocks smaller than d have no queries. All claims concern this declared space; finding it is a separate problem.

## Coverage

Query every d-subset within each block. An agreement set avoiding every query has at most min(|B_j|,d-1) elements in block j. Conversely, selecting that many coordinates independently in each block avoids every query. Therefore the query family's independence number is exactly

```
alpha = sum_j min(|B_j|,d-1).
```

If alpha<s, every agreement set of size at least s contains an independent query. Solving that query's coordinate equations recovers its parameter vector. Thus all nearby members of H appear among the query labels, for every received word.

The block-profile optimizer uses this same sum. With t queried blocks and u unqueried coordinates, alpha=t(d-1)+u. For fixed t, choose the largest feasible u subject to alpha<s and leaving at least d coordinates in each queried block. Discrete convexity of C(b,d) makes balanced queried-block sizes cheapest. Enumerating t proves optimality within this restricted construction family and the block-size cap. It does not prove that a polynomial space admits such a partition.

## Agreement from occurrence and absence

Let theta be one parameter vector. In block B, let A be its actual agreement set with the received word, and h=|A|. Each queried d-set has a unique label. That label is theta exactly when the entire query lies in A. Consequently

```
number of queries labeled theta in B = C(h,d).
```

If theta appears, h>=d. Every element of A belongs to some d-subset of A, so the union of the queries labeled theta is exactly A. This reconstructs both the agreements and the nonagreements in that block without evaluating theta at the remaining coordinates.

If theta never appears, h<=min(|B|,d-1). Absence does not mean zero agreement. The complete query transcript is essential: omitting even one query can hide a true candidate with exactly d local agreements.

For an absent block, suppose we have explicitly tested t coordinates and observed m agreements. A valid remaining upper bound is

```
u_B = min(|B|,d-1, m+|B|-t).
```

For a present block use its reconstructed exact agreement count. Sum these quantities over disjoint blocks to bound theta's total agreement. Reject theta only when the sum is less than s. If it survives testing every absent block, the sum and the assembled support are exact. The implementation also materializes and checks the emitted full word.

The producer records every candidate's positive-block counts and support union, every explicit coordinate test, and the initial and final agreement bounds. The separate reviewer reconstructs every query label by integer Cramer determinants and evaluates every candidate at every coordinate. It checks the bounds and decisions against those complete words. Integer array arithmetic is used only after establishing d(p-1)^2+(p-1)<2^63.

## Why a single test can matter

Suppose alpha=s-1 and a candidate appears in only one block, with exactly d agreements there. Its initial upper bound is s. If an unqueried singleton is a nonagreement, that bound falls to s-1 immediately. This explains the observed one-test rejections; it assumes no randomness law for the received word.

Two shortcuts would be wrong. A true candidate can occur only once, so an unconditional frequency cutoff can discard it. A true candidate can disagree at the singleton if its positive blocks supply extra agreement, so requiring every candidate to pass a fixed guard coordinate is also invalid. The saved small boundary example and the large guard-mismatch stress word check both distinctions.

There is no uniform fast-rejection theorem here. In the late-near-miss stress word at s=128, the known origin polynomial appears once but has only 127 agreements. Its agreement deficit is placed in the last tested block, and rejection takes 1,008 coordinate evaluations. True boundary candidates also require testing all absent blocks. Full output materialization is accounted for separately and can duplicate earlier evaluations.

## Below the ambient Johnson threshold

For degree<k Reed-Solomon codes, the classical agreement threshold is sqrt((k-1)n). The old n=1024,k=64,s=410 tests lie above it. This round also uses s=192 and s=128, with

```
192^2 = 36864 < 64512 = 63*1024,
128^2 = 16384 < 64512.
```

The covering and transcript arguments still apply to the supplied low-dimensional spaces. They do not supply a global affine-space discovery theorem. The original two-piece words have complete ambient oracles because any third degree<64 polynomial agrees at most 126 times. The planted stress words have complete checked lists inside the supplied space only; the two-piece ambient oracle does not apply to them.

## A remaining arithmetic packing obstruction

Zero columns must be unqueried, which the current compiler handles. More generally, a queried arc block can contain at most one coordinate from a nonzero projective class. If that class has c coordinates, a profile with t queried blocks and u unqueried coordinates requires c<=t+u. For several classes with capacities c_i and z zero coordinates it necessarily requires

```
z + sum_i max(c_i-t,0) <= u.
```

This necessary condition is stronger than the zero-column check. It does not enforce all higher-rank incidence constraints. A follow-up compiler can use it to reject impossible profiles before search and use a class-to-block flow to allocate coordinates, then verify the remaining arithmetic rank conditions. This is a concrete next tooling gap; flow and incidence capacity arguments themselves are known methods.
