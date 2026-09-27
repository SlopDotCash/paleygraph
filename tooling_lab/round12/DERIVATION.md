# Mathematical contracts

Let an affine space consist of f0 + t1 f1 + t2 f2 + t3 f3 over Fp, with basis polynomials of degree less than k and distinct evaluation coordinates x_i. Write v_i=(f1(x_i),f2(x_i),f3(x_i)). A triple determines its parameter vector exactly when its three evaluation columns have rank three. All guarantees below concern this declared space.

For agreement set A of size s, let B(A) count its independent triples. A uniformly selected unordered triple pins a given nearby member with probability B(A)/C(n,3). Origin f0 affects the coordinate equations, but not the pinning matroid.

## Multiplicities

Partition the nonzero columns by projective direction, with capacities c_i and occupancies a_i. Zero columns contribute nothing. Swapping a selected nonzero column for an unselected zero cannot increase B, so an optimum uses min(s,z) zeros, leaving t=max(0,s-z) occupied nonzero positions. Then

```
B(A) = F(a) = sum_{independent {i,j,k}} a_i a_j a_k,
0 <= a_i <= c_i,  sum a_i = t.
```

For distinct i,j, F(a+u(e_i-e_j)) has second derivative -2 partial_i partial_j F <= 0 along the feasible interval. If both occupancies are partial, its endpoints are integers, and at least one endpoint has value no larger than the starting value. Moving there removes at least one partial class without introducing any other. Repetition proves that a minimum exists with at most one partial class. This is a direct finite proof using a known rounding pattern.

The producer enumerates full-class subsets and a possible partial class. The standalone numerical verifier deliberately uses another route: recursively enumerate every feasible integer occupancy, including all allowed zero occupancies, pruning only when a nonnegative partial count is already at least the claimed target. A witness attains the target. A search limit cannot certify success.

## Pair spans and the large bound

In a simple rank-three evaluation matroid, every pair is independent. Its row span determines a unique rank-two flat L; each dependent triple belongs to exactly one such flat. The verifier groups every pair and checks that a group with m distinct endpoints contains exactly C(m,2) pairs. Therefore it has all pairs of that flat, without needing to enumerate triples. Set

```
D(A) = sum_L C(|A intersect L|,3),
d_v = sum_{L containing v} C(|L|-1,2).
```

For any s-set, 3D(A) is at most the sum of its global d_v values, hence at most the sum of the largest s degrees. Also D(A)<=D(all coordinates). These give an upper bound U on the maximum D(A); a selected witness gives a lower bound W. Consequently

```
C(s,3)-U <= min B(A) <= C(s,3)-W.
```

This uses elementary incidence counting. Both constructor and verifier require simplicity for this large method; the capacity method handles zeros and parallel columns separately.

## Sampling and complete test oracles

If B(A)/C(n,3)>=delta=a/b for every s-set, T independent uniform trials miss any fixed nearby member with probability at most (1-delta)^T. The declared space has p^3 members, so the union bound is at most p^3(1-delta)^T. The implementation chooses the smallest positive integer T satisfying

```
p^3 * (b-a)^T * 2^failure_bits <= b^T.
```

All arithmetic here is integer arithmetic. A saved trace verifies what was computed, not its randomness law. A caller-provided trace has no automatic probabilistic guarantee. Tests include one such trace that misses both true candidates.

For the large test words r, two distinct polynomials g,h of degree<k cover all coordinates: r(x) is either g(x) or h(x). Any other polynomial f of degree<k has at most k-1 roots of f-g and at most k-1 roots of f-h. Thus its agreement is at most 2(k-1). With k=64 and s=410, this is 126<410. Both supplied pieces actually agree at 512 coordinates, proving the entire ambient degree-bounded list consists of them. This special test oracle is not a theorem about arbitrary received words.
