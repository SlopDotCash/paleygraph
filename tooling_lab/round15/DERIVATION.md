# Capacity-constrained arc covers and local repair

The input is a degree-bounded affine polynomial space with direction columns v_i in F_p^d. A block of size at least d is queried on every d-subset, and each query must be independent. A partition with t queried blocks and u coordinates outside those blocks has maximum query-avoiding set size t(d-1)+u. All statements below concern this complete-block model with block size at most64, not arbitrary query hypergraphs.

## A necessary capacity constraint with a constructive relaxation

Normalize each nonzero column by its first nonzero entry. Columns with the same normalized vector form a projective class of size c_j. Let z be the zero-column count. A queried block can contain no zero column and at most one member of each projective class: otherwise a d-subset containing the offending column or pair is dependent. Thus

```
u >= z + sum_j max(c_j-t,0).
```

For fixed t, query count decreases as u increases. The largest allowed value is min(s-1-t(d-1), n-td). If this value satisfies the capacity constraint, take n-u queried coordinates with at most t per class. Such a selection exists because sum min(c_j,t) >= n-u. Assign consecutive selected members of each class to consecutive positions in a cyclic sequence of t blocks. Each class visits each block at most once, and block sizes differ by at most one. The remaining coordinates are unqueried singletons.

The forward difference of binomial(b,d) is binomial(b,d-1), which increases with b>=d. Exchanging one coordinate from a block at least two larger than another cannot increase the query cost. Therefore balanced sizes minimize cost for fixed t,u. Enumerating t gives the exact minimum in the projective-capacity relaxation. When an actual independent-query partition reaches that lower bound, it is optimal in the complete-block model itself. This is elementary counting and discrete convexity, not a new general packing theorem.

## Actual-domain obstruction and attainment

On the saved1024-point domain over F65537, let P have62 specified roots and use span(1,P,XP). Its degrees are0,62,63. At a root the evaluation column is(1,0,0); elsewhere the second and third coordinates determine x through their ratio, so all962 other projective classes are distinct. The space has rank three and no zero columns.

At agreement96, the unconstrained size optimizer chooses t47,u1 with cost70070. Its capacity for a single class is only48, so no permutation makes this profile independent. The constrained optimum is t33,u29 with cost136155. An independent enumeration of every feasible t,u checks173 profiles and gives the same minimum.

The initial cyclic allocation has two dependent triples. A bounded restart search finds a cover on attempt23. Every final query is checked by modular elimination and by integer Bareiss determinants. This attains the lower bound. All46 recorded failed-restart witnesses are independently checked, but the nonzero queries from those failed attempts are not independently replayed.

## Local exchange evidence

Swapping coordinate i of block A with coordinate j of block B preserves sizes. After checking projective capacities, only triples containing i or j can change. If R is the set of old dependencies containing an exchanged coordinate and T the set of new dependencies containing a moved coordinate, then

```
new_dependencies = (old_dependencies - R) union T.
```

The prototype accepts a swap only when |T|<|R|. Accepted swaps therefore strictly decrease the nonnegative defect count. This proves termination of accepted descent, not that a repair always exists or that local descent will find it. Unqueried singleton blocks are available for exchanges too; a polynomial control requires this case.

Two accepted swaps remove the two actual dependencies. They use1682 local query tests, giving137837 search tests including the first full scan. The observed restart search uses3131565 tests. Both then pay for final audits and certificate checks; this is an arithmetic-count comparison, not a controlled wall-clock speedup. The separate reviewer rescans both entire affected blocks, checking17110 triples, and verifies both final covers.

## Higher-rank limits and decoding scope

Every rank-r flat contains at most r coordinates from an arc block. Its size f therefore requires f<=rt+u. Projective counts capture only rank-one instances and are insufficient in general. An interpolated polynomial control on F7 has six distinct points on one line. At agreement5 the parallel relaxation allows t2,u0, but6>2*2. A five-coordinate rank-deficient witness directly obstructs every all-s-set basis cover. At agreement7, a singleton exchange produces a one-query cover. The repair tool correctly distinguishes bounded failure from a checked success; it does not yet compile general flat constraints.

The actual1024-point cover supports arbitrary received words inside the declared space. Its agreement96 is below the ambient Johnson threshold, but no list outside the supplied space is asserted. In particular, the older two-piece oracle only excludes third polynomials above agreement126 and is inapplicable here. Boundary, near-miss and guard-mismatch words are checked by independent query solutions and full evaluation of every generated candidate. They also exhibit linear verification cost for an individual candidate.

The repeated class is an opportunity as well as a packing obstruction. If m_b class coordinates require the equation theta.w=b, conditioning on that value removes one parameter and contributes exactly m_b agreements. A heavy-value branch needs s-m_b agreements outside the class; if all remaining class values have frequency at most h, their shared residual branch needs s-h outside agreements. This leads to the received-word-dependent prototype in [round16](../round16/README.md). It changes the class of query certificates, so any improvement there does not contradict this round's all-s-set complete-block optimum.
