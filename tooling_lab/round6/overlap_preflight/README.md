# Overlapping support certificates: finite success, exponential limitation

This separate prototype recovers complete scalar/codeword lists for two saved
hard pencils that cannot satisfy the disjoint partition cap. It chooses tracks
using their **actual maximal joint supports** and certifies their overlap by
checking every agreement set. The covering mechanism and interpolation proof
were already present in rounds 2 and 4. The added local capability is choosing
a smaller support cover, plus a standalone certificate readback.

No production source or cap check was changed. This is a finite exact tool,
with known mathematical foundations and unestablished historical novelty for
the particular implementation. It does not establish a prize statement or an
efficient algorithm at growing block length.

| Saved pencil | Complete candidate tracks | Selected tracks | Selected support sizes | Complete nodes | Qualifying scalars |
|---|---:|---:|---|---:|---:|
| F1009, sample 0, n16 k8 s11 | 12,870 | **118** | 118 of size 8 | 3 | 3 |
| F17, sample 3, n16 k8 s11 | 12,450 | **49** | 1 of size 10, 35 of size 9, 13 of size 8 | 18 | 12 |

Both start from all 12,870 interpolation bases and cover all 4,368 eleven-sets.
The old anchored construction selected 1,287 bases; its count is a number of
bases, which need not equal its number of distinct tracks. The new greedy
construction is not asserted to minimize the number of tracks. Its deletion
pass could remove no further selected track, which establishes only
irredundance, not optimality.

The 118 size-eight supports also form a universal coordinate-base cover for
the same n,k,s. That design can be cached and reused on another pencil by
interpolating its 118 bases, without rebuilding the full candidate catalog.
The 49-track F17 improvement uses its particular larger joint supports and
does not transfer to arbitrary inputs. Neither observation removes the
fixed-rate lower bound below.

`results.json` records construction and baseline equality. The two
`*.certificate.json` files include full inputs, selected polynomial intercepts
and slopes, maximal joint supports, a covering witness for each eleven-set,
and every scalar/codeword node with maximal agreement support. The
`*.catalog.json` files export all candidate joint supports with their first
interpolation bases; these bases and the bound input uniquely specify each
affine polynomial track. The source ledgers are shared read-only artifacts
from `../cover_barriers/`, generated barycentrically and independently replayed
using Newton differences. Only selected tracks are reconstructed again here
by coefficient Lagrange interpolation.

The complete-list contract is elementary. Let the domain consist of distinct
field elements, let 1 <= k <= s <= n, and let each selected track be
g_j(z,X)=a_j(X)+z*b_j(X), with both component degrees less than k. Its verified
joint support S_j contains exactly the coordinates where a_j agrees with u0
and b_j agrees with u1. Require

```
for every A subset [n] with |A| = s,
there is a selected j such that |A intersect S_j| >= k.
```

If f has degree less than k and agrees with u0+z*u1 on at least s coordinates,
choose any s of them as A. The condition gives k distinct roots of f-g_j(z),
so f=g_j(z). This proves that checking all selected tracks finds every
qualifying polynomial, at **every** scalar. Conversely every emitted node is
checked by its full agreement support. A track's coordinate agreement equation
is affine in z, so a coordinate contributes always, never, or at one scalar.
When its common support has size at least s, every field scalar qualifies.
Specializations of different tracks may coincide; deduplication retains the
polynomial and its entire agreement support.

Equivalently, the union of the k-subsets of the selected supports is a
k-uniform hypergraph with no independent s-set. This is a finite covering
condition, not an unproved sampling assumption. No complete candidate catalog
is needed to **verify** an exported cover: the input polynomials, all support
checks, and all s-set intersections suffice. Catalog completeness is needed
for claims about the maximum possible support size and lower bounds on the
entire candidate family.

For these two pencils the exact maximum joint supports are respectively 8
and 10. The size-only lower bounds on every disjoint partition cap are 14
and 13, both at least the threshold 11. Thus no certificate in the production
strict partition-cap family can work. The separate overlap certificates do
work. This demonstrates a representation barrier in that family, not a
mathematical impossibility for complete lists. Round 2's general basis-cover
completeness proof already allowed this possibility; its particular block
construction did not optimize these overlaps.

`verification.json` reports a second implementation that imports neither the
constructor nor its scalar-bucket routine for hard-certificate readback. It
uses powers to evaluate polynomials, literal set intersections, and direct
enumeration of every scalar/selected-track pair: 119,062 for F1009 and 833 for
F17. It checks all 8,736 agreement sets and the complete exported node lists.
It also compares 73 tiny cases against every degree-bounded codeword at every
scalar, including k=1, s=k, s=n, and a deterministic whole-field pencil.
There are 196 tiny nodes and 44 whole-field tracks. Six corruptions are
rejected: false support, false polynomial, a missing required track, a missing
node, a duplicate node, and noncanonical domain encoding.

Run with ordinary Python 3.10+ and the standard library:

```sh
python tooling_lab/round6/overlap_preflight/run_overlap.py
python tooling_lab/round6/overlap_preflight/verify_overlap.py
python tooling_lab/round6/overlap_preflight/counting_bounds.py
```

The toy interface intentionally limits prime fields to p <= 5000 and the
number of s-sets to 200,000. Generation took 25.4 and 31.2 seconds, excluding
the already completed support census; independent verification took 36.8
seconds. These measurements demonstrate certificate reduction, **not a speed
improvement** over the older tiny census. There is no extension-field or
succinct large-cover implementation in this lane.

The necessary counting bound for a support of size m uses

```
V(m) = sum_{i=k}^{min(m,s)} C(m,i) C(n-m,s-i),
```

with impossible binomial terms omitted. It covers exactly V(m) of the
C(n,s) agreement sets, so any family must satisfy sum_j V(|S_j|) >= C(n,s).
In particular, with all supports of size k, the number of tracks is at least

```
C(n,s) / C(n-k,s-k) = C(n,k) / C(s,k).
```

For n16 k8 s11 this lower bound is 78. It applies to the complete F1009
catalog. It **does not apply** to the F17 catalog: V(8)=56, V(9)=336, and
V(10)=1056. A maximum-size bound gives only 5 tracks there; accounting for
the fact that there is only one size-10 candidate improves the necessary
bound to 11. Neither lower bound is claimed tight. The 49-track certificate
shows why reducing every joint support to one arbitrary basis loses useful
coverage information.

For k=Rn+O(1), s=alpha*n+O(1), and 0<R<alpha<1, the pure-k bound is

```
2^(n * [H2(R) - alpha*H2(R/alpha)] + O(log n)).
```

The exponent is positive: the same ratio is the product of k factors
(n-i)/(s-i), each at least n/s>1. At R=1/2 and alpha=11/16 the entropy
exponent is about 0.418821 bits per coordinate. The exact count is already
at least 86,892,595 tracks at n64 and 10,165,894,361,641,766 at n128.
`counting_bounds.json` includes larger exact ratios. These are conditional
representation lower bounds for the pure-k-support regime; we have **not**
constructed a growing sequence of hard pencils with that regime, and this
is not a lower bound for arbitrary decoding algorithms. Larger joint
supports require the corrected volume formula.

The next useful question is whether actual hard-pencil support families have
enough large supports, or another algebraic structure, to avoid this cost.
A larger greedy enumeration of pure-k supports would only spend exponential
work on a known obstruction. Finding a concise verifier of a structured
cover would improve verification cost but would not by itself reduce the
number of distinct tracks that must be represented in this certificate.

The primary and local prior-art audit is in `prior_art.md`; it explicitly
records the close information-set decoding and lotto-design precedents.
