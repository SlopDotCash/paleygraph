# Complete compressed agreement-support census

This iteration turns the round1 provenance instrument into a complete census
that can skip almost all agreement subsets. It works over prime and extension
fields and retains every qualifying scalar/codeword pair with its maximal
agreement support. It is an exact algorithm for an individual input stack,
not a worst-case bound or a prize proof.

Run with Python 3.10+ and the standard library:

```sh
python3 run_experiments.py
```

`compressed_support.py` is the reusable implementation. `results.json` contains
the optimized run; `initial_results.json` preserves the first run before the
measured interpolation bottleneck was fixed. Both import field arithmetic and
read-only oracle data from `../../proximity/`. No round1 or source-repository
files were changed.

## Mechanism and completeness proof

Let `u(z)=u0+z*u1`, with `n` distinct evaluation points and RS dimension `k`.
We seek **all** pairs `(z,f)` where `deg(f)<k` and `f` agrees with `u(z)` on
at least `s>=k` coordinates.

Choose disjoint coordinate blocks `B_1,...,B_m`; coordinates outside the
blocks are ignored when choosing interpolation bases. Verify the certificate

```
outside + sum_j (k-1) < s.
```

Enumerate every `k`-subset inside each block. Every size-`s` agreement support
contains at least one of these bases: otherwise it has at most `k-1` points
in each block and at most `outside` others, contradicting the certificate.
`verify_cover` checks this universal capacity inequality and disjointness.

For each selected base, interpolate `u0` and `u1` separately to global
codewords `A` and `B`. All codewords matching a scalar word on that base lie
on the affine track `A+z*B`. At any coordinate, agreement is

```
(u0-A) + z*(u1-B) = 0.
```

It holds for all scalars, for no scalar, or for one explicitly computed
scalar. Let `t` be the number of always-agreeing coordinates. Bucket the
remaining coordinates by their matching scalar. A scalar qualifies on this
track exactly when its bucket has size at least `s-t`.

**Completeness:** any qualifying polynomial has an agreement support
containing a selected base. Uniqueness of interpolation at `k` distinct
points puts that polynomial on the corresponding track. The coordinate
histogram therefore emits it. Conversely every emitted node is directly
checked by its full support. Duplicate tracks and duplicate scalar/codeword
nodes are merged, but different codewords at one scalar are retained.

If `t>=s`, the whole field qualifies. Such tracks are stored symbolically
by `(A,B,common_support)`; fields of size at most 1000 are also materialized
for the toy verification runs. Thus a common interpolation base is never
incorrectly treated as either a single scalar or a globally correlated
stack. Only a **global** common support of size at least `s` establishes
the latter.

These are elementary correctness proofs for the algorithm, not Lean-checked
theorems. The saved small certificate proves the support family covers all
agreement supports. Checking the complete output of a large census still
requires replaying the deterministic computation; the certificate is not
a succinct proof of every negative interpolation-track result.

## Constructing the cover

For each admissible integer `m`, retain

```
t = n-s+m*(k-1)+1
```

coordinates and partition them as evenly as possible into `m` blocks. The
outside-plus-block capacity is exactly `s-1`. Select the `m` that minimizes

```
number_of_bases = sum_j binomial(|B_j|,k).
```

The single-block case is the anchored interpolation cover. The partitioned
case can be much smaller. This is a concrete elementary Turán-cover
construction; optimality among all covers is not claimed.

## Actual iterations and failures

The first candidate—enumerate every `(k+1)`-subset—has a larger combinatorial
cost than the old size-`s` census on both round1 shapes. Even enumerating
every `k`-subset was costly: on the F17 `n=16,k=8,s=11` sample, the first run
took **9.62 seconds for 12,870 bases**, whereas the full parity oracle took
**0.432 seconds for 4,368 agreement subsets**, with identical output.

The cover reduces that shape to 1,287 bases, but fewer bases alone do not
guarantee faster execution. In the initial run, the F41 cover cases took
7.7–11.7 seconds because each base reconstructed full Lagrange polynomials
with `O(n*k^2)` products. The second implementation uses prefix/suffix
products to obtain the same interpolation matrix in `O(n*k+k^2)` products.
All 221 tested matrices over prime and extension fields agree exactly with
the retained direct implementation. The optimized F41 cases took
approximately 2.4–3.1 seconds in this run. These timings include machine
load and are illustrative, not benchmark guarantees.

## Exact equality with the previous oracle

The optimized run reproduces **all 16 round1 stacks**, including every
scalar, every distinct decoded codeword, and every maximal agreement
support—not merely the scalar counts. It also reproduces all seven
extension tests in F3, F9, and F729, including whole-field correlated cases.
The cover construction was exhaustively checked on **9,216** small agreement
supports. The algorithm never samples scalars, interpolation bases within
the selected cover, or codewords.

The input stacks remain selected examples: exhaustive output for one stack
does not mean an exhaustive search over all stacks.

## Larger simplified problems

All following censuses are complete by the checked cover construction and
the interpolation argument. They were also repeated after rotating the
coordinates by seven places, giving a different complete cover; all returned
scalar/codeword/support triples agreed after undoing the rotation. This
checks a separate cover choice, not an independently implemented decoder.

| Field | n | k | s | Selected bases | Full agreement subsets |
|---|---:|---:|---:|---:|---:|
| F65537 | 32 | 8 | 20 | 4,290 | 225,792,840 |
| F65537 | 64 | 4 | 40 | 61 | 250,649,105,469,666,120 |
| F65537 | 64 | 4 | 12 | 16,815 | 3,284,214,703,056 |
| F65537 | 128 | 8 | 80 | 3,465 | 434,033,831,785,996,446,590,578,979,888,886,600 |
| F729 = F_(3^6) | 56 | 4 | 35 | 55 | 1,346,766,106,565,880 |

Each input was constructed to have two specified agreeing scalars. The
census found exactly those two scalars and two codewords in each example;
this is a complete negative check for additional codewords on those
specific stacks. The `n=64,k=4,s=12` example has `s^2<n*(k-1)`, so its
agreement threshold is below the usual Johnson agreement threshold. It is
a simplified instance in the difficult decoding range, not a proof of
worst-case behavior there. Other rows deliberately exercise less difficult
thresholds while testing the implementation's scale and completeness.

The extension example uses a subgroup of order 56 in F729, actual degree-six
arithmetic from the validated round1 adapter, and stack entries in the full
extension field. It does not reproduce the production characteristic or
domain size.

## Complexity and remaining barrier

With `C` selected bases, the current implementation uses roughly
`O(C*(n*k+k^2))` field multiplications plus inversions and track bookkeeping.
Its deduplication table can use `O(C*n)` field elements. Field operations
over polynomial-basis extensions have their own degree-dependent cost.
Whole-field tracks are represented symbolically when expanding them would
be too large.

The cover can still be exponential in dimension at fixed rate. Examples:

- `n=32,k=8,s=20`: 4,290 partition bases versus 125,970 anchored bases.
- `n=32,k=8,s=12`: no useful partition improvement; 3,108,105 bases remain.
- `n=1024,k=64,s=410`: roughly `3.57e47` bases remain. Only this arithmetic
  cost was calculated; that census was not run.

Thus the experiment closes the **tooling** gap between full agreement-subset
enumeration and a complete compressed per-stack oracle on moderate cases.
It does not provide a production-scale algorithm. A next meaningful target
is a better certified covering family, or a symbolic orbit-compression
method that certifies many interpolation tracks at once while retaining
cross-track overlap. An uncertified adaptive sample is not an adequate
replacement for the completeness certificate.

## Prior art and novelty

The local audit inspected SYZ53's subset-parity probe, the SW1 lift probe,
and searches for covering designs and interpolation covers in the probes
and knowledge base. A05's earlier note also describes full subset
interpolation through a supercode bridge; it is not a new idea to replace
decoding by interpolation candidates. No implementation of this particular
complete partition-cover-plus-track ledger was identified in that bounded
search. Some referenced historical probe paths are absent from the current
checkout, so that is not an exhaustive history audit.

Primary prior art includes [Gordon, Kuperberg and Patashnik, *New
constructions for covering designs* (1995)](https://arxiv.org/abs/math/9502238)
and [Ben-Sasson et al., *Proximity Gaps for Reed–Solomon Codes*
(2020)](https://eprint.iacr.org/2020/654). Covering designs/Turán systems,
Lagrange interpolation, and affine-pencil analysis are established tools.
The result here is a new local executable combination with exact evidence;
no claim is made that the mathematical construction or combination has
never appeared elsewhere.
