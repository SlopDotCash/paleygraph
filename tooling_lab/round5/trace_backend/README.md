# Exact prime-Paley trace inventory for the third-moment compiler

This backend computes the complete normalized cubic character-sum array with one exact cyclic convolution, retains its joint distribution with the two remaining row-edge signs, and exports the power sums needed by the generic degree-six third-moment compiler. It reaches p=1,000,033 in 6.08 seconds on the recorded shared-machine run, with 1,750 joint bins and about 40.3 MB peak child memory.

The cubic trace vector, elliptic interpretation, NTT, and trace moments are established and already appear in the local research. The contribution here is a validated exact backend and correctly normalized interface to the new global-third-moment compiler. It is not a new elliptic trace theory, a historical originality certificate, or a prize proof.

## Required schema and sign conventions

For a prime p≡1 mod4, put \(S_{xy}=\chi(x-y)\) with \(\chi(0)=0\). Normalize the three ordered rows to (0,1,t), t≠0,1. The exact edge ordering inherited from the frozen preflight is
\[
(e_{01},e_{02},e_{12})=(S_{0,1},S_{0,t},S_{1,t})
=(1,\chi(t),\chi(t-1)).
\]

The trace variable is
\[
\tau(t)=\sum_x\chi\bigl(x(x-1)(x-t)\bigr).
\]
For the Legendre curve its projective point count is p+1+τ, so the conventional Frobenius trace is **−τ**. Odd power sums retain this sign; it must not be silently reversed.

Each `inventory_p*.json` contains:

```text
{
  p, q,
  records: [{edges: [1, chi(t), chi(t-1)], tau: integer, count: integer}],
  edge_power_sums: [{edges: [...], powers: [sum tau^j for j=0,...,6]}],
  total_power_sums: [sum tau^j for j=0,...,6],
  normalized_parameter_count: p-2,
  ordered_distinct_row_weight: p*(p-1),
  n6_global_third_moment: {...},
  normalization: {...}, metadata: {...}, provenance: {...}
}
```

`count` is the number of normalized parameters t having that joint edge/trace value. Counts and power sums have **not** been multiplied by p(p−1). Records are ordered lexicographically by edges and then τ, with every nonzero count included. Power lists are indexed by the exponent, including j=0. All powers and global-moment arithmetic use Python arbitrary-precision integers. The million-prime sixth power sum requires 83 bits; a signed 64-bit accumulator would overflow.

The raw `tau_p*.bin` files contain a little-endian uint32 p followed by p little-endian signed int32 values, including the two singular parameters. They satisfy τ(0)=τ(1)=−1. Those two values are **omitted** from the normalized joint inventory and its power sums; repeated rows are added separately.

## Why the weight is p(p−1)

Given distinct ordered rows (a,b,c), write d=b−a and t=(c−a)/d. Under x=a+dy, each of the three row-sign functions becomes χ(d) times the corresponding normalized function. This follows from multiplicativity and χ(−1)=1. Their mutual signs and τ change by that same common sign. Degree-six row coefficients, and hence their third products, are invariant under this simultaneous sign complement.

For every t≠0,1 there are p choices of a and p−1 choices of nonzero d. The normalized contribution is therefore multiplied exactly once by p(p−1). There is no additional factor six: these parameters already describe ordered triples. The independent small oracle instead enumerates unordered distinct row triples and multiplies those by six, obtaining the same answer.

The output is explicitly for **prime Paley matrices**. The executable rejects p49 and other composite parameters. It neither implements the additive group of an extension field nor labels a Peisert triple correlation as a Legendre elliptic trace. The q49 twins have their separate, already verified adapters in the preflight and parent compiler.

## Exact convolution and recovery

Set \(a(x)=\chi(x)\chi(x-1)\). Because χ(−1)=1,
\[
\tau(t)=\sum_x a(x)\chi(x-t)
=\sum_x a(x)\chi(t-x)=(a*\chi)(t).
\]

Both vectors are zero-padded to the next power-of-two length L≥2p−1. The usual NTT product produces their linear convolution modulo ℓ=998244353; folding positions t and t+p gives the cyclic convolution. The backend uses exactly three transforms.

The copied NTT helper proves ℓ prime, checks its factorization ℓ−1=2²³·7·17, verifies the primitive root 3 against those prime factors, and checks every transform length. The field parameter is independently proved prime by trial division. The inverse transform uses Fermat inversion in that certified field. Runtime guards remain active in optimized builds.

Here \(\|a\|_1=p-2\), so every true convolution coefficient has magnitude at most p−2. The checked inequality \(2(p-2)<\ell\) makes the centered residue the unique integer answer. The supported transform domain is L≤2²³; oversized inputs are rejected rather than reusing an invalid transform root. Modular products stay below ℓ²<2⁶⁰ and fit signed 64-bit integers. Hasse's bound τ²≤4p is checked **after** recovery as a consistency test; recovery does not depend on assuming the Hasse interval.

The producer checks the literal edge orientation and the nonnegative integral Walsh inversion of all eight bulk triple-column types for every nonsingular t. This ensures that the exported tuple has the precise row ordering required by the generic coefficient compiler.

## n=degree=6 global third moment

For m nonzero signs of sum τ, define the standard elementary-symmetric/Krawtchouk coefficient
\[
K_6(m,\tau)=\frac{\tau^6-(15m-40)\tau^4
+(45m^2-210m+184)\tau^2-15m(m-2)(m-4)}{720}.
\]
The implementation compares this polynomial against its literal positive/negative binomial formula for every trace bin. It then exports
\[
\mathbb E[T_6(C)^3]=\frac{
-p\binom{(p-1)/2}{3}
-3p(p-1)\binom{(p-3)/2}{3}
+p(p-1)\sum_{t\ne0,1}K_6(p-3,\tau(t))
}{\binom p6}.
\]

The first term is from all-equal rows, the second from exactly two equal rows, and the last from all-distinct rows. All three exact integer numerators are retained. A second aggregation using only the power sums must agree exactly with the histogram-based sum. This n6 shortcut is not substituted for repeated-row templates at general n>6.

| p | joint bins | exact global E[T6³] at n6 | backend seconds |
|---:|---:|:---|---:|
| 13 | 7 | −173/33 | 0.00133 |
| 17 | 7 | 17/91 | 0.000195 |
| 29 | 10 | 271/117 | 0.00122 |
| 101 | 17 | 65705/14259 | 0.000515 |
| 1297 | 63 | 111595286003/23328340049 | 0.00222 |
| 65537 | 448 | 253618616559703373/51233137428508399 | 0.400 |
| 1000033 | 1750 | 41714776856822242424687/8334316710159187856293 | 6.063 |

The timings are individual measurements on a shared, concurrently used machine, not asymptotic performance claims. At the million prime L=2,097,152, the centered-recovery margin is 996,244,291 and the recorded peak child RSS is 40,304,640 bytes. The runner labels RSS as a maximum over its child-process lifetime. It does not represent an isolated allocation trace for every small case.

## Independent checks and controls

`verify_inventory.py` uses the standard library and independently constructs every character value with Euler's criterion. It imports no NTT, preflight implementation, or moment compiler.

- Every τ value is recomputed by a literal cubic character sum at p13,17,29,101,1297. At p65537 and p1,000,033 it directly checks 33 complete rows each, using 2,162,721 and 33,001,089 integer summands. Large arrays therefore have sampled independent trace validation, not a second independent full convolution.
- It rebuilds every joint histogram and every power sum from the saved complete binary array. It checks literal triple-column types for all nonsingular parameters through p101 and three selected parameters in each larger field.
- It checks τ(1−t)=τ(t) and τ(1/t)=χ(t)τ(t) at **every** nonsingular parameter in every array. Modular inverses are independently constructed. Affine/sign-complement identities are checked exhaustively through p101 and on fixed point samples at larger p; their general validity follows from the displayed prime-field algebra.
- At p13,17,29 it independently computes all distinct and repeated row-triple contributions using direct polynomial updates on the actual row products. At p13 and p17 the complete joint inventories and moments also match the frozen preflight, whose separate acceptance check used complete six-column-set enumeration.
- Guard tests reject composites, wrong congruence, missing arguments, and an oversized transform domain before creating the requested binary output.

The universal checks M0=p−2, M1=2, and M2=p²−2p−3 hold in every saved inventory. The last follows from the conference convolution identity: \(\|a*\chi\|^2=p\|a\|^2-(\sum a)^2=p(p-2)-1\), followed by removal of the two squared singular values. These low moments are controls, not new trace-moment discoveries. The individual sign-class zeroth counts also match (p−5)/4 for (++), and (p−1)/4 for each other class.

## Complexity, reproduction, and provenance

The producer uses O(p log p) arithmetic time and O(p) memory. Its shared exact convolution produces the full array once. Aggregation uses O(p) time. The sufficient compressed power inventory has four sign classes and seven integer powers per class; the optional full joint histogram is much larger but remains O(√p) bins under the standard Hasse support. The raw p-entry array is retained for replay and future adapters. No six-column subsets are enumerated by the producer.

```sh
/usr/bin/clang++ -std=c++17 -O3 -Wall -Wextra -o tooling_lab/round5/trace_backend/trace_inventory tooling_lab/round5/trace_backend/trace_inventory.cpp
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round5/trace_backend/run_inventory.py
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round5/trace_backend/verify_inventory.py
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round5/trace_backend/summarize.py
```

The runner accepts `--p 1297 65537` for a selected input list. Runtime: Apple clang and Python3.14.6, with no third-party Python dependencies. The NTT header was copied byte-for-byte from `round4/marked_counts_ablation/ntt_exact.hpp`; both hashes are recorded. Inputs, sources, executables, binary arrays, inventories, and validation files carry hashes. Earlier rounds and the preflight are unchanged.

The prior-art and ensemble-normalization audit is in `round5/preflight/README.md`. In particular, the local September4 classical experiment already implemented the elliptic trace convolution, and the September5 spectral-coupling research already examined square-restricted Legendre moments. The present data adapter retains their attribution. Pollard's finite-field FFT and centered integer-recovery mechanism are standard; the generic parent compiler's finite tau-degree argument supplies the reason that powers through six suffice. An implementation that computes these diagnostics exactly does not by itself control tails, extrema, or the full moment distribution.
