# Exact arithmetic-marked Paley relation counts

This layer supplies the coloured-cell input for the conditional moment compiler in `tooling_lab/NEXT_ITERATION.md`. For fixed marks M it counts every ordered row pair by both mark patterns and its mutual Paley sign. It reaches p=1,000,033 in 1.80 seconds, using 11 cells and producing 193 nonzero relation records. The arithmetic is exact; it never materializes the p by p sign matrix.

The new local capability is an exact, reusable arithmetic-conditioned input whose cost does not grow with the number of subsets to be averaged by the downstream compiler. Number-theoretic transforms, integer recovery, character sums, and marked/coloured refinements are established tools. This implementation does not claim historical novelty, a spectral upper bound, or a prize proof.

## Consumer contract

For a prime p congruent to 1 modulo 4, let

\[
S_{xy}=\chi(x-y),\quad \chi(0)=0,\qquad
U_\alpha=\{x:(S_{xm})_{m\in M}=\alpha\}.
\]

Marks are an **ordered** list of distinct canonical field elements, so the pattern coordinates retain their identities and their order. Marked rows are included: their patterns contain a zero. Empty cells are omitted. All pairs are ordered and include diagonal pairs.

```json
{
  "p": 13,
  "marks": [0, 1, 2],
  "cells": [{"pattern": [-1, -1, -1], "size": 1}],
  "relations": [{
    "left_pattern": [-1, -1, -1],
    "right_pattern": [-1, -1, -1],
    "sign": 0,
    "count": 1
  }]
}
```

This fragment illustrates the schema, not the entire p13 output. `cells` are lexicographically sorted. `relations` are lexicographically ordered by left pattern, right pattern, then sign −1,0,+1; **every positive count appears exactly once**. Missing relations have count zero. Additional `metadata`, `sampled_rows`, and `provenance` keys do not alter this contract.

`sampled_rows[i].cell_character_sums[j]` is the exact sum over the jth cell of \(\chi(x-y)\), where x is the recorded sampled row. These samples support independent checking; they are not the mechanism used to calculate the relation counts.

The main files are `counts_p{13,17,101,1297,65537,1000033}_marks_0_1_2.json`. At p101 the export also includes no marks, one-mark and two-mark controls, and triples (0,1,3), (0,1,4), (1,2,3), (0,4,8), and (0,2,4). `summary.json` indexes all files and their hashes.

## Exact mechanism and bounds

For each cell V, compute the cyclic integer convolution

\[
g_V(x)=\sum_{y\in V}\chi(x-y),\qquad
B_{UV}=\sum_{x\in U}g_V(x).
\]

Let \(d_{UV}=|U|\) when U=V and zero otherwise, and \(N_{UV}=|U||V|-d_{UV}\). The required counts are

\[
N_0(U,V)=d_{UV},\qquad
N_+(U,V)=\frac{N_{UV}+B_{UV}}2,\qquad
N_-(U,V)=\frac{N_{UV}-B_{UV}}2.
\]

The implementation checks parity, nonnegativity, symmetry under exchanging U,V, the zero global sum of every character convolution, and the total Paley counts \(p, p(p-1)/2,p(p-1)/2\). The independent verifier additionally checks these counts separately for each left cell.

The convolution backend adapts the existing round2 radix-2 NTT. Both inputs are padded to the next power of two L at least 2p−1. Linear coefficients at x and x+p are added modulo the NTT modulus to recover cyclic convolution. The shared character transform is computed once. Every nonsingleton proper cell costs one forward and one inverse transform. Singleton cells use an exact shifted character lookup, and a cell equal to the whole field uses its exact zero convolution.

The NTT modulus is Q=998244353=2²³·7·17+1. The backend proves Q and p prime by trial division, verifies that 3 has full multiplicative order using the prime divisors 2,7,17 of Q−1, and enforces L≤2²³. For each V it checks

\[
|g_V(x)|\le |V|,\qquad 2|V|<Q.
\]

Therefore the centered residue is the unique true integer coefficient. This is established integer-convolution recovery, not a rounding heuristic. Unsupported transform lengths are rejected; the current backend is not silently extended to larger primes.

All modular products are below Q²<2⁶⁰ and fit signed 64-bit integers. The field/length guard gives p<2²² for accepted primes, so aggregate pair counts and intermediate signed sums are below p²<2⁴⁴ in magnitude. `pair_count_bit_bound` records the bit length of p² as a conservative bound, not the largest observed individual count. There are no floating-point operations in the counts. Only timings are floating point.

## Validation and measured scale

`verify_counts.py` uses only the Python standard library. It computes characters independently using Euler's criterion, not square enumeration, and does not import the producer or use an FFT/NTT.

- It materializes independent direct sign matrices at p13, p17, and every p101 case, checking every cell and ordered-pair count. It also checks all p1297 ordered pairs without storing a dense matrix. Across the saved inputs this covers 1,794,878 ordered pairs.
- At p65537 and p1,000,033 it verifies every character value with Euler's criterion and directly recomputes 33 complete sampled rows. These compare 2,162,721 and 33,001,089 integer summands respectively. These are sampled independent large-field checks, **not** a second full independent convolution.
- Exact affine controls pass: one-mark and three-mark translation; square dilation by 4; and nonsquare dilation by 2 with every pattern sign and mutual relation sign reversed. The latter is sign covariance, not literal equality of the coloured schema. Any claim about invariance of a downstream statistic must account for its degree and this sign reversal.
- Nine rejection tests cover composite fields, wrong congruence, oversized transform domains, duplicate/noncanonical marks, and too many marks. Guards use runtime exceptions and remain active in optimized builds.

| p | cells | nonzero relations | transforms | L | backend seconds | runner wall seconds |
|---:|---:|---:|---:|---:|---:|---:|
| 13 | 9 | 113 | 9 | 32 | 0.000406 | 0.242 |
| 17 | 10 | 149 | 15 | 64 | 0.000217 | 0.00436 |
| 101 | 11 | 193 | 17 | 256 | 0.000377 | 0.00396 |
| 1297 | 11 | 193 | 17 | 4096 | 0.00230 | 0.00752 |
| 65537 | 11 | 193 | 17 | 262144 | 0.159 | 0.164 |
| 1000033 | 11 | 193 | 17 | 2097152 | 1.793 | 1.799 |

The million-prime run uses a recorded peak child RSS of 49,463,296 bytes (47.17 MiB); the runner explicitly labels this as the maximum over its child-process lifetime, not an isolated per-case measurement. It is the largest case in that run. Its counts fit the recorded 40-bit bound. The first small-case wall timing includes process-launch overhead and should not be extrapolated as algorithmic cost. Full independent validation took 17.79 seconds.

For r distinct marks the number of nonempty cells is at most 2ʳ+r: all unmarked rows have signs ±1, and each marked row supplies one possible zero-containing pattern. The implementation supports 0≤r≤6. With c cells, time is O(c L log L + pc + pr log c) and memory O(L+p+c²+rc). Fixed three marks give c≤11, hence linear memory and O(p log p) arithmetic work with this fixed-modulus backend. No collection of candidate n-sets is stored, and no memory term depends on how many sets the conditional moment compiler will average.

The optimization from 23 transforms (one shared plus two per each of 11 cells) to 17 handles the three necessarily singleton marked cells by direct shifts. This is an elementary implementation optimization, not a new transform algorithm.

## What arithmetic information is retained

At p101 the common-positive cell has size 12 for marks (0,1,2), 13 for (0,1,3), and 14 for (0,1,4). The complete profiles and relations also differ. These are ordinary marked character statistics, not a claim of a new phenomenon. Their value here is that they survive into the compiler input instead of being replaced by unmarked global pair counts.

This layer alone does not establish that the compiled conditional sixth moments differ, distinguish the Paley/Peisert examples, predict exceptional sets, or improve a worst-case inequality. Those are separate acceptance questions for the compiler. Extension fields and Peisert kernels are not accepted by this prime-field executable.

## Reproduction and prior art

```sh
/usr/bin/clang++ -std=c++17 -O3 -Wall -Wextra -o tooling_lab/round4/marked_counts/marked_counts tooling_lab/round4/marked_counts/marked_counts.cpp
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round4/marked_counts/run_counts.py --million
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round4/marked_counts/verify_counts.py
/opt/homebrew/opt/python@3.14/bin/python3.14 tooling_lab/round4/marked_counts/summarize.py
```

For a new mark list, use `run_counts.py --p 1297 --marks 0 1 4`. For the unmarked partition use `--marks` without values. The C++ executable also prints the complete schema directly, for example `marked_counts 1297 0 1 2`.

Runtime: Apple clang, Python 3.14.6; the producer and verifier require no third-party Python dependencies. Every export includes input/source/executable hashes. `validation.json` hashes every validated count file; `manifest.json` hashes final local artifacts.

The transform and centered-residue mechanism are explicitly established in Pollard, *The Fast Fourier Transform in a Finite Field*, Mathematics of Computation 25 (1971), Sections 1–3. [Primary paper scan](https://luca-giuzzi.unibs.it/corsi/Support/papers-cryptography/2004932.pdf), [publisher DOI](https://doi.org/10.1090/S0025-5718-1971-0301966-0). The publisher PDF returned a fetch error during this audit; the primary paper itself was read through the university-hosted scan. Modern public NTT implementations use the same standard modulus family; this source is adapted from the project's round2 backend, not copied from an external implementation.

The marked/coloured idea was already acknowledged in the local July11 record and in `NEXT_ITERATION.md`. This work implements the exact count interface proposed there. It does not establish that the interface or its application has never been invented elsewhere.
