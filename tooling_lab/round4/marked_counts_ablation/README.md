# Cell-profile ablation and a single-convolution replacement

**Cell sizes alone do not determine the conditional second moment once n>6.** This directory contains actual marked Paley inputs with identical complete signed/zero cell profiles and different conditional T6 second moments. It also contains an actual Paley49/Peisert49 comparison with that property. Every claimed difference is independently confirmed by enumerating every completion of the marks.

The follow-up implementation computes the one extra integer contraction retained by the root agent's three-mark reduction, using one exact character convolution rather than the full coloured relation table. This is a local mathematical and implementation specialization built from established Walsh expansions, conference identities, and NTTs. Its historical novelty has not been established; it supplies no worst-case spectral bound or prize proof.

## What the ablation holds fixed

For the sign matrix S and an ordered three-mark list M, a row x has pattern
\((S_{x,m_1},S_{x,m_2},S_{x,m_3})\). The profile records the size of **every** such cell, including the three zero-containing marked-row cells. We identify profiles under a permutation of mark coordinates, then retain explicit reordered mark lists for which the profiles match literally. We do not identify profiles under arbitrary sign flips, discard the zero cells, or merely compare an unordered list of sizes. The marked triangle is consequently held fixed too.

The target is
\[
T_6(C)=\sum_x[t^6]\prod_{y\in C}(1+S_{xy}t),
\]
averaged over all n-sets C containing M. At n=6, the compiler's cell-only shortcut is exact and all the matched-profile pairs below agree. The ablation asks if that sufficiency survives at n=7 or n=8.

The initial scan used the already independently enumerated oracle for every triple in Paley13 at n7 (286 triples) and Paley17 at n8 (680 triples). It found no matched-profile difference. These small null cases are preserved in `saved_oracle_scan.json`; they motivated the larger-prime scan rather than a sufficiency claim.

## Exact witnesses

| comparison | left marks | right marks | left E[T6²] at n7 | right E[T6²] at n7 | right minus left |
|:---|:---|:---|:---|:---|:---|
| Paley29 versus itself | (1,4,0) | (1,9,0) | 72407/575 | 188207/1495 | −256/7475 |
| Paley49 versus itself | (0,1,3) | (0,1,7) | 13991443/54395 | 13993491/54395 | 2048/54395 |
| Paley49 versus Peisert49 | (0,1,3) | (0,1,2) | 13991443/54395 | 13989907/54395 | −1536/54395 |

For the first pair, the common profile is exactly:

| pattern | size | pattern | size |
|:---|---:|:---|---:|
| (−1,−1,−1) | 4 | (1,−1,−1) | 3 |
| (−1,−1,1) | 2 | (1,−1,1) | 4 |
| (−1,0,1) | 1 | (1,1,−1) | 4 |
| (−1,1,−1) | 3 | (1,1,0) | 1 |
| (−1,1,1) | 4 | (1,1,1) | 2 |
| (0,−1,1) | 1 | | |

Both Paley29 marks have n6 second moment 7323/325. At n7 they have the same mean −213/1495, so their variances differ by the same amount as their second moments. At n8 the second-moment difference changes sign and is 256/16445. The two q49 comparisons also agree at n6, with second moment 53953/1265; their n8 gaps are respectively 75776/685377 and −18944/228459.

`witness.json`, `twins_witness.json`, and `cross_twins_witness.json` contain both exact cell/relation records, all three compiled moments, search provenance, and the compiler source hash. The corresponding `*_matrices.json` files contain the actual matrices. The q49 field encoding and Paley/Peisert construction come from the frozen `round3/pair_type_twins/results.json`; this is not a freely perturbed relation-table example. Cross-twin search scanned all 1,128 Paley49 triples containing zero before comparing Peisert49. The searches stop after a witness, so they are not complete classifications of the larger fields.

## Independent verification

`direct_completions.cpp` imports no compiler or cell-count code. It first checks each supplied sign matrix directly: symmetry, diagonal zeros, no other zeros, row sums, and \(SS^T=qI-J\). It then enumerates all C containing each M for n6, n7, and n8. It computes the pointwise sixth coefficient by polynomial updates and bit-counted positive/negative entries, and accumulates exact integer totals, squares, and full histograms.

`verify_witnesses.py` converts those totals to rational moments and compares them with the compiler. All checks pass:

- Paley29 pair: 166,660 completions across its six cases.
- Within-Paley49 pair: 3,098,238 completions.
- Cross-twin pair: 3,098,238 completions.

The combined 6,363,136 completion evaluations include repeated work when a marked graph appears in more than one comparison; they are not all distinct subsets. `direct_validation.json` preserves counts, integer sums, exact rational results, and complete conditional histograms. The oracle is deliberately bounded to q≤63 and 6≤n≤8; its totals fit signed 64-bit arithmetic in that supported domain.

## What additional information survives

The root agent's separate `marked_moments/contraction_reduction.py` reduces the three-mark conditional second moment to
\[
\mathbb E[T_6(C)^2\mid M\subset C]
=F_{q,n}(\text{cell profile})+c_{q,n}Q_M,
\qquad Q_M=h_2^T S h_3,
\]
where, **outside M**,
\[
h_2(x)=S_{x,m_1}S_{x,m_2}+S_{x,m_1}S_{x,m_3}+S_{x,m_2}S_{x,m_3},
\qquad h_3(x)=\prod_{i=1}^3S_{x,m_i},
\]
and both functions are explicitly zero on M. In particular h2 is 3 or −1 off M and h3 is ±1 there. Forgetting the h2 masking would give a different contraction.

The scalar values are:

| marked input | Q_M |
|:---|---:|
| Paley29 (1,4,0) | −42 |
| Paley29 (1,9,0) | 86 |
| Paley49 (0,1,3) | 162 |
| Paley49 (0,1,7) | −94 |
| Peisert49 (0,1,2) | 354 |

`summarize_ablation.py` recomputes each of these directly from all matrix pairs and checks the reduced formula against the original compiler at n6, n7, and n8. The exact differences equal \(c_{q,n}\Delta Q\). For p29, c is −2/7475 at n7 and 2/16445 at n8; at n6 it is zero. Thus the ablation demonstrates a specific information loss and the compressed replacement repairs it for the stated three-mark computation. It does not show that arbitrary marked moments, additional marks, higher distributions, or exceptional subsets are determined by this scalar.

## Exact single-convolution backend

`single_contraction.cpp` builds the cells, h2, and h3 in O(p) space and computes
\[
(Sh_3)(x)=\sum_y\chi(x-y)h_3(y)
\]
by one cyclic integer convolution, then takes the dot product with h2. The output is `{p, marks, cells, Q, metadata, sampled_convolution_rows}`. No relation table is constructed or consumed by that executable.

The standard NTT helper is adapted from the frozen round4 count backend. It proves the field parameter prime, enforces p≡1 mod4 and three distinct canonical marks, verifies the NTT prime/root, and enforces length L≤2²³. With modulus ℓ=998244353 and h3 supported on p−3 entries,
\[
\|h_3\|_1=p-3,\qquad |(Sh_3)(x)|\le p-3<\ell/2.
\]
These checked bounds make centered modular recovery exact. The contraction obeys \(|Q_M|\le\|h_2\|_1\|h_3\|_1\); all intermediates fit signed 64-bit integers under the field/length guards. No floating-point value is used in an arithmetic result.

The backend uses exactly **three NTT transforms**: the character vector, the signed h3 vector, and the inverse transform. The prior three-mark full-relation layer used 17 transforms. Both methods use O(p log p) time and O(p) memory for a fixed mark list; the improvement is a smaller required arithmetic input and a smaller constant, not a new asymptotic convolution algorithm.

`verify_contraction.py` compares the single-convolution result against all full exact relation inputs in the frozen count lane with three marks, plus both actual Paley29 witnesses: 13 inputs in total. Every match passes. It independently recomputes all Q sums through p1297 and 33 complete signed-convolution rows at p65537 and p1,000,033, using Euler's criterion for every character value. It also sends only the cells and Q to `reduce_from_cells`; its exact mean, second moment, and variance agree with the full relation compiler at n6/7/8 and n31 whenever valid. The reduction requires q≥17; p13 is checked for Q only.

| p; marks (0,1,2) | Q_M | NTT length | transforms |
|---:|---:|---:|---:|
| 1297 | 3330 | 4096 | 3 |
| 65537 | 196626 | 262144 | 3 |
| 1000033 | −571662 | 2097152 | 3 |

At the million prime the exact recovery margin is 996244293, the certified contraction absolute bound is 1500087001260, and the recorded maximum child RSS is 38,109,184 bytes (36.34 MiB). The RSS is a maximum over the validation runner's child-process lifetime, not an isolated allocation trace.

`contraction_benchmark.json` records three interleaved full-count/single-Q repetitions at each large prime. On the shared, concurrently loaded machine, the million-prime median wall times were 22.74 seconds and 3.96 seconds (5.75× in that run). The earlier frozen full-count run took 1.80 seconds, illustrating substantial load variation: do not combine timings from different runs or treat this ratio as guaranteed. The reproducible structural comparison is 17 versus 3 transforms. Both backends' exact outputs were checked during benchmarking.

## Reproduction and scope

```sh
python3 search_profiles.py
python3 scan_relations.py
python3 scan_relations.py --twins-only
python3 scan_relations.py --cross-twins
clang++ -std=c++17 -O3 -Wall -Wextra -o direct_completions direct_completions.cpp
python3 verify_witnesses.py
clang++ -std=c++17 -O3 -Wall -Wextra -o single_contraction single_contraction.cpp
python3 verify_contraction.py
python3 benchmark_contraction.py
python3 summarize_ablation.py
```

Run in this directory with the recorded Python runtime `/opt/homebrew/opt/python@3.14/bin/python3.14` (Python 3.14.6; NumPy 2.4.4 for the search). Set `OPENBLAS_NUM_THREADS=1` for the search. The independent pointwise and signed-convolution oracles do not use NumPy or floating transforms. Sources, inputs, outputs, and the parent compiler/reduction are hashed. Other round4 inputs and all round1–3 artifacts remain untouched.

The finite-field transform and bounded centered-residue recovery are established in Pollard, *The Fast Fourier Transform in a Finite Field* (1971), Sections 1–3. [Primary paper](https://luca-giuzzi.unibs.it/corsi/Support/papers-cryptography/2004932.pdf). The three-mark reduction uses known Walsh/Boolean Fourier and conference-matrix identities; its derivation and scope are documented in the parent lane. No exhaustive historical search has established that this particular specialization has never appeared before.

The demonstrated result is precise: matched cell profiles can lose conditional second-moment information; the extra contraction restores it for this three-mark compiler and can be obtained exactly at million-prime scale. Controlling that contraction uniformly, or using it to bound exceptional sets, remains a separate mathematical problem.
