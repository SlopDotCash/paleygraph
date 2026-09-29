# Spectral witness transport microscope

This prototype asks **what arithmetic structure is actually carried by a dangerous spectral direction?** It couples a complete near-edge eigenspace to translation tests, then exports integer witnesses for independent exact replay. It is an experimental diagnostic and hypothesis falsifier. It does not prove a Paley bound or the Proximity Prize, and no claim of historically unprecedented mathematics is made.

The concrete new work in this repository is the coupled implementation, its comparison against controls with exactly the same spectrum, its witness-extraction interface, and the documented refinement after two tempting interpretations failed. Ordinary eigenspace compression, additive energy, approximate stabilizers, subfield obstructions, and counterexample-guided refinement are established ideas.

## Why this tool was missing here

Pass 12 demonstrated that projection and anchor identities plus low-depth information can hide spectral obstructions. Pass 16 exhibited an actual square-field obstruction that preserves the genuine character structure. Pass 17 found that qualitative prime uncertainty gives an inadequate quantitative gap. Passes 24–25 showed that uniform-seed iteration misses most symmetry sectors. Pass 22 demonstrated the failure of enlarging fixed-order correlation lists alone. Those are inputs to this project, not new claims here.

The earlier scripts verify identities, moments and separately supplied structured witnesses. This tool starts with **all spectral sectors**, finds an extremal direction that optimizes a specified arithmetic feature, and examines whether that feature survives arithmetic scrambling that leaves every spectral moment unchanged. The tool never supplies a subfield indicator to the spectral optimizer. A separately supplied indicator is retained only as a baseline.

## Mathematical interface

For an actual finite field of order q, put

```
C = {x : chi(x) = chi(x-1) = 1}
H = S_C / sqrt(q) - J_C / q
Z = 2 H / sqrt(3).
```

Thus `H=2B-I` in the local two-anchor notation. Let V be the spectral subspace with `|lambda(H)| >= ||H|| - 0.025`, and U an orthonormal basis of V, extended by zero outside C. For each nonzero additive translation tau_t, compute

```
A_t = sym(U^T tau_t U)
M_V(t) = max over unit f in V of |<f, tau_t f>| = ||A_t||.
```

The equality is the real symmetric variational principle. It is basis invariant. It also uses the same vector on both sides; the generally larger singular norm of `U^T tau_t U` would allow two different vectors and answer a different question. All finite-field translations are scanned. The maximizing eigenvector of A_t supplies the witness.

The tail includes both spectral signs. Consequently it guarantees large `||Hf||`, but not a large absolute Rayleigh quotient for every mixture. **Every exported witness has its own actual Rayleigh quotient recomputed**; the implementation never infers that a mixed-sign vector has a large quotient merely because it belongs to the tail.

The revised interface rounds the selected vector to an integer vector z, then computes every translation overlap using integers:

```
a(t) = sum_x z(x) z(x+t)
n = sum_x z(x)^2
T(z) = {t : 5 |a(t)| >= 3 n}.
```

It returns `|T|`, `|T+T|`, `|T*T|`, the additive closure, and whether that closure is a proper subfield. The subfield check verifies closure under multiplication and the presence of 1, after additive closure; a finite subring of a field containing 1 is a field. Threshold 3/5 is an experimental parameter, selected for this exploratory iteration. There is no held-out classification claim.

## Iteration 1: a scalar score failed

The pilot in [pilot.json](pilot.json) used two primes and two square fields. At p=101, choosing the single numerically returned extremal eigenvector yields maximum translation overlap 0.295. Optimizing over the full near-edge subspace yields 0.797. Thus a single vector hides a substantial arithmetic direction. This also falsifies a naive scalar-score separation: the square control q=121 gives 0.801, almost the same value.

A length-31 interval in F_1009 has translation overlap 30/31 = 0.968, larger than these values, while its largest absolute H eigenvalue is only about 0.576. Translation structure alone does not imply a dangerous Paley direction. It is a separate-set stress control, not a common-neighborhood input.

The response was to retain one actual optimized vector, inspect **all of its translations together**, and test algebraic growth of its large-overlap translations.

## Iteration 2: recover an obstruction, then reject an overstrong theorem

The full sweep uses six prime fields through 4001 and six quadratic extensions through 61². It applies the same spectral-tail window, rounding, translation threshold and seeded control convention throughout.

| q | Field type | Numerical ||Z|| | Max near-edge overlap | Same-spectrum relabelled overlap | Number of large translations of rounded witness | Additive closure size |
|---:|---|---:|---:|---:|---:|---:|
| 101 | prime | 0.903722 | 0.797467 | 0.374387 | 3 | 101 |
| 257 | prime | 0.978866 | 0.476989 | 0.285625 | 1 | 1 |
| 401 | prime | 1.015557 | 0.420572 | 0.236663 | 1 | 1 |
| 1009 | prime | 0.994074 | 0.358826 | 0.143673 | 1 | 1 |
| 2017 | prime | 1.003824 | 0.268725 | 0.098547 | 1 | 1 |
| 4001 | prime | 1.010507 | 0.239538 | 0.079821 | 1 | 1 |
| 121 | square | 1.096992 | 0.800998 | 0.257819 | 5 | 11 |
| 289 | square | 1.119957 | 0.790171 | 0.215568 | 17 | 17 |
| 529 | square | 1.127001 | 0.881473 | 0.175475 | 23 | 23 |
| 961 | square | 1.132877 | 0.922132 | 0.149013 | 31 | 31 |
| 1849 | square | 1.140730 | 0.929838 | 0.080597 | 43 | 43 |
| 3721 | square | 1.142622 | 0.956498 | 0.087290 | 61 | 61 |

All six square-field witnesses recover the proper subfield. Starting at q=17², T(z) itself is exactly that subfield. At q=11², five high-overlap translations generate it. This is a **finite experimental observation about the discovered witnesses**, not a proof that this always happens, and the subfield obstruction itself is classical.

The prime p=101 false positive is resolved structurally: its three high-overlap translations generate all of F_101. At the other five prime inputs only the zero translation crosses 3/5 for the selected rounded witness.

There is a useful further failure. At prime p=401 the tool exports a mean-zero integer vector with

```
n = z^T z = 21,311,298
sum z_i = 0
A = z^T S_C z = -370,800,106
4 A^2 - 3 q n^2 = 3,602,653,237,345,732 > 0.
```

Therefore `|z^T H z|/n > sqrt(3)/2` exactly. Its only translation crossing 3/5 is zero. This falsifies the finite assertion **“every outlying witness has a nontrivial large-overlap translation.”** It does not exclude the existence of a different structured outlying witness, and it does not challenge the asymptotic Paley target. The latter allows finite outliers.

The selected translation detector therefore recognizes the persistent square-field mechanism, but is **not** a universal explanation for every finite outlier. An asymptotic inverse principle must distinguish persistent excess from isolated finite outliers; this experiment does not determine a sufficient tolerance and richer arithmetic features may be needed. This is a sharper next research question than simply asking for another moment bound.

## Controls and certification

**Exact spectrum control.** Permute coordinates inside C and fix every other field element, including the two anchors. Conjugating the full sign matrix preserves symmetry, the zero diagonal, every off-diagonal sign, `S²=qI-J`, the anchor columns, and the full and compressed spectra. Thus it preserves every spectral moment at every depth, not merely the three displayed moments. The additive-coordinate interpretation is changed. The code checks the algebraic invariants at small fields; their all-size preservation follows immediately from permutation conjugacy and the constant anchor columns on C. This is a relabelled copy of the graph, not a new graph obstruction.

The relabelled numerical overlaps above use one seeded permutation per input. They demonstrate that the feature is not determined by spectrum; they are not a null-distribution estimate or a statistical significance claim.

**Basis control.** The same subspace is tested after a seeded random orthogonal basis change. Agreement is checked within 1e-10. Floating-point eigenspace discovery is heuristic, particularly near a cutoff; it supplies candidates and approximate extrema, not a certified global upper bound.

**Exact replay.** Every rounded witness is serialized with C and its integer vector. [verify_witnesses.py](verify_witnesses.py) is an independent standard-library implementation; it imports no discovery code, NumPy or SciPy. It recomputes the field characters, every weighted difference, every translation overlap, and the signed quadratic form using Python integers. [validation.json](validation.json) records 12 passing witnesses and 2,484,083 weighted difference pairs. Exactness attaches to those witnesses and signatures, not to the floating-point maximization or asymptotics.

For a mean-zero witness the normalized-edge comparison reduces to the integer inequality above. For general witnesses the certificate records the exact algebraic expression `A/(sqrt(q)n) - (sum z)^2/(qn)`; over square fields it additionally records a rational quotient. A `false` outlier flag in validation means the particular mean-zero test did not certify it, not that the input has no outlier.

## What would make this useful at full scale?

The hypothesis worth testing is a **quantitative inverse principle for persistent excess edge**: if a sequence of actual Paley localizations has an edge excess bounded away from zero, some dangerous direction should carry a specified arithmetic signature. The present detector is only a candidate feature. No implication of that form is proved, and the exact p=401 failure prevents upgrading the current finite observations directly into a theorem.

A next prototype should optimize several translations jointly under a one-sided Rayleigh constraint, use separate positive/negative tails, and vary both the tail width and the overlap threshold on newly chosen fields and anchor depths. It should also measure multiplicative transport and sparse rational maps, with controls that preserve whichever existing features have already been imposed. The user-facing product of such work is an auditable obligation ledger: a proposed implication, the exact inputs it uses, an automatically found countermodel if one exists, and a new feature suggested by the failure.

Dense diagonalization costs O(m³) time and O(m²) memory; translation compression adds about O(q m r²) for tail rank r. These prototypes reach m=999. To move to fields of order 10^5–10^6, use FFT multiplication by the ambient additive character kernel inside an iterative extremal eigensolver, then batched FFTs for all compressed translation correlations. Residual enclosures and threshold separation would be needed before a computed tail could be treated as a certified subspace. That scalable backend is not implemented here.

The Proximity Prize would require an additional, exact bridge from a relevant failure or rank witness to such a transport representation. No such bridge is assumed or proved in this lane. This package is mathematical tooling for the spectral bottleneck, not a solution to either full problem.

## Prior-art and novelty boundary

Primary sources inspected online on 2026-09-05:

- [Kunisky, *Spectral pseudorandomness and the road to improved clique number bounds for Paley graphs*](https://arxiv.org/abs/2303.16475). Supplies the relevant localized spectral target and distinguishes limiting spectral distributions from the conjectural edge.
- [Tao, *An uncertainty principle for cyclic groups of prime order*](https://arxiv.org/abs/math/0308286). Prime support uncertainty is established mathematics; this prototype claims no new uncertainty theorem.
- [Asgarli–Yip, *The subspace structure of maximum cliques in pseudo-Paley graphs from unions of cyclotomic classes*](https://arxiv.org/abs/2110.07176). Documents the established subfield structure in square-order Paley graphs; recovery here is diagnostic rediscovery of that control.
- [Cladek–Tao, *Additive energy of regular measures in one and higher dimensions, and the fractal uncertainty principle*](https://arxiv.org/abs/2012.02747). Connects additive energy, expansion and uncertainty under regularity assumptions. Those assumptions are not established for the optimized witnesses here.
- [Cabrelli–Mosquera, *Subspaces with extra invariance nearest to observed data*](https://arxiv.org/abs/1501.03187). Shows that optimization of translation-invariant approximating subspaces is established; our finite spectral conditioning and matched-spectrum controls differ in task, not in the invention of subspace optimization.
- [Dinin–Lind, *On the eigenvalues of cyclic covers of Paley graphs*](https://arxiv.org/abs/2601.11877). Studies translation-invariant covers and prime/extension-field spectral distinctions. This is related prior art, not a result implying our diagnostic or the desired edge estimate.
- [Liu et al., CRAFT official implementation](https://github.com/yiruiliu/Craft_CounterExample). Counterexample-guided checking and repair of mathematical claims is already an explicit workflow. This project does not claim to invent that workflow.

The search included combinations of Paley eigenvectors, additive translations, spectral projectors, localization, subfields, and counterexample-guided mathematics. No source inspected described this exact coupled experiment. A bounded search cannot certify that it has never appeared, and unindexed or differently named work may overlap. The defensible status is **a new prototype in this repository with a specific experimental insight and exact finite certificates; external mathematical novelty remains unestablished**.

## Reproduce

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python3 tooling_lab/spectral/transport_microscope.py
OPENBLAS_NUM_THREADS=1 python3 tooling_lab/spectral/refine_transport.py
python3 tooling_lab/spectral/verify_witnesses.py
```

Discovery ran with NumPy 2.4.4 and SciPy 1.17.1. The fixed seed is 20260905. The first command also supports `--quick --output pilot.json`; it overwrites witness files for its four input fields, which are the same deterministic cases in the full run. Exact replay requires the full `results.json` and `refinement.json`. Source SHA-256 values are saved in each result file. No network, API key, LLM, or prize service is required to run the prototypes.

## Iteration 3: threshold, window and coordinate robustness

[robustness.py](robustness.py) reads the saved integer witnesses without replacing the pilot or main results. [robustness.json](robustness.json) records a post-hoc sensitivity audit at thresholds `1/3, 1/2, 3/5, 2/3, 3/4`.

For every one of the six square fields, the high-overlap translations generate the same proper subfield at every tested threshold. At q≥23², the high-overlap set itself is the entire subfield at all five thresholds. At q=17² it has 17, 17, 17, 13 and 3 elements respectively, always generating the same 17-element subfield. The prime p=101 vector retains nonzero large translations even at 3/4; their additive closure is the whole field. At p=257, 401 and 1009, nonzero large translations appear at 1/3 but disappear at 1/2. Thus the absence of individual translations is threshold-sensitive even where the recovered closure is stable.

This is not a learned test of whether a supplied field is prime: that information is already part of the input. The informative observation is that the discovered **spectral witness** has a large-overlap translation set confined to a particular small additive subgroup. In a quadratic extension every nontrivial proper additive subgroup is a one-dimensional line over the prime field. Hence recognizing a normalized subfield in this specific experiment is weaker than a general approximate-subfield inverse theorem.

Tail-window sensitivity is also material:

| q | Window 0.01: rank / max overlap | Window 0.025: rank / max overlap | Window 0.05: rank / max overlap |
|---:|---:|---:|---:|
| 101 | 2 / 0.429313 | 3 / 0.797467 | 3 / 0.797467 |
| 401 | 2 / 0.302851 | 3 / 0.420572 | 5 / 0.422687 |
| 289 | 1 / 0.790171 | 1 / 0.790171 | 1 / 0.790171 |
| 961 | 1 / 0.922132 | 1 / 0.922132 | 1 / 0.922132 |

The prime near-edge diagnostic can change when a nearby direction enters the tail. The separated square-field top direction is stable over the tested windows. These four sensitivity cases were selected after the main sweep, so they are not held-out confirmation.

A coordinate control found a correctable weakness in the literal “proper subfield recovered” flag. Under a residue affine change of coordinates `x -> a x+b`, move the two anchors to `b,b+a` and keep the witness values attached to their corresponding vertices. Every character entry and Rayleigh quotient is unchanged, while

```
a_new(a t) = a_old(t),  T_new = a T_old.
```

A nontrivial scalar multiple of the base subfield need not itself contain 1 or be a subfield. For example, the q=289 control uses a=18, b=3. The old literal flag switches from true to false even though the example is exactly the same arithmetic configuration. Dividing the recovered additive closure by the new anchor difference a recovers the old subfield and restores covariance.

The audit checks this identity at every field element and checks every character entry for all 12 inputs. Its repaired field-geometric signature is **the additive closure normalized by the anchor difference**, not the literal presence of 1 in an arbitrary coordinate frame. All six square cases preserve that signature; prime cases remain without such a proper recovered closure. This correction matters for applying the tool to arbitrary anchor edges, and it does not add a uniform inverse estimate.

Run the additional audit with:

```sh
OPENBLAS_NUM_THREADS=1 python3 tooling_lab/spectral/robustness.py
```
