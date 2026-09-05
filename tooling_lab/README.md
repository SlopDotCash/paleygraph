# Mathematical tooling lab: Paley and proximity

**Research and executable prototypes, September 5, 2026.** This project investigates missing mathematical information and how to detect it. It attempts neither prize proof. The first three tool families below were implemented, tested against controls, and revised after their first interpretations failed.

**Second iteration: [exact energies and agreement censuses](round2/README.md).** It adds exact all-input harmonic energies beyond field enumeration, a complete compressed agreement census, domain-preserving arithmetic lifts, and exact million-prime spectral witnesses. The new critical-size tests weaken the toy triangle's predictive interpretation; the exact toy separation remains valid.

**Further iteration: [counterfactual tests](round3/README.md).** An exact Paley/Peisert comparison demonstrates information lost by the new global-energy compiler; a phase-sensitive prototype rejects generic summaries after stronger symmetry and chirp controls. These sharpen what the next tools must retain.

**Fourth iteration: [marked contractions and certified batches](round4/README.md).** Exact arithmetic conditioning isolates one missing contraction, verifies actual equal-cell witnesses, and scales its computation to a million-point field. A norm envelope also exposes its small effect on the tested conditional averages. Coding experiments certify large interpolation families and then refine an interleaved failure with complete scalar-fiber certificates.

**Newest completed iteration: [seven trace statistics and discovered polynomial pieces](round5/README.md).** The generic global third-moment compiler now works beyond six-column sets and reduces to seven measured statistics and15 coefficient evaluations. Independent complete censuses reach450,978,066 eight-sets per order49 graph. Blind polynomial reconstruction supplies previously missing coding witnesses at n1024. The final integration audit pins113 artifacts and407 source/input bindings, with round4 preserved.

**Current work: [information limits and unanchored pencils](round6/README.md).** Exact rational envelopes quantify what ordinary trace moments discard and how small that effect is at the tested critical sizes. A new symbolic certificate handles structured pencils with no exact codeword anchor. Independent reviews of both tools passed; further work examines a universal certificate-family obstruction and overlapping support alternatives. The [live iteration plan](NEXT_ITERATION.md) supersedes the early priorities below.

The strongest result is a practical change in research method: ask **what information a proposed argument has discarded**, produce an exact witness of that loss on actual arithmetic inputs, and refine the representation before attempting an estimate. That general strategy has prior art. The potentially original work here is its precise arithmetic implementation and the resulting finite witnesses; historical novelty is not established.

![Results from the three prototype families](overview.png)

## Start here

| Tool | What it now does | Concrete result | Detail |
|---|---|---|---|
| Exact observable and realization auditor | Finds actual character sets with equal selected measurements and different target moments; extracts integer moment-preserving moves; adds incidence information. | At p=61, flat measurements leave 63 ambiguous fibers. Rooted quartic incidence leaves one. A triangle contraction of the already computed quartics resolves the remaining finite ambiguity. | [Observables](observables/README.md) |
| Spectral witness transport microscope | Searches the complete near-edge subspace, attaches an arithmetic translation signature, and exports integer witnesses. | All six square-field controls recover their subfield mechanism. An exact prime-field witness rejects an overstrong pointwise explanation of outliers. Affine controls expose and repair a coordinate dependence. | [Spectral](spectral/README.md) |
| Syzygy and realized-stack microscope | Tracks dependencies through root-motion orders, separates admissible subgroup edits, and retains actual decoded-codeword/track overlap. | First-order preservation hides second-order failure in the F41 example. Equal coset diagnostics coexist with different realized scalar counts. Extension-field checks find a bad scalar outside the prime subfield. | [Proximity](proximity/README.md) |

The [prior-art audit](novelty/README.md) and [primary-source/query ledger](novelty/literature_ledger.md) distinguish known mechanisms, local earlier attempts, and proposed combinations. The research question is broader than finding a new name for an open bound.

## What the existing record actually supports

The live Paley frontier has four distinct gaps: uniform signed moments; exceptional subgroup relations and energies; the full localized spectral edge including invisible sectors and borders; and scalar code maxima plus the separate official spot-check obligation. Previous local passes include substantial identities and finite obstructions. Their unresolved analytic inputs are not supplied by those identities.

The proximity audit uncovered a material correction to older proposal documents: the inspected record does not establish that *every* prize solution must imply the proposed Paley-style Fourier estimate. The current [necessity audit](/Users/shawwalters/proximityprize/docs/kb/deltastar-sw1-nec-2026-09-05.md) distinguishes conditional Fourier-to-incidence routes from an actual necessity theorem. Thus a direct coding route remains a separate possibility. This lab does not independently rebuild its referenced Lean cone.

Earlier local work already proposed cumulants, R-transforms, antipodal matroids, Koszul complexes, phase operators, gauge fixing and conductor gradings. Several had later negative audits. Generic use of those subjects would not meet the requested novelty standard. Nor is the classical toolbox proved exhausted.

## The most promising new research object

For a class of actual arithmetic inputs X, a target T and a measurement map f, define the finite realization envelope

```
R_f(u) = { T(x) : x in X and f(x)=u }.
```

Compare it with an explicitly larger, algebraically tractable relaxation. The **realization gap** measures what the relaxation permits that actual field inputs cannot achieve. A candidate tool should preserve enough arithmetic incidence to reduce that gap without reproducing the entire input or disguising T as a feature.

In the six-column experiment, the relaxed histogram fiber has one primitive direction `(-10,15,-6,1)` on absolute row sums `{0,2,4,6}`. Each step preserves moments through four and changes M6 by exactly 23,040. The exact p=61 census shows that 205 of 257 lower/boundary fibers have a smaller realizable upper endpoint than their nonnegative histogram relaxation. This supplies a concrete laboratory for learning what arithmetic realizability imposes beyond moment identities.

The important follow-up came from a failure: an unordered list of all quartic correlations still loses information. Attaching those same correlations to their omitted column pairs gives a weighted graph. Its incidence signatures and triangle contraction distinguish the remaining toy examples **without evaluating a new character sum**. This is the clearest actionable insight of this pass: useful information can be absent from the aggregation rather than absent from the underlying calculations.

This is a proposal for a research object and tool workflow, not a new general theorem. Integer fibers/Markov moves, separating invariants and graph contractions are known. The exact arithmetic realization tests and selected contraction are candidates for a specific contribution. A bounded literature search cannot prove that an approach has never been invented or privately tried.

## Experiments changed the designs

1. **Observable auditor, four iterations.** Start with low moments; add boundary and energy information; preserve the quartic deck; then preserve its incidence and triangle contraction. A deliberate target-leakage control is excluded from the candidate list. Across the initial suites, 569,934 normalized inputs were processed, of which 556,622 were exhaustive small-field inputs. The strong collisions at p=61 and p=257 are outside `p >= n^4`. Critical-scale samples at p=1297 and p=4099 have singleton rich-feature fibers, so their lack of collisions supplies no sufficiency evidence.
2. **Spectral microscope, three iterations.** A single eigenvector hid a high-translation direction, so the tool optimized the whole near-edge subspace. A scalar overlap score gave a prime-field false positive, so it retained joint translation structure. Threshold/window tests and affine coordinate controls then showed that the recovered additive subgroup must be normalized by the anchor difference. This was tested through p=4001 and square fields through 61². No inverse theorem for persistent spectral excess follows.
3. **Proximity microscope, three iterations.** A 15-dimensional ambient tangent space did not preserve the F41 syzygy through second order along the 24 tested straight paths. Higher jets and separate on-domain edits replaced a tangent-only diagnostic. The next experiment showed that even persistent triple identities and their edit signatures do not determine realized pencil incidence, so the tool retained codewords, maximal supports, affine tracks and overlaps. A third prototype introduced actual extension-field arithmetic and independent tiny censuses over F9 and F729.

The coding example with 20 three-point tracks but only 15 distinct covered nodes makes the loss concrete: summing the 60 memberships counts 45 repeats, and one bad scalar is outside those tracks. Counting tracks alone would discard precisely the overlap needed by an incidence argument.

## What would make these tools useful for the full problems?

| Direction | Next implementable capability | Evidence needed to continue | Remaining mathematical obligation |
|---|---|---|---|
| Character moments | Search within low-cost incidence fibers at `n approximately p^(1/4)`; retain primitive histogram moves and actual-set witnesses. | A separator that survives new primes and structured families without making nearly every fiber a singleton. | A uniform bound on the chosen contraction or a uniform characterization of realizable moment moves. |
| Spectral edge | Optimize joint additive/multiplicative transport under a one-sided Rayleigh constraint; use FFT matrix-vector products and certified residuals. | Stable arithmetic signatures for persistent excess, with matched-spectrum and affine controls. | A quantitative inverse principle or direct uniform edge estimate; the finite p=401 witness forbids an overstrong all-witness formulation. |
| Code incidence | Generate compressed agreement-support families while preserving decoded codewords, affine tracks, duplicate multiplicities and overlap. | Exact equality with full subset census on small prime and extension fields; a substantial size reduction on harder examples. | A uniform support/overlap bound at the actual code length, rate, field and soundness budget, plus all separate prize obligations. |

The main computational barriers are explicit: all-quartic enumeration, dense diagonalization and all-agreement-subset enumeration. The extension adapter removes a field-type mismatch; F729 does not approach the official field size or code length. None of the current prototypes is a production-scale prize certificate.

The next round should prioritize the **compressed agreement-support generator** and **arithmetic incidence-fiber search**. Both preserve information that the new experiments showed was being lost. Further numerical scans are worthwhile only if they test a stated mechanism or compression, not merely enlarge a table.

## Verification and reproduction

Use `/opt/miniconda3/bin/python3` on this machine, or a Python environment with [requirements.txt](requirements.txt). The proximity and independent certificate verifiers use the standard library.

```sh
python3 run_lab.py verify
python3 run_lab.py full
python3 plot_results.py
```

`verify` checks saved evidence and regenerates the tiny extension-field cross-check; `full` regenerates experiments and then verifies them. The full run writes only lab results. The small exact replay is separate from numerical discovery. Timing values vary; seeds and mathematical inputs are fixed. Numerical eigenspaces and rounded discovery vectors can vary across NumPy/SciPy/BLAS versions, so the saved integer certificates are the portable evidence.

Verification includes an independent standard-library replay of all 12 spectral witnesses (2,484,083 weighted difference pairs), all 60 first-stage moment-witness endpoints, and an independent complete p=29 realization census. Coding review independently enumerated 256 jet systems and all 729 tiny Reed–Solomon pencils, then checked all 16 saved stacks, 95 decoded codewords and 383 affine tracks. See [observable review](novelty/observable_review.md) and [proximity review](novelty/proximity_review.md). These are independent agent/code checks, not human refereeing or new Lean certification.

The original research and proximity checkout were used as read-only inputs. Source provenance, result hashes and runtime information are collected in the lab manifest. The standalone [figure](overview.png) and [PDF](overview.pdf) are generated directly from saved results.
