# Discovering the missing polynomial-piece witnesses

This prototype closes one observed interface gap: it discovers polynomial pieces from an unlabeled direction word, then sends the resulting exact partition to the frozen round4 scalar-fiber verifier. The n1024 interleaved fixtures now work without supplying their hidden labels or coefficients. The method fails on all 16 saved hard coset stacks, and does not resolve the official Proximity Prize event.

Polynomial reconstruction, finite-field factorization, Hensel lifting, polynomial root counting, and affine linear-code symmetry are established mathematics. The contribution here is the checked interface, a minimum-cover certificate, and an explicit map of where this restricted discovery method succeeds or fails. No historical novelty is claimed.

## Files and reproduction

- `piece_discovery.py`: monic relation solver, exact linear certificates, bounded branch recovery, frozen certificate integration, and export verifier.
- `run_experiments.py`: hidden-label fixtures, exhaustive tiny codeword/scalar oracles, baselines, failures, and scale cases.
- `results.json`: full certificates for inherited, generic, hard-coset and large inputs; seeded summaries for additional random hidden pieces.
- `verify_results.py`, `verification.json`: source-bound readback. This checker shares production algebra helpers; it is not an independently implemented decoder.
- `run.log`: timings and outcomes from the saved run. Independent review is maintained separately in `../review/`.

Run from the repository root, using the installed NumPy runtime:

```sh
/opt/miniconda3/bin/python tooling_lab/round5/piece_discovery/run_experiments.py
/opt/miniconda3/bin/python tooling_lab/round5/piece_discovery/verify_results.py
```

The public entry point is `discover(F, domain, k, s, word, u0=None, max_pieces=3)`. It receives no planted labels or coefficients. Coordinates must be distinct field elements and `1 <= k <= s <= n`. `word` is the direction `u1` when `u0` is supplied. All earlier rounds and the Proximity Prize repository remain unchanged.

## The exact mechanism

For a proposed piece count r, solve the affine system for

\[
P(X,Y)=Y^r+\sum_{j=0}^{r-1} A_j(X)Y^j,
\qquad \deg A_j\le(r-j)(k-1),
\qquad P(x_i,w_i)=0.
\]

There are `U=(k-1)r(r+1)/2+r` unknown coefficients. This is a restrictive monic ansatz, with ordinary value constraints and no multiplicities. It targets a cover of **every input point** by a small number of degree<k polynomials; it is not a general list decoder.

For each system the export contains selected original rows and pivot columns defining a nonsingular minor, together with a nullspace basis normalized on the free columns. The verifier checks that minor and all kernel vectors, establishing exact rank and nullity. A consistent system also has a checked particular solution. An inconsistent system has a short dual witness expressing one original row as a combination of the selected rows while its right-hand side differs. This proves that the specified monic relation cannot exist.

One deterministic solution is selected by setting free coefficients to zero. The code does **not** search the exponentially many kernel combinations. A positive-dimensional solution space can contain a useful factorable relation even when this selected solution does not; that outcome is explicitly incomplete.

At a bounded sequence of centers a, evaluate `P(a,Y)`, check squarefreeness and complete splitting with `Y^q-Y`, and find its simple roots using bounded Cantor–Zassenhaus splitting. This uses polynomial modular powering, not a scan of all q field elements. Each simple root is lifted coefficient by coefficient in `X-a` to order k. A lifted jet is accepted only if the resulting degree<k polynomial h satisfies the **full identity** `P(X,h(X))=0`. The product of all recovered `Y-h(X)` factors is then compared to P. A finite jet alone never establishes a polynomial branch.

Every recovered polynomial is evaluated on every input coordinate. If these maximal agreement supports cover all coordinates, the algorithm deterministically assigns each coordinate to one matching polynomial and submits the resulting partition to `certify_static_word`. The frozen root-count proof bounds an unlisted polynomial's agreement by

\[
\sum_j \min(|R_j|,k-1).
\]

Only a strict threshold `s` above that cap yields a complete list. If an exact code anchor `u0+z*u1=f` is found, the frozen scalar transport additionally certifies all codewords and maximal supports at every scalar symbolically. These are ordinary agreement lists; MCA nonjointness and the official prize event are separate requirements.

## What the algebra certificates tell us

If a word has an r-piece polynomial cover, the product of the r linear factors is a relation in this ansatz. Therefore inconsistent systems for all smaller r, together with a discovered r-piece cover, certify the **minimum number of degree<k polynomials needed to cover this finite input**. This is elementary polynomial reconstruction used as a diagnostic, not a new theorem. The large fixtures have certified minimum cover sizes 3 and 2. The generic n64 controls have cover size at least 4.

There is also a sufficient identifiability condition. Let `D=r(k-1)`. If there are r distinct covering polynomials and each has more than D matching coordinates, substitution gives more roots than the degree of `P(X,h(X))`, so it vanishes identically. Over `F(X)[Y]`, the r distinct monic linear factors then exhaust the degree-r polynomial, forcing their product. Thus the monic relation is unique. The measured ranks independently confirm uniqueness in the successful large cases:381/381 for three pieces,191/191 for two.

The coordinates' order is absent from this argument. It explains why interleaving destroyed the earlier coordinate-block batching advantage while leaving this reconstruction method effective. The sufficient support condition is not necessary; some smaller examples succeed without it.

## Results

All large cases use F65537, n1024, k64, s410, and the previously saved multiplicative subgroup domain. Their planting metadata is used only after the blind discovery call for comparison.

| Input | Discovery | Exact evidence | Complete scalar nodes |
|---|---|---|---:|
| Interleaved three pieces | Complete, minimum 3 | Rank 381/381; supports 342, 341, 341; cap 189 | 1 |
| Interleaved two pieces | Complete, minimum 2 | Rank 191/191; supports 512, 512; cap 126 | 131073 |
| Aligned/interleaved n48,k6,s28 | Both complete, minimum 3 | Rank 33/33; full frozen certificates | 1 each |
| Aligned n128,k8,s80 | Complete, minimum 3 | Rank 45/45 | 1 |
| Eight hidden random coefficient/assignment n192,k8 fixtures | All complete, minimum 3 | Exact planted coefficient sets recovered after blind run | 1 each |
| Four generic n64,k4,s20 words | No relation for r1,2,3 | Exact dual inconsistency witnesses | Unresolved |
| All 16 saved hard coset stacks | Discovery incomplete | No code anchor; r2 kernels of dimension 7 or 9 | Unresolved |

The n1024 runs take a few seconds each on this machine; exact timings are in `results.json`. These include relation construction, rank/kernel checks, factor recovery, anchor discovery and frozen certificate verification. The full suite takes roughly half a minute. They are measurements of this Python/NumPy implementation, not algorithmic complexity claims.

The baseline partitions consecutive coordinates into groups of at most k and interpolates each group. It is always a valid cover and is much faster. On both large interleaved cases it gives 16 pieces in about 0.27s, with cap 1008, so it cannot certify at s410. Recovered covers have caps 189 and 126. On the aligned n128 case the baseline already certifies the list with 11 pieces and cap 77<s80; reconstruction improves the cover description there but is unnecessary for certification.

There are 177 tiny complete-oracle inputs over F3, F5 and GF9. The 161 successful discovery certificates match every degree<k codeword and every scalar exactly, totaling 1091 nodes ; 143 certificates include every field scalar. The 16 other GF9 inputs are generic controls: ten recover a cover whose root-count cap is too high, and six remain incomplete. All 32 hidden two-piece GF9 fixtures succeed. Three corrupted full exports are rejected. Counts of actual full-record readbacks are reported separately in `verification.json`.

## Failure-driven iteration and remaining limits

A deliberately colliding first center initially prevents branch recovery for `(Y)(Y-X)`. Expanding the center budget recovers both branches and certifies the list. The final bounded center sampler selects distinct centers without an unbounded repeated-draw loop.

A larger budget is not a universal remedy. For example, over F3 the three distinct branches in `(Y-X)(Y+X)(Y-1)` collide at every base-field center. The simple-center method can fail even when every polynomial factor exists. Extension-center lifting with descent or a complete Roth–Ruckenstein routine would be separate, known-method extensions; they are not implemented here. Even characteristic is explicitly unsupported by this branch splitter.

The saved hard n16,k8 directions show the more relevant obstruction: the r2 system has 23 unknowns but only 16 equations and a seven-dimensional kernel. A relation exists, but the selected relation's Hensel branches are not the needed polynomial cover. This is a discovery failure, not evidence that no cover exists. Indeed, any n16 word can be covered by two polynomials by interpolating two groups of eight points. That elementary cover generally has cap 14, too large for s11; all 16 saved pencils additionally lack the anchor needed by the frozen scalar reduction. Random kernel search was deliberately not added.

A dense solve costs O(nU²) field operations; checking the selected minor adds O(U³), with O(nU+U²) storage including kernel certificates. The simple Hensel implementation has a coarse O(r²k³) operation bound. At fixed rate and fixed r, the dense algebra grows cubically with n. The generic field interface is tested over GF9, but the inherited extension-field constructor precomputes every element and is only appropriate for tiny fields. No scalable F_(p^6) backend or official-size event certificate is supplied.

The next useful structural question is whether underdetermined relation spaces contain factorable low-piece covers with sufficiently unbalanced support sizes to make the list cap strict, and whether those can be found with polynomial work. This prototype identifies that missing step and exports the exact relation spaces; it does not solve it.

## Prior-art boundary

The local repository already implements multiplicity-two Guruswami–Sudan interpolation followed by Roth–Ruckenstein polynomial-root search in `scripts/probes/probe_fsmf_predecessor_miniature_census.py` (functions at lines 181 and 201). That implementation scans field roots inside its recursion and verifies candidates by exact agreement. The current prototype uses a narrower monic complete-cover ansatz, bounded simple-root splitting/lifting, explicit algebra witnesses and the frozen scalar certificate interface. It is not a replacement for that decoder.

The local `docs/kb/hab25-2025-2110-MCA-for-RS-extracted.md`, lines 25–31, already discusses Guruswami–Sudan over a function field, discriminant-good centers, Hensel lifting and factor refinement. Existing Lean files also formalize Hensel branches. The audit therefore rules out presenting these foundations as newly invented.

Primary sources supporting the established foundations:

- [Roth and Ruckenstein, Efficient decoding of Reed–Solomon codes beyond half the minimum distance (2000)](https://cris.technion.ac.il/en/publications/efficient-decoding-of-reed-solomon-codes-beyond-half-the-minimum-/): polynomial reconstruction with low-degree univariate root finders.
- [Koetter, Ma and Vardy, The Re-Encoding Transformation in Algebraic List-Decoding of Reed-Solomon Codes (2010)](https://arxiv.org/abs/1005.5734): bivariate interpolation and factorization as established decoding stages, including their computational reduction.
- [Shoup, A Computational Introduction to Number Theory and Algebra, §20.4](https://www.shoup.net/ntb/ntb-v2_5.pdf): Cantor–Zassenhaus finite-field factorization, including its squarefree input requirement.

This bounded audit found an existing-method integration to build and test. It cannot establish that a mathematical idea has provably never been tried.
