# Next: direct inverse-coordinate cuts and orbit-realizability certificates

The previous goal turn made concrete progress. Round23 completed and independently verified compressed erasure covers, resumable refinement, a scalar-orbit intersection counter, and an integral realizability gap. The full user goal remains active: investigate missing mathematical tools beneath both Paley and proximity, implement and scale prototypes, retain failures, and keep iterating. Neither prize is proved. Worldwide historical novelty cannot be established by negative searches.

Rounds6–23 are frozen. Create round24 next. Both proof projects remain read-only; preserve the dirty checkout. No Lean builds, proof edits, commits, PRs, submissions or third-party messages. No numerical process is left running from round23. Its producers and reviews all have terminal exit0. Use the final manifest's exact artifact/binding totals.

## New finite results

On p=2013265921, N=64, g=397765732 with
f=[-1,0,-1,0,1,1,0,0,2,1,0,-1,-2,-3]+[0]*50,
any36 consecutive cyclic erasures leave at most one completion from the remaining28 digits. Round22 proved this for32 erasures. The36-erasure decoder has a27648-combination cap; uniqueness does not imply a one-candidate construction algorithm.

The complete36-erasure cover has1867 nodes,64 rational cuts and7953 exact narrowing steps, excluding3690562 nonzero visible difference representatives up to sign. At32 the new cover needs131 nodes and21 cuts for9841 representatives. At33 it needs1337 nodes and64 cuts for1328602 representatives. An initially unresolved33-erasure singleton is eliminated by reapplying later cuts, with no additional cuts or nodes. All1336 earlier finished nodes are preserved.

At40 erasures the cover is incomplete:2021 nodes,64 cuts,19 unresolved leaves containing177723199 bounded integer representatives. Pure resaturation changes none of those leaves. A rational continuous point survives.

The surviving target t=(1,0,...,0) has scalar difference345549834. An exact scalar-orbit counter proves there are zero actual pairs with that target. The tightest fixed coordinate62 has interval[978728438,1006632960], so two separate implementations exhaust27904523 candidates instead of scanning all p scalars. The first filters centered intervals; the reviewer enumerates the other scalar in the pair and compares actual centered differences. Thirty-seven coordinates are fixed and27 carry coordinates are free. This excludes one target up to sign, not all remaining40-erasure differences.

An explicit [integral gap certificate](round23/gap_certificate.json) strengthens the diagnosis: it supplies H integral with |H_i|<=p-1 and H_i=345549834*g^i modulo p. D=fH/p is integral, supported on the erased40 coordinates, has max coefficient4, and has integer coordinates in the full erasure lattice. Nevertheless no actual pair has its visible target. Restoring invisible integrality alone therefore does not remove this candidate; joint scalar-orbit realizability is missing.

Read [round23 README](round23/README.md), [derivation](round23/DERIVATION.md), [cover review](round23/cover_review.json), [orbit review](round23/orbit_review.json) and [gap review](round23/gap_review.json). The [36-erasure input](round23/erasure36_input.json), [output](round23/erasure36_output.json) and [query](round23/cover_query.py) recover scalar1234567 after a cyclic block starting at57, checking384 candidates.

## Next concrete tool: direct physical-coordinate cuts

A useful structural observation emerged from exact basis checks. Let P_j=(A_f/k)*H_j, where H_j is lattice basis row L_j embedded on the erased coordinates. For a full integer difference z, its inverse lift is

H=sum_j z_j*P_j.

The [basis information artifact](round23/lattice_information.json) and [independent matrix review](round23/lattice_information_review.json) verify P_j integrality and scalar congruences in five large geometries. Invisible directions have scalar step0 and lift to p-multiple carry vectors.

For consecutive masks of d=33,36,40, the invisible P_j are exactly signed p*e_i at i=0,...,d-14. The final binding verifier additionally checks these supports. For alternating32 there are no invisible rows. The seeded scattered32 basis has one invisible row supported on9 coordinates, with integer carry coefficients.

This suggests an exact compiler that need not discover every joint cut through LP:

1. For every physical coordinate i on which all invisible P_j vanish, compile the direct integer inequality

   |sum_(j visible) P_j[i]*t_j| <= p-1.

   Its certificate is simply the checked inverse-basis identity and the vanishing invisible coefficients. No targeted normalization or L1-optimality certificate is needed. It is valid for every actual centered difference.

2. For consecutive masks, these inequalities are the complete continuous projected body: the invisible directions are independent coordinate units, so free coordinates can be cancelled continuously. At d36 there are41 physical-coordinate cuts; at d40 there are37. At alternating32, all64 physical coordinates give direct cuts because no invisible directions exist.

3. For the scattered mask's one invisible carry row, the9 affected physical bounds impose intervals on one real carry variable. Exact pairwise interval compatibility can eliminate it, giving finitely many additional inequalities in visible t. This is ordinary Fourier–Motzkin elimination; audit its established status rather than claiming the method is new.

4. Feed these direct cuts into an exact integer-box cover, starting from the same certified individual bounds. Compare32/33/36/40 consecutive cases and the previously failed alternating/scattered32 masks. Preserve explicit cut/node budgets and unresolved regions. Do not enumerate their astronomical Cartesian boxes.

5. Avoid a type trap: round23 propagate divides coefficients using /. Pass Fraction coefficients and radius, even for integer physical cuts, so Python does not silently produce floating-point endpoints.

6. The round23 cut schema assumes a targeted norm-residual certificate. Add a new physical-coordinate certificate type in round24 and independently validate its matrix identity. Do not edit frozen round23 files to make the format fit. Exact trace and split verification must continue to prove complete coverage.

7. Once the continuous body admits a nonzero integer target, stronger consequences of the same body cannot exclude it. Use the actual orbit-realizability tool for that target, or develop certificates covering a family of such targets. The saved integral gap is the regression control.

These direct physical cuts have not yet been implemented. The observation is derived from saved, verified inverse-basis rows, not a general theorem for arbitrary relations. Degenerate p-dividing controls show why integrality and scalar congruences must be checked rather than assumed.

## Actual orbit intersection: reusable interface

For the consecutive40 basis, free inverse coordinates J={0,...,26}; I={27,...,63}. A visible target t determines fixed H(t)[i] on I and scalar delta modulo p. An actual pair with that target exists iff some a satisfies

F_a[i]-F_(a-delta)[i]=H(t)[i] for every i in I.

Equivalently F_a[i] lies in[-m,m] intersect[H(t)[i]-m,H(t)[i]+m]. Each interval has p-|H(t)[i]| integer values. The tightest interval is a bijective scalar parameter, so exhaustive filtered enumeration gives a complete pair count.

Sufficiency uses the checked invisible p*e_i basis: any actual pair meeting the fixed coordinates differs from H(t) only by integer p-multiples on J, hence is in the same visible difference class and agrees on the known digits. This is stronger than a necessary-box test.

[orbit_intersection.py](round23/orbit_intersection.py) implements complete counting with a capped witness list and explicit truncation. It currently requires at least one fixed coordinate and uses exact int64 arithmetic after an overflow check. It is an internal arithmetic interface, not a general cryptographic solver. The saved large counter runs a complete27904523-value seed interval; no time-bounded partial enumeration was called empty.

Future approaches can group targets sharing interval endpoints, compile modular interval certificates, or use lattice/hidden-number methods with exact verification. Audit against Boneh–Venkatesan1996 and known integer-programming certificate systems. A direct physical-cut compiler can improve covering costs but cannot itself resolve the saved integral gap.

## Verified coverage and controls

Round23 initial covers:5748 nodes,298 cuts,24238 exact narrowing steps across four large and seven small cases. All550 actual small ambiguity pairs survive unresolved leaves; four corrupt-cover controls are rejected. Additional refinement reviews check the33 complete and40 incomplete results, preserving previous finished nodes.

Preflight seed20260908 uses eight new scalars and eight bounded edits per mask:
-33:cap12288; all16 complete(8 singleton/8 empty),maxbox192,413 scalar candidates.
-36:cap27648; all16 complete(8 singleton/8 empty),maxbox3456,9976 candidates.
-40:cap708588;10 budget failures,6 complete(1 singleton/5 empty),maxbox139968,215424 candidates.
The separate review checks48 queries,109 cuts and109 optimality certificates,225813 scalar candidates.

Small orbit tests cover every nonzero delta and all two-lift assignments on masks{0},{0,2},{0,1,2,3} at p17,41,97:3344 intersections against246928 direct scalar-pair checks. The large orbit result is independently replayed. The integral gap has two corruption controls.

Inverse-basis review checks177 rows: all173 large rows are integral and scalar-congruent; scalar-p-squared erase1 has an integral but noncongruent row; erase3 has three nonintegral invisible rows. These controls prevent generalizing the large-case lattice properties without proof.

Eight33-erasure cyclic queries use36 scalar candidates total. CLI controls distinguish an incomplete cover, an incomplete enumeration despite certified uniqueness, a genuine17-completion small case and missing root coverage. A separate36-erasure CLI query recovers scalar1234567 at cyclic start57 using384 candidates.

The scientific figure was visually inspected. The manifest preserves round22(76), round21(87), round20(54), round19(28), round18(16), round17(22), round16(48). No jobs remain live; all tracked handles have terminal exit0.

## Broader scope retained

Round21's100000-prefix certificate is a finite natural-order state lower bound, not an arbitrary-algorithm or asymptotic lower bound. Round20's canonical diagram represents all6700417 words with36381 nodes and exact peak width2565. Compact completion representations still do not bound arithmetic norm distributions.

Round19's bounded-L1 relation obstruction still holds through budget11607941 on the large input. Current centered/orbit tools do not refute it. Round18's minimum-support words have widely varying norm defects.

Round16 proximity decoding is complete only inside a supplied affine space at n1024,k64,s96, with lists1,0,1 and revised queries71520,71520,73554. Ambient cluster discovery remains open. Return to that lane when a concrete transfer mechanism is available; these cyclotomic completion results do not automatically prove a Reed–Solomon proximity statement.

Python /opt/miniconda3/bin/python3 has SymPy1.14.0, NumPy2.4.5, fpylll0.6.4, SciPy1.15.3 and matplotlib. No packages or reset credits were used. Child agents remain unavailable after authoritative usage/authentication failures; do not retry or redeem reset credits. Reviews are separately implemented root-run checks, not human or Lean certification. Keep the goal active.
