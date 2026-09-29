# Paley graph conjecture research

**Status: open. No proof of the Paley graph conjecture or the Proximity Prize is claimed.**

The [Stepanov wave](research/stepanov-pass-summary-2026-09-27.md)
(2026-09-26/27) proves a robust form of the Hanson–Petridis inequality
and, from it, a constant-factor cancellation bound for arbitrary sets:
if `|A||B| ≥ (1/2+κ)p` then `|Σχ(a+b)| ≤ (1 − (1−(1+2κ)^{−1/2})² + O(p^{−1/2}))|A||B|`.
Chung's bound is trivial below `|A||B| = p`, and the threshold `1/2` is
sharp over all set sizes (its witness has `|A| = 2`). The key step is a Hankel determinant of the Hanson–Petridis
polynomial's derivatives, whose error part has low rank identically in
`x`. The argument was checked by the worker's verifier (1.2 million exact
checks) and by an independent exhaustive program (1.07 billion cases,
`p ≤ 29`), and an independent referee agent re-proved every step
(257.8 million further checks). It has not been reviewed by a human mathematician, and novelty is
not established beyond targeted searches. The main theorem (constant cancellation), the Hankel-minor inequality and
Hanson–Petridis itself (every proper divisor) are formalized in Lean 4 and
compile with only the standard axioms in two Mathlib versions; they are not
yet submitted to prove2me. This is
not a proof of the conjecture: it gives no power saving and does not go
below the square-root scale.

[Prove2Me project setup](PROVE2ME.md) reuses the existing private Paley
proposal and pinned Lean environment. Server proof verdicts and the
remaining formalization scope are recorded there.
The [projection algebra](research/prove2me-projection-core-2026-09-05.md)
and its [code-specific application](research/parallel26-mca-formalization-2026-09-05.md)
now have local Lean proofs with only standard axioms. The exact MCA
interleaving theorem was already present in the pinned ArkLib source;
our proof independently verifies an existing result.

The [latest assessment](research/parallel53-pass-summary-2026-09-06.md)
returns to the direct two-anchor spectral route. The ambient elliptic
kernel has inversion-odd eigenvalues at least p-32pi^2sqrt(p), by an
exact Jacobi diagonalization and a checked source discrepancy theorem.
Thus the desired2p/3 saving cannot come from an ambient norm bound.
It must use the actual common-neighbor restriction and symmetry average,
or another argument on that operator. Exact coefficient, matrix and
eigenvector checks pass; no uniform target bound improves.

The [fifty-second assessment](research/parallel52-pass-summary-2026-09-06.md)
charges new primitive collisions to the increase in cubic incidence
count at each step. An exact inequality gives a sufficient increment
bound for the energy route, but the published estimate remains too
weak. A quartic example has positive cubic excess unchanged across
one step, excluding automatic growth there. All368 distinct step checks
pass. No uniform exponent or full-goal proof is established.

The [fifty-first assessment](research/parallel51-pass-summary-2026-09-06.md)
finds an order256 quartic example with a triple collision at order128
and none in the newest primitive layer. Its full energy retains the
earlier excess. An exact compatibility equation connects lower labels
to their endpoint images, and the energy lower bound now retains its
quadratic baseline. All71 tower checks pass. These constraints give no
uniform upper bound; the conjecture and prize remain unproved.

The [fiftieth assessment](research/parallel50-pass-summary-2026-09-06.md)
recovers endpoint energy from the total first-precision triple mass,
with matching powers in its upper and lower bounds and a logarithmic
gap. A sharper recovery from an intermediate subgroup absorbs smaller
levels using the existing energy theorem. The remaining sufficient
estimate concerns only the upper part of the tower, starting near
N^(80/87) with a logarithmic correction. That estimate remains unproved.
Exact checks pass in70 towers; uniform exponents and full goals remain
unchanged.

The [forty-ninth assessment](research/parallel49-pass-summary-2026-09-06.md)
uses all coefficients of the resultant polynomial to obtain an exact
valuation formula and to show that primitive triple fibers require a
splitting prime. The target primes remain in that class. Separate
polynomial examples show why neither two resultants nor a simple sum
of collision masses always gives the exact valuation. All checks pass,
but the refined certificate equals the previous one in the checked
subgroup orders, and no uniform exponent improves. The full goals
remain unproved.

The [forty-eighth assessment](research/parallel48-pass-summary-2026-09-06.md)
certifies exactly three quartic primes at subgroup order128 with
primitive fibers of size three. This refutes the proposed uniform
maximum-fiber-two condition. Their balanced energy remains below the
quadratic bound. Complete resultant factorizations classify this
condition through order128; an aggregate of larger fibers gives a
more targeted sufficient input for the energy route. Its uniform
estimate remains unproved. All exact checks passed; the conjecture,
prize and current uniform exponents remain unchanged.

The [forty-seventh assessment](research/parallel47-pass-summary-2026-09-06.md)
uses known matrix entries to exclude both earlier uncentered interval
models. Centering and a product congruence repair those defects formally,
but all2,048 tested positions fail a field power-sum condition at one
certified quartic prime. The first n-1 power sums give an exact recognition
criterion for the subgroup's full difference profile. Finite checks
passed; this recognition result gives no new uniform coefficient bound.
The conjecture, prize and current uniform exponents remain unchanged.

The [forty-sixth assessment](research/parallel46-pass-summary-2026-09-06.md)
derives a necessary row-capacity inequality from actual incidence data.
It excludes the earlier interval with Sidon filler at sufficiently large
size. A long-interval filler preserves the difficult correlation scale
and satisfies the total-excess row condition, so no uniform exponent
improves. Exact checks passed. Local incidence allocations and the
off-diagonal multiplication identities remain to be controlled.

The [forty-fifth assessment](research/parallel45-pass-summary-2026-09-06.md)
extends the exact incidence-matrix algebra to arbitrary real weights,
including product, commutator, trace, and orthogonal projection identities.
The direct norm estimate from these identities is weaker than the current
triangle bound. Checked shifted-product results also leave the remaining
level concentration possible. All finite checks passed; no uniform
exponent or full goal improves. The next step must control the actual
weighted incidence correlation using more than the scalar moment bounds.

The [forty-fourth assessment](research/parallel44-pass-summary-2026-09-06.md)
proves a stronger mass bound when actual difference levels lie in cosets
of larger subgroups. It reaches the intermediate triangle power under
an explicit nested-coset partition hypothesis, which remains unproved
for general levels. An abstract interval family shows that the new mass
limits alone cannot improve the functional exponent. Exact matrix
identities retain further actual subgroup constraints for the next step.
Finite checks passed; no uniform exponent or full goal improves.

The [forty-third assessment](research/parallel43-pass-summary-2026-09-06.md)
improves the absolute exceptional-prime count for the intermediate
triangle target to O_c(n^(5/3)/log n), using additive energy and a
valid single-coset incidence bound. The estimate for every prime remains
open. An actual order32 subgroup now rules out diagonal domination
inside the quartic window, and exact prime-power lifts account for all
collision valuation losses. Finite checks passed; root completed the
review while the parallel lanes were unavailable at their usage limit.
The uniform energy/period estimates and full goals remain unchanged.

The [forty-second assessment](research/parallel42-pass-summary-2026-09-06.md)
bounds the nontrivial shifted-product collision excess summed over primes
at each fixed dyadic subgroup order. It gives an absolute bound on the
number of quartic-range primes that can fail the triangle criterion;
the estimate for every prime remains open. An actual subgroup with zero
diagonal rich mass and positive off-diagonal mass rules out one general
shortcut. Exact reflection smoothing also retains a recovery certificate
no better than direct sixth-moment Holder. Finite checks passed. The
parallel lanes returned findings before hitting the account usage limit;
root completed the scoped review. Energy49/20, period71/72 and the full
goals remain unchanged.

The [forty-first assessment](research/parallel41-pass-summary-2026-09-06.md)
proves a sharper uniform weighted-triangle bound,
W<<n^(86/15)(1+log n)^(8/15), retaining every incidence multiplicity.
Its power is still1/15 above the target for this route. An exact diagonal
identity shows that a nontrivial shifted-energy bound X<<n^(61/31)
would suffice to reach that triangle power; this new input is unproved.
Real-subfield recovery also reduces an arithmetic cofactor loss, while
the classical smoothed-moment estimate retains an unresolved recovery
error. Root and parallel agents reviewed the proofs and exact finite
checks. Energy49/20, period71/72 and the full goals remain unchanged.

The [fortieth assessment](research/parallel40-pass-summary-2026-09-06.md)
refutes two proposed shortcuts on actual quartic-window subgroups:
a short shell need not have a principal evaluation kernel, and distinct
edge cosets in actual difference levels can retain positive repeated
incidence fibers. An algebraic reciprocal converse has an explicit,
potentially enormous cofactor. The classical route also retains precise
structural and error-normalization restrictions. Independent exact checks
passed; the uniform estimates and full goals remain open.

The [thirty-ninth assessment](research/parallel39-pass-summary-2026-09-06.md)
constructs an exact subgroup coset refuting the literal finite amplitude
constant C=2. A cyclotomic norm cofactor forces a small shell coset;
existing sublinear cancellation prevents these small-cofactor templates
from scaling in the quartic window. A separate weighted-triangle bound
controls coincident edge cosets, leaving three distinct cosets open.
Independent agents reviewed the new constructions and the earlier
moment and lattice comparisons. No uniform exponent or full goal improves.

The [thirty-eighth assessment](research/parallel38-pass-summary-2026-09-06.md)
resolves an energy-source discrepancy through the published paper and
two parallel reviews. The old arXiv22/9 claim is not the published32/13
result; its incidence argument omits a hypothesis that fails on an exact
dyadic example. The current49/20 energy input remains in place. This
rules out an unsupported reuse of the old argument; no full target or
uniform period exponent improves.

The [thirty-sixth assessment](research/parallel36-pass-summary-2026-09-05.md)
compares centered subgroup periods with exact dual-lattice distance
deviations, with absolute constants and no exponent loss. The required
uniform annulus estimate remains unproved; the elementary cyclotomic
norm lower bound is too weak. A rational inverse-series certificate
reproduces the known maximizing coset in the order64 example without
high moments, and bounds every other coset below39.838 in absolute value.
It does not improve the existing best enclosure of the maximum. No
uniform period exponent or full conjecture bound improves. Pass39 supplies
separate-agent review of sections1-4 and clarifies the nonzero subgroup
domain; Lean review remains outstanding.

The [thirty-fifth assessment](research/parallel35-pass-summary-2026-09-05.md)
proves a Gaussian-scale lower bound for the centered opposite-free
count and compares its upper bound with the full moment at the same
degree. An established derivative-root theorem is applied to each Fourier
row before averaging. This removes the need to assume bounds at all
smaller degrees. The positive upper estimate remains unproved, so no
period exponent or full conjecture bound improves. The uniform argument
has exact finite checks and pass39 separate-agent review; Lean review
remains outstanding.

The [thirty-fourth assessment](research/parallel34-pass-summary-2026-09-05.md)
separates opposite pairs from the centered distinct-coordinate moment.
An exact positive transform and controlled inverse preserve Gaussian
hierarchy bounds up to absolute constants. The required upper bound
on the centered opposite-free aggregate remains unproved. A separate
lower bound forces at least n^6/2 distinct, opposite-free ten-term
relations in every eligible case with n>=2^35 and p<=n^4, confirming
that the principal-term subtraction cannot be dropped. No full-energy
or period exponent improves, and no full conjecture proof is claimed.

The [thirty-third assessment](research/parallel33-pass-summary-2026-09-05.md)
removes repeated coordinates from the centered subgroup moment target
with an explicit relative error tending to zero at logarithmic depth.
Each collision partition retains its own principal-term subtraction.
It also improves the recorded repeated-six-word exponent from 69/20
to 129/40. The distinct-coordinate upper bound remains unproved, so
the full sixth-energy, period, Paley, and prize targets are unchanged.

The [thirty-second assessment](research/parallel32-pass-summary-2026-09-05.md)
proves a uniform linear bound for derivations by one triangle orbit:
every generated six-term relation needs at most 8m triangle copies,
where m is the order needed to contain the normalized triangle.
A Fourier comparison with an explicitly invertible trinomial proves
the bound. It controls derivation length as the order grows; the number
of six-term outputs and all full conjecture bounds remain unproved.

The [thirty-first assessment](research/parallel31-pass-summary-2026-09-05.md)
classifies the minimum triangle derivations of every six-term remainder
orbit in the new counterexample. Twenty primitive orbits need four
triangle equations; the remaining 45 need at least 38-44, despite leaving
only six entries after cancellation. Exact quotient certificates prove
these are minimum lengths. This reveals relations missed by a few local
triangle patterns, but gives no uniform upper bound or new exponent.

The [thirtieth assessment](research/parallel30-pass-summary-2026-09-05.md)
finds D6=(5355/1024)n^3 at n=128, p=215535361, refuting the
constant-one extrapolation from the earlier census. A short zero-triangle
argument independently proves D6>n^3 in this field. General triangle-pair
counts also show that a uniform cubic remainder bound would require a
square-root-size bound on H intersect (1-H). No unspecified cubic bound,
Paley target, or prize statement is disproved or proved by this example.

The [twenty-ninth assessment](research/parallel29-pass-summary-2026-09-05.md)
proves a centered subgroup energy recurrence with the origin retained.
It removes the uniform term exactly, but the cancellation exponent does
not improve. An explicit envelope shows that deeper iteration, products
of coset moments, and feedback of the current amplitude bound still
yield at most the recorded 1/72 saving. The argument has ordinary proofs
and bounded checks; separate-author review and Lean verification remain
outstanding. No new Lean compilation was launched during this pass.

The [twenty-eighth assessment](research/parallel28-pass-summary-2026-09-05.md)
combines existing energy estimates to derive the subgroup bound
M << n^(71/72)(log n)^(1/18) throughout n^4/4<=p<=n^4. This improves
the exponent recorded in the project by 1/320. The argument retains
the origin term and has an ordinary proof; separate-author review is
outstanding. No literature novelty, current-best status, full Paley
estimate, or prize certificate is claimed. A complete finite census at
orders 4-64 also finds D6<=n^3 in all 28,774 eligible cases.

The [twenty-seventh assessment](research/parallel27-pass-summary-2026-09-05.md)
proves disjointness of the distinct balanced six-term classes and gives
a faster exact count of the remaining subgroup relations. A separate
direct enumeration agrees in all 17 cases, including six new samples at
orders 512 and 1024. The algebraic core passes Lean. No uniform upper
bound on the remainder or new cancellation exponent follows.

The [twenty-sixth assessment](research/parallel26-pass-summary-2026-09-05.md)
corrects that attribution and records the seven newly checked principal
statements. No scalar MCA bound, list bound, cancellation exponent, or
full conjecture is proved. Prove2Me has since accepted the finite
projection-counting proof, making five accepted elementary proofs. The
parallel workers last reported usage-limit errors.

The [sigma pass](research/sigma-pass-summary-2026-09-05.md) (a separate
parallel review, 2026-09-05) proves no new cancellation exponent and
confirms the conjecture is open. It machine-checks the elementary layer in
Lean 4 (shift orthogonality, second moment, Chung bound, interval-to-
non-residue link; axioms propext, Classical.choice, Quot.sound), proves that
the subgroup case B=H is exactly equivalent to the open shifted-subgroup
bound (Bourgain's Problem 5 in Chang's survey), proves that every
termwise-Weil dilation-moment bound is invariant under A↦uA and so cannot
single out t=1, proves a no-go for Burgess–Chang shift chains on
additively unstructured sets, and records that Chang's 4/9 theorem needs
small additive doubling, with no published exponent below 1/2 for two
arbitrary sets. Its referee found no counterexample to passes 9–10 in
1.4 million exact checks and classifies them as plausible with
citation-level gaps; passes 11–24 are not covered by that review. The pass-7 subgroup bound is a weaker form of
Konyagin's 2002 inequality, and its density-one class is exactly the
unproved single-prime E₃ hypothesis. A prove2me mission package is
prepared but not uploaded. A second wave then proved that the conjecture is
equivalent to domination of the dilation 2k-th moment by its square tuples
(with the dilate t=0 excluded), proved that "cancellation across cosets"
in the Gauss-period expansion of the shifted-subgroup sum is a relabelling
of that problem, and refuted the first wave's Conjecture SI with exact
witnesses (p=97, subgroup of order 8), replacing it by a signal-to-noise
version SI*. A third wave then refuted both remaining candidates (the
tuple-level Conjecture T(k) for every constant, via A=B=Q and subgroup
pairs; SI*, via designed 1×n rectangles) and obtained server-side
ACCEPTED verdicts from prove2me for the four Lean proofs, all as private
items with a private draft mission awaiting the user's launch. Thirteen
verifiers, about 2.9 million exact checks, all pass; independent human
review remains outstanding.
The [twenty-fifth assessment](research/parallel25-pass-summary-2026-09-05.md)
reduces the pinned prize combination-round count inequality exactly from
eight- and sixteen-row codes to scalar codes over the same extension
field. Its MCA portion is now identified as an existing ArkLib result.
The scalar bounds and separate spot check remain unproved; the list
argument still needs separate-author review. A uniform subgroup
bound controls repeated entries, leaving the six-distinct fully
unbalanced family and stronger error estimates open. The full spectral
decomposition shows that uniform-vector iteration cannot reach about
five-sixths of the space. Signed inversion averages retain the unknown
moment. Three workers saved their results before usage-limit errors;
root completed their checks and reviewed their proofs. No worst-case
cancellation exponent or full conjecture is proved.
The [twenty-fourth assessment](research/parallel24-pass-summary-2026-09-05.md)
gives a bounded unsigned moment for a B_h inverse of every input set;
the signed moment that transports the original input remains uncontrolled.
It evaluates the next actual spectral coefficients and proves an upper
estimate on span{u,Bu}, including its nonzero leakage. It also classifies
all quartic primes at subgroup orders 4, 8, and 16. These are limited
results: the exceptional classical inputs, growing subgroup orders,
remaining spectral operator, and official prize certificate stay open.
The four lanes have separate-author agent reviews and exact finite checks.
No worst-case cancellation exponent improves.
The [twenty-third assessment](research/parallel23-pass-summary-2026-09-05.md)
proves the sixth-moment/quartic-star upper bound for all but a quantified
vanishing proportion of Sidon inputs. Exceptional inputs remain open.
Multiplicative closure bounds the product-balanced part of subgroup
sixth relations; the unbalanced part remains open. In actual two-anchor
neighborhoods, the uniform direction's coupling tends to 1/8, so it
does not decouple from the remaining operator. A restricted code-list
bridge and exact official remainder counts are stated, but no prize
certificate follows. No worst-case cancellation exponent is improved.
The [twenty-second pass](research/parallel22-pass-summary-2026-09-05.md)
combines the completed parallel lanes. It controls the uniform vector's
Fourier mass on actual two-anchor neighborhoods, reformulates the
classical sixth-moment target as a quartic-correlation energy, and
separates the subgroup's opposite-free six-term relations. A stronger
weighted model passes all individual correlation orders and all lower
coefficient-vector moments but fails the next moment. Every required
uniform upper bound remains open; no new cancellation exponent,
all-vector spectral bound or official prize implication follows.
The [twenty-first pass](research/parallel21-pass-summary-2026-09-05.md)
reduces each fixed even-moment estimate to a one-sided upper bound on
the top squarefree character-correlation aggregate. Its negative side
already has the Gaussian order. A weighted sign model passes exact
second moments, a Gaussian fourth moment, and individual correlation
bounds through order six, but its sixth moment is too large. It is
not a prime-field character kernel or a Paley counterexample. The
actual upper estimate and every full target remain outstanding. A
separate-agent review found no mathematical correction; human review
and formal verification remain open. No new cancellation exponent follows.
The [twentieth pass](research/parallel20-pass-summary-2026-09-05.md)
uses inversion to remove short additive relations while retaining every
input element. An exact signed-moment identity and completion argument
transfer bounds from B_h sets to all sets of the same size. On the
single-size slice, h can now be as large as floor((r+1)/2) for moment
order 2r. The uniform moment bound itself is unproved; no new Paley
cancellation exponent, clique bound, subgroup bound or prize result
follows. Exact checks pass, and independent review remains outstanding.
The [nineteenth pass](research/parallel19-pass-summary-2026-09-05.md)
reduces the full classical conjecture to input sets with no nontrivial
equal h-term sums, for any fixed h. An exact partition and remainder
bound also restrict the sufficient single-size moment criterion to
such sets, with h growing more slowly than the moment order. The
required character estimate on that restricted class remains unproved.
All exact checks pass; no new arbitrary-set exponent, clique bound,
subgroup square-root bound, or official prize implication follows.
Independent mathematical review remains outstanding.
The [eighteenth pass](research/parallel18-pass-summary-2026-09-05.md)
applies an existing SL₂ expansion theorem to every polynomial-sized
Cartesian Paley family. The exact Weil trace includes an overlap
correction; controlling that centered trace remains the original
unproved cancellation problem. A concrete conjugacy average has trace
one while its operator norm tends to zero, showing why generic norm
control alone is insufficient. No new bilinear, clique, or prize bound
follows. The exact finite checks pass; independent review remains open.
The [seventeenth pass](research/parallel17-pass-summary-2026-09-05.md)
uses prime-order Fourier uncertainty and Paley determinant divisibility
to prove an explicit gap from eigenvalues 0 and 1. An exact second
moment shows why this determinant certificate decays exponentially at
the neighborhood sizes needed for the conjecture. It therefore does
not prove the required spectral edge or improve the clique bound.
The next estimate needs quantitative control of quadratic-residue
Fourier mass on the actual anchor neighborhoods. Independent review
and the full proof remain outstanding.
The [sixteenth pass](research/parallel16-pass-summary-2026-09-05.md)
tests the sharp correlations and full projection identity together.
Actual Paley graphs over square finite fields retain both, yet have
persistent localized outliers coming from a subfield clique. This
quantifies why the current route needs an additional prime-field
estimate. It is not a prime-field counterexample or a stronger clique
bound. The next decisive estimate and independent mathematical review
remain outstanding.
The [fifteenth pass](research/parallel15-pass-summary-2026-09-05.md)
extends the elliptic estimate to all translated products, with exact
repeated-shift corrections and arbitrary character twists. An abstract
sign-kernel construction preserves their square-root orders with larger
constants while producing a persistent extreme outlier. It retains the
exceptional border and symmetries, but loses the sharp constants and
ambient Paley identity; it is not a Paley counterexample. The next norm
argument must combine more of the actual arithmetic constraints. No
stronger clique bound or full proof follows, and independent review of
the cohomological argument remains outstanding.
The [fourteenth pass](research/parallel14-pass-summary-2026-09-05.md)
gives an exact elliptic-curve model of a three-anchor neighborhood,
including its exceptional points and centering term. Four symmetries
split the matrix into smaller blocks; the existing finite outlier
reduces to an 8×8 integer certificate. A uniform square-root bound
for the elliptic Fourier coefficients does not yet control the
matrix norm: the Fourier matrix still has off-diagonal terms.
No new asymptotic operator or clique bound follows from this pass.
The cohomological scalar bound needs independent mathematical review.
The [thirteenth pass](research/parallel13-pass-summary-2026-09-05.md)
improves the complete spectral bound's logarithmic depth range by separating
principal ranks from boundary errors. An exact positive weighted count gives
a smaller rank-growth rate; for two anchors the leading depth coefficient
increases by about 16.5%. Every correction remains included. The required
longer-depth estimate and any resulting improvement of the clique bound
remain unproved. This argument is source-dependent and needs independent
mathematical review.
The [twelfth pass](research/parallel12-pass-summary-2026-09-05.md)
rules out extending the short-depth estimate using only projection
identities and the selected anchor data: a rank-two modification preserves
those inputs while planting an extreme outlier. An exact actual Paley
example at p=257 separately rules out an edge interval at every finite
prime. Neither result refutes an asymptotic conjecture or improves a
clique bound. The next spectral estimate must use additional character
structure away from the anchors.
The [eleventh pass](research/parallel11-pass-summary-2026-09-05.md)
bounds the complete normalized spectral trace by (1+o(1))p in an
explicit logarithmic depth range, retaining J, anchor corrections,
exceptional fibers and the constant direction. A positive matrix-energy
recurrence isolates the still-unproved longer-depth condition. No improved
clique bound or full Paley proof follows from the currently proved range.
The argument is checked locally; independent review remains outstanding.
The [tenth pass](research/parallel10-pass-summary-2026-09-05.md)
controls arbitrary signed combinations of all word kernels in an explicit
depth range, including restriction to further adjacency conditions.
It also proves that an unrestricted isometry must fail at longer depths
by dimension. The next spectral step must exploit its specific coefficients.
These are locally checked source-dependent arguments; independent review
remains outstanding. The full spectral edge and Paley targets remain open.
The [ninth pass](research/parallel9-pass-summary-2026-09-05.md)
gives a source-dependent proof of the individual necklace estimate
at every fixed degree and length, with explicit bound
3a(2a+2)^(k−2)p^((k+1)/2) for k≥3. Rank growth and a strict
weight gap control the full convolution chain and its corrections.
The written proof has been checked locally; independent review and
formal verification have not been obtained. Its exponential constant
does not settle the growing-depth signed aggregate or spectral edge.
The [eighth pass](research/parallel8-pass-summary-2026-09-05.md)
retains the preceding all-word length-six result and its sharper constants.
The [seventh pass](research/parallel7-pass-summary-2026-09-04.md)
proves the current subgroup bound M≤(17+log N)^(1/9)N^(8/9)
on the existing density-one class of quartic splitting primes.
The uniform subgroup target, exceptional primes, stronger centered
moments and classical arbitrary-set estimates remain open.
The [sixth pass](research/parallel6-pass-summary-2026-09-04.md)
records the preceding 17/18 bound, the first adjacent-block
contraction, enlarged kernel family and typical-set fourth moment.
The [fifth pass](research/parallel5-pass-summary-2026-09-04.md) proves
the sixth-energy class, all 243 length-five words and direct anchor weights; the
[fourth pass](research/parallel4-pass-summary-2026-09-04.md) proves
the large exact fourth-energy class and the rank-pair aggregate; the
[third pass](research/parallel3-pass-summary-2026-09-04.md) proves the
arbitrary-gap family and records the conic representation; the
[second pass](research/parallel2-pass-summary-2026-09-04.md) records the
aggregate-to-clique criterion; the
[first pass](research/parallel-pass-summary-2026-09-04.md) records the
preceding 189-word length-six result and prize-side restriction.

The user's existing prize research, recovered at commit
`5b00e50c3c51b3c944201a1749a1a8132e5ce167` of
[`SlopDotCash/proximityprize`](https://github.com/SlopDotCash/proximityprize), uses
the following **primary working target**. For the dyadic multiplicative subgroup
`H=μ_n ⊂ F_p*`, bound

\[
\max_{b\ne0}\left|\sum_{h\in H}e^{2\pi i bh/p}\right|
\le C\sqrt{n\log(p/n)}
\]

uniformly in the specified thin-subgroup regime. The historical shorthand
`p≈n⁴` still needs precise quantified limits. See
[the subgroup target and audit](research/subgroup-target.md).
The full reduction to the official prize has not been independently verified.

A separate line investigates the classical prime-field quadratic-character conjecture:
for every `0 < ε < 1`, there exist `δ > 0` and `p₀` such that, for every prime
`p > p₀` and all subsets `A,B ⊆ F_p` with `|A|,|B| > p^ε`,

\[
\left|\sum_{a\in A,b\in B}\chi_p(a-b)\right|
\le p^{-\delta}|A||B|.
\]

The two-set formulation is Conjecture 7 in
[Satake (2020)](https://arxiv.org/abs/2011.02907).
It implies subpolynomial clique sizes for prime Paley graphs. A polylogarithmic
clique bound is a stronger target. A uniform `O(log p)` bound is already
ruled out by the Graham–Ringrose lower bound, as recorded by
[Magsino–Mixon–Parshall](https://arxiv.org/html/1907.05971).
For `p ≡ 3 (mod 4)`, the two-set formulation remains meaningful even though
the usual undirected Paley graph is replaced by a tournament.

The [Ethereum Proximity Prize](https://proximityprize.org/) instead states
Reed–Solomon mutual-correlated-agreement and list-decoding challenges.
**A quantitative reduction from the conjecture above to those exact challenges
has not been established here.** Shared pseudorandomness obstacles are not an
equivalence theorem. See [the source audit](research/source-audit.md).

## Results in this workspace

- [Research notes](research/moments-and-obstructions.md): proofs of the classical
  second-moment estimate, the higher-moment barrier, an exact fourth-moment
  decomposition, and obstructions to proposed stronger moment bounds. These
  are foundational results and failed-route analysis; no novelty claim is made.
- [Lean proof](research/MomentObstruction.lean): a concrete failure of a Gaussian
  eighth-moment bound over `F_1009`, with `|B| = 30 < √1009`.
- [Exact experiments](experiments/paley_exact.py): integer-only checks of
  identities, exhaustive small rectangle maxima, and adversarial examples.
- [Results](results/exact_checks.json): complete finite-check counts and witnesses.
- [Next proof obligations](research/frontier.md): the unproved estimate needed
  by the current moment approach, with its exact implication and limitations.
- [Subgroup moments](experiments/subgroup_moments.py): exact additive-energy
  computations for the primary target, with the zero-frequency contribution
  removed and multiplicative-coset repetitions accounted for.
- [Quartic-window certificate](research/SubgroupQuarticCounterexample.lean):
  exact arithmetic refuting a stronger Gaussian-moment coefficient at
  `p=6700417`, `|H|=64`. This leaves the primary asymptotic target open.
- [Finite spectral obstruction](research/finite-spectral-obstruction.md): an
  analytic certificate that the same example violates the literal amplitude
  constant `√2` with natural `ln(p/n)`. Exact moments through depth 12 give
  `43<M≤√1970`; the unspecified-constant conjecture remains open.
- [Uniform geometric-cycle bound](research/geometric-lift-and-alias.md): a
  proof of an absolute-constant moment bound for the lifted cycle, with the
  remaining prime-reduction error stated explicitly. Exact carry counts
  independently reconstruct the twelve prime energies and locate the first
  additional prime relation at order 8. The required error estimate is
  equivalent to the original moment bound up to constants and remains open.
- [Prize reduction audit](research/prize-reduction-audit.md): the primary
  MCA definition, an elementary circuit bound with exact attainment over
  sufficiently large fields, and a proof that monomial pairs need not attain
  the maximum. Explicit smooth codes at `p=2017` and `p=65537` have maximum
  56 bad scalars, while monomial pairs attain at most 40. This invalidates an
  unrestricted proposed reduction, not a Paley conjecture.
- [Coset interaction identities](research/coset-coherence.md): exact linear
  and nonlinear constraints on periods, with constructed examples showing
  that the linear constraints alone permit large spikes. The constructed
  vectors fail the arithmetic identity and are not Paley counterexamples.
  Exact checks cover five prime fields.
- [Certified maximum](research/period-polynomial-certificate.md): all signed
  moments through order 24 give rational polynomial certificates identifying
  the exact maximizing coset for `p=6700417`, `n=64`. Here `M=η₁` lies in
  `[43.802482797626304198,43.802482797626304199]`; all other cosets have
  absolute period below 42. This is a finite theorem, not an asymptotic proof.
- [Recurrence coefficient route](research/jacobi-coefficient-route.md): a
  proved implication from a stated uniform coefficient bound to the target
  estimate. The uniform hypothesis remains unproved; a new exact
  quartic-window witness with `p=67403009`, `n=128` refutes its literal
  constant `B=1` at degree two. Other absolute constants remain possible.
- [Cyclotomic bound summed over primes](research/cyclotomic-prime-average.md):
  an unconditional weighted bound on extra zero relations for a fixed
  dyadic order. It limits the count of primes with large fourth energy,
  but does not give a worst-case or logarithmic-depth estimate. Exact
  determinant audits classify every prime with extra zero quadruples
  for subgroup orders 4, 8, and 16. A sharper AM–GM version still gives
  no useful pointwise fourth-energy bound in the quartic window.
- [Official profile and trace audit](research/official-profile-and-trace.md):
  the newer official benchmark uses a base-field domain embedded in a
  sextic extension, outside the current quartic-window target. An exact
  trace-zero witness gives full-size extension-field periods, so any
  spectral bridge must account for the entire trace-zero space. This
  is a restriction on possible reductions, not a prize counterexample.
- [Subset sums and decoding lists](research/subset-sums-and-lists.md): an
  explicit character-sum bound gives subset counts, exact lists for a
  monomial word, and the bad scalars of a particular monomial pair. Root
  lifting certifies at least `2^8154` nearby codewords at radius
  `4095/8192` for the pinned official profile. This excludes that radius
  from its MCA-plus-list lower certificate, not from protocol security.
  A cyclotomic norm argument prevents extrapolating this construction
  to fixed-gap asymptotics. The general prize bridge remains open.
- [Large-spectrum and Riesz-product audit](research/riesz-tail-route.md):
  a self-contained entropy bound explains why applying dissociation to
  one multiplicative orbit cannot give square-root cancellation. An
  exact quartic-window example refutes a stronger Riesz normalization.
  The corrected sufficient estimate is equivalent to the target up to
  constants and remains unproved.
- [Known bounds and amplification](research/analytic-bounds-and-amplification.md):
  an existing theorem gives `M≤n^(2849/2880+o(1))` in the working quartic
  window. An exact exponent calculation shows the limited saving supplied
  by one proposed amplification ledger, even with ideal full-energy inputs.
  A separate centered-moment formula is screened with its quantitative
  limitations and a qualification needed in its general-function version.
- [Mixed periods and shifted energy](research/mixed-periods-and-shifted-energy.md):
  the full mixed system recovers exactly the actual periods. An exact
  correspondence identifies triple intersection counts with nontrivial
  multiplicative collisions in `(H−1)\{0}` and compares them with fourth
  additive energy. Two larger finite quartic-window subgroups are certified
  circular, with minimum fourth energy; no uniform bound follows.
- [Fourth-energy discriminant](research/kernel-discriminant.md): one integer
  polynomial classifies all primes with excess fourth energy at a fixed
  dyadic order. Complete factorizations through order 32 show minimum
  fourth energy throughout those orders' quartic windows. An exact
  valuation identity retains collisions at higher prime-power precisions.
- [Quadruple orbits and power factors](research/quadruple-orbits-and-cube.md):
  a uniform orbit argument factors the odd part of that integer into an
  explicit factor, a cube and a sixth power. Repeated-entry quadruples
  contribute only O(n²) energy; controlling the orbits of four distinct
  elements remains necessary and unproved. Exact checks cover 55 cases.
- [Dyadic descent and mixed energy](research/dyadic-descent-and-mixed-energy.md):
  an exact recurrence reduces a sufficient fourth-energy induction to
  mixed energy between the two cosets at each step. Field-degree growth
  forces all nontrivial quadruples to descend and restricts the prime
  support of the integer factors, but that growth is absent from the
  target prime-field tower. The needed mixed-energy bound remains open.
- [Positive products and high moments](research/positive-product-moments.md):
  a sufficient bound on positive products of the two child periods would
  yield the desired logarithmic-depth moment estimate. An exact centered
  balanced-count condition verifies that hypothesis in the known order-64
  example. Raw balanced counts necessarily acquire extra finite-field
  relations from order 128 onward in the working range; the uniform
  centered bound is unproved.

- [Signed quotient operators](research/signed-quotient-operators.md):
  sparse symmetric integer matrices have exactly the nonprincipal period
  and child-product spectra. Taking entrywise absolute values gives nearly
  the trivial degree throughout the quartic window, so that comparison
  cannot provide the missing cancellation. Exact checks cover 14 matrices
  and 94 trace identities; no uniform signed spectral bound follows.

- [Lists and the actual winning set](research/list-to-winning-set.md):
  a scalar-list projection argument gives a uniform unsafe suffix for the
  pinned benchmark starting at `122641/262144`, about 0.4678383, with
  score inequality 11649 centibits. The July ABF paper was recovered and
  checked. This is an ordinary proof with exact arithmetic, not a Lean
  submission, a sharp threshold, or a solution to Paley.

- [Two-anchor necklace identities](research/localized-necklace-identities.md):
  exact formulas evaluate cyclic quadratic-character sums when one of two
  singleton anchors occurs once or twice. A published Kloosterman estimate
  gives the uniform bound `k² p^(k/2)` for the latter family. Tests include
  3294 direct necklace traces and 72 literal tuple sums. Three occurrences
  lead to a weighted cubic-character sum whose needed cancellation remains
  open; no full localization or Paley bound follows.

- [Planar necklace reductions](research/planar-necklace-reductions.md):
  affine averaging evaluates every two-block word, at arbitrary length.
  Planar duality and projective changes of variables evaluate the remaining
  binary length-six patterns, giving cancellation for all 64 singleton
  words at that length. Literal graph sums check the duals, deletion terms,
  and assignments involving infinity. General longer mixed patterns and
  the full Paley targets remain open.

- [Sigma pass](research/sigma-pass-summary-2026-09-05.md): six parallel
  directions plus a barrier map ([barriers](research/sigma-barriers-2026-09-05.md),
  [sum-product](research/sigma-sumproduct-2026-09-05.md),
  [structured set](research/sigma-structured-2026-09-05.md),
  [dilation crux](research/sigma-crux-2026-09-05.md),
  [referee](research/sigma-referee-2026-09-05.md),
  [subgroup audit](research/sigma-subgroup-2026-09-05.md),
  [Lean](research/sigma-lean-2026-09-05.md)). Exact witnesses refute four
  structural hypotheses about biased rectangles, the k! tuple constant,
  container-to-subset inheritance, coset bootstrapping, and energy rescue
  of shift chains. Lean files live in `~/prove2me_workspace`. The second
  wave adds [Fourier dual](research/sigma-dual-2026-09-05.md),
  [tuple cancellation](research/sigma-tuple-2026-09-05.md) and
  [Conjecture SI](research/sigma-si-2026-09-05.md); the third wave adds
  [stress tests](research/sigma-stress-2026-09-05.md) refuting T(k) and
  SI*, the [Frobenius synthesis](research/sigma-frobenius-2026-09-05.md),
  and the [private prove2me record](research/sigma-p2m-2026-09-05.md).

## Reproduce

```sh
python3 experiments/paley_exact.py all --output results/exact_checks.json
python3 experiments/verify_lean.py
python3 experiments/subgroup_moments.py
python3 experiments/subgroup_sparse.py
python3 experiments/verify_lean.py --subgroup
python3 experiments/verify_lean.py --quartic
python3 experiments/sharp_constant_witness.py
python3 experiments/verify_lean.py --sharp
python3 experiments/quotient_moments.py
python3 experiments/dyadic_carries.py
python3 experiments/mca_circuit_certificate.py
python3 experiments/coset_coherence.py
python3 experiments/period_polynomial_bounds.py
python3 experiments/jacobi_low_order_obstruction.py
python3 experiments/cyclotomic_norm_audit.py
python3 experiments/extension_trace_audit.py
python3 experiments/subset_sum_list_certificate.py
python3 experiments/riesz_tail_audit.py
python3 experiments/analytic_bound_ledger.py
python3 experiments/mixed_period_collisions.py
python3 experiments/kernel_discriminant.py
python3 experiments/quadruple_orbits.py
python3 experiments/dyadic_energy_descent.py
python3 experiments/mixed_high_moments.py
python3 experiments/signed_quotient_operators.py
python3 experiments/list_to_winning_set.py
python3 experiments/localized_necklace_identities.py
python3 experiments/planar_necklace_reductions.py
python3 experiments/sigma_barriers_2026_09_05.py
python3 experiments/sigma_sumproduct_2026_09_05.py
python3 experiments/sigma_structured_2026_09_05.py
python3 experiments/sigma_crux_2026_09_05.py
python3 experiments/sigma_referee_2026_09_05.py
python3 experiments/sigma_subgroup_2026_09_05.py
python3 experiments/sigma_lean_2026_09_05.py
python3 experiments/sigma_dual_2026_09_05.py
python3 experiments/sigma_tuple_2026_09_05.py
python3 experiments/sigma_si_2026_09_05.py
python3 experiments/stepanov_robust_2026_09_26.py
python3 experiments/stepanov_algebra_2026_09_26.py
python3 experiments/stepanov_pencil_2026_09_26.py
python3 experiments/stepanov_adversary_2026_09_26.py
python3 experiments/stepanov_leanhp_2026_09_26.py
cc -O2 -o star experiments/stepanov_star_exhaustive_2026_09_27.c && ./star 29
python3 experiments/sigma_frobenius_2026_09_05.py
python3 experiments/sigma_stress_2026_09_05.py
python3 experiments/sigma_p2m_2026_09_05.py
```

The quotient moment, polynomial certificate, mixed high moment, localized
necklace, and planar necklace scripts use NumPy. The signed quotient script uses NumPy and SciPy; the other
Python experiments need only the standard library. Lean verification uses the
existing Lean 4.29.1 / pinned mathlib cache in `~/TheLeaningOfEverything`, read only;
`--cache-workspace PATH` selects another matching cache. The verification log is
written to `results/lean-verification.txt`.

The source PDFs are retained locally with checksums in `sources/manifest.json`.
An archived April 8 edition of the ABF prize paper was recovered locally and
its definitions and relevant statements verified. The July 6 edition was
subsequently retrieved, archived separately, and its relevant attack and
list-construction pages visually checked. Prize statements were also checked
on the official prize site. The earlier failed download remains in the manifest
as history.

Research started 2026-09-04. The main goal remains active and unachieved.
