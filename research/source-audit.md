# Source audit, 2026-09-04

Pass-53 update (2026-09-06): the
[ambient elliptic proof](parallel53-ambient-elliptic-2026-09-06.md)
uses Lu–Zheng–Zheng arXiv1305.3405v3 Theorem1.4 equation1.7, checked
on the primary page and archived HTML. Its parameters are m=1,k=1,
A_1={chi}; the delta term vanishes and the discrepancy is at most
2p^(-1/4). The varying character runs through the full nontrivial
family, with product-trivial cases excluded. No estimate for a proper
character subset or a weighted family is inferred. Root derives the
ambient Jacobi spectrum and the large odd eigenvalue from this imported
discrepancy bound; exact finite tests do not stand in for that theorem.
The existing restricted-operator identity retains all boundary terms.
Sources and dependencies are pinned in the pass audit, with the main
manifest unchanged. No uniform target improvement, independent-agent,
Lean, external review or novelty claim is made.

Pass-52 update (2026-09-06): the
[cubic increment proof](parallel52-cubic-increments-2026-09-06.md)
derives its four-cell formula and primitive-fiber cost from the existing
incidence bijection, inversion symmetry and dyadic energy recurrence.
The archived Shkredov2015 Theorem6 was reread with its product-size
condition; it gives only the existing X_N<<N^2 log N estimate here.
Its source also explicitly records the energy5/2 consequence recovered
by the current substitution. No new quantitative source theorem is
imported. Repeated cells and shifted products are counted independently;
the quartic prime and generators come from prior certificates. Source
and dependency hashes are preserved in the audit. No uniform improvement,
independent-agent, Lean, external review or novelty claim is made.

Pass-51 update (2026-09-06): the
[collision transport proof](parallel51-collision-transport-2026-09-06.md)
reuses the primitive/kernel factorization and exact energy recurrence.
Root derives their compatibility equation and retains the quadratic
baseline in the lower bound. The witness prime already appears in the
order128 certificate; a separate trial-division check and literal
order256 enumeration verify the new quartic endpoint comparison. No
new complete factorization or prime classification is asserted. A
bounded primary-source search supplied no new quantitative input.
All proof ingredients used here are derived in the local notes; source
hashes are preserved in the pass audit. No uniform improvement,
independent-agent review, Lean, external review or novelty claim is made.

Pass-50 update (2026-09-06): the
[tower recovery proof](parallel50-tower-mass-2026-09-06.md)
subtracts the two-representation baseline before Cauchy, then uses an
explicit Lipschitz recursion to compare endpoint energy with actual
first-precision triple mass. These are root derivations from the
existing primitive-fiber and energy identities. The archived MRSS
Corollary12, equation32, was reread with its M<=sqrt(p) hypothesis;
its existing49/20 energy bound yields the80/87 cutoff with the4/29
logarithmic divisor. No new external quantitative theorem is imported.
Every radical test uses rational bounds, and abstract mass profiles
are distinguished from compatible field towers. Source and dependency
hashes are recorded in the pass audit. No new uniform bound, independent
review, Lean, novelty or full-proof claim is made.

Pass-49 update (2026-09-06): the
[critical-content proof](parallel49-critical-content-2026-09-06.md)
reuses the primitive polynomial and field-degree formulas already in
the repository. The Gauss-valuation product rule is proved by reduction
of primitive polynomials; the quadratic norm argument and exact
cancellation correction are root derivations. No new external
quantitative estimate is imported. Coefficients and lifted derivative
values are checked exactly, and nonsplitting examples include small
characteristic and field degrees four and eight. Auxiliary polynomial
examples are explicitly separated from actual subgroup data. The pinned
arithmetic backend is reused. Source and dependency hashes appear in
the pass audit. No uniform improvement, Lean, separate-agent, external
review, novelty or full-proof claim is made.

Pass-48 update (2026-09-06): the
[collision eliminant note](parallel48-collision-eliminants-2026-09-06.md)
explicitly reuses the existing primitive polynomial and balanced-energy
identity. The resultant gcd, derivative valuation bound and weighted
tower estimate are root derivations. Complete fixed-order conclusions
use exact resultants, a second matrix-determinant algorithm, separate
trial-division certificates and checks at every eligible factor. No
unverified factor, finite-prefix completeness claim or uniform
maximum-fiber-two assumption is used. The
[source and tool scope](../results/parallel48_source_scope_2026_09_06.json)
records the pinned python-flint backend and official API references;
no new external mathematical estimate is imported. The main manifest
and existing primary-source hashes are preserved. No Lean, separate
agent, external review, novelty or full-proof claim is made.

Pass-47 update (2026-09-06): the
[anchored and arithmetic note](parallel47-anchored-arithmetic-2026-09-06.md)
uses the established incidence symmetries, the derivative of X^n-1,
subgroup power sums, and Newton's recurrence. Root gives the ordinary
proofs. Both earlier uncentered models fail the stronger anchored-row
condition; a centered repair satisfies the compared coarse constraints.
Its finite rejection uses certified prime arithmetic, not a probabilistic
primality assumption. No unverified shifted-character bound, new external
theorem, uniform improvement, Lean, separate-agent review or novelty
claim is imported. Recognition by n-1 moments is explicitly separated
from bounding the profile. Existing primary-source and manifest hashes
are preserved in the pass audit.

Pass-46 update (2026-09-06): the
[row inequality](parallel46-row-capacity-2026-09-06.md) is a root
derivation from symmetry, nonnegative integer entries, row mass, and
Holder. The exact autocorrelation identification uses the previously
derived matrix-square identity. The Sidon-filled family exclusion uses
only the already scoped Shkredov2015 Theorem6 consequence X<<n^2 log n
under n^2<p. The modified abstract family satisfies the total-X row
conditions; no local X_D allocation or prime-field realization is claimed.
No new external theorem, uniform exponent, Lean, or external review is
asserted. The primary-source and main-manifest hashes are rechecked in
the pass audit.

Pass-45 update (2026-09-06): root checked Warren's published2019
Corollary4 against the surviving high-level parameters. Its exponent11/9
and |A|<p^(1/4) condition yield a valid but insufficient consequence there.
Harrison–Mudgal–Schmidt arXiv2603.06483v1 Theorem1.5 is explicitly over C;
no uniform finite-field transfer is supplied. The
[source comparison](parallel45-shifted-product-scope-2026-09-06.md)
distinguishes these checked inputs from the proposed, unproved cubic
product inequality. The full weighted matrix law and trace/projection
identities are root derivations from the existing coset-convolution
algebra. Primary sources and roles are
[pinned](../results/parallel45_source_scope_2026_09_06.json).
No novelty, new uniform exponent, separate-agent, Lean, or external review
is claimed; the main source manifest is unchanged.

Pass-44 update (2026-09-06): root read Mit'kin's lemma as quoted in
Shkredov arXiv1504.04522v1, Lemma2, equation5. A singleton pair of coset
representatives gives O((|Gamma||Pi|)^(1/3)) intersection under the
explicit product-size and33^3 conditions. Using different subgroup
sizes is essential to the new coset-mass bound. The quartic nested-coset
application verifies the field-size condition; bounded small n is
absorbed by a trivial estimate. The earlier MRSS49/20 input is reused
only to handle small edges. Root-derived interval and matrix arguments
are distinguished from imported results in the
[pass summary](parallel44-pass-summary-2026-09-06.md).
The [pinned scope](../results/parallel44_source_scope_2026_09_06.json)
includes the primary HTML hash and an exact mathematical-text extract.
No claim of a proved structural partition, novelty, separate-agent,
Lean, external peer review or full proof is made. The main source
manifest remains unchanged.

Pass-43 update (2026-09-06): root visually reread published Shkredov2013
Lemma2, equation18, printed page197. For Q=H,Q1=uH,Q2=vH, the required
ratio image is exactly Q1 x Q2 and has size n^2. This valid single-coset
application yields R_max<<n^(2/3), with all incidence counts retained.
The [energy argument](parallel43-energy-prime-exceptions-2026-09-06.md)
then combines an exact positive triangle bound with the existing local
fourth-energy prime-average proof, improving the absolute prime-exception
exponent to5/3. It uses no union-of-cosets injectivity shortcut or eligible
prime-density assumption. The quartic witness and valuation identity are
exact local arithmetic, with 21 independent norm-valuation comparisons.
The [pinned scope](../results/parallel43_source_scope_2026_09_06.json)
records source and dependency hashes. Root review and finite checks were
completed; separate-agent, Lean, external peer review, novelty and full
proof claims are not made. The main source manifest is unchanged.

Pass-42 update (2026-09-06): the [primary-source review](parallel42-excess-source-audit-2026-09-06.md)
checks Ke–Kiechle2023 Proposition3.4 and Corollary3.7 for primitive
cyclotomic norms and finite circularity exceptions. Norm divisibility
implies existence of a vanishing primitive root, not a fixed-root
converse. The weighted product over every compatible tuple avoids this
quantifier issue; root derives its sharper AGM constant locally.
DoDuc–Leung–Schmidt2020 Theorems1.2–1.3 have exponential size conditions
or a prime subgroup-order restriction that do not cover growing dyadic
orders in the quartic range. Ke–Kiechle arXiv2307.05586v2 assumes
circularity in the relevant results; its tables are not a growing-order
density theorem. Macourt–Shkredov–Shparlinski arXiv1701.06192v2
Corollary4.1 retains quadratic shifted energy in the small-group range.
None of these checked statements supplies the missing uniform X saving.
The reflection proof uses the previously scoped individual Weil bound.
The [source hashes and roles](../results/parallel42_source_scope_2026_09_06.json)
are pinned; the prior main manifest is unchanged. Separate-agent review
covers the coarse norm mechanism; root completed the sharper constant
and final checks after the subagents reached their usage limit. No
novelty, external review, Lean verification or full proof is claimed.

Pass-41 update (2026-09-06): the [final triangle proof](parallel41-independent-review-2026-09-06.md)
uses the published Shkredov2013 Theorem4,d=2 tail and Shkredov2015
Theorem6 shifted multiplicative energy, under n<sqrt(p). Root visually
checked the former and read the latter's primary statement. The invalid
general incidence application from pass38 is not reused. Mixed-level
correlations and geometric summation prove W<<n^(86/15)L^(8/15).
The conditional X<<n^(61/31) criterion instead follows from exact
diagonal counts and the functional correlation lemma; it is not an
imported or proved upper estimate for X. MRSS49/20 is retained as the
energy comparison baseline, not required to prove this W bound.
The shell and unit arguments use explicit cyclotomic arithmetic; the
classical smoothing uses the individual Weil bound, checked in Volostnov
Theorem5, plus exact independent-shift moments. The ordinary arguments
and finite certificates received independent review within this task.
No novelty, external peer review, Lean proof or full-goal completion is
claimed. Sources and their roles are [pinned](../results/parallel41_source_scope_2026_09_06.json).

Pass-40 update (2026-09-06): the [shell converse](parallel40-shell-converse-2026-09-06.md)
uses the previous norm implication, integer adjugates, and an explicitly
derived elementary Gauss-sum bound at641. The counterexample requires no
class-group computation. The optional flat-Gram discussion separately
uses cyclotomic prime-ideal factorization and retains its unit-norm
condition. The [incidence correction](parallel40-distinct-edge-input-2026-09-06.md)
is an exact count on actual levels; it imports no new incidence bound.
The [classical applicability review](parallel40-classical-structure-2026-09-06.md)
checks versioned Schoen-Shkredov and Volostnov primary statements and
keeps their size/doubling hypotheses. Its relative L2 obstruction does
not apply to standard Croot-Sisask input-norm tolerance. Root reviewed
these distinctions and reproduced the finite arithmetic with independent
checks. The existing source manifest and all prior pass files are retained.
No uniform estimate, novelty claim, external peer review or Lean proof is
added by this pass.

Pass-39 update (2026-09-06): the [principal-shell construction](parallel39-shell-structural-input-2026-09-06.md)
uses elementary cyclotomic algebra and the norm-prime trinomial already
recorded in passes30-32. Its C=2 finite obstruction is checked with exact
integer/rational arithmetic; a separate agent recomputed the determinant.
The cofactor lower bound is conditional on any valid period upper bound;
using the existing sublinear input imports that input's original scope.
The [edge-coset estimate](parallel39-edge-coset-reduction-2026-09-06.md)
uses Holder and the published nonzero difference moments on PDF page10
of the pass38 source. Those moments differ from equal-sum energies.
Separate agents reviewed pass35's comparison and pass36 sections1-4;
pass36 now explicitly restricts the norm argument to A=H,a!=0.
No new uniform bound, Lean proof, external peer review or novelty claim
is made. The prior main source manifest is unchanged.

Pass-36 update (2026-09-05): the [lattice-distance comparison](parallel36-shell-inversion-2026-09-05.md)
uses NIST DLMF24.8.1 for the Bernoulli Fourier series and25.6.1 for
zeta(2),zeta(4). The primary pages were inspected. The TeX download
returned403, so the [pinned record](../results/parallel36_source_scope_2026_09_05.json)
contains equation transcriptions, explicitly not raw page archives.
The inverse, exact centering, finite truncation error and rank-one norm
argument are derived locally. The prior moment-polynomial result already
identified the same finite maximum more precisely; pass36's alternative
certificate must not be called a better enclosure. The initial comparison
only with sqrt(1970) is corrected. No uniform annulus estimate, novelty,
separate-author review, Lean proof or full-goal completion is claimed.

Pass-35 update (2026-09-05): the [single-degree comparison](parallel35-single-degree-comparison-2026-09-05.md)
imports Ravichandran, arXiv1609.04187v2, Theorem4.4: a bound on roots
of sufficiently high derivatives of a real-rooted, mean-zero polynomial
with roots in [-1,1]. Its exact statement was read in HTML and visually
on PDF page12. Each centered cosine row is scaled into those hypotheses;
the coefficient normalization and small-size branch are kept explicitly.
Both archives are [pinned separately](../results/parallel35_source_scope_2026_09_05.json),
while the prior main manifest is preserved. Averaged real-rootedness is
not assumed and is refuted by an actual subgroup example. No positive
aggregate upper bound, novelty claim, independent review, Lean proof,
new period exponent, or full conjecture proof is claimed.

Pass-34 update (2026-09-05): the [opposite-pair transform](parallel34-opposite-pair-transform-2026-09-05.md)
uses direct counting, formal power-series inversion, cycle independent
sets, and the existing pass33 centered-collision estimate. Its raw
degree-ten lower bound uses Fourier inversion, Holder, Lyapunov, and
Jensen; multiplicative closure is unnecessary for that lower bound.
The MRSS Lemma6 proof was re-read as an alternative route, without
importing a new upper estimate. The [source-scope record](../results/parallel34_source_scope_2026_09_05.json)
pins this distinction. No literature novelty, separate-author review,
Lean verification, new period exponent, or full proof is claimed.

Pass-33 update (2026-09-05): the [centered collision argument](parallel33-centered-distinct-moments-2026-09-05.md)
uses elementary partition inversion, Holder, Jensen, and the Stirling
cycle polynomial. The improved repeated-word exponent reuses the
existing MRSS energy seeds. The [recent-source audit](parallel33-source-check-2026-09-05.md)
checks the published 2026 dense uniformity normalization visually;
its diagonal terms prevent direct application to the original sparse
subgroup. The Fourier/Bohr alternative meets its size hypothesis but
does not supply the desired upper estimate. Sources and their distinct
roles are [pinned](../results/parallel33_source_scope_2026_09_05.json).
No literature novelty, separate-author review, Lean verification, or
full proof is claimed.

Pass-32 update (2026-09-05): the [uniform triangle-length proof](parallel32-linear-triangle-length-2026-09-05.md)
uses an explicit periodic inverse, a pointwise complex triangle
inequality, finite Fourier orthogonality, and coefficient norm bounds.
No non-elementary external estimate is imported. Exact finite checks
support the algebraic implementation; they do not replace the uniform
argument. A current primary-source inventory was consulted but supplies
no premise to this proof. No literature novelty, separate-author review,
Lean verification, cubic output count, or full proof is claimed.

Pass-31 update (2026-09-05): the [minimum triangle derivation proof](parallel31-short-multiples-2026-09-05.md)
uses the pass30 principal ideal and exact integer identities. No new
non-elementary theorem is needed for the finite classification. A
separate certificate checker confirms all 119 orbit representatives and
quotients. The initial source check distinguishes repeated-difference
energies from sum moments; it does not import the former as the latter.
No independent-author review, Lean verification, uniform upper bound,
or improved cancellation exponent is claimed.

Pass-30 update (2026-09-05): the [triangle-remainder argument](parallel30-triangle-remainder-2026-09-05.md)
and exact arithmetic refute D6<=n^3 at p=215535361, n=128. This corrects
only a prospective extrapolation of the finite census, not its recorded
data or an asserted asymptotic theorem. The new proof uses elementary
orbit counting and the verified integer multiplication determinant;
no new non-elementary literature input is assumed. Independent-author
review and Lean verification remain outstanding. No O(n^3) bound with
an unspecified constant, Paley conjecture, or prize statement is refuted.

Pass-29 update (2026-09-05): the [centered recurrence proof](parallel29-centered-recurrence-2026-09-05.md)
uses the already inspected published Theorem 8, equation (22), together
with an elementary exact plane-variance identity. It establishes the
needed centered convolution estimate with the origin retained, without
assuming the qualified origin-free signed-function statement. This does
not certify every step of that earlier published argument. Existing MRSS
seed bounds are reused under their original hypotheses. The
[source ledger](../results/parallel29_source_scope_2026_09_05.json) pins
the reused files. No literature novelty, independent review, full Paley
proof or prize certificate is claimed.

Pass-28 update (2026-09-05): the typeset Shkredov 2019 paper confirms
the invariant-set energy estimate (Lemma 9, equation 23) and the
arbitrary-s recurrence (equation 28). Combining them with MRSS
Corollary 7 and Konyagin's mixed-moment inequality improves the
subgroup baseline recorded here. The [proof](parallel28-moment-recurrence-2026-09-05.md)
and [source ledger](../results/parallel28_source_scope_2026_09_05.json)
retain the origin term and distinguish the new-to-project consequence
from a claim of literature novelty. Separate-author review remains open;
no general Paley or prize proof follows.

Pass-26 attribution correction (2026-09-05): the pinned ArkLib source
already contains `mcaError_interleaved_eq`. The pass-25 MCA result is an
independent rediscovery. See the [formalization and source comparison](parallel26-mca-formalization-2026-09-05.md)
for the exact assumptions, inspected proof, and limits of the new local
Lean evidence. The full Paley and prize targets remain unproved.

## Primary sources inspected

| Source | What was checked | Limitation |
|---|---|---|
| [Satake, *On the restricted isometry property of the Paley matrix*](https://arxiv.org/abs/2011.02907), v2, 2020 | Conjecture 7 and Remark 9, PDF p. 4; following clique consequence on p. 5 | The result about RIP is conditional on Paley, not a proof of Paley |
| [Fouvry, Shparlinski, Xi, *Estimates for trilinear and quadrilinear character sums*](https://doi.org/10.4171/RMI/1530), 2025 | §1.2; Appendix A, Theorem A.1 and Remark A.2; complete-sum enlargement | Structured/multilinear improvements are not arbitrary two-set cancellation at all positive exponents |
| [Official Proximity Prize](https://proximityprize.org/) | Grand MCA and list-decoding challenges, smooth domains, rates `{1/2,1/4,1/8,1/16}`, and target error | The prize statements are not stated there as the Paley conjecture |
| [Arnon, Boneh, Fenzi, *Open Problems in List Decoding and Correlated Agreement*](https://eprint.iacr.org/2026/680) | April 8 PDF: prize statements, MCA, and LD/CA implications. July 6 PDF subsequently retrieved: title and pages 21,31–33,48–50 visually checked | The July audit covers selected definitions, attacks, and list constructions, not every theorem or its dependency proof. The local April archive's June filename is not its publication date |
| [Ben-Sasson, Carmon, Haböck, Kopparty, Saraf, *Proximity gaps stop at the Johnson bound*](https://eccc.weizmann.ac.il/report/2025/169/) | Introduction and §1.4.3; different hypotheses for prime-field limitations | Their additive-subgroup sumset conjecture is not the quadratic two-set Paley conjecture |
| [Aistleitner, Frühwirth, Hauke, Manskova, *Moment generating functions and moderate deviation principles for lacunary sums*](https://arxiv.org/html/2502.20930v2), v2, 2025 | Theorem 1, sampling measure, and doubling example's cumulant expansion | Lebesgue sampling does not supply a worst-case prime-grid bound |
| [Kurlberg, *Bounds on exponential sums over small multiplicative subgroups*](https://people.kth.se/~kurlberg/eprints/short_expsum.pdf), 2007 | Theorem 1.1 on p. 2, including a rendered-page check | Exposition of a fixed power saving for fixed subgroup-size exponent; not the sought square-root bound or a claim of the best current result |
| [Official IRS benchmark source](https://github.com/proximity-prize/proximity-prize/tree/b34c0131cfa36b51111521541d7d3e35c8791082), September 4 snapshot | Exact profile and lower/upper target types; pinned CompPoly field definitions and ArkLib MCA-plus-list expression | Source-level audit, not a fresh kernel audit; the benchmark is one parameter point and is not equivalent to the grand challenges |
| [Di Benedetto et al., arXiv:2003.06165v1](https://arxiv.org/html/2003.06165), 2020 | Theorem 3.1, Lemmas 4.1–4.3, and §5 amplification exponents | Gives a quartic-window saving of 31/2880, not square-root cancellation; no best-current claim |
| [Shkredov, arXiv:1802.09066v2](https://arxiv.org/html/1802.09066v2), 2018 | Theorems 3 and 25, recurrence (54), and norm step (57) | The general-function statement as rendered needs a qualification at zero; a local counterexample does not refute the subgroup statement. Its direct exponent implication remains weak |
| [Kowalski–Untrau, arXiv:2505.22059](https://arxiv.org/html/2505.22059) | Theorem 3.8 and Lemma 3.9; v2 header July 2025, internal date August 2026 | Prime subgroup order `d=o(log p/log log p)` and a distributional conclusion; excludes the present dyadic quartic window |
| [Magsino–Mixon–Parshall, arXiv:1907.05971v1](https://arxiv.org/html/1907.05971), 2019 | Introduction's difference-graph definition and attribution of the Graham–Ringrose lower bound | Corrects the overly strong uniform logarithmic clique target; original Graham–Ringrose proof not independently audited |
| [Hoshi–Kanai, arXiv:2105.14872](https://arxiv.org/html/2105.14872) | Multiplication matrices of Gaussian periods in the introduction and §2.3 | Identifies the established framework; representation and lifting are not the desired prime-field maximum bound |
| [Garcia–Lorenz–Todd, arXiv:2112.13886](https://arxiv.org/html/2112.13886) | Symmetric period matrix, its symmetries, and circular-pair fourth-moment formula in §§2–3 | Circularity holds for fixed order outside finitely many primes, not uniformly in the working quartic window; finite counterexamples already exist here |
| [Shkredov, arXiv:1504.04522v1](https://arxiv.org/html/1504.04522), 2015 | Shifted multiplicative energy bound for subgroup size below √p, introduction and Proposition 3 | Full shifted energy is O(n² log n); subtracting the trivial matchings does not strengthen that to O(n log n) |

The accessed 2025 character-sum paper explicitly describes arbitrary small-set
bilinear cancellation as unproved. No current proof was identified by this
search. A literature search cannot certify that no new proof exists anywhere.

Satake's exact statement uses strict `|A|,|B|>p^ε` and all sufficiently large
primes. The meaningful exponent range here is `0<ε<1`; at `ε=1` the strict
size condition is vacuous. Replacing strict inequalities by non-strict ones
does not change the all-exponents qualitative conjecture, but the notes retain
the source convention when stating the target.

## Relation to prior user research

The prior research summary consulted in Codex memory reported a still-open
proximity-prize frontier involving signed/correlated cancellation, alongside
conditional Lean reductions. That is historical context, not a currently
verified reduction from this Paley formulation. The old checkout path was
not present locally during this run. No old conditional claim has been used
as an unconditional hypothesis or theorem in these notes.

Subsequently, the remote repository was recovered and inspected at
`SlopDotCash/proximityprize` commit `5b00e50c3c51b3c944201a1749a1a8132e5ce167`.
The current files distinguish the thin-subgroup single-period target from
the classical two-set conjecture and still report the main result open.
See `subgroup-target.md` for the exact provenance, the resulting change in
research priority, and a correction to an overbroad claim in an older essay.

A subsequent [reduction audit](prize-reduction-audit.md) found two material
gaps. `FarLineIncidenceEquivariance.lean` proves invariance under code
automorphisms, but `DyadicLacunaryDeltaStar.lean` cites it as proving
monomial-pair extremality. The latter conclusion does not follow, and the
new note refutes its unrestricted form. A prior LD/MCA bridge note also
overstates the radius obtained from ABF Theorem 5.1 for random RS codes:
the April paper gives the Johnson radius through that implication, not
capacity. The archived prior files are retained unchanged as sources.

## Scope of verified results

- The exact Python checks validate finite identities and explicit witnesses;
  they do not validate asymptotic statements by sampling.
- The first Lean certificate proves the explicit quadratic-character
  Gaussian-moment obstruction. Two additional certificates check finite
  subgroup arithmetic and collision counts; their normalization-to-energy
  interpretation is supplied as an ordinary mathematical proof.
- A fourth certificate proves a cosine lower envelope and applies it to the
  actual complex period at `p=6700417`, `|H|=64`, refuting the literal bound
  `M²≤2n ln(p/n)` there. The different prior envelope `M²≤2n ln p` is not
  refuted. Exact integer quotient convolution gives `M²≤1970` for this field;
  that upper computation has not been checked inside Lean.
- The asymptotic obstruction arguments are supplied as mathematical proofs
  in `moments-and-obstructions.md`; they are not yet formalized in Lean.
- A general upper bound with a logarithmic spike allowance remains a proposed
  research hypothesis. Its sufficiency for Paley is proved conditionally.
- `geometric-lift-and-alias.md` proves an auxiliary finite-cycle moment bound
  in ordinary mathematics using independent digits. The required signed
  discrepancy on reducing modulo a prime remains unproved. Exact carry
  computations locate the first additional relation in the tested example;
  they do not establish a uniform discrepancy bound. The discrepancy bound
  is equivalent to the original moment estimate up to constants.
- `prize-reduction-audit.md` proves the elementary general MCA upper bound
  `binom(n,k+1)/q`, and its attainment at the specified radii when
  `q>binom(binom(n,k+1),2)`. For length eight and dimension four, it proves
  a universal monomial-pair upper bound of `40/q`, below the attainable
  global maximum `56/q` in sufficiently large fields. Exact interpolation
  certificates verify attaining pairs at `p=2017` and `p=65537`.
  These results are not claimed as new to the literature or Lean formalized.
- `coset-coherence.md` proves exact period autocorrelations and an arithmetic
  multiplication identity. Its synthetic vectors satisfy the linear
  correlations with large spikes but fail the multiplication identity.
  They are not actual Gauss periods. Integer group-ring checks cover five
  fields; rational calculations verify four synthetic examples.
- `period-polynomial-certificate.md` proves the exact maximizing frequency
  coset at `p=6700417`, `n=64`, with `M=η₁` enclosed to width `10^-18`.
  All signed moments through order 24 are independently reconstructed by
  quotient convolution and matched to the carry counts. Rational polynomial
  positivity certificates exclude all other cosets from absolute values
  at least 42. The new result is not Lean formalized.
- `jacobi-coefficient-route.md` proves a conditional uniform implication
  using orthogonal-polynomial recurrence coefficients. The coefficient
  hypothesis remains unproved uniformly, and its literal B=1 version
  is refuted by an exact quartic-window witness at p=67403009, n=128.
  No novelty claim is made for the method.
- `cyclotomic-prime-average.md` proves a weighted bound on extra relations,
  sharpened by AM–GM. Its pointwise consequence fails to improve the
  elementary fourth-energy bound in the stated regime. At order ten,
  the principal-frequency contribution forces extra relations at every
  eligible prime once dyadic n≥1024. No worst-case centered estimate follows.
- `official-profile-and-trace.md` proves the exact trace pullback for a
  base-field subgroup embedded in an extension, and identifies the full
  trace-zero space as an additional set of full-size periods. Exact
  sextic arithmetic gives a witness in the official profile; the complete
  F_25 example checks the moment identity through order sixteen. These
  are not counterexamples to the prime-field conjecture or prize claims.
- Neither Paley nor either grand prize challenge has been resolved.
- `subset-sums-and-lists.md` gives a self-contained Fourier/coefficient
  proof of a subset-count bound, exact monomial list/MCA identities, and
  a root-lifting injection. It cites Li–Wan (arXiv:0708.2456) and Zhu–Wan
  (arXiv:1101.0289; HTML v1 inspected) for the established subset-sum
  context; no novelty is claimed. Exact checks give a list of at least
  `2^8154` codewords at radius `4095/8192` for the pinned official profile.
  This rules out that radius in its MCA-plus-list lower certificate,
  without lower-bounding actual winning-set soundness. A nonzero
  cyclotomic norm prevents this construction from extending to fixed-gap
  asymptotics. The new results are not Lean formalized.
- `riesz-tail-route.md` gives a direct dissociation entropy proof in the
  established Chang/Riesz framework; primary context is James R. Lee,
  arXiv:1508.07109v2, whose HTML was inspected. It identifies the loss
  from using one multiplicative orbit, derives an exact centered
  Riesz-product sufficient criterion, and proves that criterion is
  equivalent to the target up to constants. Exact signed-relation counts
  refute the stronger normalization Z_+(1/4)≤1 at p=67403009, n=128.
  No uniform bound, novelty, or Lean formalization is claimed.

- `analytic-bounds-and-amplification.md` records an existing uniform power
  saving in the specified window and proves an algebraic ceiling for the
  saving supplied by a stated conditional amplification ledger. Exact
  rational checks cover the complete reduced parameter range and a
  mass-at-zero test of the broad function statement in the inspected
  2018 HTML. The subgroup-centered formula is not refuted; its implications
  are analyzed conditionally. No new cancellation theorem or Lean proof
  is claimed. Four HTML snapshots and checksums are in
  `sources/analytic-bounds-2026-09-04/manifest.json`.

- `mixed-periods-and-shifted-energy.md` derives the full mixed system in
  the established multiplication-matrix framework and proves a bijection
  between its triple intersection counts and nontrivial shifted product
  collisions. The resulting fourth-energy comparison is unconditional,
  but the needed upper bound on its right side is not. Exact checks
  certify circularity and minimum fourth energy at two larger finite
  quartic-window examples. Neither a uniform estimate nor novelty or
  Lean formalization is claimed. Three primary HTML snapshots have
  checksums in `sources/mixed-periods-2026-09-04/manifest.json`.

- `kernel-discriminant.md` constructs an integer polynomial whose
  discriminant times boundary value detects every fourth-energy exception
  at a fixed dyadic order. Exact complete factorizations through n=32
  agree with the earlier independent norm enumeration where both apply.
  Prime-power pair counts verify an exact valuation identity, including
  a nonzero second-level contribution at n=16,p=17. This is a
  self-contained mathematical argument, with no new external theorem
  dependency, novelty claim, or Lean formalization.
- `quadruple-orbits-and-cube.md` proves a uniform orbit identity and
  factorization `odd(A_n)=odd(3^n−1)U_n³V_n⁶`. It treats all odd primes
  using unramified extensions and counts the three multiset types
  explicitly. Direct exact checks cover 55 cases and 109 precision
  levels, including five non-splitting cases. The calculation of A_64
  is not a complete factorization or prime classification. A literature
  screen inspected the indexed introduction of Helou's 1997 primary
  paper on Wendt's different determinant; no theorem from it is used
  as a dependency and no full-paper audit is claimed. The power
  factorization does not provide a uniform energy or spectral bound.

- `dyadic-descent-and-mixed-energy.md` proves an exact two-coset energy
  recurrence, a conditional induction with explicit constant 22, and
  a descent theorem when the field of definition doubles. Its polynomial
  factors count balanced new orbits, and its prime-support restrictions
  hold uniformly by unramified descent. These are self-contained proofs;
  they do not establish the mixed-energy bound or the Paley target.
  Exact tests cover 44 field/precision cases, ten with doubled field
  degree, and factor identities through order 64. Complete prime
  factorizations are used only through order 32. The publisher's abstract
  of Do Duc–Leung–Schmidt (2020), DOI 10.5802/alco.86, was inspected and
  archived under sources/dyadic-descent-2026-09-04/ with its checksum.
  Its reported exponential threshold does not apply to the growing-order
  target. No full-paper audit, theorem dependency, novelty claim, or
  Lean formalization is asserted.

- `positive-product-moments.md` proves an L^q recurrence retaining only
  positive products and an explicit sufficient implication to (SG) at
  one even logarithmic depth. Character orthogonality gives an exact
  stronger criterion involving centered balanced counts. A principal-term
  lower bound forces extra balanced twelve-term relations throughout
  the quartic regime once the parent order is at least 128. These are
  self-contained ordinary proofs, not new external theorem dependencies
  or a uniform cancellation result. Exact checks cover 80 signed product
  moments across nine computed levels, with dense checks in two small
  fields and a q=12 whole-tower criterion check at p=6700417,N=64.
  All arithmetic acceptance conditions are integers or fractions, and
  int64 convolutions have explicit pre-step overflow bounds. No novelty
  or Lean-formalization claim is made.


- `signed-quotient-operators.md` gives self-contained proofs of a signed
  symmetric integer representation, its exact spectrum and trace, and a
  quantitative loss when entrywise signs are discarded. It also proves
  diagonal sign conjugation under representative changes and uniformity of
  nonidentity endpoint-ratio walk counts. None of these identities supplies
  a uniform signed cancellation bound or closes the Paley/prize gap.
  Exact checks cover 14 matrices, 94 traces, 46 independent balanced counts,
  and 32 endpoint-ratio identities. Sparse integer products have explicit
  int64 row bounds, with arbitrary-precision accumulation and independent
  dense additive convolution. No novelty or Lean-formalization claim is
  made. A literature screen included the primary abstract of
  [Muzychuk–Ponomarenko, On Pseudocyclic Association Schemes](https://arxiv.org/abs/0910.0682)
  and the previously inspected Hoshi–Kanai multiplication-matrix paper.
  No theorem from this screen is adopted as an additional dependency.


- `list-to-winning-set.md` proves a single-word list projection lemma for
  the pinned fixed encoder, including a common violating instance over an
  entire admissible radius suffix. A base-field coefficient-pigeonhole list
  gives the ordinary mathematical unsafe-radius bound 122641/262144 for
  the pinned IRS profile and the exact score inequality at 11649 centibits.
  The 55-page July 6 ABF PDF is now archived separately with its hash;
  eight relevant pages were visually inspected. The proof corrects the
  local sign inconsistency in the displayed word construction in Lemma C.5
  rather than treating it as a disproof of that lemma. Existing ABF list
  projection and common-coefficient methods are acknowledged; no novelty
  or Lean-submission claim is made. Four finite cases exhaust 86086 linear
  projections and all codewords, including two with list size at least the
  field size. Production inequalities are checked with exact integers.
  Three additional pinned ArkLib definition/source files and the PDF have
  provenance in `sources/prize-attack-2026-09-04/manifest.json`.
- The screen also inspected Goyal–Guruswami–Sun–Wootters,
  [arXiv:2607.08516v1](https://arxiv.org/html/2607.08516v1), especially
  Theorem 5.6. Its improved results concern Reed–Solomon codes with random
  evaluation points and do not establish a guarantee for the fixed smooth
  domain here. The HTML is archived; its full proof chain has not been
  independently audited and no theorem from it is used as a new dependency.

- `localized-necklace-identities.md` uses the necklace framework from
  [Kunisky, arXiv:2303.16475v1](https://arxiv.org/html/2303.16475v1),
  Definitions 1.13 and Conjecture 1.14, while retaining the distinction
  between weak spectral convergence and spectral-edge control in Section 7.
  The note proves exact one- and two-occurrence formulas independently by
  matrix correlations. The analytic input is
  [Lu–Zheng–Zheng, arXiv:1305.3405v3](https://arxiv.org/html/1305.3405v3),
  Lemma 2.1 (2.3), specialized to tensor exponents (1,1), with R=1 from
  Remark 2.4. This gives `|t_j|≤(j−1)p^((j−1)/2)` after a checked Gauss
  expansion. The original lemma and its proof were read in HTML; the deep
  theorem is cited rather than claimed as independently reproved. Both
  versioned HTML files are archived and hashed. This turn used no new PDF.
  One local source correction is needed if using the displayed Jacobi
  table in Kunisky Proposition 3.10: under the stated zero-at-zero convention,
  `J(epsilon,chi)=-1`, not the displayed 0. The independent matrix proof
  explicitly keeps the trivial and quadratic directions, both with
  eigenvalue -1; it does not use that erroneous entry. This does not refute
  the separate Kloosterman-based estimate.

- `planar-necklace-reductions.md` proves the affine average of matrix powers,
  planar Fourier duality, and the even-degree/even-edge projective invariance
  used in the new reductions. It uses the already checked Lu–Zheng–Zheng
  bound as its only deep analytic input. The two graph duals are constructed
  from oriented face boundaries; literal sums independently check their
  partition functions and projective boundary contributions. These give
  arbitrary two-block identities and all binary length-six formulas, not
  the full degree-two necklace conjecture or the spectral-edge conclusion.
  No new theorem from the accompanying literature screen is adopted, and
  no novelty, Lean, or prize-submission claim is made.

Downloaded PDFs and extracted text are in `sources/`; the manifest records
original URLs, SHA-256 hashes, the earlier failed ABF download, and the
subsequent successful July-edition retrieval.
An older ABF copy was recovered locally with its provenance and hash.
Source page images
were checked for the Paley quantifiers and the Karatsuba estimate because
text extraction omitted some absolute-value signs and exponents.


## Second parallel pass: Katz Mellin input and exact spectral transfer

The primary PDF https://web.math.princeton.edu/~nmk/mellin186.pdf is archived
as sources/katz-finite-field-mellin.pdf (661647 bytes; SHA256
08e36a0fdeb8a1dc3d8c95734f44d7569fab23cd45fb1a5da8bea5576e3e1b73).
Printed pages 22, 52, and 53 were rendered and visually inspected. The
web text also supplied the good-character criterion on pages 9-10/17.
Corollary 4.2 and Theorems 15.1/15.3 give the adopted Mellin bounds for
TwLeg(1)[1] and Sym^2(Leg)(3/2)[1], respectively. Their Mellin dimensions
are two; the local monodromy at zero and infinity is unipotent. The
nontrivial characters are good. The exceptional character is handled by
exact matrix identities in research/parallel2-necklace-2026-09-04.md.
The symmetric-square stalk at t=1 requires the retained +p delta_1 term.
The analytic source is imported, not reproved or formalized here.

Kunisky's archived Section 7 was reread for its distinction between
individual necklace estimates and weighted aggregate control. The new
research/parallel2-spectral-transfer-2026-09-04.md derives its finite
projection/Chebyshev and Rayleigh identities directly, keeping every
normalization, rank term and anchor correction. It does not assert the
missing aggregate bound. The resulting conditional clique conclusion is
weaker than the full arbitrary-two-set Paley formulation and supplies no
new prize reduction.

## Third parallel pass: hypergeometric ranks and retained tensor stalks

The new primary input is Katz, G2 and Hypergeometric Sheaves,
https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf, Section 2.
It is archived as sources/katz-g2-hypergeometric.pdf, 479161 bytes, SHA256
0bf485bcc9dde2afebc8268af679566bcff65aa6d6ed8c9f4b0387d7e9e954c1.
Printed pages 3-5 were rendered and independently inspected. The source
permits repeated characters and arbitrary equal list lengths, assuming
disjointness of the two lists. It gives the adopted irreducibility, rank,
weight, trace including the singular stalk, and tame local monodromy.
No rank less than the characteristic restriction is imposed in this case.

The supporting curve machinery is Katz, Gauss Sums, Kloosterman Sums,
and Monodromy Groups, https://web.math.princeton.edu/~nmk/Katz-GKM.pdf,
archived as sources/katz-gauss-kloosterman-monodromy.pdf, 6307344 bytes,
SHA256 8711424f8edb14f38c0e61606ef49d6efc67d5a93ac9f5bd9fa3e13b7c7d9ffe.
Printed pages 32-33 and 38-41 were rendered and independently inspected.
Sections 2.3.1-2.3.2 and the proof of Section 3.6 supply Euler-Poincare,
the trace formula and the upper weight bound on compactly supported
cohomology. Text extractions and both PDF hashes are recorded in the
source manifest. These standard deep results are imported, not reproved.

research/parallel3-necklace-2026-09-04.md independently derives the
Gauss-factor trace identity, the fixed constant twist, the rank-mismatch
vanishing, and the actual tensor-stalk quotient. This retains dimension
j-1 at t=1 instead of replacing the tensor by its middle extension.
The resulting H_c^1 dimension is j+1 for j!=2. The equal-rank exception
uses the previous exactly corrected symmetric-square operator.

The subgroup lane also checked Do Duc, Leung and Schmidt, Main Theorem 1,
https://arxiv.org/html/1903.07314. Its constant cyclotomic-number bound
requires characteristic exponential in the subgroup order here, outside
the quartic window for the growing dyadic orders. It is not used to
claim a new uniform bound. The classical lane uses the previously
checked rank-two Kloosterman bound in Lu-Zheng-Zheng, Section 2; its new
Fourier identity and interval obstruction are derived directly.

## Fourth parallel pass: arithmetic duality and the growing prime class

The kernel and necklace arguments reuse the archived Katz G2 Section 2
and Katz-GKM trace/weight inputs. An independent reviewer rendered and
visually inspected Katz G2 printed pages 4-5 again because extracted
text drops conjugation bars. In the equal trivial/quadratic lists,
the additive-character switch has trivial translation and constant
twist. With the checked Gauss normalization this gives the arithmetic
self-duality F_r^vee=F_r(r-1). The identity endomorphism supplies the
one-dimensional quotient between the tensor of actual invariant
stalks and the invariant stalk of the tensor. Nontrivial Kummer
self-twists are excluded by inertia at zero. These details underlie
research/parallel4-kernel-aggregate-2026-09-04.md and were independently
audited rather than inferred from small-field tests.

The new necklace family derives RawHyp[(chi psi,1);(chi,psi)](u)
=p I_psi(u) for every nontrivial psi, including u=1, and over every
finite extension with its cardinality in place of p. The lists are
disjoint for precisely those characters. The trivial character is
evaluated separately. The rank-two equality exception has no invariant
because of its zero-monodromy characters. The literal raw sums and
zero-coordinate corrections are also checked in exact arithmetic.

The primary Thorner-Zaman source is archived as
sources/thorner-zaman-2108.10878v2.html, from
https://arxiv.org/html/2108.10878v2, 428325 bytes, SHA256
ecc1ba8e04d9ce868f150e6542e312152b58241b2e32d2493b6f93c93edd938d.
Corollary 3.1, Section 3.1, and equation (3.2) were read in the
versioned HTML. For q=N=2^m, rad(q)=2, a=1, x=N^4, and h=3N^4/4,
fix epsilon_0=1/12. The exceptional zero is absent for all sufficiently
large N by the source's small-radical criterion, and h/phi(N)
>=x^(7/12+epsilon_0). Its effective error tends to zero, yielding
|P_N|~3N^3/(8 log N). Constants are absolute along these dyadic
moduli; no numerical onset is claimed. Both the primary agent and an
independent reviewer checked the specialization and the log(2)/36
proportion. The source theorem is imported, not reproved.

The primary Alsetri-Shao source is archived as
sources/alsetri-shao-2509.07765v1.html, from
https://arxiv.org/html/2509.07765v1, 371582 bytes, SHA256
4b7ad66cccedcdaef2dbd4740948562c3dc645e74ac2a97afd9f52a30f08203a.
Its introduction states the classical Burgess bound; Theorem 1.1
supplies the proper rank-two progression bound. Both were read and
their uniform translation/dilation scope checked. The derived
directional estimate uses these as established inputs and does not
assign an unsupported numerical value to delta(epsilon).

The two versioned HTML files are pinned in sources/manifest.json.
The fourth-pass integration audit hashes the new files and all
recorded local proof inputs. No new deep theorem is claimed to have
been independently reproved by the finite verifiers.

## Fifth parallel pass: direct twists and a checked amplification theorem

The new anchor and necklace arguments reuse the archived Katz Section 2
and GKM cohomology inputs. Printed page 5 of the G2 PDF was rendered
and visually checked again for the tame infinity characters. The
rank-j kernels have only the quadratic character there. This excludes
the invariant in both new arguments: a single translated anchor
contradicts a possible equal-rank self-twist, and the three-factor
Mobius tensor has three quadratic inertia factors at infinity while
the pulled rank-two factor is regular. Actual stalks, omitted fibers,
and all field-zero values are treated explicitly in the proof notes.
An independent reviewer checked the anchor proof and its adjacency-cell
correction; the primary agent checked the full necklace proof.

The subgroup proof extends the previously proved normalized-tuple norm
budget to six terms, using T6(s)=15s^3-45s^2+40s and the checked
Thorner-Zaman quartic prime count. That yields an unconditional
density-one class with near-quadratic fourth and near-cubic sixth
energy. The analytic maximum consequence uses the fixed selection
orders from Di Benedetto et al., arXiv:2003.06165v1, Section 5 and
Lemma 4.1. Its existing archive is
sources/analytic-bounds-2026-09-04/di-benedetto-et-al-2003.06165v1.html,
287859 bytes, SHA256
1df3718ccfe5d44f2844cfe09e641fe48deec2c5a2f6b0f6978f98f8c37031f8.
The actual proof was read; the new note retains A and B in E2<=An^2,
E3<=Bn^3. Its dyadic thresholds lose constants, while the three
support estimates contribute log^5 after multiplication. The
three zero deletions and phase weights were independently audited.

The original weighted trilinear input was also checked in
Petridis-Shparlinski, arXiv:1604.08469v4, Theorem 1.1 and its definitions,
https://arxiv.org/html/1604.08469v4#S1.SS3. The source orders the three
set sizes. Permuting sets and their factorwise weights gives the
stronger expression p^(1/4)(XYZ)^(3/4)min(X,Y,Z)^(1/8), which implies
the labeled bound used here. A nonzero frequency is absorbed into
one set. There is no p^(2/3) condition in this trilinear theorem;
the nearby quadrilinear theorem has that separate restriction.
The primary HTML is archived as
sources/petridis-shparlinski-1604.08469v4.html, 652698 bytes, SHA256
94d8be554227b46460f64ac07c348fb8791181ef60801c5bbc35556ca3726168.

Two independent proof reviews accepted the subgroup normalization,
all-level exceptional proportion, selection/deletion conditions,
Delta^72 and log^5 powers, and the final 23/24, log^(7/72) bound.
These are ordinary source-dependent mathematical arguments, not a
claim of formal verification or literature novelty. The exact
high-moment target and the exceptional primes remain uncontrolled.

The classical single-size reduction and lower-additive-energy
obstruction are proved with subset averaging, a digit/power-sum
construction, Newton identities and an elementary row double count.
The construction includes its polynomial prime-size control and
does not depend on an unquantified prime-in-progression assertion.

## Sixth parallel pass: elementary amplification and translated kernels

The stronger subgroup consequence, exponent 17/18 on the sixth-energy
class, uses an elementary weighted bilinear inequality proved directly
by additive orthogonality. No new external analytic theorem is used.
The all-level sixth-energy prime budget and the Thorner-Zaman prime
count are inherited from passes five and four with their existing
archives. The source link is pinned to arXiv:2108.10878v2. Two proof
reviews independently checked the subgroup invariance identity,
Jensen exponent, origin deletion, negative nonzero-coordinate constant,
and the explicit constants in the 17/18 and conditional 15/16 bounds.
The higher centered eighth-energy budget remains an unproved input.

The multiple-seed proof reuses Katz G2 Section 2 and Katz-GKM curve
cohomology. Translating the singularity from 1 to b excludes a
distinct-seed invariant at b; a nonempty anchor set disjoint from
both seeds instead supplies scalar quadratic inertia. Restoring the
actual tensor stalks gives dimensions mrs+r+s for distinct seeds
and mrs+r+s-1 for coincident seeds. The equal-seed diagonal uses
the prior arithmetic-duality correction translated to b. Independent
review checked these dimensions, the complex Gram row/column bounds,
the optional zero-anchor twist and the full-space limitation.

The new necklace operator is an exact principal compression of an
augmented single-anchor operator. Its only imported analytic norm
input is the already checked TwLeg Mellin theorem from Katz,
Sato-Tate theorems for finite-field Mellin transforms, Theorem 15.1
and Corollary 4.2. The archive remains
sources/katz-finite-field-mellin.pdf, 661647 bytes, SHA256
08e36a0fdeb8a1dc3d8c95734f44d7569fab23cd45fb1a5da8bea5576e3e1b73.
The previous proof's good-character conditions and exceptional
trivial mode are retained. Root checked the new projective pole,
three added eigenmodes, graph normalization and coincident-anchor
term independently; no new purity assertion is needed.

The classical fixed-size mean/variance formulas use only elementary
two-point character identities, set partitions and inclusion
probabilities. The companion symbolic calculation clears denominators
and proves zero polynomial residuals for both signs of chi(-1).
The typical-set corollary uses Chebyshev. It is not a claimed
uniform moment theorem, literature-novelty result or Paley proof.


## Seventh parallel pass: signed multipliers and middle convolution

The subgroup multiplier theorem uses only elementary Fourier inversion,
Jensen against a nonnegative sum distribution, and bilinear orthogonality.
Its density-one consequence inherits the fifth-pass sixth-energy budget
and the archived Thorner-Zaman quartic prime count. The subgroup worker
independently accepted the root proof before its account limit; root
checked the signed-multiplier obstruction and continuous-order argument.
Finite rational-interval certificates also include complex periods.

The necklace crossing theorem imports Katz, Rigid Local Systems,
corrected manuscript, archived as sources/katz-rigid-local-systems.pdf,
from https://web.math.princeton.edu/~nmk/wholebookRLScorr.pdf,
1,115,958 bytes, SHA256
ca6eb8d5d9e21076dda4e154e83dfa2821f586d6ccbe8c8f1b371b4352f753d1.
Root read the relevant complete sections and visually inspected full
PDF pages 63, 67, 104, 107, 109 and 139. Lemma 2.9.4 supplies the
constant infinity correction; Theorems 2.9.7 and 3.3.3 preserve the
applicable irreducible objects; Corollaries 3.3.6 and 3.3.7 supply tame
local monodromy and rank; Section 5.5.5.10 identifies the pure middle
image of weight w+1. These are one-based PDF page indices.

The two inputs are irreducible rank-two middle extensions, so neither
is an excluded punctual, rank-one Kummer or Artin-Schreier object.
Finite invariant codimensions and the twisted infinity invariants give
rank four on both sides. Their common local type at the crossing is
quadratic plus three trivial characters. Twisting one side reverses
these multiplicities, excluding the geometric invariant. Root checked
all zero masks, the sign and constancy of the infinity correction,
coincident-anchor terms, exceptional fibers and label-transfer errors.
The optional exact infinity trace uses a split semistable line-conic
fiber with two rational intersection points; the uniform estimate
needs only its weight bound. No rank<p restriction was assumed.

The primary proof is supplemented by exact original matrix sums and
48 traces over F5, F25, F125 and F625. The resulting 12 degree-four
polynomials have exactly certified root modulus 5; this supports the
arithmetic normalization but does not numerically prove irreducibility
or the uniform theorem. All source hashes are pinned in the manifest.
The classical countermodel argument is elementary and retains its
artificial-function scope. No general prize reduction is supplied.


## Eighth pass: a second middle convolution

The new proof reuses the archived Katz Rigid Local Systems input
and the curve-cohomology estimates already audited. The primary
corrected PDF URL was rechecked. No new outside theorem is imported.
Root checked Corollaries 3.3.6 and 3.3.7 against the full local-type
tables for both singleton twists, including the intermediate
infinity representation 1+chi^3. The next finite codimensions are
4,2,1, while the twisted infinity invariants have dimension one.
This gives rank six in all four ordered anchor pairs. The inputs
remain irreducible of rank at least two, outside the excluded
punctual, Kummer and Artin-Schreier objects. Purity raises weight
two to three in the second middle image.

The rank-two/rank-six tensor has no invariant. Its rank is twelve,
with four tame punctures, giving H_c^1 dimension 24 and a p^(5/2)
inner bound. Two constant infinity corrections and the finite
spike are kept symbolically. Only their weight bounds are needed;
the verifier permits every complex correction within those bounds
and does not guess the second scalar. Two final necklace reductions
are direct traces; the nested reduction inherits the fully checked
adjacent-C graph normalization with its coincident-anchor term.

Root developed and audited this extension locally. No second worker
review is claimed: the parallel agents remain stopped by the account
usage limit. Exact computations check Jordan-block bookkeeping,
original character matrices, all 42 new words in five fields and
literal sums at p=5. They supplement the source-dependent proof;
they are not a formal verification of cohomology or of the full
Paley conjecture.


## Ninth pass: exact compact convolution and a strict weight gap

The new all-degree necklace proof reuses the Katz Rigid Local Systems
archive and the prior curve-cohomology source. In addition to the
previously inspected local-monodromy, irreducibility, inversion,
infinity-kernel and purity pages, root visually inspected the complete
PDF pages 45 and 50. Sections 2.3.3 and 2.6.1 give the exactness
needed for the raw perverse chain. Only compact convolution with
the fixed property-P Kummer object is treated as exact; middle
convolution itself is not assumed exact on arbitrary objects.

The additional source is Beilinson-Bernstein-Deligne, Faisceaux pervers,
from the primary author archive:
https://www.math.tau.ac.il/~bernstei/Publication_list/publication_texts/faisceaux_pervers.pdf
Archived at sources/bbd-faisceaux-pervers.pdf, 6,284,485 bytes,
SHA256 6e976a78e7607c3026d006d319eea06df44ef75f80f35d9e475b3a8fbc1288d9.
Complete PDF spreads 65,66,67,69,70,71 were rendered and visually
inspected, corresponding to printed pages 124-129 and 132-137.
Sections 5.1.8-5.1.9 fix upper weights by ordinary stalk cohomology;
5.1.13-5.1.14 give compact-image and product bounds; 5.3.1-5.3.2
give stability under perverse subquotients and middle-extension purity.

The constructed error is perverse of weights at most i, while
the principal middle extension has perverse weight i+1. On the
lisse open these correspond to usual weights at most i-1 and i.
Finite punctual kernels introduced by newly trivial monodromy
have usual weights at most i-1 before convolution; the infinity
kernel is a constant sheaf shifted by one. Both remain below the
new principal weight. The concrete {alpha},{beta},{alpha} sequence
has a new one-dimensional finite invariant, so the raw chain
cannot be replaced by the pure middle chain alone.

Root checked arbitrary-anchor rank growth, both moving-point parity
types, the geometric inverse, exceptional constant/Kummer/punctual
constituents, additive mass bounds, finite stalk dimensions and the
final H_c^2 exclusion. No restriction rank<p is needed. The primary
Kunisky text was rechecked for Definition 1.13, Conjecture 1.14 and
Theorem 1.17. Its weak-convergence implication is explicitly imported;
no extreme-eigenvalue or literature-priority claim is made.

The exact checker verifies bounded local-type enumeration and original
character sums, not the imported sheaf theory. The proof has no
independent worker review: the current agent-status check confirms
all three workers remain stopped by the account usage limit.


## Tenth pass: word recovery, duality and simultaneous correlations

The proof imports the ninth-pass rank chamber and weight-gap theorem.
The additional primary-source check is Katz, Rigid Local Systems,
Section 2.5.1, complete PDF page 49. Root rendered and visually
inspected the full page; it states the duality exchanging compact
and ordinary convolutions, compatibly with arithmetic Galois action.
Together with the middle-image definition and duality of intermediate
extension, this gives F_w^dual=F_w(|w|). The scalar invariant in
the diagonal tensor consequently has the exact H_c^2 eigenvalue
p^(|w|+1); no guessed arithmetic normalization is used.

Geometric inequivalence is proved by recovering the last nonempty
label after inverse convolution and repeating. It does not assume
that local types classify arbitrary sheaves. For the raw generic
rank recurrence, finite conductors are additive stalk Euler
characteristics. The proof checks the constant, Kummer and punctual
exceptions as well as ordinary middle extensions. It does not treat
ordinary stalk dimensions as separately additive.

The full-word Gram and fresh-anchor bounds use the already audited
tame Euler characteristic and weight estimates. Root checked all
Tate shifts, diagonal invariants, zero masks and half-valued cell
boundary corrections. The explicit raw-rank second-moment recurrence
is derived by summing subset updates; its growth rate is the Perron
root of the displayed 2-by-2 matrix.

The finite checker uses integer circular convolution with an explicit
overflow guard, arbitrary-precision integers where needed, rational
upper square roots, and exact fraction-free positive-definiteness
certificates. Its numerical constituent inventory records only local
types, not an unjustified classification of sheaves by those types.
The separate null-vector certificate uses only original integer
character matrices and exact rational elimination. It establishes a
dimension obstruction to unrestricted long-depth isometry, not a
counterexample to the spectral aggregate or Paley conjecture.
No independent mathematical review or formal verification is claimed.


## Eleventh pass: the full spectral energy, including corrections

The new scalar energy recurrence, curvature identity and two-sided
trace/energy comparisons use only finite-dimensional linear algebra
and the previously checked Paley projection. The logarithmic-depth
bound imports the ninth- and tenth-pass sheaf arguments. No additional
external theorem is imported; the versioned Kunisky primary text's
discussion of the long-depth spectral requirement was rechecked.

The new finite-point estimate uses the already established raw
conductor c_v and weight gap. Ordinary H^(-1) stalk dimension is
at most the generic raw rank r_w, by the long exact sequence and
simple-constituent dimensions. The conductor identity then gives
dim H^0≤c_v≤r_w at the finite singularities. Both ordinary trace
weights are at most j because the principal object has no H^0 and
the error has perverse weights at most j. Thus the actual path is
bounded by 2r_w p^(j/2), including at the moving point and anchors.

The proof keeps the specific all-one word coefficient vector in the
soft-mask expansion, restores the constant direction annihilated by S,
and handles the J and anchor terms inside the full powers through a
Frobenius perturbation bound. Root checked all p and N normalizations,
the exceptional columns, both trace/energy inequalities, and the
logarithmic range. The positive energy criterion is explicitly shown
to retain the original long-depth difficulty up to controlled factors.

The checker uses exact integer object matrices, rational pairs u+v√p,
algebraic sign comparisons and rational upper square roots. No floating
eigenvalue approximation is used for acceptance. The independent
two-dimensional Jordan case is an abstract projection block, explicitly
not a Paley graph. Finite verification does not prove the remaining
long-depth row-energy estimate or any prize reduction.


## Twelfth pass: an elementary countermodel and exact finite outlier

No new cohomological source input is asserted. The rank-two projection
construction, preservation of anchor columns, trace/norm perturbation
bounds, and planted eigenvalue are proved directly by linear algebra in
research/parallel12-bootstrap-obstruction-2026-09-05.md. The inference
that short-depth control survives takes the previous source-dependent
estimate as an input; independent review of that earlier proof remains
outstanding. None of its sheaf descriptions is claimed for the modified
matrix away from anchors.

Kunisky's archived v1 HTML was reopened on the official arXiv page.
Section 7 states an asymptotic spectral criterion with depth much larger
than log dimension. This supports the distinction between that target
and a zero-error statement at each finite prime. No new PDF pages were
used or implicitly treated as newly reviewed in this pass.

The exact p=257 Rayleigh certificate uses only Euler's criterion,
integer quadratic forms and positive integer comparisons. The associated
transfer spectral-radius conclusion uses the already written elementary
two-projection calculation. Numerical eigenvectors were used only for
discovery; the retained small integer vector is independently accepted by
exact arithmetic. The modified projections are not Paley matrices, and
the actual finite outlier is not an asymptotic counterexample. The final
results and artifact audit pin the inputs and preserve pass eleven.


## Thirteenth pass: principal rank sums and error inertia

Root reopened the primary Katz Rigid Local Systems PDF and visually
rechecked complete PDF pages 49,50,63,64,67,107,109,139. The review
covers arithmetic convolution duality, exact compact convolution with
a property-P Kummer object, the constant infinity kernel, applicable
geometric inversion, tame local quotient/rank rules, and purity of the
middle image. This is a targeted root review, not independent validation
of all preceding arguments. No new BBD page was inspected in this pass.

The new error-inertia observation uses the already proved additive
finite conductor: c_y(K_w)=c_y(F_w[1])=1 implies c_y(E_w)=0.
Nonnegative simple-constituent conductors show that the error is lisse
at y. Its tensor with the geometrically irreducible principal constituent
has no invariant or coinvariant, giving a cross-correlation bound via
the existing tame Euler-characteristic and weight inputs. This argument
does not guess arithmetic values for boundary constituents.

The summed principal-rank recurrence and positive weighted count are
new algebraic derivations from the local rules. Coefficient identities
reduce to lambda^2-h(a+1/a)lambda+N=0. Their proof applies at arbitrary
depth; finite enumeration and quadratic-field arithmetic are additional
checks, not a replacement for that derivation or for the sheaf inputs.
The complete-energy improvement restores every term retained in pass
eleven and remains restricted to a bounded multiple of log p. Sources
and previous artifacts are pinned by the final audit.

## Fourteenth pass: elliptic halving and Lang-character twists

The primary versioned HTML of Bekker–Zarhin, arXiv:1702.02255v2,
was reviewed and archived as sources/bekker-zarhin-1702.02255v2.html.
Theorem 2.1 supplies the rational-halving criterion on a split elliptic
curve in characteristic other than 2; the three differences being
squares also handles rational halves of 2-torsion. This is an imported
standard criterion, not a novelty claim. No PDF from this source was
used. The archived HTML is 1077009 bytes, SHA256
49a902d503a66c7702e599a51322e362d13fbd3fb8a721aefb7573a742686096.
Its manifest entry was appended without changing the preceding entries.

Root rendered and visually reviewed the complete Katz GKM PDF spreads
21,25,34,35: printed pp.32–33,40–41,58–61. The relevant inputs are
Euler–Poincare and the trace formula, the weight bound on H_c^1, the
Lang isogeny and its character-sheaf trace convention for any smooth
connected commutative algebraic group, and tame Kummer sheaves.
The Lang statement therefore applies to the elliptic curve; it is not
restricted to the additive or multiplicative group examples afterward.

For the quadratic Kummer sheaf of y on E minus E[2], the four odd
orders in div(y) give four nontrivial tame local monodromies. Tensoring
with a Lang-character sheaf leaves these unchanged and introduces no
puncture or Swan conductor. Hence H_c^2=0 and dim H_c^1=4. Its
trace bound is 4sqrt(p); dividing the four-to-one rational G[2]
quotient gives the stated scalar bound sqrt(p). This application was
checked by root but has not received independent mathematical review.

The exact kernel, its exceptional block, centered compression, Klein
four projectors and Fourier-basis formula are elementary derivations.
The retained finite tests use integers and cyclotomic polynomial
reduction. Off-diagonal witnesses explicitly prevent reading the
four-product formula as a diagonalization. No theorem about additive
quadratic-phase Weil gates has been assumed to control the different
quadratic-character multipliers here. The norm and full goal remain
unproved. The final audit pins the new inputs and preserves pass thirteen.


## Fifteenth pass: parity of translated divisors and a modified sign kernel

The current correlation argument reuses the unchanged primary Katz GKM
inputs reviewed in pass fourteen: complete PDF spreads 21,25,34,35,
printed32–33,40–41,58–61. No new PDF page was rendered or represented
as newly reviewed in this pass. The source hash is checked again.
For s distinct quotient shifts, their E[2] cosets are disjoint, giving
exactly 4s odd-order punctures. The Lang twist stays unramified on E,
so it cannot cancel that local monodromy. Tame Euler characteristic
and the weight bound give the sum estimate. Positive even factors are
removed by an exact finite-function identity retaining their zeros;
no stalk value is guessed. This application needs independent review.

The archived primary Bekker–Zarhin v2 HTML was reopened at Section4,
Theorem4.1, to verify Hasse's interval and the resulting size bound
L>=16ceil(sqrt(p)) for p>=10000. The ordinary two-factor structure of
the finite elliptic group follows from prime-to-p torsion rank two;
p cannot divide its order here, since the order is divisible by16
and lies below2p. The rectangle construction thereafter uses finite
abelian groups, sumsets, L1 telescoping and an explicit Rayleigh vector.
The prior source-manifest entries and source files remain unchanged.

The modified function need not retain the original four-puncture
Kummer realization. Its bound has a larger fixed constant50, so it
does not refute the original constants1 ands. It also fails the ambient
Paley spectral identity: the planted clique with the three anchors
has Rayleigh quotient greater than sqrt(p). Accordingly this is an
obstruction for the orders of abstract sign-kernel correlation bounds,
not a Paley or prize counterexample. It preserves a different set of
inputs from the pass-twelve projection modification. The final audit
keeps that distinction and preserves all five pass-fourteen artifacts.

An additional exact audit recomputes the original constant Fourier
coefficient directly as one quarter of the full curve sum. Since only
negative signs are changed, the modified coefficient adds twice the
changed-sign count. At p10009 this is 3+384=387, with square margin
139760 above p; at p65537 it is −1+2828=2827, with margin7926392.
This directly certifies loss of the original sharp scalar bound and
prevents presenting these examples as retaining its constant1.

## Sixteenth pass: the field of definition is essential

The primary Kunisky versioned HTML, arXiv:2303.16475v1, was reopened
at Section 2. Its asymptotic parameter is explicitly restricted to
primes p congruent to 1 modulo 4. The archived HTML is unchanged and
its hash is verified again. Thus square-order examples are not described
as counterexamples to its prime-order conjecture.

The primary Asgarli–Yip arXiv abstract, 2110.07176, was checked only
for background attribution of the classical subfield-clique phenomenon.
No PDF of that paper was opened, and no classification theorem is
imported. The subfield clique and exact eigenvector in the new note
are proved directly by elementary finite-field character identities.

The same unchanged Katz GKM PDF supports the all-order twisted sums.
Its previously reviewed complete spreads 21,25,34,35 (printed32–33,
40–41,58–61) give the finite-field trace and weight bounds, Lang
isogeny and tame Kummer inputs. No new PDF rendering or independent
review is claimed in this pass. Replacing Frob_p with Frob_q leaves
the 4s nontrivial tame punctures and quotient factor four unchanged.
This source-dependent application still requires independent review.

Actual quadratic-field Paley matrices retain the sharp constants and
ambient identity simultaneously; the script never plants signs or
modifies projections. Positive-semidefinite certificates allow a zero
pivot only when the remaining row is zero. Scalar sharp certificates
are singular in all three checked elliptic examples, so a strict
positive-definiteness test would incorrectly reject valid equality.
The final audit pins the inputs, preserves all five pass-fifteen
artifacts, and checks that the source manifest is unchanged.

## Seventeenth pass: prime Fourier minors and exact determinant gaps

The primary Tao arXiv v6 paper, math/0308286v6, was archived as
sources/tao-uncertainty-math0308286v6.pdf (130591 bytes, SHA256
244f4e79e667ee831e65d0cb0d8d42de4e354e83d576c65acaf391b0da519880)
and the matching versioned HTML (161191 bytes, SHA256
268bcdcd3a1bcb487c1ea03d6927d4f63e72ba4e0312e2bd3ecc58f1e33f20f8).
Root rendered and visually reviewed complete PDF pages 2,3,4,5.
Theorem1.1, Lemmas1.2–1.3 and Corollary1.4 provide the prime-order
support inequality and nonvanishing Fourier minors. Their hypotheses
and proof were checked against the HTML. The publisher PDF endpoint
returned403; the primary arXiv files were retrieved successfully.
The Cauchy–Davenport discussion is not an input to this pass.

These are imported classical results, not new discoveries. The new
note applies the minor theorem to the two rank-(p-1)/2 frequency
projections. Its determinant divisibility proof explicitly uses the
integer norm in Z[(1+sqrt(p))/2]; omitting this step would lose a
factor4^m. The exact generalized compression, determinant product,
second-moment expansion and certificate decay were derived locally.
No cohomological argument is newly imported or visually reviewed.

The exponential upper bound applies to the derived gap CERTIFICATE,
not to the true Paley gap. The separate interval-frequency construction
does not preserve Paley signs or anchor neighborhoods. Its small p5
case coincides with a Paley projection, while exact cyclotomic tests
at p13,29 verify the loss of off-diagonal signs. The asymptotic loss
is proved by the geometric-series formula. No complete Paley
counterexample or official prize implication is claimed.

The exact verifier checks both real embeddings, independent integer
determinants, rational second-moment bounds, shifted positive minors,
and two singular cases beyond the rank cutoff. The initial run reached
all Paley checks but failed while formatting an unnecessarily large
binomial fraction; the final bounded certificate set completed with
exit0. No mathematical assertion was weakened to pass a test.
Two source-manifest entries were appended, preserving all twenty
previous entries, and all five pass-sixteen artifacts are unchanged.
Independent mathematical review and formal verification remain open.

## Eighteenth pass: primary SL₂ and Weil-representation inputs

[Thomas, arXiv math/0610644v3](https://arxiv.org/html/math/0610644v3)
was archived as sources/thomas-weil-math0610644v3.html, 661651 bytes,
SHA256 5abc992c438e373bfcf243ec61247d9ac3152dd5973f0a2d64a3342a3984d7b5.
Section 2 gives the finite SL₂ character formulas, including the
unipotent exception; Section 3.2 fixes the Weil-index convention.
Root derived the normalized kernels and character values directly.
The exact verifier checks both generator multiplication identities for
all group elements at p=5,13 and every character value there.

[Lyamkin, Sbornik: Mathematics 213:10](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9707&wshow=paper),
DOI 10.4213/sm9707e, supplies the imported expansion input. Root
reviewed Section 1.1's unnormalized Fourier operator convention and
Section 1.5's subgroup classification (Theorem 10/Lemma 5) and
Theorem 11. The latter assumes uniform escape from all proper
subgroup cosets; the Cartesian transporter argument verifies that
hypothesis for all sufficiently large primes. The resulting exponent
is unspecified. The classification and expansion inputs are classical.

Primary live HTML retrieval succeeded. Direct Lyamkin HTML and PDF
archival requests returned 403, so its manifest entry records only
metadata and review scope, with no local file or invented content hash.
The 2022 journal citation and April 7, 2023 webpage publication date
are kept distinct. No PDF was newly rendered or visually reviewed
for the proof inputs in this pass.

The proof retains the exact overlap correction in the Weil trace.
The trace-one conjugacy average is a probability average in the same
representation but is not claimed to preserve the Cartesian input
form. It therefore supplies a limitation of the generic norm argument,
not a Paley counterexample. The verifier uses exact arithmetic;
an integer square root replaced the initial floating primality cutoff
before the final successful run. All five pass-seventeen artifacts
and all twenty-two preceding manifest entries are preserved.
Independent mathematical review and formal verification remain open.

## Nineteenth pass: elementary B_h extraction and quantifier audit

[Shkredov, arXiv:2103.14670v1](https://arxiv.org/html/2103.14670v1)
was reviewed in primary HTML for the Introduction's Sidon definitions
and extraction context. It is archived as
sources/shkredov-sidon-2103.14670v1.html, 394457 bytes, SHA256
bb9624db0131540fe4feb5f414814f0631f5c05a59a8dfd32d9b5de030175f6e.
No sharp extraction result from that paper is imported. The extension
rule, partition, bounded-kernel transfer and sufficient moment
criterion are proved directly in the new note. The background search
also found other higher-energy papers; their theorems were not used.
No new PDF was rendered or visually reviewed in this pass.

The proof includes repeated summands in B_h and uses p>h to divide
by multiplicities. It retains both remainders in the two-sided transfer
and the one-sided remainder in the moment argument. Every chosen
moment and relation order is fixed before taking p sufficiently large.
The restriction h_j/r_j→0 is sufficient for the displayed argument;
no universal necessity or optimality claim is made. The all-exponents
Paley equivalence is to cancellation on B_h pairs, not to the moment
criterion SS-B, for which no converse is established.

The final verifier independently enumerates multiset sums, tests
exact minimal energies, raises the moment transfer to an integer
power, and multiplies actual SL₂ quotient matrices. The high orders
in the rational exponent ledger do not represent tested character
moments at those orders. All exact checks pass. The five pass-18
artifacts and all twenty-four preceding source-manifest entries are
preserved. Independent mathematical review and formal verification
remain outstanding.

## Twentieth pass: inversion, repeated relations and sign completion

No new external source or theorem is a proof input. The pole count
uses nonzero rational functions and the elementary polynomial root
bound. The completion uses the unchanged pass-19 forbidden-point
argument, and the final sampling step uses the unchanged pass-5
single-size transfer. The new inversion formula is derived directly
for both congruence classes of odd primes. No PDF was newly rendered
or visually reviewed, and the source manifest remains unchanged.

Root checked repeated multiset entries, cancellation before forming
the numerator, its nonzero value at every surviving pole, and the
degree reduction from zero total coefficient. The hypothesis p>h
ensures that no surviving multiplicity vanishes in the field. The
moment identity keeps both the removed row and the added zero row.
The signs are coefficients, while only the auxiliary subset choices
are probability averages. No positivity of signed weights is assumed.

The exact-cardinality completion estimates only k-set moments; the
2k-set is used to supply fillers. Both finite inequalities needed
for good poles and completion are explicit. Root verified the
critical leading coefficient (h-1)/(h!)² and the sufficient threshold
K_h. The orders and constants stay fixed before p grows. The new
moment equivalence is up to an explicit factor; it does not prove
either uniform estimate or a new Paley exponent.

The verifier checks numerator polynomials independently of direct
inverse-set sum collisions, pointwise signed kernels, corrected
moments, completions and all sign patterns in the designated cases.
Rational shifted-polynomial certificates check all k above K_h for
the bounded h ledger; the prose proof covers arbitrary h. The final
run passes on current proof and script inputs. The five pass-19
artifacts and the complete source manifest are preserved. Independent
mathematical review and formal verification remain outstanding.

## Twenty-first pass: squarefree recurrence and weighted sign measure

[McDonald, Sahay and Wyman, arXiv:2210.03789v2](https://arxiv.org/html/2210.03789v2),
The VC dimension of quadratic residues in finite fields, was reviewed
in primary HTML, Section 2 and especially Lemma 2.1. Only the standard
bound (d-1)sqrt(p) for a quadratic character applied to a product of
d distinct linear factors is used for comparison. No VC-dimension
result is imported. The local HTML has 339531 bytes and SHA256
3ae1bc2ca60ef412b271e3a03056df70fe1465fda26a89727e1929cfe9e1f622.
The manifest appends this entry after all 25 previous entries.
No new PDF was used, rendered or visually reviewed.

The elementary-symmetric recurrence is derived from the generating
polynomial, and its root bound from a symmetric tridiagonal matrix.
Root checked the N<2r and N=0 cases and the one-zero character rows.
The moment comparison requires only a one-sided top aggregate upper
bound; the opposite side is bounded directly. All constants depend
on the fixed order r. No uniform growing-order estimate is inferred.

The weighted model is an explicit finite positive rational measure,
not a list of p unit-weight rows. Each zero-column fibre has total
mass one but may contain multiple rows. Permutation/sign symmetries
supply the exact first moments and Gram matrix. Distinct-product
correlations are derived by binomial coefficient averaging; the
zero-row coefficient is k-d, not k. The interval containing the
bulk second moment and all factorial denominator bounds are retained.

The construction uses k=floor(p^(1/4)) for every sufficiently large
prime; it does not require a prime in every interval between adjacent
fourth powers. Finite examples choose primes in five such intervals
by exact primality checks, but do not infer a prime-gap theorem.
The large fixture computations are compressed weighted models, not
actual character kernels over billion/trillion-sized fields. Sidon
column labels are independently checked but impose no field-difference
relation on the sign rows. The counterexample therefore concerns
only deductions from the stated relaxed constraints.

The final verifier compares the recurrence with independent binomial
expansion, direct field-product sums with row-polynomial sums, and
compressed weighted correlations with explicit small row measures.
It retains the sixth-moment exceptional-row correction. All arithmetic
acceptance is exact; decimal fixture ratios were displayed separately
for inspection only. The initial patch-writing call failed due to a
missing diff marker and made no file change; the corrected file write
and the first full mathematical verifier run succeeded. No theorem
was weakened to pass a check. Independent mathematical review and
formal verification remain outstanding pending any separately recorded
review result.

The subsequent separate-agent review is recorded in
[its pinned report](parallel21-independent-review-2026-09-05.md).
It found no mathematical correction, checked the algebra and source
code without rerunning the verifier, and explicitly limited its scope.
The report pins the unchanged proof, verifier, transfer and comparison
source. It supersedes the proof snapshot's original pending-review
annotations. This is not human review, external expert certification
or formal verification. Those remain outstanding. The completed
review file is unchanged; that worker then began the spectral lane.

## Twenty-second pass: parallel results, cross-reviews and normalization repair

The root all-orders weighted construction imports no new external
mathematical theorem. Positivity is checked by exact factorizations;
Walsh orthogonality proves every distinct-product correlation, and
pairing even-multiplicity words proves the arbitrary-real-coefficient
moment bound. Every r is fixed before primep grows. Integer model
parameters in the independent verifier are labeled nonprime; the
root verifier separately checks 17 genuine prime-size slices. Both
measures are rationally weighted and do not obey actual character
translate identities. B_h labels do not change this limitation.

The quartic-star lane retains the full p-column kernel identity and
all zero/overlap corrections. Root and a separate reviewer checked
the exact equivalence constants. Its squarefree quartic summands do
not remove the degree-six problem in their mutual Gram entries.
The ordinary Weil comparison uses the unchanged McDonald-Sahay-Wyman
primary HTML source. No VC-dimension theorem is imported.

The subgroup identity follows by unique multiset deletion and exact
extension weights. Root and a separate reviewer checked the U/J
corrections, the small-n absent multiplicity classes and Konyagin's
conditional payoff. The symmetric Sidon obstruction lacks multiplicative
closure. Its unbounded dyadic prime-field family uses the already
imported pass-4 progression nonemptiness. Root freshly reviewed
[Thorner-Zaman, Section3.1, Corollary3.1 and equation3.2](https://arxiv.org/html/2108.10878#S3.SS1)
for the powerful-modulus hypothesis and removal of the exceptional
real zero at fixed radical2. At q=N,x=N^4,h=3N^4/4, the required
h/phi(q)>=(x)^(7/12+1/12) holds and the relative error tends tozero.
The existing v2 HTML archive remains unchanged; no numerical eventual
threshold is claimed. The subgroup agent's trilinear background search
imports no new theorem into its relation decomposition.

The spectral lane derives the Jacobi-square and anchor identities
explicitly. Its original assertion that every needed identity was
derived locally omitted the sign of the quadratic Gauss sum needed
to name the residue frequencies. The separate reviewer identified
this source-scope gap and supplied the classical normalization from
[NIST DLMF20.11.1-2](https://dlmf.nist.gov/20.11#E2). Root read the
official HTML and its coprimality/parity hypotheses. Specializing to
(m,n)=(2,p), p1mod4, gives the positive sign. This changes no spectral
formula; it makes the omitted imported input explicit. The projection
identities alone give the positive eigenspace, not its frequency label.
The finite surd checks do not independently verify that Gauss sign.

Root amended only the spectral note's normalization/source-scope prose,
preserving its original reviewed note/result bytes in
results/parallel22_spectral_reviewed_snapshot_2026_09_05.json. The
unchanged verifier was rerun on the amended note for final provenance;
the primary-source argument, rather than this finite rerun, supplies
the sign normalization. The independent review's original hashes are
resolved against the preserved snapshot where the current note/result
changed. Other review inputs remain unchanged.

NIST's official HTML displayed version1.2.7, release2026-06-15.
Live primary review succeeded; direct archival retrieval returned403.
The source manifest therefore appends a metadata-only entry after all
26 preceding entries, with no fabricated byte count or content hash.
The old manifest bytes are preserved at /tmp/paley22-manifest-before.json.
No PDF was newly opened, rendered or visually reviewed in this pass.
The spectral agent's Bhowmik-Barman abstract lookup is background only;
its full text was not reviewed or imported.

All four new mathematical outputs received separate-agent review;
one reviewer also wrote an independent exact weighted-model program.
This is not human refereeing or formal proof certification. One
spectral cross-review assignment failed immediately at model capacity;
another available worker completed that review. No model override,
usage reset, purchase, external message, proof submission, automation
or memory write was performed. The full goal remains unproved.

## Twenty-third pass: actual partial upper bounds and spectral coupling

The [pass assessment](parallel23-pass-summary-2026-09-05.md) distinguishes
four limited deductions from the full open targets. The classical almost-all
Sidon result uses only locally derived sampling, slice-variance and complete
quadratic-character identities. Its coefficient 216/5 is an upper-bound
constant, not an exact variance or a worst-case exponent. One result-label
issue found in review was corrected: the conditional Sidon certificate is
null outside n>=6 and n^4<=p. The verifier was rerun with identical check
counts; no proof formula changed.

The subgroup product-ratio bijection and balanced-part estimate use the
existing [Shkredov primary HTML](https://arxiv.org/html/1504.04522v1), Theorem 6.
Root freshly checked the nonzero shifts, size product below p, and
zero-inclusive energy convention. This proves a bound for a subset of
opposite-free relations, not their total count. The existing local HTML
remains the source input.

The spectral coupling limit imports Grove's published square-parameter
second moment from [the primary publisher article](https://link.springer.com/article/10.1007/s11139-026-01372-y),
Theorem 1.3 at moment 2. The introductory point-count normalization,
the boundary values, and Theorem 2.3(2) at k=2 were checked. The HTML archive
has 450232 bytes and SHA-256
1174074d2d3104414f03446b060a245d935c0fe2ca80beb2c14a01bbe10ec440.
The transfer to C and block-norm computation are derived in the local note.
The finite verifier does not independently compute modular coefficients.

The prize analysis uses the archived September-4 contract. Root additionally
retrieved four ArkLib definition dependencies at its already pinned commit
e65197892890b8fd9b0dc05b8980273cf1d595cc, preserving source-version scope.
[The new source ledger](../results/parallel23_prize_additional_sources_2026_09_05.json)
records the exact bytes and URLs. Lambda is a maximum of distinct lists;
IsMCA is an existential projected-code witness with a nonexplainable input.
These definitions support the note's distinct remainder and ratio counts.
Neither production maximum nor the required prize inequality is evaluated.

The manifest preserves all 27 previous entries and adds Grove plus the
four exact-pinned Lean files, for 32 entries. No PDF was newly used.
Separate-agent reviews and final hash reconciliation are linked from the
assessment; they are not human refereeing or formal certification. The full
Paley, subgroup, spectral and prize targets remain unproved.

## Twenty-fourth pass: reused primary inputs and exact finite arithmetic

The [pass-24 assessment](parallel24-pass-summary-2026-09-05.md) records
four bounded lanes and separate-author reviews. Root read the versioned
[McDonald–Sahay–Wyman HTML](https://arxiv.org/html/2210.03789v2), Lemma
2.1, live for the quadratic-character Weil bound. The new inversion
note reduces repeated multiplicities explicitly to a monic squarefree
polynomial and missing-root errors. The squared complete-sum inequality
was already derived in the sigma sum-product note; only its combination
with the good-pole count is new within this project. Neither a uniform
signed bound nor literature novelty is claimed.

The spectral author and separate reviewer inspected the pinned primary
Katz sources: *G₂ and hypergeometric sheaves*, section 2, and *Gauss
sums, Kloosterman sums, and monodromy groups*, the tame Euler
characteristic, trace formula, and compact-cohomology weight bound.
The disjoint lists, fixed-base-field Gauss twist, rank mismatch, and
different local Jordan forms are the inputs for the cross-pullback
estimate. The result controls the actual second Krylov space and its
nonzero leakage. It is not a full spectral-edge theorem. Root reviewed
the derived argument; the primary PDF inspection in this pass belongs
to those worker lanes. No PDF inspection by root is claimed here.

The subgroup norm criterion, dyadic irreducibility and orthogonality,
intrinsic product-ratio ledger, and finite order-16 classification are
proved or exhaustively certified locally. They import no new source
estimate for the unbalanced aggregate. The local exceptional-set lane
retains the high moments needed by its conditional variance bound and
does not turn the previous global average into an adversarial estimate.

All 32 source-manifest entries remain unchanged. Existing archive hashes
are rechecked by the [final audit](../results/parallel24_pass_audit_2026_09_05.json).
No official contract, field profile, or prize status was refreshed.
The original full goal remains unproved; finite exact checks and agent
reviews do not constitute human refereeing or formal certification.

## Twenty-fifth pass: scalar projection and subgroup energy inputs

The [pass-25 assessment](parallel25-pass-summary-2026-09-05.md) records
four deductions and their limits. The new random-projection proofs for
MCA and list maxima are self-contained. The primary
[Gopalan–Guruswami–Raghavendra abstract](https://arxiv.org/abs/0811.4395)
was read only to establish provenance for classical interleaving and
list-decoding work; no theorem, finite constant, or formula is imported
from it. The full arXiv HTML request returned 404. No PDF was used.
The [source ledger](../results/parallel25_prize_projection_source_2026_09_05.json)
pins the actual abstract HTML. The official combination-round condition
continues to use the same archived contract and ArkLib commits. The
scalar reduction does not refresh or settle the official prize.

Root read [MRSS Theorem 3 and Section 2](https://arxiv.org/html/1712.00410v1)
in primary HTML. Its subgroup additive-energy bound is
E(H)≪|H|^(49/20)log(|H|)^(1/5) for |H|≤sqrt(p). Its E(H)
equals this project's E2(H). Its separately defined cubic moment of
difference multiplicities does not equal this project's equal-triple
energy E3(H); no such cubic estimate is imported. The
[new HTML ledger](../results/parallel25_subgroup_source_2026_09_05.json)
records this exact scope. Earlier source archives remain intact.

Root also rechecked [Shkredov Theorem 6, Section 4](https://arxiv.org/html/1504.04522v1)
in primary HTML: taking both subgroups H, nonzero shifts −1, and n²<p
gives the shifted multiplicative-energy bound used in the balanced
remainder. These conditions hold in the stated quartic window.
The two energy inputs bound a proper part of the six-term relation
count; six distinct fully unbalanced relations remain uncontrolled.

Signed inversion and the full spectral decomposition use explicit
identities and previously stated inputs. The spectral note compares
its sectors to the pass-13 depth allowance; it does not supply a new
proof of that allowance or bypass its earlier review limitations.
The [root review](parallel25-root-review-2026-09-05.md) checks the three
worker deductions and current result hashes. The prize note is
root-authored and its separate-author review did not finish before
worker usage-limit errors. No human or formal certification is claimed.

The source manifest retains all 32 previous entries and adds the two
primary HTML archives above, for 34 entries. The
[final audit](../results/parallel25_pass_audit_2026_09_05.json) verifies
their bytes and hashes alongside prior artifact preservation. There is
no new worst-case cancellation exponent or full goal proof.


## Twenty-sixth pass: existing MCA theorem and independent Lean proof

Root read the complete `mcaError_interleaved_le` and
`mcaError_interleaved_eq` proofs in the previously archived Errors file,
and checked the live primary file at the same pinned ArkLib commit.
This exact result was already present before pass 25. Its proof avoids
at most |F| proper subspaces simultaneously, preserving all original
bad scalars. The earlier MCA progress framing is corrected to independent
rediscovery. No new scalar cancellation theorem is claimed.

Two additional files, InterleavedCode.lean and Basic/LinearCode.lean,
were archived from that same commit. They identify projected code
membership with codeword agreement and rowwise membership. The
[ledger](../results/parallel26_primary_sources_2026_09_05.json) pins their
bytes. The main manifest preserves all 34 previous entries and adds
these two, for 36 entries. The current
[formalization note](parallel26-mca-formalization-2026-09-05.md) distinguishes
the seven compiled Mathlib statements from the source-level ArkLib
comparison; the latter dependency closure was not rebuilt.

No matching finite list collision-count transfer was found in the
inspected ListDecodability and new code files. This is a limited search,
not a claim of novelty or absence elsewhere. The existing project
coefficient-pigeonhole threshold was also recognized and not counted
again as progress. All full mathematical targets remain open.


## Twenty-seventh pass: source inputs bound computation cost

The [new exact-count note](parallel27-six-distinct-2026-09-05.md) uses
MRSS Theorem3 and Shkredov's shifted-energy Theorem6 only to bound its
candidate-enumeration cost. Root rechecked both exact-version statements
in primary HTML, with n<=sqrt(p) and n²<p respectively. Both hold in the
quartic window. The resulting O(n^(49/20)polylog n) cost is not an upper
bound on the six-distinct count D6. The separate cubic difference-energy
notation is not substituted for the equal-triple energy.

The [source scope record](../results/parallel27_source_scope_2026_09_05.json)
pins the existing archives. All36 source-manifest entries are unchanged.
A limited fresh search did not yield an additional bound used in this
argument; it does not establish literature exhaustiveness or optimality.
The unique-split and exact unmarking arguments are self-contained project
derivations, without a literature novelty claim. Only their displayed
algebraic core is checked in Lean; the full counting formula is not.

## Thirty-eighth pass: published energy theorem and missing hypothesis

The [publisher-version audit](parallel38-published-energy-2026-09-06.md)
resolves the discrepancy left in pass37. Shkredov's published paper,
Moscow Journal of Combinatorics and Number Theory3(2013),issues3-4,
pp189-239, has Theorem8 on printed216 (PDF29), with the small-subgroup
energy term n^(32/13)log^(41/65)n. The older arXiv1208.2344v3 Theorem34
instead advertises22/9. The published32/13 result is compatible with
MRSS49/20; the old statement is not imported into this project's bounds.

Published Lemma2 on printed197 (PDF10) adds an equality for the size of
the normalized pair image. The [independent audit](parallel38-independent-incidence-2026-09-06.md)
proves that it requires injectivity of the corresponding coset-triple
map. This is absent from the preprint's Lemma7 and fails for diagonal
multi-coset triples in its proof. Exact examples and a conditional
weighted repair are recorded in the [pass summary](parallel38-pass-summary-2026-09-06.md).
The audit does not claim the numerical22/9 bound itself has been refuted.

Raw publisher bytes, provenance HTML, and the two inspected page images
are in the separate [source manifest](../sources/parallel38-energy-source/manifest.json).
The publisher PDF SHA256 is
79dbfd8841ae220880ad8c692684bcd0889e417b700a96e2ed64f3ac7983c206.
Root and both parallel reviewers inspected the decisive pages. The main
36-entry source manifest remains unchanged. No full published proof,
Lean formalization, uniform period improvement, or full target is claimed.
