# Root review of the three worker deductions

The signed inversion, growing subgroup, and full spectral notes have been
read and checked by root, a different author from each worker. No
mathematical correction was found. This review includes a separately
written verifier that imports none of the worker functions. It is an
agent review with finite checks, not human refereeing or formal proof.

The root-authored prize projection note has only its author's proof
audit and finite verification. Its assigned separate-author review did
not finish: all three workers ended with usage-limit errors. None is
currently running. The spectral author saved its proof and verifier;
root completed the verifier run and generated its result after the error.
No claim that all four lanes received separate-author reviews is made.

## Signed inversion

The two actual sign halves retain the pole row correctly: the formal
value at the pole and the missing transformed value differ by a sign,
which disappears in the sixth moment. Their joint identity has a
positive copy of the original moment. In the completed average,
T1 vanishes and T2=p(nM4−A4)−M6 includes the negative diagonal term.
The odd mixed sums control the difference of the half moments, not
their common magnitude. The joint good-pole upper bound still contains
the original M6 and cannot close the desired estimate.

The general fractional-linear change multiplies each input weight by
chi(cu+d). Keeping the missing row gives exactly
M6(new)+|sum(new weights)|^6=M6(old)+|sum(old weights)|^6.
The two-inversion calculation composes these weights; it supplies no
independent random signs. These points agree with the stated scope.

Independent direct character computations at p=31 and 41 check 130
actual half moments, 65 joint identities, both completed mixed identities,
and six augmented moment identities with mixed starting weights. The
author's larger integer verifier also passes. These computations check
the identities; they do not turn the remaining asymptotic hypothesis
into a proved bound.

## Growing subgroup orders

Marking an equal positional pair in a zero-sum six-word gives exactly
15n r4(2). Intrinsic words contribute 15n(6n−8), so the difference
is a nonnegative weighted count of nonintrinsic repetitions. Each
opposite-free repeated word has positive weight, establishing its stated
upper bound. Cauchy gives r4(2)≤E2. Summing product-ratio fibers costs
no further factor of n because those fibers partition the fixed split.

The root inspected the primary HTML of MRSS Theorem 3 and its energy
definitions. The imported bound is E2(H)≪n^(49/20)log(n)^(1/5),
under n≤sqrt(p). The paper's third difference-multiplicity moment is
not the six-term energy used here and is not imported. Shkredov's
Theorem 6 was also rechecked in primary HTML with the two nonzero
shifts −1 and n²<p. These hypotheses hold in the stated quartic
window. The resulting error in E3=T6+D6+error is bounded independently
of D6. No estimate on the six-distinct fully unbalanced D6 follows.

The norm argument correctly distinguishes rational rank from reduction
modulo p. A nonzero polynomial of degree below deg(Phi_n) has nonzero
norm. Splitting at distinct primitive roots makes the rank defect count
the vanishing roots; divisibility only supplies at least one p factor
per such root. Rational independence of all cyclic shifts does not
supply separate factors. The AM–GM bound uses squared coefficient mass
Q and exponent n/4, with Q=6 for six distinct entries.

Independent weighted enumeration checks the incidence subtraction,
partition, and inequalities at (p,n)=(97,8),(33713,16),(37201,16).
The first is explicitly outside the quartic window and tests only the
elementary identities. A separate integer Bareiss determinant calculation
at n=256 confirms norm 4419283227438373820201815271409083396 and
valuation one at p=1073748737. Finite cases do not bound D6 uniformly.

## Full spectral operator

The reflection and inversion permutations preserve the actual two-anchor
set and sign matrix. Their S3 representation gives all orthogonal
sectors, with the standard representation appearing as two equivalent
blocks. The uniform vector belongs only to the trivial sector, so all
its polynomial iterates stay there. The exact orbit count gives that
sector asymptotic dimension m/6, leaving about 5m/6 dimensions unseen
regardless of iteration depth. The small-prime extremizers in other
sectors are examples, not asymptotic counterexamples.

The full square identity retains the inversion correction −J. The
leakage identity is obtained from the actual full-field involution and
keeps both rational and sqrt(p) coefficients, including the rank-two
border on the trivial sector. On the other sectors the desired edge
requires an upper bound (2/3+o(1))p on the elliptic conjugacy average;
the elementary bound remains p. The rank-two border does not vanish
on the trivial sector. Finite block decomposition does not remove the
growing-depth difficulty because each nonempty asymptotic sector still
has dimension proportional to p.

Independent computations at p=13,17,61,269 build the full field matrix
first, derive the elliptic kernel by direct matrix multiplication, and
calculate off-block leakage directly. They check both its rational and
surd coefficients, all sector projections, dimensions and invariance.
The p=269 case is beyond the author's range. Integer accumulation has
an explicit overflow guard. The author's 24-prime verifier also passes.

## Evidence

The [independent verifier](../experiments/parallel25_root_review_2026_09_05.py)
and its [result](../results/parallel25_root_review_2026_09_05.json) pin the
three reviewed proof notes, their results, and the newly archived MRSS
HTML. The final pass audit checks this ledger against current bytes and
also pins all author scripts. No worker note or prior proof was edited
to make these checks pass. The [pass assessment](parallel25-pass-summary-2026-09-05.md)
records the remaining proof obligations and the interrupted fourth review.
