# Current proof frontier

Stepanov wave (2026-09-26/27, [summary](stepanov-pass-summary-2026-09-27.md)):
a Hankel-minor refinement of the Hanson–Petridis polynomial proves
`N₋ ≥ (e+1)[n − (d−e+r)/(m−2e)]` non-residue sums for every `e ≤ (m−1)/2`,
hence robust Hanson–Petridis with `f(η) = (1−√(2η))^{−2}` and a constant
bias saving `≍ κ²` whenever `|A||B| ≥ (1/2+κ)p`. The threshold `1/2` is
sharp. This closes the constant-saving question at the square-root scale;
the conjecture needs a power saving below it, and the method's degree
budget `(e+1)(d−e)` still caps it at `|A||B| ≍ p`. Next targets: the
optimal `η(κ)`, Lean formalization, and any way to shrink the degree budget.


Latest: the [fifty-third pass](parallel53-pass-summary-2026-09-06.md)
returns to the direct two-anchor operator. The elliptic kernel on all
quadratic residues diagonalizes with eigenvalues Re(J(psi,chi)^2).
Lu–Zheng–Zheng Theorem1.4 equation1.7, with one fixed quadratic
character and one full character family, gives inversion-odd eigenvalues
at least p-32pi^2sqrt(p). They exceed2p/3 for every eligible p>=2^20
and have ratio tending to1. This rules out an ambient norm saving
asymptotically; it does not refute the compressed S_3 average required
by pass25. The exact p89 eigenvalue73 has an integer odd vector with
nonzero support outside C. All4,096 cyclic coefficients and3,198
restricted-identity entries pass. The next direct spectral input must
control the joint C-support restriction and S_3 average. All uniform
exponents and full goals remain unchanged.

The [fifty-second pass](parallel52-pass-summary-2026-09-06.md)
uses four-child incidence coarsening to prove Delta_s=X_s-X_(s/2)>=0,
Delta_s>=36Y_s and Y_s^2<=sDelta_s/24. It also controls the increase
in normalized energy excess and characterizes zero cubic increment.
Holder gives sum_(M<s<=N)sqrt(Y_s)<<N^(1/4)(X_N-X_M)^(1/4), so
X_N-X_M<<N^(5/3) at the pass50 cutoff would suffice. That increment
estimate remains unproved. The existing X_N<<N^2 log N returns only
energy5/2, weaker than49/20. At p2144280833,N256 the cubic excess
is1260 at both128 and256; no fixed-factor growth at every step holds.
All368 distinct steps,71 towers and397 cutoff checks pass. Uniform
exponents and full goals remain unchanged. The next task is a genuine
quantitative saving, beyond the new accounting identity.

The [fifty-first pass](parallel51-pass-summary-2026-09-06.md)
tracks primitive labels into the full endpoint kernel. Its exact
compatibility equation is U+2C+H=A_tot/2+2sum_s T_s/s, where U counts
extra mergers under the power maps, C counts coincidences across levels,
and H counts the distinguished label. The energy lower bound sharpens
to E_N>=3N^2-3N+6Nsum_s Y_s. At the actual quartic endpoint
p2144280833,N256, Y128=3 but Y256=0, while D128=D256=48.
Thus a pointwise top-layer domination for every tower is false; this
finite example does not rule out an eventual or additive-error variant.
All71 towers,397 levels and2,731,208 literal pairs pass. No uniform
exponent improves. The next task remains the upper-tower estimate,
using the full compatible label data when needed.

The [fiftieth pass](parallel50-pass-summary-2026-09-06.md)
proves 3N+6N sum_s Y_s<=E_N and
sqrt(E_N)<=sqrt(14N^2-25N)+2sqrt(2N)sum_s sqrt(Y_s).
The resulting upper bound28N^2+16N(log_2 N-1)sum_s Y_s compares
energy with unweighted actual triple mass up to a logarithm. Recovery
from any intermediate subgroup improves the earlier sufficient weighted
criterion. Using the existing49/20 energy theorem handles all levels
through the largest dyadic M<=N^(80/87)/(1+log N)^(4/29). The remaining
sufficient input is sum_(M<s<=N)sqrt(Y_s)<<N^(2/3), still unproved.
Exact checks cover70 towers,390 levels and460 cutoffs. Abstract
equal-mass profiles only establish the limitation of the scalar Cauchy
step; no field realization is claimed. Uniform exponents and all full
goals remain unchanged.

The [forty-ninth pass](parallel49-pass-summary-2026-09-06.md)
replaces the two-resultant gcd by Ccrit_n, the gcd of all coefficients
of Res(R_n,R_n'+T R_n''). Its valuation is the sum of the minima of
the two derivative valuations at the lifted roots. It dominates all
triple masses over p-adic precisions, with an explicit nonnegative
reciprocal-sum correction. Generic polynomial examples prevent assuming
either equality with the coarse gcd or equality with that collision sum.
Norm geometry and field-degree descent show that every odd prime
dividing Ccrit_n is 1 modulo n. The target primes satisfy this condition,
so no target estimate improves. Through order128 the refined integer
equals G_n, and all70 checked splitting cases have zero higher or
cancellation contribution. These finite facts are not uniform identities.
The next task remains a bound for the splitting-prime aggregate.
Uniform exponents and full goals remain unchanged.

The [forty-eighth pass](parallel48-pass-summary-2026-09-06.md)
uses the gcd of two derivative resultants to classify primitive triple
fibers through order128. Precisely three primes in its quartic interval
have maximum fiber3; all three have B5120<8192. The proposed uniform
maximum-fiber-two condition is false, while the aggregate quadratic
bound survives throughout this fixed-order interval. For general dyadic
n, Y_n=sum_(c>=3)c(c-2)<=v_p(G_n) and B_n<=n^2/2+nY_n.
The weighted tower sum V_N gives E_N<=44N^2+22N^2 V_N. A uniform
V_N<<N^(1/3) would reach the intermediate energy7/3 and triangle17/3
targets, but remains unproved. Complete factorizations, prime
certificates and literal counts passed. The next task is an aggregate
estimate, with no reliance on the refuted pointwise condition. Uniform
exponents and full goals remain unchanged.

The [forty-seventh pass](parallel47-pass-summary-2026-09-06.md)
derives h_i>=a_i(a_i-1)+a_(-i)(a_(-i)-1) for i!=0 from known matrix
entries. Both earlier uncentered interval families fail this condition.
A centered repair satisfies it, the formal product congruence, and every
total-X row condition at large n, while retaining the difficult norm
power61/15. Actual field power sums impose stronger constraints: the
first n-1 characterize the complete profile by Newton reconstruction.
All2,048 tested positions of the repair fail the first moment at one
certified quartic prime. This finite exclusion and recognition theorem
give no uniform coefficient or excess bound. That remains the next task,
using field arithmetic or the full matrix law. Exact checks passed;
uniform W86/15, energy49/20, period71/72, absolute exception5/3 and
the full targets remain unchanged.

The [forty-sixth pass](parallel46-pass-summary-2026-09-06.md)
proves H_D<=nd max(1,sqrt(R_D))+6^(1/3)d^(2/3)X_D^(2/3), where h_i
is the row's repeated-incidence mass, H_D its sum, and R_D the maximum
outside mass. For actual C, h is determined by a's autocorrelation.
The Sidon-filled interval forces X>=c n^(11/5), excluding its realization
at sufficiently large size under the source bound. A long-interval filler
retains the correlation power while satisfying every total-X consequence.
It has no assigned local X_D or actual matrix. The next task is to use
those local allocations or the off-diagonal multiplication identities.
Finite checks passed; no uniform exponent or full target improves.

The [forty-fifth pass](parallel45-pass-summary-2026-09-06.md)
derives L(u)L(v)=n<u,v>I-nuv^T+L(L(u)v) for arbitrary real weights.
Exact rank corrections, cubic and quartic traces, and the projection onto
the shifted-matrix span are retained. For u=a^2, the remaining centered
correlation contains W. Substituting only the current scalar moments into
the resulting Cauchy bound gives power7, weaker than86/15. Warren's
checked shifted-product consequence is also insufficient at the surviving
interval parameters; a specified stronger cubic product input is unproved.
The checked2026 complex-group theorem is not a finite-field substitute.
All exact finite checks passed. The next task is targeted control of the
actual incidence correlation using the full multiplication law, rather
than a scalar norm substitution. Uniform exponents and full goals remain
unchanged.

The [forty-fourth pass](parallel44-pass-summary-2026-09-06.md)
proves sum_D a<<n^(2/3)|D|^(1/3) for every quotient subgroup coset D.
Whole level cosets therefore satisfy T^3|D|^2<<n^2. Nested level-coset
pieces contribute at most the target triangle power17/3, and a stated
polylogarithmic partition of levels above n^(23/60) would suffice up to
logs. That partition is not established. Abstract weighted intervals
still attain correlation power61/15 while satisfying all these coset
mass caps and the other compared scalar bounds. They are not actual
subgroup data. Actual C satisfies C^2=nI-ne0e0^T+L(a) and an exact
shifted Frobenius identity, now derived and checked. The next task is
to use this additional arithmetic against general high-level
concentration. Uniform W86/15, energy49/20, period71/72, the absolute
prime-exception exponent5/3, and all full targets remain unchanged.

The [forty-third pass](parallel43-pass-summary-2026-09-06.md)
improves the absolute exception count for W<<n^(17/3) in a quartic
prime interval to O_c(n^(5/3)/log n). The exact inequality
W<=R_max E_*^3/n^2 and the published single-coset bound R_max<<n^(2/3)
make energy E<<n^(7/3) sufficient; the existing fourth-energy prime
average controls how many primes can fail that criterion. This leaves
individual exceptional primes possible and does not assert a proportion.
At p278177,n32, minimum additive energy coexists with X_dist=X72 and
zero diagonal rich mass, ruling out all constant-factor diagonal
domination even inside the quartic window. Compatible root lifts give
v_p(mathcal_P_n)=sum_r X_r, with exact finite norm comparisons.
The uniform W, energy49/20, period71/72 and full targets do not improve.
The next task is individual exceptional-prime control or stronger actual
incidence structure beyond the ruled-out diagonal comparison.

The [forty-second pass](parallel42-pass-summary-2026-09-06.md)
proves sum_(p=1 mod n) X_p log p <=(M/2)log(S/M), with
M=(n-1)^4-2(n-1)^2+(n-1) and S=8n^2(n-1)^2-2n^4.
For each fixed dyadic n>=4, this gives an absolute exception count
O_c(n^(63/31)/log n) for X_p>n^(61/31) in a quartic prime interval.
The pass41 criterion gives W<<n^(17/3) outside those exceptions.
No bound for every prime or eligible-prime proportion follows.
The actual subgroup p353,n16 has X_dist=X72 and zero diagonal rich mass,
refuting every constant-factor diagonal domination in the general sparse
range; it lies outside the quartic window. Matched reflection projections
give a positive smoothed sixth moment but an exact Holder recovery
certificate at least as large as the direct one. Norm and reflection
checks passed. The full goal, uniform energy49/20, and period71/72
remain unchanged. The next task is individual exceptional-prime control
or additional restrictions on the actual joint incidence distribution.

The [forty-first pass](parallel41-pass-summary-2026-09-06.md)
proves W<<n^(86/15)(1+log n)^(8/15) for the full nonzero weighted
triangle sum, retaining all multiplicities and removing the initial
dyadic logarithm loss. The power gap to17/3 is1/15. The exact dependence
of the difference tail on nontrivial shifted-energy excess X shows that
X<<n^(61/31), up to logarithms, would suffice to close this triangle gap.
That input remains unproved. Real-subfield recovery changes the known
scalar cofactor1217^31 to1217; centering and full-space unit balancing
still do not supply the uniform shell estimate. A positive smoothed sixth
moment has a quantified unsmoothing gap. Independent checks passed;
energy49/20, period71/72 and all full-conjecture targets are unchanged.

The [fortieth pass](parallel40-pass-summary-2026-09-06.md)
rules out reversing the short-shell construction into a norm-p generator:
the actual order64 example has V<1 and a nonprincipal kernel. Reciprocal
recovery gives cofactor tau^(N-1), with exact scalar-family loss1217^31
in that example; other certificates can have smaller cofactors.
The actual quartic order128 difference levels also have positive repeated
fibers with three distinct edge cosets. Correct incidence52736 becomes
51968 after discarding multiplicity. The collision correction is now an
explicit sum of multiplicative level intersections, still unbounded at
the needed scale. Classical Sidon translation/container bounds and the
distinction between relative and input-norm approximation are recorded.
Root independently reviewed and checked the results. No uniform exponent
or full target improves.

The [thirty-ninth pass](parallel39-pass-summary-2026-09-06.md)
proves a cyclotomic cofactor restriction and an exact finite obstruction
to the literal constant C=2: p215535361,n128,a38468180 has V=1 and
eta>2 sqrt(n log(p/n)). The bound with an unspecified constant remains
open. Norm(f)=kp with p not dividing k forces V<=k^2 for some nonzero
coset; any existing M=o(n) implies k>=(1-o(1))sqrt(n)/(2 pi), preventing
small-cofactor constructions from scaling in the quartic window.
The coincident-edge triangle contribution is <=3 F3* F4*/n; its
three-distinct-coset complement still needs a bound. Independent agents
reviewed these constructions and the prior moment/shell comparisons.
No full target or uniform exponent improves.

The [thirty-eighth pass](parallel38-pass-summary-2026-09-06.md)
resolves the22/9 versus49/20 energy discrepancy. The published source
has exponent32/13 and adds an injectivity condition to its incidence
lemma. That condition fails for diagonal triples of sets containing
multiple subgroup cosets. Exact enumeration verifies failure on an old
proof level: p1153,n8,|S|24, ratio image1600 instead of1728. A separate
p97 example counts48 incidences versus16 after multiplicity is dropped.
These are hypothesis and counting failures, not asymptotic energy
counterexamples. Two independent reviews support the source comparison
and multiplicity identity. The current49/20 energy input,129/40 repeated-six
bound, and71/72 uniform period exponent remain unchanged; the full goals
and positive moment upper bound are unproved.

The [thirty-sixth pass](parallel36-pass-summary-2026-09-05.md)
gives an exact Bernoulli-series and Dirichlet-inverse comparison between
centered periods and the shortest squared norms of dual-lattice cosets.
For maximum deviation D and centered period maximum M_f,
12p^2/(p^2-1)D<=M_f<=36p^2/(p^2+1)D; the left constant becomes24
when multiplication by2 preserves the set. The same holds in every
normalized l^r norm, r>=1. The desired uniform thin-annulus estimate
is open. A rank-one cyclotomic norm bound gives only V>=p^(-2/N),
far from the required mean near N/12. Truncated inversion yields an
exact finite certificate without high moments and recovers the known
order64 maximizing coset. It does not improve the prior tight enclosure
of M. No uniform exponent or full target is proved. Pass39 independently
reviews sections1-4 and clarifies A=H,a!=0; Lean review remains open.

The [thirty-fifth pass](parallel35-pass-summary-2026-09-05.md)
proves, for every symmetric nonzero set and every degree s>=1,
B_s>=-(256sn)^s and T_s<=16^s B_s+2(4096sn)^s, together with
|B_s|<=4^s[T_s+(64sn)^s]. The negative side is controlled, and the
remaining sufficient input is a one-sided upper bound at one selected
degree B_s<=(Ksn)^s. It no longer requires a hierarchy or the pass33
epsilon condition. The proof applies Ravichandran's derivative-root
theorem before averaging and covers all sizes. An actual quartic
p1153,n8 example refutes real-rootedness of the averaged generating
polynomial. The positive upper estimate, period improvement, and full
Paley/prize targets remain unproved. Pass39 supplies separate-agent
review of the comparison; Lean review remains outstanding.

The [thirty-fourth pass](parallel34-pass-summary-2026-09-05.md)
separates opposite pairs with an exact positive transform and an inverse
whose coefficients count independent sets in cycles. For n=2N and
2s<=N, the centered distinct counts Q_s and centered opposite-free
counts B_s satisfy Gaussian hierarchy bounds together, up to absolute
changes in constants. Combined with pass33, this holds at logarithmic
depth in the quartic window. The missing input is a one-sided uniform
bound B_t<=K^t(2t-1)!!n^t; it is not proved. For every eligible symmetric
set with n>=2^35 and p<=n^4, O_10>=n^6/2, so the raw high-degree count
cannot replace B_t. Opposite-free six-words are not merely D6. No full
energy or period exponent improves; all full targets remain open.

The [thirty-third pass](parallel33-pass-summary-2026-09-05.md)
proves a centered distinct-coordinate reduction with explicit relative
error product_(j=1)^(2s-1)(1+j/sqrt(alpha))-1, where
alpha=n(p-n)/(p-1). In the quartic window this is O(s^2/sqrt(n))
and tends to zero at logarithmic depth. Distinct counts subtract
(n)_(2s)/p, while full moments subtract n^(2s)/p; every collision
partition is centered before Holder. For repeated six-words, subgroup
averaging gives R6<=15 sqrt(E2 E3), improving the recorded component
exponent to 129/40 with log power3/5. The remaining distinct-coordinate
upper bound is unproved and includes opposite-pair configurations;
it is not simply D6. No full-energy or period exponent improves.

The [thirty-second pass](parallel32-pass-summary-2026-09-05.md)
proves that every six-term relation generated by a fixed triangle orbit
has a derivation using at most 8m copies, with m the least subgroup
order containing the normalized triangle. The inverse trinomial has
coefficient norm at most 4m/3, by comparison with the explicit inverse
of 1+X+X^2 and Parseval. A primitive generated target lies in one
coset of that subgroup and has minimum cycle rank at most 4m-2.
Integrality of the quotient remains essential. This is a uniform
length estimate, not a bound on the number of outputs or a stronger
cancellation exponent; separate-author and Lean review remain open.

The [thirty-first pass](parallel31-pass-summary-2026-09-05.md)
classifies all 119 remainder orbits at p=215535361, n=128 by their
minimum number of scaled zero-triangle terms. The unique integer
quotient by 1+X+X^19 proves the minima: 54 triangular orbits need 2,
20 primitive orbits need 4, and 45 primitive orbits need 38-44 terms.
Every minimum cancellation graph of a primitive target is connected
with cycle rank L/2-2, so the latter 45 require 17-20 cycles. These are
exact finite certificates, not a uniform growth theorem or a new
bound on the signed spectral quantity. The full goal remains unproved.

The [thirtieth pass](parallel30-pass-summary-2026-09-05.md)
refutes the literal D6<=n^3 extrapolation with the actual quartic-window
field p=215535361, n=128: D6=10967040. A zero triangle gives a short
analytic lower certificate; sparse and direct counts agree on the full
answer. The triangular part has 54 scaling orbits and the primitive
part 65; each exceeds n^3 in ordered-word count. A general two-sided
bound for triangle pairs shows that D6=O(n^3) would imply
|H intersect (1-H)|=O(sqrt(n)). All relations at this example lie in
the same principal ideal generated by 1+X+X^19, whose norm is p.
Uniform upper bounds and all full targets remain open. The prior
finite census is unchanged; it never covered order128 exhaustively.

The [twenty-ninth pass](parallel29-pass-summary-2026-09-05.md)
proves T_(2s)<=D n^(2s-1/2)T_s for exact-centered subgroup energies,
with an absolute constant uniform in s. It keeps signed weights and
the origin under control and removes the principal term. The recorded
seeds, interpolation, this recurrence, products of coset moments, and
feedback of the mixed-moment amplitude bound have an explicit closed
power envelope with saving 1/72. No stronger cancellation exponent or
full target follows. The proof is ordinary mathematics from a pinned
incidence theorem; independent review and Lean verification remain open.
Prime-subfield examples prevent an automatic extension-field transfer.

The [twenty-eighth pass](parallel28-pass-summary-2026-09-05.md)
derives M << n^(71/72)(log n)^(1/18) for multiplicative subgroups in
the quartic window by combining existing energy inputs with the mixed
(3,6) Konyagin inequality. The nonnegative weighted lifting controls
the origin separately. This is an ordinary proof with separate-author
review outstanding, and improves the prior project exponent by 1/320.
No novelty or current-best claim is made. Iteration of these energy
inputs supplies a maximum saving of 1/72 in the recorded exponent
ledger; it does not approach the square-root target. The D6 census now
covers all 28,774 eligible cases at orders 4-64, with 75 independent
direct comparisons. Growing-order D6, classical signed correlations,
full spectral bounds, and scalar prize bounds remain unproved.

The [twenty-seventh pass](parallel27-pass-summary-2026-09-05.md)
proves that a distinct opposite-free six-word is balanced at at most one
split. Exact unmarking of repeated words and the existing coset formula
then compute D6 in O(n^(49/20)polylog n) field/comparison operations,
with the subgroup supplied and O(n²) storage. This is a computational
bound, not an upper bound on D6. A separate direct program matches all
17 cases, including six new samples at orders512 and1024. Three algebraic
statements pass Lean; the full counting formula has an ordinary proof.
The six new samples have E3=T6+D6, but provide no uniform or worst-prime
estimate. The signed, full spectral, and scalar prize bounds remain open.
The previously queued private incidence lemma has now received an
ACCEPTED Prove2Me verdict, for five accepted elementary proofs in total.
This hosted theorem is not the code-specific reduction or full conjecture.

The [twenty-sixth pass](parallel26-pass-summary-2026-09-05.md)
formalizes one-projection preservation of every MCA-bad scalar and exact
equivalence of uniform count bounds, including linear-code agreement
witnesses. Seven principal statements pass Lean with standard axioms only.
The pinned ArkLib source already proves exact MCA interleaving invariance;
the pass-25 MCA deduction is a rediscovery. This corrects its attribution
without changing its validity. Scalar MCA/list bounds and every full
target remain unproved. The list transfer still has only an author audit
and finite checks. The private Prove2Me statement job was pending then;
the three workers are stopped by usage limits.

Sigma pass (2026-09-05, [summary](sigma-pass-summary-2026-09-05.md)):
the exact form of the obstruction is now written down. Every
termwise-Weil dilation-moment bound depends on (A,B) only through
dilation-invariant ratio-class sums, so it cannot distinguish t=1 from
the other p−1 dilates; the needed input must depend on r_{A+B}. The
subgroup case B=H is equivalent to the open uniform shifted-subgroup
bound M_H≤p^(−δ)|H| (Bourgain's Problem 5), so even one structured set
is out of reach. No published exponent below 1/2 exists for two arbitrary
sets; Chang's 4/9 requires small additive doubling, and shift chains
cannot save more than a constant on additively unstructured sets. The
subgroup exponent 8/9 is Konyagin 2002 at k=l=3 and becomes unconditional
exactly when E₃(μ_N)≤N^(3+o(1)) holds at every quartic prime; even ideal
fourth energies cap that method at 7/8. Passes 9–10 survive 1.4 million
exact checks with citation-level gaps only. Conjecture SI (biased dilate
sets have size ≤8 log₂p) was refuted by the second wave (p=97, A=2H,
B=H, |H|=8, |D|=88) and replaced by SI* with a signal-to-noise
hypothesis S²≥8R ln p, supported on 371,000 rectangles and unproved. The
conjecture is exactly the statement that the non-square Weil terms of the
dilation 2k-th moment over t∈F_p^* sum to at most p^(−2kδ')(mn)^(2k);
the only dilation-breaking data there are the signed ratio-class sums.
The third wave refuted the candidate sufficient condition T(k) for every
constant (A=B=Q gives ratio (2k−1)!!√p; subgroup pairs give ratios
growing like |H|^0.65) and refuted SI* (designed 1×n rectangles with
|D|/log₂p growing to 6.2 at p≈10^6). The four elementary Lean theorems
are now ACCEPTED by prove2me's server in a private mission draft.

The [twenty-fifth pass](parallel25-pass-summary-2026-09-05.md)
removes interleaving from the pinned combination-round count condition:
B8+L16<=R iff B1+L1<=R, over the same full extension field. MCA bad-scalar
maxima are exactly invariant under interleaving for any linear code;
the list projection condition follows from the verified official budget.
The scalar quantities and separate spot check remain unbounded. This
root-authored reduction has finite checks but its separate-author review
was interrupted. For growing subgroup orders, repeated opposite-free
six-words are O(n^(69/20)(1+log n)^(1/5)); the remaining six-distinct
fully unbalanced count is unbounded and the controlled error still
exceeds cubic scale. The actual spectral operator splits into all S3
sectors; a uniform-vector iteration sees only about one-sixth of the
space at any depth. Its other sectors need an unproved 2/3 elliptic
operator bound, and the trivial sector retains a nonzero rank-two border.
Signed inversion averages control asymmetry but retain the unknown M6.
All four author verifiers pass; root separately reviewed the three worker
proofs after the workers reached usage limits. No full target or new
worst-case cancellation exponent follows.

The [twenty-fourth pass](parallel24-pass-summary-2026-09-05.md)
proves a bounded unsigned B_h representative in every inversion family.
At k=floor(p^1/4), every input has a Sidon inverse with
M6<=(20+o(1))pk^3. The transported weights chi(d) remain uncontrolled;
this does not supply the uniform completion input in SS-B*. Exact local
drift and shell estimates explain the loss in one attempted way of
eliminating exceptions. In the spectral lane, the next coefficients
are a=1/2+O(p^(-1/2)) and c^2=3/64+O(p^(-1/2)); upper estimates
hold on the actual span{u,Bu}, retaining its leakage. The full operator
still needs growing-depth and remaining-sector control. For subgroup
orders 4,8,16, all quartic primes are classified: R6=480 only at the
three order-16 primes 33713,37201,41521 and zero otherwise. This finite
result does not extend automatically to growing orders. All four lanes
have a reviewer distinct from the author; no worst-case exponent or
official prize result follows, and the full goal remains unproved.

The [twenty-third pass](parallel23-pass-summary-2026-09-05.md)
proves an actual upper bound for almost all Sidon inputs at n=floor(p^1/4):
M6<=15pn^3, hence L_out<=(2/3)p^2n^3, with exceptional proportion
at most (216/5+o(1))n^2/p. The universal exceptional-set bound remains
unproved. For actual subgroups, product-balanced opposite-free sixth
relations are O(n^3(1+log n)); the unbalanced remainder remains open.
The two-anchor uniform direction has projection coupling tending to
1/8, so deleting its cross terms does not have vanishing norm error.
The full spectrum remains uncontrolled. A higher-coefficient restricted
list bridge is proved conditionally, and a maximum polynomial-phase
extension is obstructed at fixed gap. Exact arbitrary-word remainder
counts needed by the pinned official certificate are unbounded here.
No new worst-case cancellation exponent or full goal is established.

The [twenty-second pass](parallel22-pass-summary-2026-09-05.md)
audits the completed parallel results. The uniform vector on every
actual two-anchor prime neighborhood has residue Fourier mass
3/8+O(p^(-1/2)), separated from the endpoints for p>=73. The
other directions and coupling are not controlled. In the classical
lane, the missing positive T6 bound is equivalent to an off-set
quartic-star energy estimate. In the subgroup lane, an exact relation
decomposition isolates opposite-free six-term counts. A weighted
model passes all individual orders and all lower even moments for
arbitrary coefficients but fails the next moment. It lacks actual
field-translate structure. All decisive upper estimates remain open.

The [twenty-first pass](parallel21-pass-summary-2026-09-05.md)
shows that the moment target is equivalent, up to fixed-order constants,
to a one-sided upper bound on its top squarefree correlation aggregate.
The aggregate's negative side already has Gaussian scale. A rationally
weighted sign kernel preserves the exact Gram matrix, the Gaussian
fourth-moment order and individual bounds through degree six, but its
sixth moment is too large on k=floor(p^(1/4)). This is a relaxed-model
obstruction, not a character-kernel counterexample. The actual positive
upper bound remains outstanding. A separate-agent review found no
mathematical correction; human review and formal verification remain open.

The [twentieth pass](parallel20-pass-summary-2026-09-05.md)
removes the partition loss in the restricted moment transfer. A good
inversion pole makes the entire input B_h; an exact missing-row term
and a size-2k completion handle the introduced signs. Restricted and
unrestricted size-k moment bounds are equivalent up to a factor
3^(2r), under explicit finite conditions. On k=floor(p^(1/(r+1))),
these hold eventually for 2h≤r+1, including the boundary. The new
criterion SS-B* permits h roughly r/2, but its uniform upper bound
remains unproved. No new cancellation exponent or broader target is
established; independent mathematical review remains outstanding.

The [nineteenth pass](parallel19-pass-summary-2026-09-05.md)
proves that the full all-exponents classical conjecture is equivalent
to its restriction to B_h input pairs, for each fixed h. A greedy
partition with an explicit remainder also shows that SS need only be
proved on B_(h_j) sets when h_j/r_j→0. This SS-B hypothesis remains
unproved; minimal additive energies alone do not establish it. The
extraction loss prevents using a Sidon-only fourth-moment bound to
claim the earlier ε>1/3 consequence. No new cancellation exponent,
clique bound or official prize bridge follows. Independent review
remains outstanding.

The [eighteenth pass](parallel18-pass-summary-2026-09-05.md)
proves uniform proper-coset escape for Cartesian Bruhat families and
applies an existing SL₂ expansion theorem. The Weil trace represents
the Paley sum plus an exact overlap correction. A genuine conjugacy
average has trace one and norm about 2/p; it is not asserted to be
Cartesian and is not a Paley counterexample. Generic operator control
does not give the centered trace estimate. No new cancellation or
clique bound follows, and the source-dependent application needs
independent review. The three workers remain stopped at usage limits;
root completed this pass locally.

The [seventeenth pass](parallel17-pass-summary-2026-09-05.md)
proves an explicit prime-field endpoint gap using Fourier minor
nonvanishing and the integer p^(m+1)*det(X)*det(Y). The exact Paley
second moment bounds the stronger determinant certificate by an
exponentially small quantity when m is proportional to p. This bounds
the certificate, not the actual eigenvalue gap. It cannot give the
positive constant needed for the spectral edge. The next estimate
must use quantitative quadratic-residue Fourier information on the
actual anchor neighborhoods. No clique bound or full target is proved.

The [sixteenth pass](parallel16-pass-summary-2026-09-05.md)
shows that actual Paley graphs over square finite fields retain both
the sharp elliptic correlations and full projection identity but still
have persistent localized outliers. Their subfield clique gives a norm
limit 2^a/(2*sqrt(2^a-1)) at fixed anchor depth a, and exponential
long-power energy and even-trace obstructions. These are not prime-field
counterexamples. The next argument needs a quantitative input exploiting
prime field order. No stronger clique bound follows; the prior depth
range and other open targets remain unchanged.

The [fifteenth pass](parallel15-pass-summary-2026-09-05.md)
bounds every translated elliptic product with arbitrary character twists:
s odd-multiplicity shifts and e positive even-multiplicity shifts give
s*sqrt(p)+e, with an exact formula when s=0. An abstract sign-kernel
construction retains these square-root orders with larger constants,
the full exceptional border and four symmetries, but has persistent
normalized outliers. It loses the sharp constants and ambient Paley
identity and is not a Paley counterexample. The next norm argument
needs more of the combined arithmetic constraints. The previous
logarithmic energy range and clique bounds are unchanged. Independent
review remains outstanding; the subgroup exponent stays 8/9 on its
existing density-one quartic prime class. Historical sections retain
earlier limitations.

The primary prize-linked target is now documented in
[`subgroup-target.md`](subgroup-target.md), after recovering the existing
repository. Its central obligation is (SG), a centered additive-energy bound
at logarithmic moment depth. The notes below track that target and the
distinct classical quadratic-character conjecture; their quantifiers and
connections to the prize remain separate.

The latest [geometric-cycle analysis](geometric-lift-and-alias.md) proves a
uniform moment estimate for an auxiliary lift of any even-order subgroup.
The prime-field estimate still requires an unproved upper bound on the
signed reduction error (D); treating it as automatically nonpositive is false.
That error estimate is equivalent to (SG) up to constants.

The [coset-coherence analysis](coset-coherence.md) identifies exact nonlinear
relations between periods and the cosets of `1+H`. Linear autocorrelations
and Gauss-sum magnitudes by themselves admit synthetic almost-maximal spikes;
those vectors are excluded by the nonlinear identity. No resulting bound
on actual periods is proved. A necessary lower-order consequence of the
target is `E_2=O(n² log n)` in the quartic window, also unproved here.

The [polynomial certificate](period-polynomial-certificate.md) now uses all
signed moments through degree 24 to identify the exact maximizing coset in
the resonant test field. The separate
[recurrence coefficient criterion](jacobi-coefficient-route.md) states a
uniform sufficient estimate with its implication proved. Its hypothesis
remains open; the verified single field is not evidence of uniform closure.

The proposed literal coefficient constant `B=1` is now refuted at
`p=67403009`, `n=128`, already for β₂. Uniform coefficient bounds also
require `E₂=O(n²)` in the quartic window. The new
[cyclotomic norm argument](cyclotomic-prime-average.md) proves a weighted
bound on extra relations summed over primes, and bounds the number of
fourth-energy exceptions. It supplies neither worst-case control nor
the logarithmic-depth centered estimate; its higher-order numerator is
too large. The recent exact-decomposition polynomial theorems do not
apply to these collision counts without a further bridge.

The refined norm estimate (N') still gives a pointwise fourth-energy
upper-bound expression larger than n³ for dyadic n≥16 in this window.
For dyadic n≥1024, extra zero ten-tuples are forced at every eligible
prime p≤n⁴ by the principal-frequency term. Thus high-order extra
relations cannot be treated merely as rare prime exceptions.

The [official profile audit](official-profile-and-trace.md) has now pinned
the newer benchmark and its relevant dependencies. Its domain lies in
the KoalaBear prime field inside a sextic coefficient field; neither field
size is in the present quartic window. Nonzero trace-zero frequencies
give periods of size n. A correct extension-field spectral argument must
account for the full trace-zero space, and a bridge to the actual
MCA-plus-list certificate remains unproved.

The [subset-sum analysis](subset-sums-and-lists.md) now proves a limited
bridge: prime-field character cancellation controls distinct subset sums,
which give exact lists and bad scalars for a specific monomial pair one
step below capacity. Root lifting gives a finite list lower bound of at
least `2^8154` at radius `4095/8192` in the official profile, excluding
that radius from its MCA-plus-list lower certificate. It does not bound
winning-set soundness or arbitrary-pair MCA. For fixed dyadic image order,
odd zero-sum subsets disappear at sufficiently large primes by a nonzero
cyclotomic norm bound; thus this construction cannot reach a fixed-gap
asymptotic counterexample. General vanishing polynomials at larger gaps
impose several coefficient constraints, beyond the linear subset-sum
condition. Those constraints and the uniform (SG) estimate remain
uncontrolled.

The [large-spectrum/Riesz audit](riesz-tail-route.md) proves that the direct
dissociation entropy estimate remains compatible with M=n/2 even when a
dissociated subset of the maximizing orbit has the largest possible size.
Thus this particular argument cannot yield M=o(n). A whole-subgroup Riesz
product avoids that dimensional loss but requires a nonprincipal mean
bound. The stronger proposed mean bound Z_+(t)≤1 is refuted exactly at
p=67403009, n=128, t=1/4. The corrected critical-scale condition
Z_±(t)≤exp(Cnt²) is equivalent to the desired maximum bound up to constants,
so it is not a weaker remaining obligation. Arithmetic dependence still
needs to be bounded, not replaced by independence.

The [analytic bound audit](analytic-bounds-and-amplification.md) establishes
an applicable literature baseline: the 2020 Di Benedetto et al. theorem
gives `M≤n^(2849/2880+o(1))` throughout the working quartic window. This is
an existing bound, not a new or claimed current-best result. For a stated
abstract extension of its three-selection argument, even granting optimal
full-energy exponents makes the saving in that ledger at most 1/16.
The claimed existence of general selections is not used. This calculation
limits the formula's output, not all amplification methods.

A 2018 centered-moment formula would give a best fixed-depth saving of
11/1280 with the recorded E₂ input, still weaker than the baseline. Its
unqualified general-function statement as rendered in v2 HTML admits a
mass-at-zero counterexample, so that statement is not adopted as an
independently checked dependency; the subgroup formula itself is not
refuted. The 2025/2026 Gaussian-equidistribution paper treats prime orders
far smaller than ours and does not control the maximum here. Arithmetic
control beyond these direct applications remains necessary.

The [mixed-period analysis](mixed-periods-and-shifted-energy.md) now
retains every product x_j x_(j+t), with cyclotomic matrix C. The normalized
mixed system has exactly the actual shifted period vectors as solutions;
its integer multiplication matrix has precisely the nonprincipal periods
as eigenvalues. This is an exact representation, not an estimate.

An explicit bijection proves
`X(H)=E×((H−1)\{0})−(2n²−5n+3)=Σ_(t,s)(C_ts)_3`.
It counts all ordered triples in intersection cells, with the two trivial
multiplicative matchings removed. Parity of row zero then gives
`E₂(H)≤3n²−n+nX(H)/3`. The known bound `X≪n²log n` is too weak in
this comparison; uniform `X=O(n log n)` would give the necessary fourth
energy bound, but is neither proved nor shown necessary for the spectral
target. Concentration in row zero may be a weaker quantity to control.
Exact tests certify X=0 and minimum fourth energy for n=256 and 512 in
two quartic-window primes. Existing resonant examples have X=114 and 720,
so uniform circularity is false. None of these results gives (SG).

The [kernel discriminant](kernel-discriminant.md) now gives an exact
criterion for fourth-energy excess: `p | A_n`, where
`A_n=disc(P_n)P_n(2^n)` and the roots of P_n are `(1+ζ_n^j)^n`, one per
inverse pair excluding ±1. Complete integer factorizations classify all
exceptions through n=32; none lies in those orders' quartic windows.
The order-32 conclusion is complete, not a prime-prefix scan. The exact
identity `4v_p(A_n)=Σ_r D_r` includes higher p-adic collisions, which
cannot always be omitted. The elementary height bound is too large.

The [quadruple orbit argument](quadruple-orbits-and-cube.md) further proves
`D_r=4ε_r+12u_r+24v_r`. Here ε records three equal entries, u counts
orbits with one repeated entry, and v counts orbits of four distinct
elements, after removing opposite pairs. The first two types contribute
only O(n²) energy. Thus `E₂=O(n²log n)` is equivalent to
`v₁=O(n log n)` in the working window. This estimate remains unproved.
The same argument over unramified extensions proves the uniform integer
factorization `odd(A_n)=odd(3^n−1) U_n³ V_n⁶`, with
`v_p(U_n)=Σ_r u_r` and `v_p(V_n)=Σ_r v_r`. It supplies structural
information but no useful bound on the four-distinct contribution.
Factoring out these powers changes the height estimate by constants,
not the missing power of n. Exact checks include non-splitting primes
and the unfactored 40003-bit A_64. Neither the fourth-energy estimate nor
the higher centered-moment estimate (SG) follows from these identities.

The [dyadic descent analysis](dyadic-descent-and-mixed-energy.md) now
proves `E₂(H_(2k))=2E₂(H_k)+6B+8T` and
`T²≤(E₂(H_k)−k²)B`, with B the mixed additive energy of the two cosets.
Uniform `B≤Ck²Λ` at every level of a target tower would imply
`E₂(H_s)≤22Cs²Λ`; taking Λ=max(1,log p) would give the required
fourth-energy consequence. This is a whole-tower sufficient hypothesis,
including smaller orders outside the fixed quartic window, and remains
unproved. The term B requires control only of the new balanced
four-distinct orbits; the other new orbits could be handled by Cauchy–Schwarz.

When `ord_(2k)(p)=2ord_k(p)`, a field basis separates the cosets and
forces `B=k²,T=0`, with no new nontrivial orbits at any p-adic precision.
This cannot be applied in the target, where all those orders equal one.
It yields the new-prime support restriction
`p≡1 or n/2−1 mod n` for U_n/U_(n/2),V_n/V_(n/2), plus stabilization
at fixed p. None of these restrictions excludes target primes, which
already satisfy p≡1 mod n. The norm cutoff `p^ord_n(p)≤2^(n/2)` for
four-distinct exceptions likewise ceases to help in the quartic window
once n≥64. Exact checks cover 44 field/precision cases, including ten
with doubled field degree; the mixed-energy induction is conditional.

The [positive-product moment criterion](positive-product-moments.md)
extends this recurrence to the required high depth. With normalized
nonprincipal norms, `A_s(q)=||η_s||_(2q)^2` satisfies
`A_(2k)≤2^(1/q)A_k+3||(η_k(b)η_k(bg))_+||_q`.
Uniform positive-product bounds `≤Cqk`, at one even q comparable to
log((p−1)/N) and along the needed tower, would give
`Q_q(H_N)≤m(4CqN)^q` and the target maximum bound. The low steps k≤q
already satisfy the hypothesis trivially. No uniform estimate for the
remaining steps is proved.

The stronger absolute-product condition is exactly
`p Z_(q,q)(K,gK)−k^(2q)≤(p−1)(Cqk)^q`, for even q. It is checked at
q=12 with C=1 throughout the p=6700417,N=64 tower, using two new
integer quotient computations and three trivial lower-level bounds.
This is finite verification of the criterion, not uniform closure.
The raw complex count `T_q(k)^2` cannot replace the centered expression:
already at q=6, every eligible prime p≤(2k)^4 with dyadic k≥64 has
`Z_(6,6)>T_6(k)^2`, forced by the principal term. Cauchy–Schwarz on the
centered off-diagonal correlation is still too weak for the induction.
The next obligation is genuine arithmetic control of positive product
moments or that centered balanced correlation at logarithmic depth.


The [signed quotient construction](signed-quotient-operators.md) gives
sparse symmetric integer matrices with exactly the nonprincipal period
spectrum, and similarly for the two-child product. The principal frequency
is absent directly. For integer H-invariant weights of mass d and squared
mass E, an exact Frobenius identity bounds how many entries can cancel.
Consequently the entrywise absolute period matrix has spectral radius at
least `n−2/n` throughout the quartic window; the corresponding product
matrix has radius at least `k²−k/4`, where n=2k. These are nearly the trivial
degrees, so discarding signs cannot provide the needed estimate.

Changing coset representatives merely conjugates the matrix by diagonal
signs and leaves every closed-walk product unchanged. Exact counts also
show all nonidentity multiplicative endpoint ratios have equal numbers of
closed quotient walks. The remaining signed trace is precisely the centered
count already needing a bound. The identities, checked in 14 matrix cases
with 94 trace comparisons, do not establish cancellation at high depth.


The [list-to-winning-set analysis](list-to-winning-set.md) now gives a
verified ordinary-mathematics consequence for the pinned prize benchmark.
For an injective linear encoder, a single-word list of size L gives a
winning-set lower bound L/(q+L−1) by taking the second word zero and claims
(0,1); if binom(L,2)<q, the bound improves to L/q. The same violating
instance works throughout the larger-radius suffix below minimum distance.
This does not assume the maximum list is smaller than q.

A base-field coefficient-pigeonhole construction gives a list of
85677801616821870413970774 at radius 122641/262144 for the pinned profile.
Its projection can be injective, and its winning density exceeds 2^-128.
Exact integer comparisons verify the suffix and the score inequality with
11649 centibits. This is optimal only within that particular pigeonhole
lower bound, not an optimal code or prize threshold. The proof is not Lean
formalized or submitted. The July 6 ABF edition has now been recovered and
its relevant pages checked; the full Paley and prize targets remain open.

The alternate [localization analysis](localized-necklace-identities.md)
now evaluates an infinite subfamily of classical quadratic-character
necklaces. For two singleton anchors, one occurrence gives `N=-t_k`;
two occurrences separated by ℓ edges give an exact expression in the
monochromatic t_j, including both exceptional directions and zero-coordinate
corrections. The Lu–Zheng–Zheng twisted Kloosterman estimate implies
`|N|≤k²p^(k/2)` for every k and prime p≡1 mod 4. This covers all binary
singleton patterns through length five, not all degree-two label patterns.

At three occurrences the exact averaging formula is a cubic character
correlation weighted by three explicit kernels. Termwise Weil and Frobenius
bounds give only `O(p^(k/2+1))`, with no saving at the required scale.
Exact checks cover 3294 necklaces, 420 kernel trace products, 72 literal
tuple sums, and six cubic averaging identities. The missing signed bound,
the general localization estimates, and extreme-eigenvalue control are
not established. No thin-subgroup or prize reduction is inferred.

The subsequent [planar reductions](planar-necklace-reductions.md) resolve
the length-six singleton examples that first appeared in that cubic sum.
Affine averaging proves `N(0^r1^s)=p t_r t_s-t_(r+s)` for arbitrary
positive block lengths. The other balanced length-six classes satisfy
`N(001011)=p(t_5+2t_3)-t_6` and
`N(010101)=p(p-5)t_4+p²(p-2)-t_6`.
These follow from planar Fourier duality and projective invariance, with
all infinity/deletion terms retained. Together with the earlier formulas,
they prove `|N(w)|≤36p^(7/2)` for all 64 binary singleton words of length
six. The direct absolute-value barrier is thereby bypassed at this depth.
General longer mixed patterns, words using the doubleton label, and the
spectral-edge/growing-depth requirements remain unresolved.

The subsequent [parallel pass](parallel-pass-summary-2026-09-04.md) extends
this scope. An inversion transfer controls the change when a singleton and
doubleton label are exchanged, at a cost at most k p^(k/2). Thus 189
length-six words using at most two of the three nonempty labels now satisfy
the uniform bound 42p^(7/2); 540 words using all three labels are not covered.
The exact identity N({0},{1},{0,1})=t_4 also completes all degree-two words
through length three. General growing-depth cancellation remains missing.

In the subgroup lane, projection onto the index-two symmetry improves the
fourth-energy propagation constant from 22 to 7+4 sqrt(3), with a centered
version for higher convolution counts. It supplies no uniform balanced-count
bound. In the classical lane, biased row patterns sharpen the necessary
logarithmic allowance but do not prove or refute (LM). On the prize side,
base-field-valued received pairs have MCA-bad challenges only in the base
field, hence error at most p/q=p^-5 for the pinned sextic profile. This is
not an upper bound for arbitrary extension-valued pairs or for protocol
winning sets. The linked pass assessment gives the separate proofs.

**Unresolved:** cancellation for arbitrary two sets of size `p^ε`, especially
when both are smaller than `√p`. No mathematical obstruction encountered so
far blocks further research; the goal is neither complete nor marked blocked.

## A sufficient moment estimate with necessary spike growth

The preceding analysis eliminates constant-coefficient Gaussian corrections.
A possible replacement to investigate is the following **unproved hypothesis**:
for every fixed integer `r≥2`, some constant `C_r` satisfies

\[
M_{2r}(B)\le C_r\left(p|B|^r+|B|^{2r}\log p\right)
\tag{LM}
\]

for every odd prime `p` and every nonempty `B⊆F_p`.
Here `log` is natural logarithm. The common-neighborhood lower bound explains
why a logarithmic term is necessary in some ranges; it does **not** establish
the upper bound (LM). This is a candidate research hypothesis, not an assertion
that it is a standard equivalent formulation or that it is true.

### Conditional implication, with exponents

Suppose (LM) holds. Fix `0<ε<1` and choose an integer `r≥2` with `rε≥1`.
For `m,n>p^ε`, Hölder gives

\[
\left(\frac{|S(A,B)|}{mn}\right)^{2r}
\le C_r\left(\frac p{mn^r}+\frac{\log p}m\right)
\le C_r p^{-\varepsilon}(1+\log p).
\]

For sufficiently large `p`, the last expression is at most `p^{-ε/2}`.
Therefore the Paley conjecture follows with `δ=ε/(4r)>0`.
**All the new difficulty is in proving (LM).** This implication is not an
unconditional result about Paley. The fourth-moment identity is one starting
point for investigating (LM), but one fixed moment cannot cover all `ε>0`.

## An equivalent tail formulation that avoids Gaussian assumptions

Let `T_B(τ)={x: |F_B(x)|>τ|B|}`. A precise qualitative reformulation of
the full Paley conjecture is:

For every `0<ε<1`, there is `η>0` such that, for every sufficiently large
prime `p` and every `|B|>p^ε`,

\[
|T_B(p^{-\eta})|\le p^{\varepsilon/2}.
\tag{Tail}
\]

**(Tail) implies Paley.** For `|A|>p^ε`, split the sum over `A` into
the exceptional rows and the remainder. The normalized absolute sum is at
most `p^{-η}+p^{-ε/2}`. For sufficiently large `p` this is bounded by
`p^{-δ}` for any fixed `0<δ<min(η,ε/2)`.

**Paley implies (Tail).** Apply Paley at exponent `ε/4`, obtaining `β>0`,
and use threshold `p^{-β}`. If the tail has more than `p^{ε/2}` elements,
one sign class has more than `p^{ε/2}/2 > p^{ε/4}` elements for sufficiently
large `p`. Summing over that sign class and `B` contradicts the Paley bound.
Thus (Tail) holds with `η=β`.

This is only a reformulation, not independent progress. A bound for complete
moments that assumes (Tail) would be circular as a proof of Paley.

## Second parallel pass: the current aggregate obligation

The [new growing-family proof](parallel2-necklace-2026-09-04.md) uses Katz's
twisted Legendre and symmetric-square Mellin bounds to control A^rBC and
A^rBAC for every r>=1. Both normalized bounds tend to zero uniformly at
logarithmic length. At that stage the general two-gap identity was proved,
but gaps of at least three remained outside the result; the third pass
below closes that family. Arbitrary words remain outside the result. The
quadratic Mellin mode of size p^2 is subtracted explicitly; ignoring it
would invalidate the smaller operator estimate.

The [exact spectral transfer](parallel2-spectral-transfer-2026-09-04.md)
defines W=(sqrt(p)S-J)(2^a E-I), where E indicates the common neighborhood
of an a-clique. It proves an exact two-projection/Chebyshev trace identity,
including all rank and anchor terms. A bound on
abs(tr(W^k))/(p sqrt(2^a-1))^k with subexponential growth at even
k much larger than log p gives the explicit clique bound in that note.
At a=floor(sigma log_2 p), this would imply O(log p+p^((1-sigma)/2))
for each fixed sigma<1. The aggregate hypothesis remains unproved. Even
uniform individual necklace bounds of the usual square-root type lose
an exponential factor when combined by the triangle inequality.

In the [subgroup lane](parallel2-subgroup-2026-09-04.md), the primitive
polynomial's fiber multiplicities give B=2k sum(c_i^2). A defect bound
or absence of triple roots is sufficient for near-quadratic energy, but
neither is established uniformly in the quartic family. The new finite
witness p=17189277697,n=512 has B=67584>65536 and two double fibers.
This refutes collision-freeness in that order, not the subgroup target.

The [classical fourth-moment audit](parallel2-classical-2026-09-04.md)
retains an exact signed elliptic cross-ratio pairing and subtracts its
random-set mean. The proposed sufficient complete-L2 bound for the
remaining coefficient vector is false on an unbounded family of intervals.
Their fourth moments already satisfy the ordinary Weil estimate. This
rules out the L2 route without refuting (LM); directional cancellation
in the signed pairing is still needed.

## Third parallel pass: arbitrary gaps and the same aggregate obligation

The [arbitrary-gap proof](parallel3-necklace-2026-09-04.md) now establishes
abs(N(B A^(j-1) C A^(m-1))) <= (min(j,m)+1)p^((j+m+1)/2)+1 for every
positive j,m and every prime p=1 mod 4. Katz's equal-rank hypergeometric
sheaf has the exact product-convolution trace, including t=1. Pairing a
rank-j kernel with the rank-two kernel leaves a cohomology space of
dimension j+1 whenever j!=2. The previous explicit correction treats
j=2. The source has no rank<p restriction. Label permutations give the
same normalized vanishing for lengths o(sqrt(p)). Both exceptional labels
still occur just once; arbitrary word patterns and their aggregate remain
uncontrolled.

The [conic representation](parallel3-conic-2026-09-04.md) changes the
two-anchor problem into the kernel f(t/u)f(tu), keeping an explicit
coupled boundary block and the centering term. Its Mellin matrix has
Jacobi-product entries and nonzero off-diagonal couplings. The six
anharmonic automorphisms give exact smaller blocks, but no improved
uniform bound on their extreme eigenvalues.

The [subgroup index identity](parallel3-subgroup-2026-09-04.md) gives
B<=k^2+4k v_p(I_(2k)), where I is the index of one explicit cyclotomic
order in another. A pointwise O(k log k) index valuation would suffice
for near-quadratic mixed energy at the relevant endpoints, but is unproved
and can be stronger than the first-precision collision requirement.
Quartic examples retain repeated image values despite a reduced ambient
algebra, nonzero derivatives, and a global cyclic Galois group.

The [classical Fourier calculation](parallel3-classical-2026-09-04.md)
turns the elliptic pairing into an exact Kloosterman-weighted Fourier sum.
Its unrestricted Fourier-L1 majorant is false even after any fixed number
of exceptional cross-ratios are removed. The counterexample primes have
no polynomial upper bound in the set size. Consequently this does not
rule out a norm bound restricted to n>=p^epsilon for fixed epsilon,
which remains relevant to Paley. No moment or Paley conjecture is refuted.

The [third-pass assessment](parallel3-pass-summary-2026-09-04.md) separates
these proved facts from the remaining uniform estimates. The signed
aggregate criterion of the second pass is still the main clique obligation.

## Fourth parallel pass: growing classes and an aggregate across ranks

The [fourth-pass assessment](parallel4-pass-summary-2026-09-04.md)
records three new proved inputs and a structured-set corollary. They
change which subproblems are open, without closing the main obligations.

The [subgroup result](parallel4-subgroup-2026-09-04.md) proves a lower
limiting proportion at least 1-log(2)/36 of primes p=1 mod N in
N^4/4<=p<=N^4 with E_2(H_s)=3s^2-3s at every dyadic divisor s>=2.
It combines the normalized cyclotomic norm budget, a negligible
repeated-entry prime set, 24N divisibility of the remaining extra
quadruples, and the checked Thorner-Zaman prime count for dyadic moduli.
A density-one theorem also controls the whole tower with vanishing
relative fourth-energy error. Worst-case primes and centered high
moments are still unbounded; fourth energy does not determine the
maximum period. A useful next step is to test whether the norm/orbit
budget extends to centered higher relations without counting their
necessary principal-frequency contribution as an error.

The [necklace result](parallel4-necklace-2026-09-04.md) proves
abs(N(A^r B^h C)) <= (s+1)p^((k+1)/2)+abs(t_s-(-1)^s),
s=min(r,h), k=r+h+1. This permits growing repeated B labels. At length
five it covers 30 of the 90 words of multiplicities (2,2,1). The
remaining representatives AABCB and ABABC have an exact two-variable
quartic-kernel reduction with two retained zero-coordinate terms.
That kernel is not a function of the ratio, and both natural graphs
have K_3,3 minors. A new two-variable trace or operator estimate is
needed; the direct planar proof does not apply.

The [kernel aggregate](parallel4-kernel-aggregate-2026-09-04.md)
controls every Mellin twist of |sum_(r<=R) a_r h_r(t)|^2, with error
[(3R^2-R-2)/(2sqrt(p))+2/p]||a||_2^2. It is useful for R=o(p^(1/4)).
The diagonal rank correction includes both +p^(r-1) delta_1 and
-p^(r-1); dropping either changes the invariant contribution.
This is an aggregate across lengths, not across arbitrary label words.
For the next-anchor weight chi(1-t), Mellin L1 costs order sqrt(p),
removing the present saving. A direct weighted bound or cancellation
among Mellin coefficients is the next obligation. No lower bound on
the actual weighted sum is implied by failure of that upper estimate.

The [classical result](parallel4-classical-2026-09-04.md) gives a
power saving in the signed elliptic pairing for polynomial-size
progressions and controlled perturbations, using Burgess, and for a
range of proper rank-two progressions, using Alsetri-Shao. Its lower
bound 0<=R<=3pn^2 in the exact M4 decomposition permits the transfer.
It supplies no arbitrary-set (LM) estimate or refutation. Repeating a
structured-set theorem alone will not address the remaining quantifier.

## Fifth parallel pass: the translated weight is controlled, higher moments remain

The [fifth assessment](parallel5-pass-summary-2026-09-04.md) records an
actual stronger maximum bound on a density-one quartic class:
M<=C N^(23/24)(log N)^(7/72). The
[subgroup proof](parallel5-subgroup-2026-09-04.md) gives all-level
sixth energy E3<=(15+log N)s^3 outside O(1/log N) of eligible primes.
Combining it with fourth energy and the checked fixed-order trilinear
amplification gives the new exponent. The zero deletions, factorwise
phase weights, set-ordering condition and log^5 selection loss are
explicit. The conclusion is not uniform over its exceptional primes.

The next arithmetic budget is already unproved at order eight:
sum_(p in P_N)(E4(H_N)-N^8/p)log p << N^7(log N)^A.
The raw norm method gives scale N^8 and subtracting the principal term
does not remove that leading order. The new maximum gives centered
E4 << N^(59/12)(log N)^(43/36), still above the desired N^4 scale.
All higher centered bounds from interpolation retain the exponent
23/24, far above the target 1/2. The fixed quartic example p262657,
N32 disproves exact-fourth-energy implying exact-sixth-energy.

The [three-gap theorem](parallel5-necklace-2026-09-04.md) resolves the
two (2,2,1) patterns left open above, and all gaps a,b,c>=1. The inner
Mobius pullback has scalar quadratic inertia at infinity, so its
tensor has no invariant for any ranks. The full-field values at zero
are kept distinct from the zero-extended multiplicative convolution.
This proves all 243 individual words of length five, and 639 of 729
at length six. The next finite family is the 90 balanced (2,2,2)
words. None of these counts gives the signed aggregate or clique bound.

The [direct anchor theorem](parallel5-anchor-aggregate-2026-09-04.md)
also closes the preceding kernel note's next-anchor estimate: direct
Kummer twisting avoids the Mellin L1 cost. Products of m distinct
anchor weights have a quadratic kernel-span bound O((m+1)R^3/sqrt(p)).
It gives the expected mass in every adjacency cell with
a<=(1/2-epsilon)log_2 p, uniformly for R=O(log p). This controls
only the indicated R-dimensional span. Showing that it captures or
approximates extremal eigenvectors of the restricted Paley matrix is
unproved. Substituting R~p makes the estimate useless.

The [classical result](parallel5-classical-2026-09-04.md) reduces a
sufficient moment input to one set cardinality per fixed rank, along
an unbounded sequence. The single-size upper estimates are unproved.
An elementary construction rules out bounding a higher Gaussian
character moment solely from exact minimal lower additive energies,
even for polynomial-size sets at every sufficiently large prime.
It does not refute LM with its logarithmic spike term or the new
single-size sufficient hypotheses.

## Sixth parallel pass: a better maximum, and explicit limits of the kernel span

The [sixth assessment](parallel6-pass-summary-2026-09-04.md) improves
the subgroup exponent from 23/24 to 17/18 on the sixth-energy class
alone, with no additional fourth-energy intersection:
M<=2^(1/6)(15+log N)^(1/18)N^(17/18). The
[proof](parallel6-subgroup-2026-09-04.md) gives a centered mixed-order
weighted bilinear inequality valid for all subgroups and all r,s.
Its r=s=3 specialization converts E3<=Bn^3 and p<=n^4 directly into
Delta^18<=8B/n. The proof retains both the zero-sum atoms and the
negative constant term on nonzero coordinates. It uses no trilinear
selection theorem. Root and an independent reviewer accepted the
argument; the exceptional-prime proportion is still O(1/log N).

The new interpolation bound is centered E4<<N^(44/9)(log N)^(10/9).
It does not prove the averaged centered eighth budget stated above.
Feeding the available interpolated moments into the same coarse
mixed-order gate is proved to have maximum saving 1/18 over every
fixed pair r,s. A sharper eighth-energy bound E4-n^8/p<=Dn^4 would
instead give M<=3^(1/4)D^(1/16)n^(15/16) by the proved gate at
r=1,s=4. The averaged eighth budget thus has a concrete spectral
payoff, but no new arithmetic estimate establishes it.

The [necklace proof](parallel6-necklace-2026-09-04.md) gives
||W_C+S||<=2p, where W_C=S entrywise-times (S D_C S). A projective
map identifies this with a principal compression of an augmented
single-anchor kernel, including the pole. Its three added modes
are evaluated exactly. Normalizing adjacent C vertices then proves
abs(N(AABBCC))<=7p^(7/2), with an exact coincident-anchor correction.
Its twelve-word orbit raises the length-six catalog to 651/729.
The remaining 78 lie in the AABCBC, AABCCB, ABACBC and ABCABC classes.
The first two have exact contractions recorded in the note, but the
new norm gives only the order-p^4 scale for them. The signed full
aggregate remains unbounded.

The [seeded-kernel result](parallel6-seeded-kernels-2026-09-04.md)
settles one preceding uncertainty negatively: all old kernels are
invariant under Jf(t)=chi(t)f(t^(-1)), and at p=13 the omitted sector
contains the largest eigenvector on the common neighborhood {4,10}.
Thus the old span cannot universally capture the extremal eigenvector.
Multiple seeds repair that example and admit proved correlation and
cell estimates. With L<=p^sigma, R=O(log p), the relative cell error
tends to zero for a<=(1/2-sigma-epsilon)log2p and sigma+epsilon<1/2.
Seeds must be disjoint from nonzero cell anchors; a collision can
create an order-p rank-one correlation. Taking all seeds gives an
invertible matrix and a dimension obstruction to full-space cell
isometry. This does not refute a spectral bound for the restricted
Paley matrix; that still needs an independent argument.

The [classical result](parallel6-classical-2026-09-04.md) proves exact
fixed-cardinality fourth-moment means and variances, with variance
asymptotic to 24pn^4. At n=floor(p^(1/3)), M4<=3pn^2 holds for a
proportion 1-(6+o(1))p^(-1/3) of all sets, and the centered signed
pairing is similarly controlled for typical sets. This supplies no
uniform single-size estimate. A fixed p=1009 example refutes only
coefficient 3 on that slice. Global mean and variance allow rare
outliers and cannot control adversarial subsets by themselves.

## Other continuing work

1. Try to refute (LM) before investing in a proof: optimize higher moments for
   intervals, multiplicative subgroups, unions, and multirow sign constraints.
   Track `M_{2r}/(pn^r+n^{2r} log p)` across primes and fixed `r`. Finite large
   ratios do not refute an unspecified constant; an unbounded family is needed.
2. Analyze the signed quartic aggregate `T_4(B)` using its elliptic-curve
   interpretation. Any bound must retain the dependence among the four labels;
   summing absolute pointwise Weil bounds has already been accounted for.
3. Investigate tail bounds directly if the moment hypothesis proves too strong.
   Allow a growing exceptional set: a bounded number of removed rows is known
   to be insufficient.
4. Complete the exact prize reduction before claiming a connection. The April
   ABF source and the July revision's relevant attack/list pages have now
   been recovered and checked. The unrestricted monomial-pair reduction is refuted in
   `prize-reduction-audit.md`. State the field, evaluation domain,
   rate, dimension, error parameter, and every intermediate hypothesis.

No theorem currently supplies (LM) or (Tail) for arbitrary small subsets.


## Seventh parallel pass: exponent 8/9 and the next arithmetic input

The [signed multiplier theorem](parallel7-multiplier-2026-09-04.md)
gives M<=(17+log N)^(1/9)N^(8/9) on the existing sixth-energy class.
The centered eighth-energy feedback exponent improves to 43/9;
it is still above the N^4 scale. A new input centered E4<=DN^4 would
give M<=2^(1/4)D^(1/8)N^(7/8), but that input is unproved. The same
averaged centered eighth-energy budget remains a concrete next target.

The [subgroup audit](parallel7-subgroup-2026-09-04.md) proves the
current fixed real/integer order substitutions cannot exceed the
1/9 saving. It also certifies negative multiplier coordinates, so
an unrestricted second Jensen step with that signed function is
invalid. Neither statement excludes a method using additional
arithmetic correlations or cancellation in the bilinear expression.

The [crossing class](parallel7-necklace-2026-09-04.md) adds 36 words
with bound 58p^(7/2). The remaining fixed length-six classes are
AABCCB, ABACBC, ABCABC, containing 18,18,6 words. Middle convolution
and a local multiplicity mismatch now handle one equal-rank crossing.
The nested ABBA contraction and the two classes with no adjacent
equal labels remain open, as do growing-depth signed aggregation
and the full restricted operator outside the seeded kernel span.

The [classical countermodel](parallel7-classical-2026-09-04.md)
shows why mean, variance, affine invariance, quartic degree,
nonnegativity and square-root coefficient bounds alone cannot
supply a worst-case fourth moment on the p^(1/3) slice. The required
extra input must use actual character-sum structure. This is not a
counterexample to the desired moment statement or the conjecture.


## Eighth pass: complete length six, then growing depth

The [two-convolution theorem](parallel8-necklace-2026-09-05.md)
handles the remaining nested and nonadjacent classes. The intermediate
rank-four sheaf has types J3+1, J2+J2, chi+1^3 and 1+chi^3 at the
two anchors, moving point and infinity. Either singleton twist then
middle convolution has rank six: finite invariant codimensions
4+2+1, minus one twisted infinity invariant. Pairing it with the
rank-two quartic Legendre sheaf excludes an invariant by rank.
All correction terms are retained in the resulting bound.

All 729 degree-two words at length six now satisfy O(p^(7/2)).
The last three class constants are 48,47,45. These are fixed-length
bounds, not an aggregate estimate. The next structural question is
whether longer-chain rank and monodromy transformations can be
controlled uniformly, including changes in finite stalks when a
Kummer twist creates new invariants. Any such extension must also
retain arithmetic constants and quantify the weighted aggregate;
counting or individually bounding a finite catalog does not do so.

The seventh-pass subgroup exponent 8/9, its exceptional-prime gap,
centered eighth-energy budget, classical exceptional-set obstacle,
full restricted operator and official prize bridge are unchanged.


## Ninth pass: all individual necklaces, then signed aggregation

The [general chain proof](parallel9-all-degrees-2026-09-05.md) replaces
fixed-degree catalog enumeration with a rank induction. At every finite
anchor, the difference between the numbers of trivial and quadratic
Jordan blocks is at least the preceding rank increment; the reversed
difference holds at infinity. A nonempty twist flips at least two
of these places, so the next rank increment cannot decrease and is
at least one. The principal factor can never return to rank one.

The complete raw chain is an exact extension of that principal
middle sheaf by a perverse error of strictly lower weight. Its mass
is at most (2a+2)^j after j steps. This proves the bound
3a(2a+2)^(k−2)p^((k+1)/2) for k≥3 and every nonempty-label word
on a anchors. It supplies Kunisky's fixed-degree necklace conjecture
and, through the published implication, weak spectral convergence.
It is a locally checked source-dependent proof, not an independently
reviewed or formally verified result.

The next obligation is cancellation in the signed word aggregate.
The absolute-value substitution alone grows exponentially like
√p((2a+2)√(2^a−1))^k after the spectral normalization, already
for its S-only portion. This failure is a limit of the bound,
not a lower bound or counterexample for the true aggregate. The
full trace also contains rank-one J terms and anchor corrections.
Possible next work must control their joint sum or bypass this
trace criterion; larger finite word catalogs alone are now redundant.
The subgroup and classical obligations remain as in the seventh pass.


## Tenth pass: all-word correlations, with an actual depth obstruction

The [word-aggregate proof](parallel10-word-aggregate-2026-09-05.md)
shows that geometric inverse convolution and local block-count signs
recover each word uniquely. Arithmetic self-duality and Schur's lemma
give a principal Gram matrix equal to identity up to operator error
(aΣd_w²+1)/√p. A raw-rank recurrence and the weight gap control
the error in replacing the original paths by principal traces. Further
disjoint anchor weights are controlled by a local inertia mismatch.

For the complete word family, the squared raw-rank budget grows as
Λ_a^m, with explicit Λ_a; Λ_2=(9+√65)/2. The resulting useful
range is Λ_a^m=o(√p), or 2^s(a+s)Λ_a^m=o(√p) with s fresh
adjacency conditions. The result applies to all coefficient vectors
in this range, rather than only chosen examples or a single-rank span.

This unrestricted isometry cannot be the final spectral argument:
when (2^a−1)^m>p−a−1 its evaluation map has a nonzero kernel.
The exact p=13 raw-path null vector certifies this on an actual
27-word family. That is an obstruction to the proposed universal
norm comparison, not to the Paley conjecture or its particular
cyclic trace. Next: retain the actual spectral coefficient vector,
or identify and use relations in a quotient of word space, while
handling the full J and exceptional-fiber terms. Do not simply
extrapolate the all-coefficient estimate to m/log p→∞.


## Eleventh pass: full powers and the remaining row-energy estimate

The [complete-energy proof](parallel11-full-energy-2026-09-05.md)
uses the actual coefficients rather than an unrestricted isometry.
The squared rank-sum budget is F_j=L_j²/(2^a−1)^j, with
L_j=[a(2^(a−1)(a+1)−1)^j−(2^a−1)^j]/(a−1) for a≥2.
It grows at rate Θ_a=[2^(a−1)(a+1)−1]²/(2^a−1), giving
Θ_2=25/3. Raw finite-point dimensions sharpen the boundary bound
to 2r_w p^(j/2) for each path. The proof then restores exceptional
columns, the missing constant direction, J, and anchor corrections
by exact identities and a Frobenius-norm perturbation estimate.

The resulting complete aggregate is at most (1+o(1))p when
j≤(1/2−ε)log p/log Θ_a, for fixed a. This is an advance beyond
the open-curve word family, but still does not imply a new clique
bound. The spectral criterion needs j/log p→∞.

For the exact normalized matrix T, define H_j=||T^j||_F² and
ρ_j=2^a||ET^j||_F²/H_j−1. The next sufficient target is
Σ_(l<j)log(1+(2^a−2)ρ_l/(2^a−1))≤o(j) at those longer depths,
uniformly in the required anchor cliques. No assumption about the
signs of ρ_l is available. The converse bound
H_j≤p+[(N−1)²/N]j²(|tr T^(2j)|/2+p), N=2^a−1,
shows this is a positive reformulation with essentially the same
strength as the missing trace condition. The arbitrary-two-set,
subgroup and prize obligations remain separate and unproved.


## Twelfth pass: an obstruction to a projection-only extension

The [projection countermodel](parallel12-bootstrap-obstruction-2026-09-05.md)
replaces one direction of the Paley projection while keeping its rank,
constant null direction, all chosen anchor columns, common-neighborhood
mask, and S^2=pI-J. The change has rank at most two and plants a
localized eigenvalue at zero or one. Its other entries no longer form
an exact character matrix, which is the condition the countermodel loses.

For the complete normalized matrix, H'_j <= (sqrt(H_j) +
2 sqrt(2) j N^(j/2))^2, while H'_j >= N^j. The prior short-depth
range satisfies j^2 N^j=o(p), so its (1+o(1))p bound survives the
modification. The required long-depth bound does not. A successful
extension must use information not captured by these preserved inputs,
such as the full arithmetic of the character entries away from anchors.

For the actual p=257 localization at anchors {0,1,62}, an integer
vector has normalized Rayleigh quotient -812/(19 sqrt(1799))<-1.
This rules out a zero-error edge interval for every prime. It does not
rule out a vanishing error as p tends to infinity. The next target
must respect those asymptotic quantifiers. No clique or subgroup bound
improves in this pass, and the exact official prize bridge is still open.


## Thirteenth pass: the principal rank rate and a separated error term

The [principal-budget proof](parallel13-principal-budget-2026-09-05.md)
keeps the moving-endpoint singularity separate from the lower-weight
error. Additivity gives c_y(E_w)=c_y(K_w)-c_y(F_w[1])=0. Every
error constituent is lisse there, so it cannot match a principal word
constituent. The resulting principal/error cross correlation is bounded
by a*d_u*e_v/p after normalization, rather than a pointwise norm estimate.

Let D_j be the sum of principal ranks and L_j the raw sum. Summing local
Jordan types over all words yields an exact linear recursion. A positive
weighted count with weights linear in k plus a multiple of a^(-k)
proves D_j is comparable to lambda_a^j, where lambda_a is the larger
root of lambda^2-h(a+1/a)lambda+N=0, h=2^(a-1), N=2^a-1.
It is strictly below the raw rate h(a+1)-1.

The generic-column budget is now 1+(a D_j^2/N^j+1)/sqrt(p)+
[2a D_j(L_j-D_j)+(L_j-D_j)^2]/(p N^j). All finite rows,
exceptional columns, the constant direction, J and anchor corrections
are restored as before. With Psi_a=lambda_a^2/N and Theta_a the raw
budget rate, the useful depth is the minimum of
(1/2-epsilon)log(p)/log(Psi_a) and
(1-epsilon)log(p)/log(Theta_a). At a=2, Psi=(19+5sqrt(13))/6,
and the leading logarithmic depth coefficient improves by about 16.5%.

The remaining dominant term is the principal/principal Frobenius
correlation budget. Neither the new rank identity nor the error
separation controls that term at j/log(p) tending to infinity. The
full spectral edge, clique improvement, arbitrary two-set and uniform
subgroup bounds, and the exact prize bridge remain unproved.

## Fourteenth pass: an elliptic kernel with four character blocks

The [elliptic derivation](parallel14-elliptic-model-2026-09-05.md)
normalizes three anchors to {0,1,r} and uses E:y²=x(x−1)(x−r).
Under rho(P)=x(2P), ordinary common-neighbor fibers have eight
points and the four exceptional fibers contain all sixteen rational
4-torsion points. The complete lifted matrix is exactly
M(P,Q)=f(P+Q)f(P−Q), f(x,y)=chi(y), f(O)=0. The full quotient
has block [[8S_C,8sqrt(2)1],[8sqrt(2)1^T,12]] plus −4I_3.
Its border coupling is of square-root order and is not discarded.

On G₀=G minus G[4], the operator
(M_G₀−J_G₀/sqrt(p))/(2sqrt(7p)) has exactly the old normalized
localized spectrum, with additional zeros. The four fractional-linear
symmetries r/x, (x−r)/(x−1), r(x−1)/(x−r), and identity give
four character blocks, including fixed-point corrections. The existing
p=257, r=62 outlier lies in character (1,−1,1,−1). Its reduced
8×8 block has a vector of squared norm 38 and quadratic form −406,
reproducing the prior exact normalized Rayleigh quotient below −1.

The function f descends to H=G/G[2]. Every Fourier coefficient is
bounded by sqrt(p), using a rank-one tame sheaf with four punctures
and unramified Lang-character twists. The Fourier matrix entry is
|H|^(-1) sum_{gamma²=alpha/beta} F_gamma F_(gamma beta), with four
terms or none. These entries are not diagonal. Their correlations
and the punctured centered restriction still require a sharp norm
estimate; scalar coefficient bounds alone have not supplied it.

The verifier checks 27 curves, 136704 kernel entries, 976 inverse
root triples, 5632 periodicity identities, 108 character projectors,
352 positive leading minors and 4768 exact cyclotomic Fourier entries.
No new asymptotic operator or clique bound is claimed. The pass-thirteen
logarithmic range, subgroup/classical/prize gaps and need for independent
review remain unchanged.


## Fifteenth pass: complete translated products and an order-level obstruction

The [higher-correlation proof](parallel15-elliptic-correlations-2026-09-05.md)
works on the quotient elliptic group H. For a multiset of shifts, let O
be its odd-multiplicity set and B its positive even-multiplicity set.
With g_O(h)=product_(t in O) f(h+t), the actual product equals

g_O(h) - sum_(b in B) g_O(-b) 1_(h=-b).

For s=|O|>0, every character-twisted sum of g_O is at most s*sqrt(p).
There are 4s disjoint tame punctures upstairs on E, with dim H_c^1=4s;
dividing the four-to-one quotient gives the constant s. The full
product is bounded by s*sqrt(p)+|B|. For s=0 it is exactly character
orthogonality minus the missing points. This also controls the stated
joint Fourier convolution and every four-point row correlation of M².

For every eligible p>=10000, the same H supports a modified even sign
function with zero set {0}, full exceptional border and four symmetries.
A rectangle A of size between ceil(sqrt(p)) and twice that size has a
sum/difference region of size less than 24ceil(sqrt(p)). Setting the
nonzero signs in that region to +1 produces a clique block. The exact
centered normalized Rayleigh quotient is

4(|A|-1-|A|/sqrt(p))/sqrt(7p) >= (4/sqrt(7))(1-2/sqrt(p)).

The modified translated-product bounds are at most 50s*sqrt(p)+|B|
at every order and twist, by an L1 perturbation estimate. Thus their
square-root orders do not force edge convergence. This does NOT
preserve the original sharp constants, original Kummer realization,
or ambient projection. It is not a Paley counterexample; the previous
projection obstruction and this sign-kernel obstruction retain different
constraints and cannot be combined into one purported example.

Exact checks cover 1047 parity patterns, 162 twisted supports, 6864
positive leading minors, 5168 M² entries and 100 cyclotomic convolution
identities. The p10009 and p65537 examples cover cyclic and balanced
quotients and preserve 75552 exceptional entries. Their rational
normalized Rayleigh lower bounds exceed one. The next step needs
sharp constants and/or additional arithmetic relations together with
the ambient projection constraints. No asymptotic Paley norm or clique
bound improved; all broader obligations remain open.

## Sixteenth pass: square fields retain both sharp inputs

The [square-field derivation](parallel16-square-field-barrier-2026-09-05.md)
uses actual Paley matrices on F_q, q=ell^2, and the classical clique
F_ell. Its centered indicator is an exact positive eigenvector, with
squared-mass fraction 1-1/ell on that subfield. For a fixed set of a
anchors in the subfield, the localized projection has Rayleigh quotient
1-(a+2)/(2*ell)+a/(2*ell^2). The centered normalized norm therefore
converges to 2^a/(2*sqrt(2^a-1)), greater than one for a>=2.

The SAME examples retain the full ambient Paley identity, all actual
sign entries, full elliptic model and exceptional border, four symmetries,
original Kummer realization, and sharp twisted-product bounds s*sqrt(q)+e.
The finite-field Lang and tame-cohomology inputs permit q=ell^2.
Hence merely combining the two constraints lost by earlier abstract
examples is insufficient. The characteristic also tends to infinity.

For three anchors and ell>=29, the exact transfer spectral radius is
greater than 3*sqrt(7)/4. Energy exceeds (63/16)^j and even trace exceeds
(63/16)^j-q. The needed long-depth estimates fail on this family.
The primary Kunisky conjecture requires prime order, so this is not a
counterexample to it. A quantitative prime-field input is still needed.

Five exact quadratic fields pass both ambient identities, 25 localizations,
and the subfield eigenvector checks; three also pass the complete elliptic
model and sharp all-character certificates for selected shift supports.
Singular certificates explicitly accommodate equality. The uniform proof
and finite checks remain separate, and independent review is outstanding.
The prior prime-field depth range, clique and subgroup bounds, and exact
official prize bridge are unchanged.

## Seventeenth pass: prime-field nondegeneracy and determinant decay

The [prime uncertainty derivation](parallel17-prime-uncertainty-2026-09-05.md)
uses the classical Chebotarev Fourier-minor theorem in Tao's primary
paper. For m<=r=(p-1)/2, both compressed Paley projections X,Y are
positive definite, with X+Y=D=I-J/p. In particular 0<X<I. This is
a true prime-field input; it excludes exact endpoint eigenvalues.

There is a positive integer k(C)=p^(m+1)*det(X)*det(Y). The rank-one
J term extracts (sqrt(p))^(m-1) from det(2pX). The integral norm in
Z[(1+sqrt(p))/2] removes the power of two, proving divisibility by
p^(m-1). Positivity gives k>=1. With Delta=1-m/p and
theta=4^m*k/(p^(m+1)*Delta^2), every eigenvalue of X lies between
g=Delta*(1-sqrt(1-theta))/2 and 1-g. In particular the gap is at
least k*4^(m-1)/(p^m*(p-m)), and hence at least the same bound
with k=1.

For A=D^(-1/2)*S_C*D^(-1/2), the exact identity is
tr(A^2)=m*(m-1)+2*||S_C*1||^2/(p-m)+(1^T*S_C*1)^2/(p-m)^2.
Thus theta<=exp(-m*(m-1)/p) and the certificate g is at most
Delta*exp(-m*(m-1)/p)/2. At m~p/2^a this certificate tends to
zero exponentially. This is NOT an upper bound on the actual gap.
It proves that this single determinant-product certificate cannot
establish the required constant edge gap for a>=2.

A separate prime cyclic interval-frequency projection has an explicitly
localized binomial vector with exponentially small Fourier leakage.
It retains real circulant structure, zero transformed diagonal and
the ambient square identity, but loses Paley off-diagonal signs for
large p. Its coordinate interval is not an anchor-neighborhood claim.
It illustrates the limits of qualitative Fourier uncertainty and is
not a counterexample satisfying the complete Paley hypotheses.

Exact checks pass 4132 Paley compressions, 41420 positive conjugate
leading minors, 4132 independent integer determinant identities and
second-moment bounds, and 518 shifted-gap minors. The rank cutoff
is checked with two singular larger cases. All prior bounds and
broader obligations remain unchanged; independent review is open.

## Eighteenth pass: Cartesian expansion and the exact trace gap

The [Möbius analysis](parallel18-mobius-trace-2026-09-05.md) uses
G(U,V)={u_s w u_t}. Point transporters over F_p and its quadratic
extension imply |G∩hH|≤max(2max(m,n),120) for every proper subgroup
H, by the classical subgroup classification. For m,n≥p^ε this meets
Lyamkin's Theorem 11, yielding a normalized nontrivial-representation
operator bound O_ε(p^−κ(ε)). The exponent is unspecified. The exact
energy n²E₊(U)+m²E₊(V)−m²n² is also proved.

For the Weil average T over G(A+2,−B),
tr(T)=[S(A,B)+√p|A∩B|]/(mn). Its (0,0) entry is 1/√p.
Thus the generic dimension-times-norm inequality cannot supply the
desired trace saving. The exact factorized norm and diagonal
Cauchy–Schwarz identities also retain the familiar threshold.

The genuine Weil conjugacy average of [[3,−1],[1,0]] equals
(I+χ(5)R)/(p+χ(5)); its trace is one while its norm is 2/(p+χ(5)).
Every higher-power trace is exactly ((p+χ(5))/2)^(1−j). This proves
a specific limitation of recovering the first trace from generic
operator or higher-power bounds. It is not a Cartesian-family example
and does not refute Paley. The needed centered trace estimate is
exactly the original cancellation problem, not a new proved reduction
that removes its difficulty. Independent review remains outstanding.

## Nineteenth pass: relation-free input sets and restricted moments

The [relation-free reduction](parallel19-relation-free-reduction-2026-09-05.md)
proves an exact extension criterion for B_h sets, p>h: the forbidden
new points are S union the sets j^(-1)(hS-(h-j)S), 1≤j≤h.
This partitions an arbitrary input into B_h blocks of size k with
remainder at most L_h(k)=(k-1)+Σ_j(k-1)^(2h-j).

Cancellation on all polynomial-sized B_h pairs, for any fixed h,
therefore implies the full classical two-set conjecture. The converse
is immediate. This equivalence is unconditional as a reduction; the
restricted cancellation statement is still unproved.

For the existing slice k=⌊p^(1/(r+1))⌋, a bound
M_(2r)(C)≤C_r p^(1+β_r)k^r is now sufficient when imposed only
on B_(h_r) sets along an unbounded sequence of fixed r with β_r→0
and h_r/r→0. This hypothesis is called SS-B. It has not been
established, and no converse from Paley to SS-B is claimed.
The normalized transfer includes the explicit remainder term
L_h(k)/|B|. The sufficient fixed-order threshold is
ε>max(1/(r+1)+β,(2h-1)/(r+1)). For h=r=2 it is vacuous in
the meaningful ε<1 range; a Sidon-only fourth moment cannot be
substituted for the unrestricted fourth-moment input without this cost.

The previous forced-row obstruction at size p^(1/r) does not refute
SS-B at the smaller slice. The simplest mixed Weil product instead
recovers the existing norm identity and adds no estimate. The new
partition identities and exact finite checks pass. Independent review,
the restricted uniform character estimate, the full classical and
subgroup targets, and the quantitative official prize bridge remain open.

## Twentieth pass: inversion and equal-size moment transfer

The [inversion reduction](parallel20-inversion-moments-2026-09-05.md)
strengthens the preceding sufficient criterion. For an arbitrary
k-set C and p>h, all but at most
Q_h(k)=k+(2h-2)binom(binom(k+h-1,h),2) inversion poles produce
a B_h set D={1/(c-z)}. Distinct h-term multisets give nonzero
rational functions with numerator degree at most2h-2, even when
summands repeat or cancel. This proves the pole count directly.

Weights w_(1/(c-z))=chi(c-z) obey
M_(2r)(D,w)=M_(2r)(C)-|F_C(z)|^(2r)+k^(2r)>=M_(2r)(C).
Completing D to a B_h set D union E of size2k expresses w1_D as
two averages of unsigned k-set indicators plus s/k times1_E,
where |s|<=k. This costs at most3^(2r) in the moment bound.
All sets whose moments are estimated still have size exactly k.

If p>Q_h(k) and p>f_h(2k-1), uniform restricted and unrestricted
moment bounds are equivalent up to3^(2r). The explicit sufficient
threshold K_h=max(2,4h(h-1),(h+1)2^(2h-1)+1) gives both
conditions for k>=K_h,p>=k^(2h). Therefore on the SS slice
k=floor(p^(1/(r+1))), h can be any integer2<=h<=floor((r+1)/2).
At2h=r+1, the bad-pole leading coefficient is (h-1)/(h!)²<=1/4;
the completion obstruction has smaller degree. This includes the
boundary without assuming an unspecified small constant.

The new SS-B* criterion uses this range and beta_j->0 along
unbounded fixed r_j. The old h_j/r_j->0 requirement is unnecessary
for this stronger transfer. No actual uniform SS-B* estimate is
proved. A Sidon sixth moment on the p^(1/4) slice would transfer;
the Sidon fourth moment on the p^(1/3) slice is outside the proved
range. Exact checks and the final artifact audit pass. All full
targets and independent mathematical review remain outstanding.

## Twenty-first pass: the one-sided squarefree target

For C of size n and T_(2r)=sum_(Q subset C, |Q|=2r) K(Q), the
pointwise elementary-symmetric recurrence proves
T_(2r)>=-(8r)^r p n^r/(2r)! and
M_(2r)<=2^r(2r)!T_(2r)+2(16r)^r p n^r.
Conversely |T_(2r)|<=[M_(2r)+(8r)^r p n^r]/(2r)!.
Thus one-sided T upper bounds replace the SS-B* moment hypotheses
with no change to the permitted relation order or size slice.
The upper bound itself remains open. No estimate uniform in growing
r is claimed; every r and its constants are fixed before p grows.

At r=3 the weighted construction on k=floor(p^(1/4)) has totalmassp,
column zero mass1, mean0 and Gram pI-J; odd correlations0,
C2=-1, |C4|<=2sqrt(p), |C6|<=5sqrt(p), M4<=3pk², and
M6/(pk³)>=k/2. Additive Sidon labels exist, but there is no
connection between those labels and the sign rows. The measure is
rationally weighted, not the uniform p-row character measure. This
rules out using just the retained relaxed constraints to bootstrap
M6. It does not rule out a character-sensitive argument.

Next sufficient input: prove a positive upper bound for the actual
squarefree aggregate on the SS-B* slice, or another sufficient
character estimate. Exact verification is complete for the stated
finite cases; mathematical review and all broader targets remain open.

## Twenty-second pass: three actual-kernel interfaces and all-order obstruction

Classical: for n>=6,n^4<=p, U_C(t)=sum_(|Q|=3)K(Q union {t}),
L_out=sum_(t outside C)U_C(t)^2. Exact full-translate identity and
zero-row corrections imply L_out<=20pT6+p²n³ and
20pT6<=L_out+3p²n³. The off-C quartic sums are squarefree.
Their disjoint-triple Gram entries contain pK(Q union R)-K(Q)K(R),
so ordinary termwise Weil returns the original degree-six problem.
No uniform L_out estimate is proved. Sidon transfer starts eventually;
its convenient n25 threshold is not certified by n6..10 fixtures.

Subgroup: let F4=E2-(3n²-3n), R6 count zero-sum six-tuples with
no opposite pair, K=r2(2), epsilon3=1_(3 in H). Then
E3=(15n³-45n²+40n)+(15n-60)F4+60n(K-1)-30n epsilon3+R6.
In particular intrinsic E2 implies E3=intrinsic6+R6. With
F4<=An²L and R6<=Dn³L, Konyagin gives
M<=[15+(15A+D)L]^(1/9)n^(8/9) in p<=n^4. Neither uniform
input is proved. Symmetric non-subgroups can have intrinsic E2 but
E3>=n^4/20 in the dyadic quartic regime; multiplicative closure matters.

Spectral: for p1mod4,C={chi(x)=chi(x-1)=1},m=(p-5)/4,
R_C=(J(eta,chi)^2+conj(J(eta,chi))^2-6p+36)/16.
The uniform direction has total QR Fourier fraction
q_C=(3p+5)/(8p)+R_C/(2m sqrtp)=3/8+O(p^-1/2), with
q_C in[1/4,1/2) atp>=73. This is one direction only.
Its residual variance is (1/(16m))sum_C(L(x)-6)^2-(R_C/m)^2,
L(x)=sum_y chi(y(y-1)(y-x)); this coupling and nonconstant
vectors are still uncontrolled. The all-vector gap/depth does not improve.

Weighted obstruction: for each fixedr>=3,k=floorp^(1/(r+1)),
k>=2(r+1), extreme massk², zero-fibre mass1 percolumn and bulk
W=p-k-k² with cube density1-((k²+1)/W)e2 in[1/2,3/2].
All odd products0,C2=-1,all evenCd=k² for d>=4; Gram pI-J.
All real-coefficient moments at2s<=2r-2 have constant
(3/2)(2s-1)!!+1 timesp||a||2^(2s), yet M2r/(pk^r)>=k/2.
B_h labels exist for all2h<=r+1. This is a weighted measure,
not uniform field translates, and therefore no Paley counterexample.

The missing positive estimates must use more actual arithmetic
structure. Separate-agent reviews and exact checks support the stated
limited results; human review and formal verification remain open.
