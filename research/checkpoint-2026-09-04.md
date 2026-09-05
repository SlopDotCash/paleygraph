# Checkpoint: main goal active and unproved

Primary target: square-root additive-character cancellation for thin dyadic
multiplicative subgroups, recovered from the existing prize repository.
See `subgroup-target.md` for the precise working statement and its limits.
The classical quadratic-character Paley formulation is a separate secondary
line; it must not be substituted for the prize-linked target.

Completed in this run:

- Retrieved current primary literature and prior repository files pinned
  to `5b00e50c3c51b3c944201a1749a1a8132e5ce167`.
- Reconstructed exact second/fourth quadratic-character moment identities.
- Proved obstructions to uniform Gaussian and constant-corrected moment bounds
  for arbitrary quadratic-character input sets, including a logarithmic lower
  allowance from common-neighborhood counting.
- Checked a quadratic-character eighth-moment counterexample in Lean.
- Recovered the exact centered additive-energy/coset moment identity for the
  actual subgroup target, with an explicit conditional logarithmic-depth bound.
- Found a structural obstruction to the exact real-Gaussian coefficient when
  the subgroup contains 2. Supplied a mathematical proof of `E₂(H)≥3n²+n` for
  `-1,2∈H` in characteristic greater than 3.
- Independently verified the finite case `p=6700417`, `H=⟨2⟩`, `n=64`,
  `E₂=12864`, in the concrete window `n⁴/4≤p≤n⁴`.
- Lean checked the finite subgroup arithmetic and normalized collision count;
  the bijection converting that count into additive energy is proved in prose.
- Verified all three Lean logs against their current source hashes. Their axiom
  censuses contain only standard Lean axioms; no holes or new axioms were used.

Further progress in the logarithmic-depth investigation:

- Built exact additive convolution on the multiplicative quotient, reducing
  the target field to 104695 stored values including zero. Checked orbit labels,
  mass, array overflow bounds, independent dense small-field convolutions,
  and agreement with the previous sparse energies.
- Computed `E_r` and the centered `Q_r` exactly through `r=12=ceil(ln104694)`.
  The last moment gives `M²≤1970`. The weaker proposed (SG) with `K=2` holds
  at all twelve tested orders; the exact Gaussian coefficient fails at
  every order from 2 through 12.
- Found an analytic finite counterexample to the literal amplitude constant
  `√2` in `M≤√(2n ln(p/n))`, using the same quartic-window subgroup. A rational
  cosine envelope proves `Re η₁>43`, while `ln(p/64)<12`. See
  `finite-spectral-obstruction.md` and `SharpConstantCounterexample.lean`.
- Distinguished that failed constant from the pinned repository's unspecified
  squared coefficient and its different envelope `M²≤2n ln p`. The exact
  upper bound shows the present example satisfies the latter envelope.
- Saved reproducible rational witness data and all twelve exact moments.
  The quotient computation is not a Lean proof; the analytic lower certificate
  uses the actual complex exponential sum and standard Lean axioms.
- The completed fourth Lean certificate passes with six printed axiom censuses,
  all containing only `propext`, `Classical.choice`, and `Quot.sound`.
  All four verification logs match their current proof-source hashes.

What was NOT achieved:

- No proof or disproof of either Paley formulation.
- No asymptotic counterexample to a bound with an unspecified constant.
- No independently audited theorem chain establishing equivalence to the
  official Proximity Prize statements.
- No prize submission, external messages, or changes to the prior repository.

Next concrete research step: estimate centered convolution growth in (SG)
uniformly across fields, with enough slack to allow the
arithmetic resonance `2∈H`, while keeping the constant uniform as the moment
order grows like `log((p-1)/n)`. The exact Gaussian coefficient is not a valid
universal substitute; neither is the literal amplitude constant `√2` with
natural `ln(p/n)` at all finite sizes in the chosen quartic window.
Any proposed restriction on exceptional primes must be
proved compatible with the sponsor's worst-case quantifiers.

Further progress from the geometric-cycle investigation:

- Strengthened the general resonance count to `E₂≥3n²+9n` for `-1,2∈H`
  in characteristic greater than 5, by counting all twelve permutations of
  `(u,u,2u,-4u)`. The current example attains this lower bound exactly.
- Proved in ordinary mathematics a uniform bound `(160rn)^r` for the average
  nonprincipal `2r`-th period moment on the full cycle modulo `a^(n/2)+1`.
  The proof uses independent base-a digits, a uniform bounded-difference
  estimate, and exact removal of the principal frequency. It applies to every
  integer `a≥2` and every moment order. It has not been formalized in Lean.
- Derived the exact signed prime-reduction error `D_r` and proved that a
  uniform upper estimate `D_r≤(Krn)^r` at logarithmic depth would imply (SG).
  This error bound remains unproved; the auxiliary theorem does not close it.
- Implemented exact carry counting modulo `a^N+1`, cross-checked against
  dense convolution at all 222 residues of seven small moduli through order 8.
- Independently reconstructed all twelve prime energies through order 24.
  For generator 2, additional prime relations first occur at order 8; one is
  `1-2^7-2^9+2^14-2^17+2^19-2^21+2^23=p`.
- The computed `D_r` is negative through depth 12 for generator 2. The same
  subgroup with generator 8 has positive `D₂>765`, disproving a proposed
  unqualified nonpositive-error shortcut. The lift is generator-dependent.

Updated next step: estimate the signed arithmetic error (D) from
`geometric-lift-and-alias.md` uniformly at logarithmic depth, or find a
different route to the prime-field estimate. The full-cycle concentration
bound cannot be transferred by conditioning without an exponential loss.

Further progress from the prize-reduction audit:

- Proved the converse direction for the signed lift error: (SG) implies
  `D_r≤(Krn)^r` directly from nonnegativity of the auxiliary centered moment.
  Together with the previous implication, these bounds are equivalent up to
  constants. The lift has not made the remaining theorem weaker.
- Recovered the original April 8 ABF paper from a local archive, with hash
  recorded; visually checked pages 5, 17, and 23. The July revision remains
  unavailable, so its additional lower bound has not been audited.
- Verified the exact MCA event, including the same-witness and
  no-joint-explanation conditions. No Paley reduction was found in the
  inspected primary paper.
- Proved an elementary circuit theorem: at every radius, MCA error is at
  most `min(1,binom(n,k+1)/q)`. At `δ≥(n-k-1)/n`, this is attained when
  `q>binom(binom(n,k+1),2)`, by choosing a degree-k direction and an offset
  outside finitely many collision hyperplanes. No novelty claim is made.
- Proved monomial pairs have at most 40 bad scalars for dimension-four RS
  on `μ_8`, whereas the global maximum is 56 in sufficiently large fields.
  The proof groups 24 complementary triples containing an opposite pair;
  their scalars lie in an eight-element set. This is an ordinary proof.
- Produced explicit attaining pairs `f=X^5+2X^6+3X^7`, `g=X^4` at primes
  2017 and 65537. The first lies in the quartic window for n=8.
  Exact modular arithmetic checks all 56 circuits and 64 monomial pairs
  in each field, independently verifies functionals by polynomial division,
  and saves all 112 cubic interpolation witnesses. Both monomial maxima
  are exactly 40; the unrestricted maximum is exactly 56.
- Identified the unsupported extremality link in the prior chain:
  `FarLineIncidenceEquivariance.lean` proves automorphism invariance, not
  reduction to monomial pairs. Existing prior files already record other
  monomial-domination failures. The new result does not rule out a monomial
  direction with arbitrary offset; its own direction is monomial.
- Corrected a prior note's claim that the cited LD-to-MCA theorem gives
  random RS MCA up to capacity; the verified April source gives Johnson
  through that implication because of the square-root loss in radius.

Next step: distinguish an actual theorem linking arbitrary MCA pairs to a
Paley estimate from special-family calculations. Any proposed bridge must
survive the circuit obstruction and preserve the sponsor's full parameters.
The prime-field moment estimate remains a valid separate research target,
with its unproved obligation stated without claiming the lift has reduced it.

Further progress from the multiplicative-symmetry investigation:

- Rechecked a primary exposition of BGK (Kurlberg, Theorem 1.1, 2007).
  It gives a fixed power saving for a fixed subgroup-size exponent, contrary
  to an older note's `n^(1-o(1))` description in the fixed quartic regime.
  This is a baseline theorem, not a new or optimal bound claimed here.
- Derived the exact coset correlations `Σ_j x_j x_(j+s)=p1_(s=0)-n`.
  Equivalently all nontrivial multiplicative Fourier coefficients have
  magnitude `sqrt(p)` and the trivial coefficient is -1.
- Constructed explicit synthetic vectors meeting those correlations,
  the sum -1, and the coordinate bound n, with a spike `n-(n+1)/m`.
  They are not actual periods. Four exact examples have prime p and prime
  quotient m≡3 mod4 in the quartic window. No unbounded family of such
  simultaneous primes, or actual period counterexample, is claimed.
- Derived the exact additional identity `x_j²=n+Σ_s k_s x_(j+s)`, where
  k_s counts members of 1+H in the corresponding coset. Inversion in H
  proves every kernel coefficient is even except the coefficient at 2H.
  Cubic and quartic moment identities follow; the latter recovers the
  earlier energy calculations rather than improving them.
- Checked the complete group-ring identity in five prime fields, counting
  all ordered pairs independently. In the original resonant example the
  nonzero kernel coefficients are 28 copies of 2, one 3, and one 4.
  Each synthetic vector fails the actual subgroup's nonlinear equation;
  the nonzero residual is verified in Q(sqrt(D)) without floating point.
- Added the necessary energy consequence of the target maximum:
  `E_2≤n^4/p+C²n²(1-n/p)ln(p/n)`, hence `O(n² log n)` in the quartic
  regime. This necessary estimate is itself unproved here and is not a
  substitute for the full requested objective.

Updated next step: investigate the arithmetic kernel of 1+H and the
nonlinear period relation for information beyond low moments. Linear
coset correlations alone do not exclude almost-maximal constructed spikes.
It is not proved that even the displayed nonlinear relations alone suffice
for the desired bound. The main goal remains active and unachieved.

Further progress from signed moments and polynomial certificates:

- Independently reconstructed every signed period moment through order 24.
  For `f_r=1_H^{*r}`, the odd zero counts use `Σ f_(r-1)f_r`; even counts
  use `Σ f_r²`. Quotient convolutions through r=12 suffice, with checked
  int64 bounds and arbitrary-size integer products. All 25 zero counts
  match the separate carry computation, including the previously unused
  odd orders.
- Built degree-twelve rational orthogonal polynomials and directly audited
  the full Gram matrix. The associated sum-of-squares bounds prove every
  period is in `(-26,43.81)` and at most one coset has period at least 42.
  The latter uses positive Bernstein coefficients over the entire interval
  `[42,43.81]`, not sampled point evaluations.
- Used Machin's identity, alternating arctangent/cosine series, and outward
  rational rounding to certify
  `43.802482797626304198≤η₁≤43.802482797626304199`.
  Combining this with the polynomial filters proves `M=η₁` and shows all
  maximizing frequencies are exactly H. All other cosets have absolute
  period below 42. This is an exact finite result, not a Paley proof.
- Proved the general conditional recurrence criterion (J): bounds
  `|α_j|≤A√(n(j+1))`, `β_j≤Bnj` through d imply
  `M≤[A+(2+m^(1/(2d)))√B]√(nd)`. At d=ceil ln m, absolute A,B would
  supply the target bound. Uniform control of these coefficients is not
  established or assumed; the criterion is not claimed easier than (SG).
- The test field passes the coefficient inequalities through degree 12
  with A=B=1. Its central fourth moment nevertheless exceeds the Gaussian
  coefficient: the exact decomposition includes
  `β₁(α₁-α₀)²`, while β₂ itself is below 2β₁. This accounts for how
  skewness and the Gaussian comparison can differ.
- Saved all rational filter coefficients, positivity margins, moment counts,
  recurrence coefficients, and source hashes in
  `results/period_polynomial_bounds.json`. No new Lean proof is claimed.

Updated next step: derive or refute uniform control of the recurrence
coefficients from the actual subgroup arithmetic. The finite maximum is
fully certified; further numerical precision on that instance is unnecessary.
The full Paley estimate and the exact prize reduction remain unresolved.

Further progress from testing the recurrence hypothesis:

- Found and independently certified a counterexample to the literal B=1
  coefficient bound inside the quartic window: p=67403009, n=128,
  H=⟨64701253⟩, β₂=76801470694776/277291762225>256. The kernel has
  κ_H=0 and squared sum 277. The energy is exactly 51840.
- Enumerated all zero quadruples. Besides the degenerate 3n²-3n
  ordered tuples, there is precisely one H-orbit of n unordered
  nondegenerate quadruples, each with 24 orderings. A representative is
  (1,1074964,1550267,64777777), whose entries sum to p. There are no
  zero triples. Pair sums and direct Gram–Schmidt independently verify
  the energy and coefficient calculation. No new Lean proof is claimed.
- Proved that (J) through degree two already entails
  E₂≤n⁴/p+(A⁴+8A²B+3B²)n². Thus the proposed route requires a
  uniform O(n²) energy estimate in the quartic window. A weaker
  degree-scale coefficient condition is equivalent to the desired
  maximum bound up to constants, rather than an established shortcut.
- Inspected Theorem 1.1 and Proposition 2.7 of Yip–Yoo 2608.02568v1
  and Theorem 1.1 of Kim–Yip–Yoo 2309.09124v4. Exact decompositions
  and full product-rectangle containment are not supplied by the
  collision counts; no application to the required bound was proved.

Further progress from a canonical cyclotomic lift:

- Proved the weighted prime-sum bound (N): for fixed dyadic n and k,
  Σ_{p≡1 mod n}(Z_k(p,n)-T_k(n))ln p≤(n^k-T_k(n))ln k, where T
  counts zero tuples of complex n-th roots and Z the finite-field
  tuples. The differences are nonnegative. The proof uses integer
  multiplication determinants modulo X^(n/2)+1 and their nullities
  over each split prime; no assertion of novelty or formalization.
- Derived an absolute O(n²/(t ln n)) bound on the number of primes
  in the quartic window with E₂>3n²-3n+tn². This is not a density
  theorem without a prime-count lower bound, and not worst-case control.
- Completely audited all normalized quadruples for n=4,8,16, factored
  every nonzero norm, and checked nullities by primitive-root evaluation
  and energies by pair sums. This classifies all primes with extra
  zero quadruples for those orders, not only primes in a scanned range.
- Identified the remaining limitation: the numerator n^(2r)ln(2r)
  grows too fast at r≈ln((p-1)/n), and the estimate does not remove
  the principal-frequency term. It does not prove (SG).

Updated next step: investigate whether the cyclotomic structure controls
an individual exceptional prime, or whether a sharper centered estimate
can avoid the high-moment growth in (N). Do not extend B=1, treat an
absolute count as a density statement, or claim the decomposition
theorems automatically apply. The full Paley and prize goals remain open.

Further progress from testing the individual-prime norm route:

- Sharpened (N) by AM–GM to (N') with right side
  (n^k-T_k)/2 * ln(k*n^k/(n^k-T_k)). The exact determinant suite
  verifies its denominator-cleared power inequality with integers.
- Proved its direct pointwise fourth-energy upper-bound expression is
  larger than n³ for dyadic n≥16 and p≤n⁴, hence weaker than the
  elementary bound. This is a limitation of the displayed consequence,
  not a general prohibition on cyclotomic methods.
- Proved Z_10-T_10≥n⁵(n-945)>0 for every dyadic n≥1024 and every
  prime p≡1 mod n, p≤n⁴. Extra relations are therefore universal at
  this order in that range, not rare exceptions. The principal-frequency
  contribution explains this without contradicting the centered target.

Further progress from the live official-profile audit:

- Located the Ethereum Foundation's August 20 better.codes announcement,
  followed its public repository link, and verified main/contract commit
  b34c0131cfa36b51111521541d7d3e35c8791082 (September 3 commit).
- Archived IRSProfile, both protected target files, README, and the
  dependency manifest at that exact revision. Archived the relevant
  CompPoly files at 641694629e4557520a1539b272ec338c9f3044c7 and
  ArkLib files at e65197892890b8fd9b0dc05b8980273cf1d595cc. All eleven
  files are checksummed in sources/official-prize-2026-09-04/manifest.json.
- The current profile is F_(p⁶), p=2130706433, with the size-2¹⁸
  domain in F_p, dimension 2¹⁷ per row, eight rows, and 128 repetitions.
  Neither p nor p⁶ lies in the present quartic window. This source audit
  is not a new Lean build or a proof of the full dependency closure.
- The protected lower statement bounds certifiedGammaError by 2^-128;
  the pinned source expands this as MCA plus the two-interleaved list
  term divided by q. It separately scores (1-δ)^128. The upper target
  concerns a worst-case winning-set quantity on an entire radius suffix.
  No source-level equality to our subgroup maximum was established.
- Proved η_q(b)=η_p(Tr b) for H⊂F_p, and hence maximum n for every
  extension degree d>1. All nonzero trace-zero frequencies maximize
  when n>1. The centered coset moments obey
  Q_q=p^(d-1)Q_p+(p^(d-1)-1)n^(2r-1). Removing only b=0 is insufficient;
  deleting the whole trace-zero space recovers the prime-field moment.
- Independently checked the official base prime, subgroup generator order,
  sextic irreducibility via Rabin, trace on the power basis by Frobenius
  and multiplication matrices, and the trace-pairing discriminant -3⁹.
  The explicit nonzero witness θ has trace zero, giving period n=262144.
  Exhausted F_25 with H=F_5*, including moment checks through order 16.
  Saved experiments/extension_trace_audit.py and its exact JSON output.

Updated next step: seek a centered prime-field estimate beyond the norm
average, while checking any proposed prize bridge against the actual field,
domain, MCA, and list parameters. The new official benchmark is context for
the bridge, not a replacement objective. Neither the trace obstruction nor
the high-moment barrier disproves the prime-field Paley conjecture.

Keep the main goal active. The current obstacle is unsolved mathematics, not
a missing user permission or an external access condition preventing progress.

Further progress from the subset-sum and list-decoding investigation:

- Proved the prime-field bound
  |N_s(c)-binom(n,s)/p| ≤ (p-1)/p * binom(U+t-1,t), with
  t=min(s,n-s) and integer U≥max(1,M(D)). Fourier inversion and a formal
  exponential coefficient majorant prove it because n<p. Complementation
  gives the smaller t. This is established subset-sum methodology, not a
  novelty claim; primary context is Li–Wan and Zhu–Wan.
- At radius 1-(k+1)/n, the list of X^(k+1) is exactly the number of
  size-(k+1) zero-sum subsets of the domain. The bad scalars for the pair
  X^(k+1),X^k are exactly the negatives of all such subset sums. The
  no-joint-explanation condition follows from the degree-k direction.
  This concerns one pair, not worst-case MCA over arbitrary pairs.
- If all sums occur for D⊂F_p⊂F_q, that pair has exactly p bad scalars,
  hence probability p/q. For the official q=p^6 this is p^-5<2^-128,
  although the list is enormous. An extension-field period estimate
  that ignores the base-field support would be invalid.
- Proved a root-lifting injection: for a|n,k, zero-sum subsets T of
  H^a of size k/a+1 give degree-<k polynomials
  X^(k+a)-product_{z∈T}(X^a-z) agreeing on k+a points. This also embeds
  in any interleaving with the same column agreement.
- At official p=2130706433, n=262144, k=131072, choose a=32. The image
  order is 8192, its exact fourth energy is 203464704, centered coset
  moment Q2=52370599862533, and M≤2691. The kernel has one entry 1,
  4029 entries 2, and 33 entries 4. An independent 8192² membership
  count gives 24837 normalized quadruples and the same energy.
- The exact binomial lower bound has 8155 bits, hence at least 2^8154
  codewords at radius 4095/8192. It exceeds q, so the official finite
  MCA-plus-list certifiedGammaError exceeds 1 for every larger radius.
  Thus its protected lower certificate cannot use these radii. This
  does not prove a winning-set lower bound or an upper-track claim.
- Proved an obstruction to extrapolating root lifting: an odd subset
  of complex dyadic Nth roots cannot sum to zero; its nonzero norm is
  at most t^(N/2). A prime-field zero sum requires p to divide the norm.
  For fixed N≥4 and t=N/2-1, no such subset remains when p>t^(N/2).
  The image order N stays fixed at fixed gap a/n=1/N, so this method
  cannot yield a fixed-gap asymptotic counterexample.
- Exact DP checks cover 3490 small-field size/residue combinations.
  Exhaustive enumeration verifies all 11440 subsets in the p=17 case;
  114 reconstructed polynomial witnesses check all small-field bad
  scalars. The p=97 root-lifting test verifies all 64 injected codewords
  for a=2 and the associated complement norms divisible by 97.
  All eleven pinned official source hashes and both helper hashes are
  recorded. experiments/subset_sum_list_certificate.py passes and saves
  results/subset_sum_list_certificate.json. These are not Lean proofs.

Updated next step: pursue the uniform centered prime-field bound (SG),
or a legitimate arbitrary-word/pair reduction with the required radius.
For a degree-(k+a) monomial offset, top-coefficient cancellation imposes
a elementary-symmetric constraints; the single linear subset-sum bound
only handles the first. Polynomial-phase character sums would appear in
a Fourier treatment of the further constraints, and no sufficient
uniform bound is established. Do not extend the particular-pair identity,
the finite list obstruction, or the fixed-order norm argument to broader
quantifiers. The full goal remains active and unachieved.

Further progress from the large-spectrum and Riesz-product audit:

- Returned to the uniform prime-field goal after the finite decoding
  bridge. Rechecked the subgroup target, nonlinear period identities,
  recurrence route, and their current limitations from workspace files.
- Proved a direct entropy estimate: if a dissociated set Λ of D
  frequencies has Fourier magnitude at least δn on A⊂F_p, |A|=n,
  then D I(δ)≤ln(p/n), with
  I(δ)=((1+δ)ln(1+δ)+(1-δ)ln(1-δ))/2. A normalized Riesz product,
  Jensen, and the endpoint chord of ln(1+ty) prove it. This is standard
  Chang-type methodology, not a novelty claim.
- Applying it to one maximizing multiplicative orbit cannot force
  M=o(n): D≤log_2 p by distinct subset sums, and even that maximal D
  leaves the displayed inequality compatible with δ=1/2 in the quartic
  window. This is a limitation of that argument, not an actual family
  of subgroup counterexamples or a prohibition on all spectrum methods.
- Defined R_σ(b;t)=product_{h∈H_+}(1+σt cos(2πbh/p)), with one member
  from each opposite pair, and its mean Z_σ(t) over nonzero b. The product
  is constant on multiplicative cosets. The same chord estimate gives
  an explicit pointwise bound in terms of ln(m Z_σ(t)).
- At t=sqrt(ln m/n)≤1/2, Z_±(t)≤exp(Cnt²) implies
  M≤(2C+8/3)sqrt(n ln m). The converse holds with C=A/2 if the
  target has constant A. This is an equivalent formulation, not a
  reduction to an easier unproved estimate. Small n uses the trivial bound.
- Proved the exact normalization
  Z_+(t)=(p*a_0(t)-(1+t)^(n/2))/(p-1), where a_0 counts signed zero
  subsets with no repeated opposite pair. The principal deletion must
  be retained; no independence of the subgroup is justified.
- Independently enumerated all signed zero triples and quadruples in
  H=<64701253>⊂F_67403009, n=128. There are no triples and exactly
  128 squarefree signed quadruples, matching the known H-orbit.
  Hence a_0(t)≥1+8t⁴ for t>0 and
  Z_+(1/4)≥((33/32)p-(5/4)^64)/(p-1)>201/200. This refutes only the
  stronger proposed Z≤1 normalization, not an unspecified-C bound.
- The standalone standard-library script experiments/riesz_tail_audit.py
  passes. It checks exact positive/negative group-ring products, coset
  invariance and principal deletion in three small fields, with exhaustive
  dissociation checks at dimensions 5, 4 and 8. It rechecks prime/order
  data and saves all 128 signed witnesses plus the strict rational
  certificate in results/riesz_tail_audit.json. No Lean proof is claimed.
- Read James R. Lee's primary arXiv:1508.07109v2 HTML for the established
  entropy/Riesz context. The candidate proof in this note is self-contained;
  no current-best energy or exponential-sum bound is asserted from the
  broader literature search.

Updated next step: obtain genuine arithmetic control of the centered
relation counts or of the nonlinear period system. Reusing an independence
bound, the direct Chang estimate, or the equivalent Riesz criterion would
be circular or quantitatively insufficient. The conjecture and the full
prize reduction remain unproved; the main goal stays active.

Further progress from the primary-source and amplification audit:

- Read Di Benedetto et al., arXiv:2003.06165v1, Theorem 3.1 and its §5
  proof. Its existing bound M≲n^(2689/2880)p^(1/72) applies to the entire
  working window n⁴/4≤p≤n⁴, dyadic n≥4. It yields
  M≤C_η n^(2849/2880+η) for every η>0. No novelty, explicit useful finite
  constant, or best-current status is asserted.
- Checked the actual proof's exponent elimination:
  p≳n^(191/40)Δ^72, giving saving 31/2880. Its inputs are E₃≪n⁴log n
  and E₂≪n^(49/20)(log n)^(1/5); it does not depend on the questioned
  general-function theorem below.
- Proved an algebraic bound for a precisely stated conditional extension
  with three selection orders r,s,ℓ. The saving supplied by that ledger is
  (2r+2s+2ℓ−e_r−e_s−e_(2ℓ)/2−4)/(8ℓrs). Optimal full-energy powers
  cannot lie below max(j,2j−4), so the formula is at most 1/16. A clamp
  reduces the global optimization to 32 triples, all checked exactly;
  equality at (1,4,1),(4,1,1). With r,s≥2 the maximum is 3/64.
  At the actual source orders (3,3,1), ideal energies give only 1/24.
  Arbitrary selections, zero-sum deletions, and their thresholds have not
  been proved; this is a bound on the formula, not all amplification.
- Read Shkredov, arXiv:1802.09066v2, Theorems 3 and 25 and the latter's
  proof. The displayed general-function hypotheses admit
  f=1_{0}−1/p and H=F_p*. This mean-zero invariant function satisfies
  f*f=f, ||f||₁=2(1−1/p), and T_j(f)=1−1/p. At k=2 the statement
  would force 1≤65536 C_*(log p)^4(1−1/p)^4/sqrt(p−1), impossible for
  large p. Equation (57) visibly omits the origin. Record this as a
  qualification needed in the inspected HTML's broad function version;
  it does not refute subgroup Theorem 3 or establish a repaired proof.
- Granting the separate subgroup-centered formula (C), its direct
  fixed-k saving is (k+1−2e₂)/(4·2^k). The best is 11/1280 at k=5 for
  e₂=49/20, or 1/64 at k=4,5 for hypothetical e₂=2. An exact difference
  formula proves those global maxima. At r≈log n the saving shrinks to
  O(log log n/log n); the formula does not give (SG).
- Read Kowalski–Untrau Theorem 3.8 and Lemma 3.9. Their Gaussian theorem
  concerns prime order d=o(log p/log log p), not dyadic quartic subgroups
  or maximum control. Retained the arXiv header/internal-date discrepancy
  without relabeling the version.
- Corrected README's clique wording: uniform O(log p) is false because
  Graham–Ringrose gives c log p log log log p infinitely often, as recorded
  by the inspected Magsino–Mixon–Parshall introduction. Polylogarithmic
  remains the relevant stronger candidate relative to subpolynomial.
- Added research/analytic-bounds-and-amplification.md and
  experiments/analytic_bound_ledger.py. The standard-library script passes
  exact exponent elimination, reduced/global optimization, four rational
  convolution origin tests, and hashes of all four archived HTML files.
  Results are in results/analytic_bound_ledger.json. No Lean proof added.

Updated next step: seek arithmetic information that improves centered
high moments or uses the nonlinear period identities beyond the full-
energy ledger. The primary target, classical Paley conjecture, and full
prize reduction remain open. Do not repeat this optimization as if it
were an untried route or use the unqualified general-function theorem.

Further progress from mixed period identities and shifted energy:

- Classified the preceding turn as progress: it changed the authoritative
  analytic baseline and ruled out a quantitative use of one amplification
  ledger. Re-read current source, frontier, coset, and recurrence notes.
- Derived all mixed identities x_j x_(j+t)=n*1_(t=0)+Σ_s C_ts x_(j+s),
  where C_ts counts z∈g^tH with 1+z∈g^sH. C is symmetric, obeys
  C_ts=C_(-t),(s-t), and has row sums n−1_(t=0). The symmetric matrix
  with blocks 0, sqrt(n)e₀ᵀ, sqrt(n)e₀,C has spectrum n plus every
  nonprincipal period once. The principal projection is explicit.
- The integer matrix L=C−n e₀ 1ᵀ has precisely the periods as eigenvalues.
  Simplicity of those eigenvalues follows from the cyclotomic minimal
  polynomial. Hence any vector of sum −1 satisfying the m mixed equations
  with base coordinate y₀ is exactly a shifted actual period vector.
  This establishes completeness of that arithmetic system; it is not a
  new bound or easier equivalent target.
- Proved direct counting identities ΣC=p−2 and
  ΣC(C−1)=(n−1)(n−2). Most importantly, with R=(H−1)\{0},
  X(H)=E×(R)−(2n²−5n+3)=ΣC(C−1)(C−2).
  The inverse bijection from nontrivial shifted-product equality
  (a−1)(d−1)=(b−1)(c−1) is
  x=(b−1)/(a−b), y=ax, z=cx. All denominator and distinctness conditions
  are proved; the two trivial pair matchings are exactly the exclusions.
- X=0 iff all C entries≤2, equivalently all distinct affine copies aH+b
  have pairwise intersection size≤2 (circularity). This gives a finite
  certificate using O(n²) product counts rather than all ambient copies.
- Row-zero parity gives D=Σκ²−(2n−3)≥0 and
  D≤Σ_s(κ_s)_3/3+2≤X/3+2. Therefore
  E₂(H)≤3n²−n+nX/3; when X=0 the exact value is 3n²−3n.
  Uniform X=O(n log n) would give necessary E₂=O(n² log n), but this
  stronger excess estimate is not established or shown necessary for
  the desired spectral bound. The known X≪n²log n is too weak here.
- Read primary Hoshi–Kanai multiplication matrices, Garcia–Lorenz–Todd
  moment/circularity framework, and Shkredov's shifted-energy theorem.
  Archived all three HTML sources under sources/mixed-periods-2026-09-04/
  with checksums. The constructions are in established frameworks; no
  novelty or best-current claim is made.
- Added experiments/mixed_period_collisions.py and
  research/mixed-periods-and-shifted-energy.md. Exact checks pass in nine
  fields: all mixed group-ring identities and independent cell enumeration
  in five small fields, traces of L through order eight in two of them,
  and forward/inverse collision bijections plus additive/product counts
  throughout. All arithmetic acceptance conditions are integers.
- Resonant p=6700417,n=64 has X=114; p=67403009,n=128 has X=720.
  All their cells of size at least three are recovered and saved.
  New quartic-window p=1073748737,n=256,g=1064280392 and
  p=17179869697,n=512,g=13395504394 have X=0. Thus circularity and
  exact E₂=195840 and 784896 are certified there. Primality is checked
  by trial division; generator order by dyadic modular powers. These
  are finite properties, not uniform cancellation or a conjecture proof.

Updated next step: control the nontrivial shifted-energy excess, its
concentration in the distinguished row, or signed higher powers of the
full arithmetic matrix. The global comparison (9) alone is too weak
with known inputs, and low-order control still does not give (SG).
Do not assume circularity throughout the window or extrapolate the new
finite examples. The full goal remains active and unachieved.

Further progress from the kernel discriminant and quadruple orbits:

- Classified the preceding turn as progress: it derived exact mixed-period
  and shifted-product collision identities and changed the proof frontier.
  Independently reviewed the new polynomial proof and its finite certificates.
- Constructed P_n(Y)=product_(1<=j<n/2)(Y-(1+zeta_n^j)^n) in Z[Y].
  Its roots are distinct and real, and A_n=disc(P_n)P_n(2^n)>0.
  Power sums give a binomial/Newton recurrence with checked integer
  divisions. For splitting primes, excess fourth energy occurs exactly
  when p divides A_n. This is a self-contained argument, not a conjecture
  proof or novelty claim.
- Proved 4v_p(A_n)=sum_r D_r over compatible prime-power residue rings.
  The exact sequence at n=16,p=17 is D=(196,12,0), so first-level energy
  does not always exhaust the valuation. Direct pair counts at every
  precision independently verify the kernel-label calculation.
- Factored A_n completely for n=4,8,16,32. There are respectively
  1,2,6,37 exceptional splitting primes. The first three classifications
  agree exactly with the prior independent normalized-quadruple norm
  enumeration. None of these orders has an exception in n^4/4<=p<=n^4.
  At n=32, every exceptional prime is <=194977 except 21523361, above
  the window [262144,1048576]. Thus E2=2976 throughout that entire window.
  This is complete finite classification, not a scan or spectral bound.
- The existing n=64 and128 resonances have v_p(A_n)=3 and6, with D_2=0.
  The two circular n=256,512 examples have valuation zero. No attempt
  was made to refine the already certified n=64 spectral maximum.
- Proved the general orbit identity D_r=4epsilon_r+12u_r+24v_r.
  Here epsilon records membership of 3 in the lifted subgroup, u counts
  nontrivial zero-sum multiset orbits of type (2,1,1), and v counts those
  with four distinct entries. Opposite-pair solutions are removed first.
  The dyadic group acts freely on the four-distinct multisets because
  any nontrivial stabilizer contains -1 and would force opposite pairs.
- Normalizing the repeated entry gives t_r=1+2epsilon_r+2u_r, where
  t_r counts pairs in H_r with sum -2. The boundary-root multiplicity
  equals epsilon_r+u_r. Over an unramified extension this works at every
  odd prime, giving the uniform integer factorization
  odd(F_n)=T_n U_n, odd(disc(P_n))=U_n^2 V_n^6,
  odd(A_n)=T_n U_n^3 V_n^6, with F_n=P_n(2^n),
  T_n=odd(3^n-1), and positive odd integers U_n,V_n.
  Their p-valuations are sum_r u_r and sum_r v_r respectively.
- The repeated-entry types contribute only O(n^2) energy, without any
  arithmetic hypothesis. In particular
  3n^2-3n+24n v_1 <= E2 <= 9n^2-15n+24n v_1.
  Hence the needed E2=O(n^2 log n) is equivalent to v_1=O(n log n).
  This orbit bound remains unproved. A corresponding bound on v_p(V_n)
  would suffice but is potentially stronger due to higher precisions.
  The elementary height estimate after removing powers still has order
  n^3/log n and does not improve the trivial orbit bound.
- Added research/kernel-discriminant.md and experiments/kernel_discriminant.py,
  plus research/quadruple-orbits-and-cube.md and experiments/quadruple_orbits.py.
  The latter checks direct additive counts, independent multiset orbit
  enumeration, polynomial reductions, and compatible lifts in 55 cases
  spanning 109 precision levels. Five cases use explicitly verified
  unramified quadratic extension rings, so non-splitting primes are tested.
  The power lift is independent of the earlier Taylor correction.
- The exact determinant A_64 has 40003 bits and satisfies the power
  factorization, but it is not fully factored. Its U_64,V_64 have
  1730 and5107 bits. No complete order-64 prime classification is claimed.
  At the known n=64 resonance, (u_1,v_1)=(1,0); at n=128, (0,1).
- A literature screen found related power factorizations for the different
  Wendt resultant in Helou (1997), whose indexed primary introduction
  was inspected. It is not an external dependency of the new proof;
  no full-paper audit, identification with A_n, or novelty claim is made.

Updated next step: bound the number v_1 of four-distinct zero-sum multiset
orbits in the growing-order quartic regime, or find arithmetic control of
higher centered moments directly. The orbit/power identities are exact
but do not provide that bound. Do not repeat fixed-order factorizations
or remove higher p-adic levels as though either established uniformity.
The primary subgroup target, classical Paley conjecture, and complete
reduction to the official prize remain open. The goal stays active.

Further progress from dyadic descent and mixed energy:

- Classified the previous turn as progress: it proved the uniform power
  factorization and isolated four-distinct zero-sum orbits. Read the live
  orbit, subgroup-target, and geometric-lift notes before continuing.
- Proved the exact recurrence for H=K union gK, |K|=k:
  E2(H)=2E2(K)+6B+8T, where B counts two entries from each coset and
  T counts three from K and one from gK. The mixed energy B is also
  sum_x r_(K+gK)(x)^2. Cauchy-Schwarz gives
  k^2<=B<=E2(K), T^2<=(E2(K)-k^2)B. Simply using B<=E2(K) does not
  close a quadratic-energy induction.
- Proved a conditional whole-tower criterion: B<=C k^2 Lambda at each
  step, C,Lambda>=1, implies E2(H_s)<=22 C s^2 Lambda at all levels.
  The exact constant inequality uses sqrt(22)<19/4. Lambda=max(1,log p)
  would give the necessary fourth-energy bound. This hypothesis includes
  smaller levels outside the fixed quartic window; it is not established
  or claimed equivalent to an endpoint bound alone. No independent upper
  estimate on T is required if B is controlled.
- Split new orbits into balanced (2+2 coset entries) and unbalanced
  (3+1). Exact identities are B=k^2+4k u22+8k v22 and
  T=k(epsilon_H-epsilon_K+3u31+6v31). Thus only balanced four-distinct
  orbits need separate control for the conditional induction. Their
  uniform bound remains open.
- With roots lambda_i of P_k, proved that the balanced repeated-entry
  count is #{lambda_i=-2^k}, and the balanced four-distinct count is
  #{i<j:lambda_i=-lambda_j}. Set C_k=|P_k(-2^k)| and
  J_k=|product_(i<j)(lambda_i+lambda_j)|. Both are nonzero integers,
  and disc(product(Y-lambda_i^2))=disc(P_k) J_k^2. Summing precision
  counts gives U_k odd(C_k)|U_(2k), V_k odd(J_k)|V_(2k); their quotient
  valuations count the unbalanced new orbits.
- Proved descent when ord_(2k)(p)=2 ord_k(p): {1,g} is a basis over
  the smaller unramified field, forcing B=k^2,T=0 and preserving all
  nontrivial orbit counts at every precision. The target has p=1 mod N,
  hence ord_s(p)=1 at every level. Membership g outside K is not field
  linear independence; this descent cannot be used along the target tower.
- Proved no nontrivial zero quadruples when p=-1 mod n: Frobenius
  inverts all roots, giving an even root polynomial and opposite pairs.
  Combining this with dyadic order formulas shows that new odd prime
  divisors of U_n/U_(n/2),V_n/V_(n/2) satisfy p=1 or n/2-1 mod n.
  Fixed-p valuations stabilize at m=2^v2(p-1) if p=1 mod4, and at
  m=2^(v2(p+1)+1) otherwise.
- Combined that theorem with the complete lower-order factorizations:
  for every dyadic n>=32, V_n=7*17^5*47*79*Z_n, all prime factors of
  Z_n congruent to +/-1 mod32. The fixed valuations are exact. No
  complete factorization of V_64 is asserted. This support restriction
  does not exclude primes already satisfying p=1 mod n.
- Proved the norm cutoff p^ord_n(p)<=2^(n/2) for p dividing V_n.
  A four-distinct sum without opposite pairs has exactly four +/-1
  coefficients in the dyadic cyclotomic basis, so its mean squared
  conjugate magnitude is four. AM-GM and residue degree give the bound.
  It removes four-distinct orbits in the quartic windows through n=32,
  but its threshold is above n^4 for every dyadic n>=64. It permits
  starting the mixed-energy induction at a level comparable to log p;
  the remaining larger levels still need (MB).
- Added research/dyadic-descent-and-mixed-energy.md and
  experiments/dyadic_energy_descent.py. Exact pair counts, independent
  orbit enumeration, and root labels pass in 44 field/precision cases,
  including ten field-doubling cases. Quadratic-over-quadratic tests
  verify the minimal base-field degree and the generator orders before
  using a binomial extension. Factor identities are checked through
  doubled order 64, complete prime support through order 32, and selected
  stabilized valuations at orders 32 and64.
- The two known quartic resonances separate the mixed terms: p=6700417,
  n=64 has B=1024,T=96 and one new unbalanced repeated-entry orbit;
  p=67403009,n=128 has B=4608,T=0 and one new balanced four-distinct
  orbit. Thus the literal claims B=k^2 and T=0 both fail in the target
  window. Larger same-prime and extension examples used for the descent
  audit are explicitly outside that window where applicable.
- Inspected and archived only the publisher abstract of Do Duc, Leung,
  Schmidt (2020), DOI 10.5802/alco.86. Its reported condition
  p>(sqrt14)^(n/ord_n(p)) is outside the growing dyadic quartic target.
  No full-paper audit, theorem dependency, novelty claim, or Lean proof
  is asserted. Archive and checksum: sources/dyadic-descent-2026-09-04/.

Updated next step: seek an actual upper bound on mixed energy B, or on
opposite-root collisions in P_k, along the split prime-field tower between
orders comparable to log p and the quartic endpoint. The recursion and
integer factorization are exact but do not bound these collisions. Do not
assume field-degree growth where all roots already lie in F_p, and do not
repeat prime-support or fixed-order checks as though they gave uniformity.
Even an adequate fourth-energy bound leaves the high centered moments,
the classical two-set Paley conjecture, and the full prize reduction open.
The goal remains active and unachieved.

Further progress from positive products and logarithmic-depth moments:

- Classified the preceding turn as progress: it established exact dyadic
  recurrences, conditional fourth-energy induction, and prime-support
  restrictions. Re-read the live dyadic, target, recurrence-coefficient,
  Riesz, and analytic-bound notes before extending the high-moment route.
- For the two real child periods X=eta_k(b), Y=eta_k(bg), proved
  (X+Y)^2<=max(X^2,Y^2)+3(XY)_+. With norms over nonzero frequencies,
  A_s(q)=||eta_s||_(2q)^2 therefore satisfies
  A_(2k)(q)<=2^(1/q) A_k(q)+3||(XY)_+||_q.
  The positive part suffices; negative products may cancel the parent.
- Proved that, for a fixed q>=4 and C>=1, the whole-tower hypothesis
  ||(XY)_+||_q<=Cqk implies A_s(q)<=4Cqs. The coefficient check is
  2^(1/q)<5/4 and 4*(5/4)+3=8. Taking the smallest even q>=log m,
  with minimum four, gives Q_q(H_N)<=m(4CqN)^q and
  M_N^2<=8e C N log m when m=(p-1)/N>=e^2. This reaches (SG) at the
  needed logarithmic depth, conditionally. The uniform hypothesis remains
  open, and it is not claimed necessary from an endpoint bound alone.
- Initial steps k<=q satisfy the hypothesis by |XY|<=k^2<=qk. The
  remaining condition begins at child orders larger than log p in the
  target regime. Also proved the unrolled weighted-sum version of the
  recurrence, allowing nonuniform bounds across those larger levels.
- For even q, defined Z_(q,q) to count q entries from each child coset
  with total zero. Exactly ||XY||_q^q=(p Z_(q,q)-k^(2q))/(p-1).
  The corresponding arithmetic condition (CM) is sufficient for the
  positive-part hypothesis. Odd q instead gives a signed moment and
  must not be treated as an absolute moment.
- The stronger absolute-product condition can overcount cancellations:
  exact moments at p=17,|K|=8 prove XY=-4 at every nonzero frequency
  (mean -4, second moment 16, zero variance). Its positive part is zero.
  A second constant-negative example has p=17,|K|=4 and product -1.
- Rewrote the arithmetic condition as a centered correlation:
  Z_(q,q)-k^(2q)/p=sum_x f_q(x) f_q(gx), where f_q is the q-term
  additive-sum count from K minus its mean. Cauchy-Schwarz gives only
  ||XY||_q<=A_k(q), which leaves the recurrence factor 3+2^(1/q)>2.
  It therefore does not close the required induction. The exact paired
  difference identity retains the squared imbalance discarded there.
- Over the complex roots of unity the balanced count is T_q(k)^2,
  by quadratic field separation. In the finite field the principal
  frequency alone gives Z_(6,6)>=k^12/p>=k^8/16 for p<=(2k)^4.
  Since T_6(k)=15k^3-45k^2+40k<=15k^3, every dyadic k>=64 forces
  Z_(6,6)>T_6(k)^2 at every eligible prime in that range. This includes
  all working quartic windows with parent order N>=128; no prime-count
  existence hypothesis or actual enumeration of the large count is used.
- At k=64, the complex balanced count is 14065500160000 and the
  principal lower bound at the largest allowed p is 17592186044416,
  forcing at least 3526685884416 extra ordered balanced relations.
  The known quartic witness p=67403009,N=128 is separately checked by
  primality, exact parameter inequalities, and its stronger integer
  lower bound. These facts refute raw complex-count replacement, not
  the centered mixed condition or the subgroup conjecture.
- Added research/positive-product-moments.md and
  experiments/mixed_high_moments.py. Exact quotient convolution and a
  fixed-point-free involution on K-cosets compute the centered product
  moments. Every int64 update has a checked overflow bound; dot products
  and principal subtraction use Python integers. Dense convolution
  independently checks every computed level in F_17 and F_97.
- Computed 80 signed product moments across nine levels. For the known
  p=6700417,N=64 example, rational exponential bounds certify q=12 as
  the smallest even integer >=log(104694). The three lower child orders
  2,4,8 satisfy (CM), C=1, trivially. Exact calculations verify (CM), C=1,
  at child orders 16 and32. Thus the full critical-depth hypothesis is
  certified in this single field. The archived parent q-th moment also
  satisfies the resulting coarse conclusion. No uniformity or new
  spectral-maximum refinement is asserted.
- New exact counts retain the mixed relations: at p=6700417, k=32,
  Z_(3,3)=26112 while the complex count is zero; Z_(4,4)=11677184
  versus complex 8856576. At k=16, equality still holds at q=4, while
  Z_(6,6)=2562534400 exceeds complex 2556313600. The critical-depth
  calculations retain these correlations without an independence model.
- Updated README, target, frontier, prior dyadic note, and source audit.
  The new proof is self-contained; the literature screen supplied no
  adopted external theorem. No novelty or Lean-formalization claim.

Updated next step: bound the positive product moments, or the stronger
centered balanced correlation (CM), at one logarithmic depth along the
split prime-field tower above its trivial initial levels. The off-diagonal
correlation needs new arithmetic control; Cauchy-Schwarz and complex-field
separation do not supply it. Do not ignore the principal-frequency term,
use odd signed moments as absolute ones, or repeat this finite tower
certificate as evidence of uniform closure. The entire subgroup target,
classical Paley conjecture, and full prize reduction remain open. The
persistent goal stays active and unachieved.


Further progress from signed quotient operators:

- Constructed an m-by-m symmetric integer matrix A(w) for every nonnegative
  integer H-invariant weight w with w(0)=0, using the nontrivial quadratic
  character of the dyadic group H. Its spectrum is exactly the nonprincipal
  Fourier values, once per H-coset, without a dense principal projection.
  This follows by restricting convolution to functions f(hx)=psi(h)f(x).
- For w=1_H, its spectrum is the subgroup periods. For n=2k and
  w=1_K*1_(gK), it is the child products eta_K(b)eta_K(gb). Thus the
  positive-product condition becomes a positive-part matrix trace bound.
  At every integer depth r, tr(A(w)^r)=(p*(w^r)(0)-d^r)/n, with
  convolution powers understood and d=sum(w).
- If B(w) is the unsigned quotient, proved
  sum(B_ij^2-A_ij^2)=d^2-2*sum(w^2). Its row sum total is md-d/n.
  Integer entry weights give total absolute-entry loss at most
  (d^2-2*sum(w^2))/2. A constant-vector Rayleigh quotient then proves
  rho(|A(w)|)>=d-[d/n+d^2/2-sum(w^2)]/m.
- Throughout the quartic window this yields rho(|S_H|)>=n-2/n for the
  subgroup matrix, and rho(|R_H|)>=k^2-k/4 for the product matrix.
  Both are almost the trivial degree. Hence taking entrywise absolute
  values cannot supply the intended cancellation through this comparison.
  This is a limit of that comparison, not a disproof of the signed target.
- Changing representatives is exactly A'=DAD for a diagonal sign matrix D.
  It changes no eigenvalue, absolute entry, or product around a closed walk.
  Random representatives therefore do not create independent edge signs.
- For weighted closed quotient walks of length r, the count with identity
  multiplicative endpoint ratio is m*c_r; the count for every other ratio
  u in H is exactly (d^r-c_r)/n. Here c_r is the weighted zero-sum count.
  Taking the quadratic sign gives the centered trace formula again.
  This exact balance does not bound the remaining identity-count excess.
- Added research/signed-quotient-operators.md,
  experiments/signed_quotient_operators.py, and its result JSON. Checked
  14 matrices at seven field/order pairs, 94 exact trace identities,
  46 independently computed balanced counts, and 32 endpoint-ratio
  identities. Matrix dimensions reach 1087. All acceptance tests use
  integers; sparse multiplication runs only under d^r<=2^63-1, and
  accumulation/subtraction uses Python integers. The n=32 product case
  uses depth six to stay within the checked bound. No new maximum
  computation at the completed n=64 example was performed.
- Updated README, frontier, target, positive-product note, and source audit.
  No external theorem from the literature screen was adopted. These are
  ordinary proofs and exact computations, with no novelty or Lean claim.

Updated next step: investigate arithmetic cancellation among signed closed
walks at logarithmic depth, or the positive spectrum of the product matrix.
Do not replace the signed matrix by its absolute value or treat diagonal
sign switching as independent randomness. The uniform centered estimate,
classical Paley conjecture, and complete prize reduction are still open.
The persistent goal remains active and unachieved.

Further progress from the July prize audit and list-to-winning-set proof:

- Classified the preceding turn as progress: it constructed the signed
  quotient matrices and proved the near-trivial absolute-matrix bounds.
  Direct closed-walk rearrangement in this turn supplied no additional
  cancellation theorem. The full spectral obligations remain open.
- Successfully recovered the July 6, 2026 ABF paper at its public PDF URL,
  after the previous HTTP failures. It has 55 pages and SHA-256
  ccb3eea9966485d9dd312a0eb46b3219d9a92cbe4fd191a88ad9976115ed892a.
  The title and PDF metadata agree on the date. Rendered and visually
  inspected pages 1,21,31,32,33,48,49,50. The April archive remains intact.
- Checked the new Lemma 4.16, Definition 6.11, the list-attack Lemma 6.12
  and its L<q hypothesis, Remark 6.14 separating CA from MCA, and Appendix
  C's common-coefficient list construction. The display in Lemma C.5 has
  a local sign inconsistency; the independent construction uses u=f-Q,
  so f-u=Q. This corrects the formula without refuting the list-size lemma.
- Archived the exact pinned ArkLib Definitions.lean, SoundnessBounds.lean,
  and SimplifiedIOR.lean. Verified the fixed-encoder relaxed relation,
  winning set, supremum over violating instances, and the upper target's
  required whole suffix below minimum distance. No converse equality to
  minimal game error or fresh Lean dependency audit is claimed.
- Proved for any injective linear encoder: a single-word list of L
  messages at radius delta0<d_min gives winning density >=L/(q+L-1),
  throughout [delta0,d_min). The witness has second word zero and claims
  (0,1). A suitable linear functional has image at least qL/(q+L-1), by
  the exact pair-collision average and Cauchy-Schwarz. The same instance
  stays invalid because any second message with claimed inner product 1
  has a nonzero codeword of weight at least d_min.
- If binom(L,2)<q, a union bound gives an injective functional on the list
  and the stronger winning lower bound L/q. Unlike the cited maximum-
  two-word-list result, the general single-word lemma permits L>=q.
  This is ordinary mathematics in an established list/projection framework,
  with no novelty claim.
- Proved a scalar list lower bound ceil(binom(n,r)/p^(r-k)) for any
  n-point domain in F_p and coefficient field containing F_p: partition
  r-subsets by the top r-k nonleading coefficients of their vanishing
  polynomials, and subtract each polynomial from their shared high part.
  Each selected codeword agrees at exactly r points. The full list at
  that center is precisely the selected coefficient fiber. Embedding in
  one interleaved row preserves distance. No subgroup cancellation is used.
- At the pinned p=2130706433, q=p^6, n=262144, scalar k=131072, eight rows,
  the exact r=139503 list lower bound is 85677801616821870413970774.
  Its unordered-pair count is smaller than q and its ratio to q exceeds
  2^-128. Thus winningSetDensity >2^-128 for every
  122641/262144<=delta<131073/262144. The lower endpoint is about
  0.4678382874. The same existential center and functional prove the
  entire suffix; no general monotonicity theorem is assumed.
- At the next agreement count r=139504, the pigeonhole guarantee is only
  35350350170772326, below the needed 274980728111395088. Exact successive
  ratios prove optimality only within this particular bound, not a safe
  lower radius or a sharp prize threshold. An independent exploratory
  binary search agreed with the recorded exact recurrence result.
- Exact integer inequalities 2^11649*139503^12800>=262144^12800 and
  2^11648*139503^12800<262144^12800 give the least centibit score for
  this radius, 11649. This is an upper-track score inequality, not a
  lower-track security guarantee. No production word or functional was
  explicitly enumerated, no Lean certificate was built, and nothing was
  submitted or sent externally.
- The prior scalar list of at least 2^8154 at radius 4095/8192 now gives
  actual winning density >1-2^(-7968) throughout that suffix, using q<2^186.
  This new conclusion comes from the single-word lemma, not from Gamma>1.
  A simple July antipodal family also gives a concrete center
  X^133120-X^131072 at radius 63/128 and binom(63,32) codewords; it has
  common subset sum 1, not the zero sum ruled out in an earlier scope.
- Added research/list-to-winning-set.md, experiments/list_to_winning_set.py,
  and results/list_to_winning_set.json. Four small-field tests exhaust
  86086 projection vectors and all codewords, check the exact collision
  average, list fibers, invalidity via minimum distance, and winning sets.
  Two cases have lists larger than the field and attain every challenge.
  Production arithmetic uses only integers and checks the selected and
  next binomials independently. Source/input hashes are recorded.
- The literature screen also found Goyal-Guruswami-Sun-Wootters,
  arXiv:2607.08516v1. Its Theorem 5.6 concerns random evaluation points;
  it does not supply the fixed smooth-domain guarantee. The HTML is
  archived, but its full proof chain is not adopted or audited here.

Updated next step: independently audit and formalize the single-word
list-to-winning-set lemma against the pinned fixed-encoder interface, then
its coefficient-pigeonhole specialization if pursuing the benchmark proof.
Keep the original Paley objective separate: no uniform signed-walk or
centered-moment cancellation has been proved, and the benchmark result does
not establish equivalence to either Paley formulation. The main goal stays
active and unachieved. There is no external blocker and no live process left
from these completed computations.

Final readback for this checkpoint verified all seven result-source hashes,
all five successful entries in the new source manifest, the root July-PDF
manifest entry, and the recorded exact list, projection, radius, and score
inequalities. The finite-case projection counts sum to 86086. Checked 60
local document links and corrected two display-math closing delimiters.
The fixed-encoder definitions and upper target were reread directly. The
pair (w,0) has a joint MCA explanation (u,0) for every folded agreement
witness; the attack instead violates the added inner-product claim. The
proof note now explicitly records why its winning-density lower bound is
not a lower bound for the grand MCA challenge's error.

Further progress on the alternate classical localization route:

- Classified the previous turn as progress: it proved and checked the
  fixed-encoder winning-set lower bound. No live process carried over.
  Rechecked the current frontier, preserved the goal, and returned to the
  classical Paley side rather than redefining the goal around the benchmark.
- Screened the localized-spectrum route in Kunisky arXiv:2303.16475v1.
  Read the necklace definitions, constant-label result, Appendix B.1, and
  Section 7's distinction between weak spectral limits and extreme values.
  This is a separate route; no implication for the dyadic additive-period
  target or sponsor equivalence has been established.
- Proved an exact identity for two singleton anchors: if one occurs once,
  the necklace equals -t_k, where t_j=tr(D_0 S)^j/(p-1) and t_0=1.
  Translation and scaling make every nonzero second anchor equivalent;
  summing that anchor then vanishes. This is valid for all lengths.
- Proved the two-occurrence identity, with r=k-ell, s=min(ell,r), d=|ell-r|:
  N_(k,ell)=p t_ell t_r-p^s t_d-t_k+2(-1)^k sum_(j=1)^(s-1) p^j.
  The proof uses full kernels K_j=S(D_0 S)^(j-1), with zero coordinates
  retained, and the elementary two-point correlation p*1_(x=y)-1.
- The nonzero-coordinate matrix C=D_* S_* is normal. Its trivial and
  quadratic directions have eigenvalue -1; on their orthogonal complement,
  CC^T=pI. This gives the exact cross trace
  tr(K_ell K_r)=p^s(p-1)t_d+2(-1)^k(p-p^s).
  Independently established the diagonal and boundary entries and
  ||K_j||_F^2=(p-3)p^j+2p. No generic eigenvalue assumption is used.
- Retrieved and read Lu-Zheng-Zheng arXiv:1305.3405v3 Lemma 2.1 (2.3)
  and Remark 2.4. With tensor exponents (1,1), R=1, and the exact
  identity t_j=p^(-j/2) sum_(a!=0) chi(a)|Kl_j(a)|^2, it gives
  |t_j|<=(j-1)p^((j-1)/2). This is an existing deep theorem used with
  its checked hypotheses, not an independent proof of the sheaf result.
- Therefore |N_(k,ell)|<=k^2 p^(k/2), uniformly for k>=2 and primes
  p=1 mod 4. Normalization by p^(k/2+1) gives <=k^2/p, even when
  k=o(sqrt(p)). This covers binary singleton patterns with a minority
  anchor occurring at most twice, hence all such patterns through length
  five. It does not cover general words involving the doubleton label.
- At three occurrences, proved the exact average over the cubic correlation
  C_3(x,y,z), weighted by K_ell1(x,y)K_ell2(y,z)K_ell3(z,x).
  Repeated indices have elementary values; distinct indices give a signed
  elliptic-curve trace. Termwise Weil plus Frobenius gives only the unsaved
  O(p^(k/2+1)) scale after anchor averaging. No estimate closing this gap
  is asserted. The length-six patterns 000111 and 010101 remain examples
  outside the proved family, with exact finite data recorded.
- Added research/localized-necklace-identities.md,
  experiments/localized_necklace_identities.py, and
  results/localized_necklace_identities.json. Seven prime fields, depths
  ten/twelve: 516 one-occurrence traces, 2778 two-occurrence traces,
  420 cross traces, 72 literal coordinate-tuple sums, and six independently
  contracted cubic-averaging identities passed. A separate multiplicative
  convolution reconstructs every monochromatic trace. All calculations use
  Python integers in NumPy object arrays; no floating-point threshold.
- Archived both primary HTML sources, with root manifest hashes. Noted the
  inconsistent trivial-character entry of Kunisky's displayed Jacobi table
  under its zero-at-zero convention; the matrix derivation keeps both -1
  directions and does not rely on the erroneous entry. No broader theorem
  refutation or novelty claim is made. No PDF or new Lean proof was used.

Updated next action on this route: seek cancellation in the weighted cubic
sum for three marked occurrences, or a simultaneous bound for the mixed
necklaces that can survive at growing depth. Keep the full two-set Paley
quantifiers and dyadic uniform target intact. Neither has been proved;
even full fixed-degree necklace convergence would leave the edge estimate
unproved. The goal remains active, with no external blocker. All matrix
verification runs completed with exit code zero; no process remains live.
The final run explicitly includes square and nonsquare second anchors in
every field. Final readback verified all four result-source hashes, both
new source-manifest entries, 49 local document links, the recorded check
counts and cubic contractions, Python syntax, and display-math delimiters.

Further progress from affine averaging and planar/projective reductions:

- Classified the previous turn as progress: it proved the one-/two-occurrence
  necklace formulas and exposed the weighted cubic sum. Rechecked that proof
  and the live goal, which remains active. No live process carried over.
- Proved the exact affine matrix average sum_b C_b^j=t_j(pI-J), where
  C_b=D_b S and t_j=tr(C_0^j)/(p-1). Translation/scaling give one diagonal
  and one off-diagonal value; C_b kills the constant vector. It follows
  that N(0^r1^s)=p t_r t_s-t_(r+s) for arbitrary positive block lengths.
  The existing analytic input gives |N|<=(r+s)^2 p^((r+s)/2), even at
  logarithmic length. This does not control other mixed words at that depth.
- Proved planar character-partition duality
  Z_p(G)=p^(v-e/2-1) Z_p(G*) directly by the quadratic Gauss expansion,
  vertex orthogonality, and the p-to-one dual-potential parametrization of
  divergence-free flows. For the two-anchor frame of a length-k binary
  necklace, this gives N(w)=Z_p(G_w*)/(p-1)-t_k.
- Proved the projective invariance needed for even-degree graphs with an
  even edge count. Assigning exactly an independent set I to infinity
  contributes the affine Z_p(G-I). The sum conditioned at any one vertex
  being infinity is independent of the chosen vertex. All uses retain these
  boundary terms; projective and affine sums are not interchangeable.
- For 010101, the planar frame is a cube and its dual is an octahedron.
  Normalize three distinct projective points to (0,infinity,1), and retain
  the case where the first nonedge pair coincides. The remaining twisted
  elliptic second moment W=sum_t chi(t)E(t)^2 is exactly t_4 by the
  K_2 kernel identity. Subtracting six singleton and three paired infinity
  contributions gives N(010101)=p(p-5)t_4+p^2(p-2)-t_6.
- For 001011, the dual H is K_6 with path edges da,ae,eb removed. Its
  explicit edge list is recorded. Complete it with a vertex g adjacent to
  its odd-degree vertices a,c,e,f. Condition first at the now-universal c,
  obtaining a five-rim wheel; condition at g, obtaining H and two deletions.
  Each deletion is minus the K_4 partition, by summing a degree-two vertex
  whose neighbors are adjacent. Hence N(001011)=p(t_5+2t_3)-t_6.
- Balanced binary length-six words have three symmetry classes, represented
  by 000111,001011,010101. The new identities and the previous minority-at-
  most-two formulas cover all 64 singleton words. A uniform ordinary proof
  gives |N(w)|<=36p^(7/2), hence o(p^4), for primes p=1 mod 4. This
  resolves the previously named length-six examples without proving the
  general three-occurrence weighted bound or the full degree-two conjecture.
- Added research/planar-necklace-reductions.md,
  experiments/planar_necklace_reductions.py, and its result JSON. Updated
  the preceding note/data to point to the subsequent formulas instead of
  describing those individual examples as unresolved. Historical checkpoint
  entries are retained. No new deep analytic dependency, Lean certificate,
  novelty claim, prize submission, or message to another person was made.

Next research obligation on this route: extend the graph reductions to
general longer mixed patterns and doubleton labels, or prove a uniform
signed bound for the resulting graph partitions. Fixed-length convergence
still does not control the spectral edge. The original arbitrary-small-set
Paley quantifiers, the uniform dyadic target, and the prize bridge remain
unproved; none has been replaced by this narrower intermediate theorem.

Verification completed: 462 two-block traces, all 64 binary length-six
words in seven fields (448 checks), and 18 complete affine matrix averages.
All 62 nonconstant binary length-six words have checked oriented sphere
embeddings. Literal affine/projective graph sums in F_5 and F_13 verified
four dual pairs, the two completion conditions, both deletion terms, and
the octahedron's signed elliptic moment and infinity contributions.
The earlier necklace verifier also passed after its labels were updated.
Final readback checked both result-source hash maps (four entries each),
the exact reported identities/counts, 55 local links, Python syntax, and
display-math delimiters. Every verification session returned exit code zero.
The live goal is active and unachieved; there is no external blocker and no
remaining process from this checkpoint.


User-requested parallel pass and status assessment:

- The user asked whether the work was closer, then explicitly authorized
  parallel agents. Reported that additional cases were proved but there
  was no complete proof route or justified distance-to-completion estimate.
  Spawned three bounded agents: subgroup moments, necklace extension, and
  classical arbitrary-set moments. The primary agent worked on the prize
  connection and independently reviewed all three proof notes.
- Necklace inversion swaps the singleton {1} and doubleton {0,1}, fixing
  {0}, on nonzero coordinates. The full-trace error is at most
  (n_B+n_C)p^(k/2), by a rank-one telescoping argument. This transfers all
  known binary estimates to the other two-label alphabets. All 189
  length-six words using at most two labels satisfy 42p^(7/2). The other
  540 words are not covered. The exact wheel identity N(ABC)=t_4 also
  completes the full degree-two family through length three.
- The subgroup symmetry projection gives 2T^2 <= B(E+B-2k^2), with a
  centered version at every convolution depth. The conditional energy
  induction constant improves from 22 to 7+4 sqrt(3). Balanced counts at
  the required scale remain unbounded; there is no exponent improvement.
- Biased row patterns strengthen the classical moment obstruction to
  M_(2r)(B) >= (D_r+o(1)) n^(2r) log(p/n) for some |B|=n whenever
  p/n tends to infinity, where D_r > 1/log 2. This is a necessary
  allowance for proposed bounds, not a proof or refutation of LM.
- On a base-field RS domain, base-field-valued pairs have exactly the same
  MCA-bad challenge set after extension of the coefficient field. Non-base
  challenges admit separate F-linear projections of the folded explaining
  codeword. Thus every such pair in the pinned sextic profile has MCA error
  <= p^-5 < 2^-128. More generally, alphabet subspaces U,V restrict nonzero
  bad challenges to those with U intersect gamma V nonzero, with at most
  (|U|-1)(|V|-1)/(p-1) such challenges. This does not control arbitrary
  extension-valued pairs or the protocol winning-set quantity.
- All four lanes have separate parallel-2026-09-04 notes, scripts, and JSON
  results. The readable assessment is research/parallel-pass-summary-2026-09-04.md.
  README and frontier point to it. No central file was edited by an agent.

Verification: the necklace lane checked 35 operator conjugations, 3516 each
of restricted inversion/rank-one/transfer/reflection statements, 945
length-six values, 735 two-block values, 30 ABC permutations, 156 literal
coordinate sums, and two literal wheel sums. The subgroup lane checked 52
fourth energies (six nonminimal), 608 higher-depth statements, and 76 parity
identities. The classical lane checked 52074 subsets, 33 double counts,
99 moment inequalities, and 1484 hypergeometric bounds. The prize lane
checked 2187 exhaustive base-pair bad-set comparisons, 256 further
comparisons using enumerated degree-one codewords, and 288 subspace-pair
thresholds. Finite checks supplement the ordinary proofs; no new Lean
certificate or external submission was made.

The root's final artifact audit verified all 13 recorded input hashes,
parsed all four scripts, and checked 18 local note links. The audit file
also hashes all five new notes, all four scripts, and their four results.
All agents and verification processes finished. The live goal read at the
start of this user-requested pass reported paused after the interruption;
no goal status change or scheduled background task was created. The full
Paley objective and the general prize reduction remain unachieved.


Second authorized parallel pass:

- Classified the preceding user-requested pass as progress. Re-read the
  authoritative notes and live goal; the goal is active and unachieved.
  All three prior agents were terminal. Reused them for a second bounded
  pass on subgroups, necklaces and classical moments; root worked on the
  quantitative spectral transfer. All four lanes are now terminal.
- Two growing three-label families are proved using checked Katz input:
  |N(A^rBC)|<=2p^((r+3)/2)+1 and, k=r+3,
  |N(A^rBAC)|<=2p^((k+1)/2)+p^(k/2)+p-1, for every r>=1,
  p=1 mod 4. The exceptional quadratic mode of V is exactly p^2-2p-2;
  H=V+pI-p vv^T replaces it by -2. The t=1 stalk correction is retained.
  General separating gaps at least three and weighted sums remain open.
- Root proved an exact two-projection trace identity for
  W=(sqrt(p)S-J)(2^a E-I), with E the common neighborhood of an a-clique.
  Its Chebyshev form retains the residual term
  ((p-r-m)+(-1)^k(r-m))/(2^a-1)^(k/2). A subexponential normalized
  aggregate bound at even k much larger than log p would give explicit
  clique bounds for growing a, without a separate localization-size or
  degree estimate. That aggregate is unbounded here; individual-word
  triangle bounds still lose an exponential factor. Degree a=1 is an
  unconditional orthogonal-involution calibration at square-root scale.
- The subgroup lane isolated a primitive-root polynomial R_n with
  B=2k sum(c_i^2). Its squarefree defect d gives
  k^2+4kd<=B<=k^2+2kd(d+1); no triple root gives B<=2k^2. The condition
  is an equivalent pointwise collision restriction, not a uniform theorem.
  A new quartic witness p=17189277697,n=512 has B=67584>65536,
  E_2(K)=195840,T=0,E_2(H)=797184, and two double fibers.
- The classical lane proved an exact signed cross-ratio representation
  with the random-set mean removed. A proposed complete-L2 sufficient
  estimate is false on intervals in an unbounded prime family, even in
  a range where the ordinary Weil bound establishes the fourth-moment
  estimate. Thus this refutes that L2 route, not LM or Paley.
- Archived the primary Katz PDF and extracted text, pinned its hash in
  the source manifest, and visually checked printed pages 22,52,53.
  No new Lean certificate, external submission or scheduling occurred.

Exact verification: necklace checks covered 192 growing-family cases,
216 gap identities, 288 label permutations, 10 literal sums, and 18 matrix
positivity certificates with 480 leading principal minors. The subgroup
pass checked 6144 quartic primes in fixed deterministic prefixes, 46
historical exception profiles, and the prior two quartic witnesses; this
is finite evidence, not a density theorem. The classical pass checked
645264 fractional-linear identities and 8440 fourth-moment/centering
cases. Root checked 344 exact trace recurrences in 43 projections across
five fields, 28 anchor expansions and 122 clique Rayleigh numerators.

The readable assessment is research/parallel2-pass-summary-2026-09-04.md.
The full two-set Paley conjecture, uniform dyadic estimate and general
prize bridge remain open. The next substantive obligation is signed
aggregate control, or another route that yields an actual uniform upper
bound. No external blocker has been identified and the goal remains active.

Final second-pass readback verified 13 recorded input/source hashes, all four Python scripts, 14 local note links, and the Katz manifest entry. The separate result audit hashes all five notes, four scripts, four results, and both archived source formats. Read-only PDF rendering intermediates were removed after inspection. The live goal was rechecked as active; no completion or blocked status was set.


Third authorized parallel pass:

- Continued the three existing agents on necklaces, subgroup mixed energy,
  and classical fourth moments. Root derived the two-anchor conic model
  and audited the analytic theorem inputs. All lanes and independent
  reviews are now terminal; no external blocker is present.
- Proved the arbitrary-gap family for A={0}, B={1}, C={0,1}:
  abs(N(B A^(j-1) C A^(m-1))) <= (min(j,m)+1)p^((j+m+1)/2)+1,
  every positive j,m and every prime p=1 mod 4. Exact hypergeometric
  identification gives the rank-j kernel including t=1. The tensor of
  the actual middle extensions has stalk rank j-1 at 1; rank mismatch
  removes H_c^2 for j!=2 and gives H_c^1 dimension j+1. The equal-rank
  case uses the prior exact correction. No rank<p restriction is needed.
  Label permutations retain normalized vanishing for k=o(sqrt(p)).
- Archived and independently inspected Katz G2 Section 2 (printed pp.3-5)
  and Katz GKM Euler-Poincare/trace/weight input (printed pp.32-33,38-41).
  PDFs and text extractions are pinned in sources/manifest.json. Deep
  source theorems are imported, not reproved or formalized.
- Root's conic map x(t)=((t+t^-1)/2)^2 has four-to-one ordinary fibers
  over R and two-to-one fibers over each anchor. Its lifted matrix has
  kernel f(t/u)f(tu), f(t)=chi(t^2-1), and an exact quotient consisting
  of [[4S_R,4*1],[4*1^T,2]] plus eigenvalue -2. Mellin entries are two
  Jacobi products and have off-diagonal couplings. The S3 action supplies
  block dimensions but no uniform norm improvement. The retained border
  coupling has size 4sqrt(|R|), relevant to extreme eigenvalues.
- The subgroup lane proved disc(R_n)=I_n^2 2^(d-1)d^d and
  B<=k^2+4k v_p(I_n), with I an explicit order index. The pointwise
  valuation bound needed for near-quadratic energy remains unproved.
  Quartic examples demonstrate failure of collision-free reasoning from
  unramifiedness, nonzero derivative, or Galois-stable fibers. Root fixed
  the d=1 endpoint of a height inequality to <=; final hashes refreshed.
- The classical lane proved Ehat(t)=tau chi(t)e(-t/2)Kl(t^2/16) for
  nonzero t. Its unrestricted Fourier-L1 sufficient bound fails, even
  after any bounded number of cross-ratio coordinates are removed:
  intervals have Omega_h(n^2 log n) configurations at each -h^2.
  Independent review confirmed the construction and clarified that no
  polynomial upper bound on the chosen p is supplied. FL1 restricted to
  n>=p^epsilon for a fixed epsilon is not refuted. LM and Paley remain open.

Exact verification: 320 gap identities/bounds, 192 full averages, 46
operator certificates with 864 positive principal minors, 48 literal
cyclotomic hypergeometric fibers, and 10 literal necklaces; five integer
index/discriminant certificates and three modular index certificates;
106 classical Fourier identities, 900 symmetries, 54 directional/affine
checks and 20 integer configuration records with 12 independent counts;
93312 conic entries in 14 fields, 2800 cyclotomic DFT entries, 84 quotient
and centering identities, and 84 automorphisms. Independent conic review
prompted the explicit centering assertion and reproduced the stored cases.

The readable assessment is research/parallel3-pass-summary-2026-09-04.md.
results/parallel3_pass_audit_2026_09_04.json verifies 14 recorded input
hashes, four script parses, 25 local links, two new source manifest entries,
and the 279-of-729 length-six coverage count for the specific proved
families. It hashes five notes, four scripts, four results and four source
files. The coverage count is not a measure of progress toward full Paley.

The full two-set target, thin dyadic target, and general prize reduction
remain unachieved. The signed weighted aggregate is still unproved, so
the new individual estimates do not yield the desired clique bound.
Next substantive work should retain this aggregate or tackle word
patterns with repeated exceptional labels, beginning with the uncovered
(2,2,1) multiplicities. The rank-mismatch method is available, but any
new tensor identification and invariant contribution must be checked.
No new Lean proof, prize submission, external message, or automation was
made. No completion or blocked goal status is justified.


Fourth authorized parallel pass:

- Continued the three existing agents on necklaces, subgroup arithmetic,
  and classical moments. Root proved pairwise kernel estimates and a
  signed aggregate across ranks. All four verification lanes and the
  independent reviews are terminal; no external blocker is present.
- The subgroup lane proves a lower limiting proportion at least
  1-log(2)/36=0.9807459116... of quartic splitting primes with exactly
  intrinsic fourth energy at every dyadic level. This is an ordinary
  asymptotic theorem, using the existing cyclotomic norm budget and the
  checked Thorner-Zaman PNT for dyadic moduli. It is not a finite-scan
  density estimate or a percentage of Paley proved. A boundary product
  excludes repeated-entry extra quadruples on a negligible prime set;
  only outside it does the 24N orbit quantum apply. The known order-64
  resonance has 12N excess and lies in that boundary support. A separate
  density-one all-level theorem has tolerance eta with an explicit
  bad-prime count; eta=N^(-1/2) gives vanishing relative energy error.
  Worst-case primes and high centered moments remain uncontrolled.
- The necklace lane proves the growing family A^r B^h C, with bound
  (s+1)p^((k+1)/2)+abs(t_s-(-1)^s), s=min(r,h), k=r+h+1.
  Its raw rank-two hypergeometric identification retains u=1 and all
  nontrivial Mellin characters; the trivial character is handled exactly.
  The full trace has a retained zero-coordinate correction. The two
  remaining length-five (2,2,1) bracelet classes have a full quartic
  reduction and explicit K_3,3 minors. Thus this pass covers 30 of 90
  such words, not all 90 or their full signed aggregate.
- Root proves unequal-rank bounds for every Mellin character and the
  corrected equal-rank formula q_r=k_r^2+p^(r-1)delta_1-p^(r-1),
  with Mellin norm at most (2r-2)p^(r-1/2). Arithmetic self-duality,
  the missing identity tensor in the singular stalk, and exclusion of
  Kummer self-twists were independently audited against Katz's PDF.
  The resulting quadratic aggregate across ranks R has normalized
  error ((3R^2-R-2)/(2sqrt(p))+2/p)||a||^2 for every complex a and
  every Mellin character. The next-anchor Mellin triangle inequality
  still loses sqrt(p); the full necklace aggregate remains unproved.
- The classical lane proves 0<=R<=3pn^2 in the exact fourth-moment
  decomposition. Burgess yields a directional saving for progressions
  and growing perturbations when p^(13/30+beta)<=n<=sqrt(p).
  Alsetri-Shao supplies a proper rank-two analogue. The bounds are
  consequences of established structured-set inputs and do not supply
  the arbitrary-set moment hypothesis. The relevant primary HTML files
  were independently checked, archived, and hashed.

Exact verification: 16 subgroup prime certificates across 68 dyadic
levels; 192 necklace block checks, 200 Mellin factorizations, 144 raw
hypergeometric fibers, 32 norm certificates with 608 positive leading
minors, 89 three-gap identities, and 9 literal necklaces; 42208 kernel
entries, 83 norm certificates with 1452 positive leading minors, 162
aggregate cases, and stalk ranks 1 through 8; 3869 classical exhaustive
cases, 20 progression cases including 8 perturbations, 99 Burgess
parameters, and 4 rational exponent choices. Small-field tests do not
establish the analytic asymptotics or replace the ordinary proofs.

The readable assessment is research/parallel4-pass-summary-2026-09-04.md.
results/parallel4_pass_audit_2026_09_04.json records input hash checks,
syntax checks, local links, new primary-source manifest entries, word
coverage, and final artifact hashes. No new Lean certificate, prize
submission, external message, or automation was made. The primary
goal remains active and unachieved; completion or blocked status is
not justified. Next work should address centered higher relations,
the two-variable necklace kernel, or the translated-anchor aggregate.


Fifth authorized parallel pass:

- Classified the fourth pass as progress and rechecked the live active
  goal and current notes. Continued the three existing agents on
  subgroup higher moments, the remaining necklace patterns, and
  arbitrary-set classical moments. Root handled translated-anchor
  quadratic aggregates. All verification runs and independent reviews
  completed; no external blocker is present.
- Proved a centered sixth-moment class at every dyadic level. Its
  explicit bad-prime proportion is at most
  (8 log6)/(21w)*(1+o(1)). With w=log N and the prior fourth-energy
  class, this gives E2<=4N^2 and E3<=(15+log N)N^3 on a density
  1-O(1/log N) of quartic splitting primes. Principal terms are
  retained exactly in the centered formulas.
- The fixed three-selection amplification, with the actual source
  hypotheses checked, gives M<=C N^(23/24)(log N)^(7/72) on that
  class. Root and an independent agent checked the support deletions,
  thresholds with no log loss, phase weights and exponent calculation.
  The original Petridis-Shparlinski trilinear theorem's size ordering
  is handled by permuting sets and weights. It has no additional
  quadrilinear size restriction. The source HTML is archived and hashed.
  The weaker proposed fixed linear saving was not used as an advance
  over the existing uniform sublinear baseline.
- The resulting centered eighth energy is bounded by
  N^(59/12)(log N)^(43/36) up to a constant; higher interpolated
  moments still encode exponent23/24, not the target1/2. The next
  proposed averaged eighth budget sum(E4-N8/p)log p << N7 polylog N
  remains unproved. A fixed p262657,N32 witness has exact E2 but
  additional six-term relations, excluding exact-fourth=>exact-sixth.
- Proved the full three-gap necklace family:
  abs(N(B A^(a-1) B A^(b-1) C A^(c-1)))
  <=(3a+1)min(b,c)p^((k+1)/2)+5p^(k/2)+2, k=a+b+c.
  The Mobius tensor has scalar quadratic inertia at infinity and
  actual-stalk Betti dimension(3a+1)b. The master identity keeps
  full-field kappa_j(0)=(-1)^(j-1), U=U_G+2(-1)^k, both other
  boundary fibers and the quartic constant. Root independently
  audited these signs and the five remainder bounds.
  All243 length-five words now satisfy abs(N)<=15p^3. At length6
  coverage is639/729, leaving90 balanced222 words. The full signed
  spectral aggregate remains open; these are finite family counts.
- Root removed the preceding translated-anchor Mellin-L1 loss by
  direct Kummer twisting. Pair rank bounds extend to any product of
  distinct nonzero anchors and every Mellin character. This gives
  O((m+1)R^3/sqrt(p)) for all complex quadratic combinations in the
  R-kernel span, and expected L2 mass in all adjacency-sign cells with
  a<=(.5-epsilon)log2p and R=O(logp), after half-valued anchor terms
  are removed. An independent review passed. This does not control
  every vector supported on the cell or approximate an extremal
  eigenvector; the dimension gap is explicit.
- The classical lane proves a one-cardinality sufficient moment
  reduction for a sparse sequence of fixed ranks, without claiming
  those upper bounds. It constructs polynomial-size sets at every
  sufficiently large prime with exact minimal additive energies
  through h<r but divergent M2r/(p n^r). This excludes the specified
  lower-energy Gaussian bridge, not LM with its logarithmic term.

Verification: subgroup12 fixed cases/47 moment levels, bounded six/eight
norm certificates and rational selection/exponent audits; necklace1580
inner checks,472 master/full-boundary checks,1215 length-five cases,
9 literal necklaces,95 Mobius maps and22 Parseval checks; root8848
kernel entries,340 pointwise checks,636 pair sums,294 all-Mellin
semidefinite certificates (3726 positive and90 zero pivots),648 signed
quadratic and1512 cell cases (978 nonzero boundary corrections);
classical12384 subset checks,92410 rectangle checks,8 structured
witnesses,38 energy and24 moment checks,99 exponent audits.
Finite arithmetic supports the proofs without establishing the
analytic asymptotics or imported cohomology/trilinear theorems.

The readable assessment is research/parallel5-pass-summary-2026-09-04.md.
The integration audit results/parallel5_pass_audit_2026_09_04.json
records final input/source hashes, scripts, links, and finite coverage.
The original goal stays active and unachieved. No new Lean proof,
prize submission, external message or automation was made. Continue
with centered eighth moments, balanced222 necklaces/full aggregate,
or the single-size upper estimates. No complete/blocked status is
justified by the current evidence.

### Sixth parallel pass completed

The original proof goal remains active and unachieved. Three authorized
agents completed the subgroup, necklace and classical lanes; root
completed the seeded-kernel lane and reviewed all three agent proofs.
An independent reviewer accepted both the subgroup and root proofs.

- The subgroup bound is now
  M_N<=2^(1/6)(15+log N)^(1/18)N^(17/18) on the existing sixth-energy
  class, of relative size 1-O(1/log N) among quartic splitting primes.
  It needs no additional fourth-energy intersection. The new centered
  mixed-order inequality is proved for every r,s by subgroup invariance,
  Jensen and elementary bilinear orthogonality. Origin atoms and the
  negative nonzero-coordinate constant are retained exactly.
- Higher-moment feedback gives centered E4<<N^(44/9) times
  (log N)^(10/9). All fixed-order substitutions into the coarse gate
  are proved unable to improve its saving beyond 1/18 with the present
  energy inputs. A genuinely N^4-scale centered E4 bound would give
  exponent 15/16. The averaged centered eighth budget remains open.
- The necklace lane proves ||W_C+S||<=2p by a pole-inclusive projective
  compression of an augmented single-anchor operator. Its three extra
  eigenmodes are -(p-1),0,p+1. Normalizing the adjacent C vertices gives
  N(AABBCC)=f0^T W_C S W_C f1-p t4+2 and a bound 7p^(7/2).
  All twelve words in its orbit agree exactly. Length-six coverage is
  now 651/729; the remaining 78 lie in four balanced classes. The next
  two exact contractions still receive only the order-p^4 bound here.
- Root proved every old kernel has inversion symmetry and identified
  an exact missed largest eigenvector on the p=13 common neighborhood.
  Multiple seeds capture that vector and have proved correlation,
  Gram and cell estimates. Seeds must avoid cell anchors; collisions
  can give order-p rank-one correlations. For L<=p^sigma and R=O(log p),
  relative cell error tends to zero up to
  a<=(1/2-sigma-epsilon)log2p, with sigma+epsilon<1/2. All seeds at
  rank one form an invertible matrix, so an isometry on a proper cell
  cannot hold on the full coefficient space. No full spectral bound
  is inferred from this enlargement.
- Classical exact fourth-moment variance is asymptotic to 24pn^4.
  On n=floor(p^(1/3)), the coefficient-3 Gaussian bound holds for
  1-(6+o(1))p^(-1/3) of all sets. A signed-pairing bound also holds
  for typical sets. A fixed p=1009 example refutes only coefficient 3;
  the note records why variance alone cannot control exceptional sets.
  No worst-case classical improvement was obtained.

Terminal exact checks: subgroup 8 fixed cases in 5 fields, 32 moments,
216 mixed-order certificates covering 1440 frequency/order pairs,
64 centering coefficient identities, 4096 exponent pairs; three sharp
equalities use algebraic certificates. Necklace 3005 compression
entries, 15 extra eigenvectors, 95 zero-sum basis vectors, 105 pole
entries, 20 anchor covariance pairs, 10 norm certificates, 60 block
words, 15 contraction identities, 6 literal paths, 3 full graph
partition checks and 5 original necklace sums. Seeded kernels:
18 inversion checks, 912 pointwise bounds, 2800 full-seed Gram entries,
2064 pair sums, 768 all-Mellin PSD certificates, 480 quadratic
combinations and 1440 adjacency cells. Classical: two exact symbolic
variance identities and 37 mean/variance classes covering 24678 subsets.
These checks supplement proofs; they do not establish asymptotics
or reprove imported analytic and cohomology theorems.

Read research/parallel6-pass-summary-2026-09-04.md and the four lane
notes for the current state. results/parallel6_pass_audit_2026_09_04.json
pins final input/source hashes, local links, counts and reviews.
No Lean proof, prize submission, external message or automation was
made. Continue with the centered eighth budget, remaining balanced
classes/full aggregate, or uniform single-size moments; do not count
additional kernel-span estimates alone as a full-operator advance.


## Seventh-pass checkpoint — completed locally September 5, 2026

The full goal remains active and unproved. Read
research/parallel7-pass-summary-2026-09-04.md for the latest assessment.
Three parallel workers hit the account usage limit. Root preserved
their drafts, completed the missing subgroup verifier, audited the
necklace primary source and boundaries, and ran all four verifiers
successfully. This is an external worker interruption, not a proof
that further mathematical progress is impossible.

Accepted results:
- Signed Fourier multiplier plus one positive-measure Jensen step:
  M<=(17+log N)^(1/9)N^(8/9) on the fifth-pass sixth-energy class,
  replacing exponent 17/18. No parity or zero-triple exclusion.
- Centered eighth-energy feedback exponent 43/9. A hypothetical
  N^4-scale centered E4 would imply exponent 7/8 with explicit constant.
- Independent subgroup review accepted the multiplier theorem;
  actual negative multiplier cosets and fixed-real-order optimization
  establish the limits of the current substitutions.
- Rank-four middle convolution and a local multiplicity mismatch give
  |N(ABABCC)|<=56p^(7/2) and <=58p^(7/2) on its 36-word orbit.
  Fixed length-six coverage is 687/729; remaining classes have 18,18,6 words.
- Classical conditional fourth-moment identity, local robustness, and
  an artificial affine quartic countermodel with the true mean and
  variance. No new uniform arbitrary-set upper bound.

All four final verifier outputs are in results/parallel7_*_2026_09_04.json.
The pass audit pins source, input, artifact and central-file hashes.
The Katz Rigid Local Systems primary PDF is archived and manifested;
complete relevant pages were visually checked. No new Lean proof,
prize submission, external communication or automation was made.
The uniform subgroup target, exceptional primes, logarithmic-depth
centered moments, signed aggregate, full operator and official prize
bridge remain open. Do not present the fixed-word fraction as a
percentage of completion of the conjecture.


## Eighth-pass checkpoint — September 5, 2026

Previous goal turn: progress, with the seventh pass fully integrated.
This turn: root proved the three remaining length-six classes via a
second middle convolution, giving rank six against rank two.
Read research/parallel8-necklace-2026-09-05.md and
research/parallel8-pass-summary-2026-09-05.md.

The additional class sizes are 18,18,6; their uniform constants are
48,47,45 times p^(7/2). All 729 words of length six are now covered.
The proof keeps both infinity constants, finite spikes and every
exceptional fiber. The nested graph correction is -p*t4+2; the two
nonadjacent cases are direct matrix identities.

The verifier experiments/parallel8_necklace_2026_09_05.py completed
successfully on p=5,13,17,29,41, with integer and exact surd comparisons.
Results and final file/source hashes are pinned in
results/parallel8_necklace_2026_09_05.json and the pass audit.
No new parallel worker was started while the existing workers remain
interrupted by their account usage limit. Root performed this work.

The full goal remains active and unproved. Next: growing-depth
rank/monodromy transformations with actual middle-extension stalk
corrections, then the weighted aggregate. A finite all-word theorem
alone cannot be substituted for the needed spectral conclusion.
The subgroup estimate and its open obligations remain as in pass seven.
No Lean formalization, prize submission or external communication
was made in this pass.


## Ninth-pass checkpoint — September 5, 2026

Previous turn classification: progress, with fixed length six complete.
Current pass: a source-dependent proof for every individual necklace
at every fixed degree and length. Read
research/parallel9-all-degrees-2026-09-05.md and
research/parallel9-pass-summary-2026-09-05.md.

For k>=3 the bound is 3a(2a+2)^(k-2)p^((k+1)/2), uniformly in
distinct anchor positions. Every nonempty twist followed by quadratic
middle convolution strictly raises the principal rank. Exact compact
convolution preserves a strict weight gap for all boundary errors;
their constituent mass is at most (2a+2)^j. Rank-one pairing then
gives the final character cancellation. The proof establishes the
fixed-degree statement of Kunisky Conjecture 1.14 and hence the weak
spectral convergence supplied by published Theorem 1.17. It has been
checked locally; independent mathematical review remains outstanding.

The verifier enumerates 114510 transitions/inverses and 1691778 local
block identities for a=1..5 at respective depths 40,8,5,4,3. Original
integer sums check 1170 endpoint identities, 980530 trace-sign entries,
340 literal paths, 6821 mask rows, 420 low-order identities and 410
literal necklaces. Results and final hashes are recorded in
results/parallel9_all_degrees_2026_09_05.json and
results/parallel9_pass_audit_2026_09_05.json. The small-field bound
checks are explicitly not asymptotic evidence.

The BBD primary author PDF is archived and manifested. Its complete
relevant scanned pages and the additional Katz exact-functor pages
were visually inspected. All pass-eight artifacts are preserved.
The three existing workers are terminal with account usage errors;
root continued locally without resetting or buying credits.

The full goal remains active and unproved. Next necklace work must
control signed aggregation: the new exponential constant still fails
the edge criterion. The subgroup square-root target, exceptional
primes, centered eighth-energy budget, full restricted operator,
classical arbitrary-set bound and exact official prize bridge remain
open. No Lean proof, external communication, prize submission or
automation was made in this pass.


## Tenth-pass checkpoint — September 5, 2026

Previous goal turn: progress, with the all-degree individual necklace
argument integrated. Current pass: simultaneous quadratic bounds for
the whole word family and a proved dimension obstruction to pushing
an unrestricted isometry to arbitrary depth. Read
research/parallel10-word-aggregate-2026-09-05.md and
research/parallel10-pass-summary-2026-09-05.md.

Inverse middle convolution recovers the last label by signs of local
Jordan-block differences, then recovers the whole word. Principal
word sheaves are geometrically distinct and arithmetically self-dual
with twist |w|. Their Gram error is at most (a*sum(d_w^2)+1)/sqrt(p).
The ninth-pass weight gap transfers this to actual raw character paths.
Fresh adjacency conditions are handled by scalar quadratic inertia
at a new anchor and exact half-valued boundary corrections.

An exact overwrite recurrence tracks raw generic ranks with every
boundary constituent retained. Its sum-of-squares budget grows as
Lambda_a^m; Lambda_2=(9+sqrt(65))/2 and Q starts 1,17,199,2001.
The useful all-word isometry range is Lambda_a^m=o(sqrt(p)), with
2^s(a+s) as the additional relative-error cost for s fresh anchors.

When (2^a-1)^m>p-a-1, a nonzero coefficient vector evaluates to zero.
The checker supplies an exact integer null vector for p=13, A={0,1},
y=2, m=3 (27 columns and 10 field points). Hence a uniform relative
error below one for every coefficient vector is impossible there.
Do not extrapolate this result to the needed m/log(p)->infinity.
The next spectral step must exploit its specific coefficients or a
quotient of word space, while retaining full cyclic J/anchor terms.

The exact verifier records 7541 word recoveries/distinct signatures,
773 constituent-inventory comparisons, 114 convolutions, 573 literal
entries, 172 masks and 70 positive-definiteness certificates for
35 Gram/cell bounds; six bounds are nonvacuous. The final results and
audit pin proof, script and sources and preserve all pass-nine artifacts.
Katz PDF page 49 was newly rendered and visually inspected for duality.
The proof is checked locally; independent mathematical review is open.

The full goal stays active and unproved. Subgroup/classical/prize
obligations are unchanged. The parallel workers remain terminal at
their account limit. No credit reset, purchase, external message,
submission, Lean formalization or automation was made.


## Eleventh-pass checkpoint — September 5, 2026

Previous goal turn: progress, with the all-word quadratic estimate and
dimension obstruction integrated. Current pass: a complete aggregate
bound at an explicit logarithmic range, retaining all corrections.
Read research/parallel11-full-energy-2026-09-05.md and
research/parallel11-pass-summary-2026-09-05.md.

For T=(S/sqrt(p)-J/p)(bE-I)/sqrt(N), b=2^a, N=b-1, define
H_j=||T^j||_F^2. Then H_(j+1)=H_j/N+(N-1/N)||ET^j||_F^2,
|tr T^(2j)|<=H_j, and H_(j+1)-2H_j+H_(j-1)=((N-1)^2/N)tr T^(2j).
The converse H_j<=p+((N-1)^2/N)j^2(|tr T^(2j)|/2+p) prevents
claiming that positive energy has weakened the original difficulty.

The actual coefficient-vector budget grows as Theta_a^j, where
Theta_a=[2^(a-1)(a+1)-1]^2/(2^a-1), and Theta_2=25/3.
The proof restores all finite rows, exceptional anchor columns,
the constant direction, J and the anchor correction inside powers.
For fixed a>=2 and epsilon>0, j<=(1/2-epsilon)log(p)/log(Theta_a)
implies H_j<=(1+o(1))p and the same bound for the full normalized
even aggregate. This range does not yield an improved clique bound.

The next sufficient target is an o(j) upper bound for the exact
cumulative logarithmic row-energy bias at j/log(p)->infinity,
uniformly in the required anchor cliques. No positivity, independence
or monotonicity of individual biases is assumed. The full Paley,
uniform subgroup and exact prize targets remain open.

The verifier checks 28 clique cases, 168 energy recurrences, 140
curvature identities, 168 comparisons in each trace/energy direction,
168 complete upper bounds, 168 constant-direction identities, 672
boundary checks, 30 exact all-word matrix sums and 117120 pointwise
bounds. It also checks an abstract edge Jordan block and a p=10009
first-step bound below 2p. Results and the final artifact audit pin
all inputs and preserve pass ten. The proof remains source-dependent
and locally checked; independent mathematical review is outstanding.
No reset, purchase, new worker, external message, prize submission,
Lean formalization or automation was made. The full goal stays active.


## Twelfth-pass checkpoint — September 5, 2026

Previous goal turn: progress, with the complete aggregate bound integrated.
Current pass: a proved obstruction to bootstrapping from the projection
identities and short-depth estimates, plus an exact finite Paley outlier.
Read research/parallel12-bootstrap-obstruction-2026-09-05.md and
research/parallel12-pass-summary-2026-09-05.md.

For v=e_x-e_y supported on two common neighbors, s=2, w=Pv and
q=v^TPv=1-chi(x-y)/sqrt(p), put u=sw-qv. The update
P^-=P-ww^T/q+uu^T/[sq(s-q)] is a projection of the same rank,
annihilates v and 1, and preserves every anchor column. The alternate
P^+=P-ww^T/q+vv^T/s fixes v. S'=sqrt(p)(2P'-I+J/p) retains
S'^2=pI-J and the anchor characters, but its non-anchor entries fail
the exact zero-diagonal/plus-or-minus-one condition. This is a
countermodel to the projection-only bootstrap, not to Paley.

Fixed-word traces change by at most 4k p^(k/2). Full transfer energies
satisfy H'_j <= (sqrt(H_j)+2sqrt(2)jN^(j/2))^2 and H'_j >= N^j.
Thus the pass-eleven short-depth bound survives while the required
long-depth condition fails. Any further extension needs arithmetic
information not preserved by the update. This is no new asymptotic
upper bound.

For the actual p=257 graph with clique {0,1,62}, the 32-coordinate
integer vector in the proof has squared norm 152, coordinate sum zero,
and character quadratic form -1624. Its normalized Rayleigh quotient
is -812/(19sqrt(1799))<-1, certified by square margin 9905. A second
margin 17702209 proves transfer spectral radius >9/8, hence energy
H_j>(81/64)^j at this fixed prime. It excludes a zero-error finite
edge interval, not asymptotic convergence or the Paley conjecture.

The exact verifier checks six projection constructions at p=13,17,41,
366 preserved anchor entries, 78 word perturbations, 44 power trace
bounds, and 22 each of power norm bounds, planted powers and energy
recurrences. It independently reconstructs the actual p=257 cell and
integer certificate. The final audit pins new artifacts and preserves
all four pass-eleven artifacts. The official arXiv HTML was revisited;
no new PDF input, Lean formalization, external message, submission,
reset, purchase or automation was made. Existing workers remain
terminal at the previously observed account limits; root worked locally.

The goal remains active and unproved. Next: exploit exact character
structure to control the actual rare outliers, allowing asymptotically
vanishing edge error. The source-dependent earlier arguments need
independent review, and subgroup/classical/prize obligations remain open.


## Thirteenth-pass checkpoint — September 5, 2026

Previous goal turn: progress, with the projection obstruction and actual
finite outlier integrated. Current pass: a source review followed by a
smaller principal-rank budget and stronger principal/error correlation.
Read research/parallel13-principal-budget-2026-09-05.md and
research/parallel13-pass-summary-2026-09-05.md.

For each raw error E_w, c_y(E_w)=1-1=0. Nonnegative conductors of
simple constituents imply it is lisse at y, unlike the irreducible
principal F_w. For equal-length words the normalized principal/error
correlation is at most a*d_u*e_v/p. It remains source-dependent on
geometric irreducibility, the weight gap, and tame curve cohomology.

For a>=2, h=2^(a-1), N=2^a-1, the total principal rank D_m has
rate lambda_a, the larger root of lambda^2-h(a+1/a)lambda+N=0.
Let u=(a-1)N/(lambda-N), A*=a(u-1)/(a-1), B*=(u+a)/(a-1),
f_+(k)=uk-A*(1-a^(-k)), f_-(k)=uk-B*(1-a^(-k)). These are
positive and comparable to k. Applied to the summed local Jordan
counts, their weighted sum Q_m obeys Q_(m+1)=lambda Q_m+N^(m+1),
Q_0=1. Thus D_m is comparable to lambda^m with explicit constants.
This rate is strictly below the raw rate h(a+1)-1.

For full transfer energy, replace the generic budget by
V_j=1+(a D_j^2/N^j+1)/sqrt(p)+[2a D_j E_j+E_j^2]/(p N^j),
E_j=L_j-D_j. Keep U_j=pV_j+4(a+1)L_j^2/N^j+(a+1)N^j and
H_j <= (sqrt(U_j)+j d0 N^((j-1)/2))^2 with the prior d0. The
new range is the minimum of (1/2-epsilon)log(p)/log(Psi_a) and
(1-epsilon)log(p)/log(Theta_a), where Psi=lambda^2/N and Theta
is the prior raw rate. For a=2, Psi=(19+5sqrt(13))/6, and the
leading depth coefficient rises from 0.2358197 to 0.2747391.
This is about 16.5%, still far short of j/log(p) tending to infinity.

The exact verifier passes: 8516 enumerated words, 5053 complete error
inventories unramified at y, 32 aggregate levels, 729 weighted
recurrences, 1458 coefficient identities in exact quadratic fields,
and 168 new full upper bounds; all 168 are below their prior versions.
A p=10009 case is nonvacuous. Complete Katz PDF pages 49,50,63,64,
67,107,109,139 were visually rechecked. No new BBD page, independent
review or Lean formalization is claimed. The final audit preserves
pass twelve and pins all inputs. Existing parallel workers remain
unavailable at their previously observed account usage limits.

The goal remains active and unproved. The remaining principal/principal
Frobenius correlations still have an exponential rank budget. Next work
must control them at much longer depth or use a different arithmetic
argument for rare extreme eigenvalues. The projection countermodel,
finite p=257 outlier, subgroup/classical obligations and exact official
prize gap remain in force. No account reset, purchase, external message,
prize submission or automation was made.

## Fourteenth-pass checkpoint — September 5, 2026

Previous goal turn: progress, improving the logarithmic full-energy
range. Current pass: an exact elliptic-group model using the actual
character entries, rather than another rank-budget refinement.
Read research/parallel14-elliptic-model-2026-09-05.md and
research/parallel14-pass-summary-2026-09-05.md. Goal active, unproved.

For three anchors {0,1,r}, E:y²=x(x−1)(x−r), G=E(F_p), the map
rho(P)=x(2P) covers the common neighborhood C eight-to-one. The
exceptional fibers over 0,1,r,infinity have four points each and are
exactly G[4]. All 16 are rational. The root-triple inverse has
D=(1−r)a+rb−c=(a−b)(a−c)(b−c), t=r+r(r−1)(a−b)/D,
y=r(r−1)/D; equivalently t=(a+b)(a+c), y=(a+b)(a+c)(b+c).
This gives a direct inverse for every a²=x,b²=x−1,c²=x−r.

With f(O)=0, f(x,y)=chi(y), f is even and G[2]-periodic. The full
projective character kernel is M(P,Q)=f(P+Q)f(P−Q), including every
exceptional pair. On normalized fiber indicators the border is
[[8S_C,8sqrt(2)1],[8sqrt(2)1^T,12]] plus −4I_3. On G₀=G minus
G[4], the centered kernel (M_G₀−J_G₀/sqrt(p))/(2sqrt(7p)) has
the normalized localized spectrum plus zeros. Do not drop the border
from the full matrix or the J term from the punctured restriction.

The maps x, r/x, (x−r)/(x−1), r(x−1)/(x−r) form V4 and preserve
S_C. Character dimensions are (m+epsilon0*t0+epsilon1*t1+
epsilonr*tr)/4 with fixed counts ti<=2. The existing p257,r62
outlier is in character (1,−1,1,−1), with eight free orbits represented
by 2,16,18,26,30,60,123,141. The written 8×8 block and vector
w=(-1,1,-3,3,-3,1,2,2) give norm²38, quadratic form −406 and the
same Z Rayleigh −812/(19sqrt1799)<−1, square margin9905. This
is only the previous finite outlier in smaller coordinates.

On H=G/G[2], order L=2m+4, all fbar Fourier coefficients F_eta
satisfy |F_eta|<=sqrt(p). The proof twists the rank-one Kummer
sheaf chi(y) by an unramified Lang-character sheaf. Four nontrivial
tame punctures give dim H_c^1=4 and full sum bound4sqrt(p); the
quotient sum is one quarter of the full sum. The Fourier-basis entry
of Mbar is L^-1 sum_{gamma²=alpha/beta} F_gamma F_(gamma beta).
There are four roots or none; there are nonzero off-diagonal entries.
No operator bound follows merely by treating the F_eta as eigenvalues
of Mbar. No additive-phase Weil-representation norm theorem applies
without a new identification, which was not obtained here.

Exact verification: 27 curves, 136704 kernel entries, 976 inverse
triples, 5632 periodicity checks, 108 character projectors, 352
positive leading minors for the finite Fourier bounds, and 4768
cyclotomic DFT entries with three off-diagonal witnesses. The complete
border, centered compression and reduced outlier are checked. Source
inputs: primary Bekker–Zarhin v2 HTML, Theorem2.1; complete Katz
GKM PDF spreads21,25,34,35 visually reviewed for printed32–33,
40–41,58–61. No independent review or Lean proof is claimed.

The final audit preserves all five pass-thirteen artifacts and all
preceding source-manifest entries. The three existing workers were
rechecked and remain terminal at account usage limits (reported
retry September11); this pass ran locally. No reset, purchase,
external message, submission, automation or memory write was made.

Next: exploit correlations of the elliptic Fourier coefficients for
the exact centered norm, or resume another arithmetic long-depth
approach. This new representation does not improve the asymptotic
edge/clique bound yet and only treats three-anchor localization.
The previous logarithmic range, uniform subgroup and exceptional-prime
obligations, classical arbitrary-set conjecture and exact official
Reed–Solomon prize bridge remain open. Do not mark the goal complete.


## Fifteenth-pass checkpoint — September 5, 2026

Previous goal turn: progress, the complete elliptic model and scalar
Fourier estimate. Current pass: all translated correlations, followed
by an abstract obstruction that preserves sign-kernel structure and
square-root orders with larger constants. Read research/parallel15-elliptic-correlations-2026-09-05.md and
research/parallel15-pass-summary-2026-09-05.md. The full goal is active
and unproved; no new Paley clique or asymptotic norm bound.

For shifts t_i in H=G/G[2], let O contain odd-multiplicity shifts and
B positive even-multiplicity shifts. Product_i f(h+t_i)=g_O(h) minus
sum_(b in B)g_O(-b)1_(h=-b). For s=|O|>0, every character-twisted
g_O sum has absolute value <=s sqrt(p); full products are bounded
by s sqrt(p)+|B|. If O is empty, the sum is exactly L times the
trivial-character indicator minus the missing-point character values.
The proof has4s disjoint E[2]-translated tame punctures upstairs,
no geometric invariants after the unramified Lang twist, dimH_c^1=4s,
and divides the full sum by4. Higher Fourier convolution formula is
L^(1-k) sum_(alpha1...alphak=eta) product_j F_alphaj alphaj(t_j).
It is not an arbitrary contraction of the Fourier indices.

M²_uv is the four-shift product sum for {u,-u,v,-v}. Its strata are
L-1 on identical H[2] rows; L-2 on distinct H[2] rows or ordinary
v=+/-u; <=2sqrt(p)+1 when exactly one row lies in H[2]; <=4sqrt(p)
when both rows are ordinary and distinct modulo sign. Every boundary
correction is explicit. Entrywise accumulation still loses the norm.

For all eligible p>=10000, H=Z/n x Z/d with d|n both even has
L>=16N, N=ceil sqrt(p), by Hasse. Choose k=min(N,floor(n/4)),
ell=ceil(N/k), A=[1,k]x[0,ell-1], a=kell in [N,2N). It avoids
H[2] and its negative. D=(A+A) union(A-A) union(-A-A) has |D|<24N.
Set f'=1 on D except0, retain the other values. It is even and has
exactly the same zero. Delta changed signs give product L1 error
<=2sDelta, so all full twisted bounds become <=50s sqrt(p)+e.
The exact all-even formula is unchanged. This retains square-root
ORDERS, not the original sharp constants or Kummer realization.

M'=f'(u+v)f'(u-v) preserves all exceptional couplings, the J4-I4
border, plus/minus fibers and H[2] translations. On A union(-A)
it is four copies of J_a-I_a. The centered normalized quotient has
Rayleigh4(a-1-a/sqrt(p))/sqrt(7p), at least
(4/sqrt7)(1-2/sqrt(p)), persistently >1. It cannot retain the ambient
Paley identity: the clique A plus the three anchors has Rayleigh
a+2>sqrt(p). Do not combine this with the projection-only obstruction
as though one example retained both sets of hypotheses. It is NOT a
Paley counterexample, and does not rule out using sharp constants.

The exact verifier passes1047 parity patterns,162 twisted supports,
6864 positive leading minors,5168 M² entries,100 cyclotomic product
identities. Actual group reconstruction at p10009,r2 gives H factors
2,1252, A size101,405 modification-region points and192 changedsigns;
Rayleigh rational lower bound9899/7575. At p65537,r2 the factors
are128,128, A dimensions32,9 andsize288,2672 region points and1414
changed signs; lower bound2287/1542. Both preserve all75552 exceptional
entries in total,8192 sampled translation identities, and48 exact
L1 perturbation certificates. The large group coordinates are certified
by both generating translations and a bijection, with no floatingpoint
acceptance. General theorem proofs are separate from finite tests.

Sources: earlier complete Katz GKM review reused; primary Bekker–Zarhin
HTML Theorem4.1 additionally checked for Hasse. No new source or
source-manifest mutation, no independent review or Lean proof. Preserve
all five pass-fourteen artifacts. The previous live worker check showed
all three terminal at account usage limits; they were not restarted.
No reset, purchase, external message, submission, automation or memory
write. Next: combine the sharp arithmetic correlation information and
ambient projection rather than rely on unspecified square-root orders.
The prior logarithmic range, uniform subgroup and exceptional-prime
problems, classical arbitrary-set cancellation and exact official prize
bridge remain open. Do not mark the goal complete.

Pass-fifteen final audit additionally checks the sharp-constant failure:
the zero-frequency coefficient changes from3 to387 at p10009 and
from−1 to2827 at p65537. Square margins above p are139760 and
7926392. These witnesses are stored in the final audit; they preserve
the distinction between bounds of the same order and the original
sharp bound. Previous artifacts and the source manifest are unchanged.

## Sixteenth pass — completed field-of-definition obstruction

Previous turn classification: progress. This pass is progress, not a
proof or disproof of the prime-field goal. New files are
research/parallel16-square-field-barrier-2026-09-05.md,
research/parallel16-pass-summary-2026-09-05.md,
experiments/parallel16_square_fields_2026_09_05.py,
results/parallel16_square_fields_2026_09_05.json, and
results/parallel16_pass_audit_2026_09_05.json.

For ell an odd prime, q=ell^2, the classical subfield clique K=F_ell
has exact eigenvector w=1_K-(1/ell)1 with S*w=ell*w and P*w=w.
Its squared mass on K is 1-1/ell. For a fixed anchor set A in K,
B=K\A gives localized projection Rayleigh
1-(a+2)/(2*ell)+a/(2*ell^2). With b=2^a,N=b-1 and the existing
Z normalization, ||Z|| tends to b/(2*sqrt(N)), above one for a>=2.
At a3 the limit is 4/sqrt(7). Characteristic tends to infinity too.

These SAME actual examples retain all actual Paley entries, S^2=qI-J,
full elliptic model and exceptional border, four symmetries, original
Kummer realization and SHARP all-order twisted bounds s*sqrt(q)+e.
No constant50 or planted modification. The Lang twist over F_q uses
1-Frob_q; the same 4s disjoint tame punctures and factor4 prove the
estimate. Independent review of that application is still needed.
This shows the proposed combination of sharp bounds and ambient
projection cannot yield a field-uniform edge. Prime cardinality is
an explicit hypothesis of Kunisky; it excludes this family.

The exact transfer radius tends to sqrt(N). For a3,ell>=29, it is
greater than 3*sqrt(7)/4; energy >(63/16)^j and even trace
>(63/16)^j-q. The explicit starting margin at ell29 is665.
This rules out long-depth polynomial-q energy and trace bounds on
the square-field family. It is not a prime-field counterexample.

Final verifier exited0 using fields49,121,289,361,841. Totals:
938165 entries of each ambient identity,25 localizations,4983
subfield-eigenvector entries,12 power traces,266496 complete elliptic
kernel entries,16656 squared entries,12 sharp all-character certificates.
Elliptic quotient sizes16,80,100 have scalar certificate nullities2,4,1.
All arithmetic is exact; singular PSD pivots require remaining rowzero.
Unit cases explicitly reject indefinite matrices. The three-anchor
Rayleigh bounds exceed one at q289,361,841, and the analytic family
proves persistence separately from these finite computations.

Sources: Kunisky v1 primary HTML reopened for prime-order scope;
Asgarli–Yip abstract used for background only; existing complete Katz
GKM review reused, no newly rendered PDF pages. Five pass-fifteen
artifacts and sources/manifest.json remain unchanged.
The final live check confirmed all three workers remain terminal at
account usage limits. This pass ran locally without restarting them or consuming
a reset. No external messages, submission, automation or memory write.

Next: derive a quantitative estimate using prime field order, such as
quantitative additive Fourier uncertainty, that excludes concentration
strongly enough to reach the desired spectral edge. The qualitative
absence of a proper subfield does not supply that bound. No such
decisive estimate was proved here. The previous logarithmic range,
classical arbitrary-two-set cancellation, uniform subgroup target and
exceptional primes, and exact official prize bridge are still open.
Keep the full goal active; it is unachieved.

## Seventeenth pass — prime uncertainty and determinant limitation

Previous turn classification: progress. This pass is progress: it adds
a prime-field nondegeneracy theorem with an explicit endpoint gap and
proves why the associated determinant certificate cannot close the
spectral route. The full conjecture and official prize remain unproved.

Files: research/parallel17-prime-uncertainty-2026-09-05.md,
research/parallel17-pass-summary-2026-09-05.md,
experiments/parallel17_prime_uncertainty_2026_09_05.py,
results/parallel17_prime_uncertainty_2026_09_05.json,
results/parallel17_pass_audit_2026_09_05.json.

Let p=1mod4 prime, r=(p-1)/2, C size1<=m<=r, X=P_C,Y=Q_C,
D=X+Y=I-J/p, Delta=1-m/p. Imported Chebotarev nonzero Fourier
minors imply X,Y>0, so0<X<I. The source is Tao arXiv math/0308286v6,
Lemma1.3/Corollary1.4; this is classical, not a novelty claim.

New exact invariant k(C)=p^(m+1)det(X)det(Y) is a positive integer.
For B=2pX=pI-J+sqrt(p)S_C, rank-one J extracts(sqrt(p))^(m-1)
from detB, so its conjugate norm is divisible by p^(m-1).
The entries of pX lie in Z[omega],omega=(1+sqrt(p))/2, whose
conjugate norm is integral; gcd(p,4)=1 removes the factor4^m.
This proves the divisibility and positivity of k.

R=D^(-1/2)XD^(-1/2), theta=4^m*k/(p^(m+1)*Delta^2) in(0,1].
Every R eigenvalue mu has4mu(1-mu)>=theta, hence X and I-X
are bounded below by g(C)=Delta*(1-sqrt(1-theta))/2.
Rational lower bound g>=k*delta>=delta, where
delta=4^(m-1)/(p^m*(p-m)). All inequalities retain the J term.

For A=D^(-1/2)S_C D^(-1/2), exact
trA²=m(m-1)+2||S_C1||²/(p-m)+(1^TS_C1)²/(p-m)².
Theta=product_i(1-a_i²/p) <=(1-trA²/(mp))^m
<=exp(-m(m-1)/p). Thus g(C)<=Delta*exp(-m(m-1)/p)/2.
THIS UPPER BOUND IS ON THE CERTIFICATE g, NOT THE ACTUAL GAP.
At fixed-anchor m~p/2^a the certificate is exponentially small;
the desired edge needs gap1/2-sqrt(2^a-1)/2^a>0 for a>=2.
This specific determinant argument is insufficient, even retaining
the exact k. It does not rule out more information from minors/ratios.

Separate prime Fourier example: p=4h+1, frequency interval[h+1,3h]
gives a real rank-r circulant projection Pi with Pi1=0. For coordinate
interval0..n,m=n+1<=r, v_j=(-1)^j binomial(n,j), Fourier leakage
<=2^n/binomial(2n,n)<=(2n+1)2^-n. S'=sqrt(p)(2Pi-I+J/p)
has zero diagonal and (S')²=pI-J, but S'01~ -2sqrt(p)/pi,
so loses Paley signs at large p. Smallp5 is a coincident Paley case;
p13,29 sign failures are checked exactly. C is not asserted to be
an anchor neighborhood; no complete Paley counterexample claimed.

Final verifier exit0:4132 compressions,41420 positive leading minors
across both embeddings,4132 independent integer determinant identities,
4132 exact second-moment/rational decay bounds,518 shifted-gap minors,
226 independent permutation determinant checks,47 complete cyclotomic
convolution rows across three interval projections. Exhaustive sets
containing0: p5,size<=2 count5; p13,size<=6 count1586; p17,size<=5
count2517. These are supplemented by24 actual localizations throughp257.
At p29,m15 and p41,m21, beyond rankr, determinants are exactlyzero.
Binomial leakage examples p257,1021,4093,12289 are all below1e-6.
At p257,anchors0,1,62,m32, certificate upper bound<1/100 while
the target gap for a3 exceeds1/6; both comparisons exact.

Primary Tao sources archived: PDF130591 bytes/hash
244f4e79e667ee831e65d0cb0d8d42de4e354e83d576c65acaf391b0da519880,
HTML161191 bytes/hash
268bcdcd3a1bcb487c1ea03d6927d4f63e72ba4e0312e2bd3ecc58f1e33f20f8.
Root rendered and visually reviewed complete PDFpages2,3,4,5,
checked theorem hypotheses/proof against versionedHTML. No other
PDF is newly reviewed. Cauchy-Davenport is not used. Manifest has
two appended entries, preserving the precedingtwenty. Five pass16
artifacts are unchanged; /tmp/paley17-before.json pins them and the
old manifest, /tmp/paley17-manifest-before.json preserves its contents.
New deductions have no independent review or formal proof yet.

Three workers were terminal at account limits at the last live check
(pass16). They were not restarted or represented as running. Root
worked locally. No reset, purchase, external message, submission,
automation or memory write. Keep the full goal active and unachieved.
Next: quantitative control of quadratic-residue Fourier mass on
actual anchor neighborhoods, or another uniform long-depth argument;
qualitative support nonvanishing and this determinant product are
insufficient. The prior logarithmic range and clique bounds, classical
arbitrary-two-set cancellation, uniform subgroup target/exceptional
primes, and exact official prize bridge remain open.

## Eighteenth pass checkpoint: Cartesian expansion and centered Weil trace

Files: research/parallel18-mobius-trace-2026-09-05.md,
research/parallel18-pass-summary-2026-09-05.md,
experiments/parallel18_mobius_trace_2026_09_05.py,
results/parallel18_mobius_trace_2026_09_05.json,
results/parallel18_pass_audit_2026_09_05.json.

p=1mod4 prime, A,B nonempty, sizes m,n. Set
g(s,t)=[[s,st-1],[1,t]]=u_s w u_t. Its action is s-1/(z+t).
Rational point transporters have at most max(m,n) elements;
finite-finite and nonrational quadratic-extension transporters at most
min(m,n). By Dickson classification every proper-subgroup coset meets
G(U,V) in at most max(2max(m,n),120). For m,n>=p^epsilon,
the normalized mass is <=p^(-epsilon/2) for sufficiently large p.
Lyamkin Theorem 11 then gives normalized operator norm
O_epsilon(p^(-kappa(epsilon))) on every nontrivial irrep, kappa>0
unspecified. Existing expansion theorem, not a new one. Exact energy
E(G)=n²E+(U)+m²E+(V)-m²n², proved by computing g^-1 g'.

Normalized Weil model: U_t f(x)=e_p(t*x²)f(x),
F_xy=p^(-1/2)e_p(-2xy). For c!=0, kernel
rho([[a,b],[c,d]])_xy=chi(c)/sqrt(p)*e_p((a*x²-2xy+d*y²)/c).
For c=0 it is chi(a)e_p(ab*x²) times 1_(y=ax).
Right-generator identities and the Gauss formula prove the group law.
Character equals chi(tr(g)-2) off the trace-two locus; for trace-two,
c!=0 it is sqrt(p)*chi(c); for nonidentity upper unipotent it is
sqrt(p)*chi(b); identity trace is p. Never omit this exception.

T_AB averages rho(g(a+2,-b)). Exact
tr(T_AB)=[S(A,B)+sqrt(p)*|A intersect B|]/(mn).
Kernel factorizes as p^-1/2 e_p(2x²-2xy) alpha_A(x²)alpha_B(-y²).
Its00entry=p^-1/2, so norm>=p^-1/2; dimension-times-norm cannot
give a trace saving. Exact squared Frobenius norm equals
(1/p)*(p/m+sqrt(p)S(A,A)/m²)*(p/n+sqrt(p)S(B,B)/n²).
Diagonal Cauchy gives |S(A,B)+sqrt(p)overlap|²
<=(sqrt(p)m+S(A,A))*(sqrt(p)n+S(B,B)). No new threshold follows.

For p>=13 and g0=[[3,-1],[1,0]], conjugacy average
V=(I+chi(5)R)/(p+chi(5)), Rf(x)=f(-x). Commutant span{I,R}
proved directly by commuting with U_t and F. Parity modules are
irreducible/nontrivial. TraceV=1, normV=2/(p+chi(5)),
tr(V^j)=((p+chi(5))/2)^(1-j). This is an actual positive probability
average and a PSD scalar parity projection, but not asserted to be
G(A+2,-B). It is a generic trace-recovery limitation, NOT a Paley
counterexample. The missing centered-trace power saving is EXACTLY
the original Paley target; no new independent sufficient lemma proved.

Final verifier exit0 after integer-primality cutoff correction:
1324 quadratic Gauss identities;744192 right-generator kernel checks;
2304 all-group character checks at p5,13;1028 Cartesian cases (all961
nonempty subset pairs p5 plus67 p13);6704 rational transporter rows;
29672 nonrational quadratic-extension rows;1028 exact energies,
centered traces and factorized norm identities each. Complete class
averages p13,17,29 check1299 matrix entries, summing836642 entries.
Cyclotomic identities, integer finite-field/group counts and rational
norms/powers use no floating-point acceptance. These finite checks
do not establish the imported asymptotic theorem or cancellation.

Thomas primary HTML arXiv math/0610644v3 archived661651bytes,
hash5abc992c438e373bfcf243ec61247d9ac3152dd5973f0a2d64a3342a3984d7b5.
Sections2 and3.2 reviewed for trace/unipotent/normalization inputs.
Lyamkin primary live MathNet HTML sm9707e Section1.1 definitions and
Section1.5 Theorem10/Lemma5/Theorem11 reviewed. Direct HTML/PDF
archival attempts returned403. Metadata-only source manifest entry;
no local archive or content hash claimed. No PDF newly rendered or
visually reviewed in this pass. Manifest22old+2new=24entries.
/tmp/paley18-before.json pins five pass17 files and the old manifest;
/tmp/paley18-manifest-before.json preserves the old manifest bytes.

Live worker check in pass18: all three existing workers remain
terminal at account usage limits, retry September11. No restart,
reset, purchase, external message, submission, automation or memory
write. Root worked locally. Independent review of the new arguments
and the full goal remain outstanding. Next useful work must exploit
additional Cartesian character structure, retaining the overlap
correction, or return to a different uniform estimate. Do not treat
operator expansion as proof of Paley. Keep goal active/unachieved.

## Nineteenth pass checkpoint: B_h reduction and SS-B

Files: research/parallel19-relation-free-reduction-2026-09-05.md,
research/parallel19-pass-summary-2026-09-05.md,
experiments/parallel19_relation_free_2026_09_05.py,
results/parallel19_relation_free_2026_09_05.json,
results/parallel19_pass_audit_2026_09_05.json.

Previous turn classified as progress after checking pass18 artifact
hashes against its audit. Root continued locally; workers last checked
pass18 were terminal at usage limits, retry September11. No restart,
reset, purchase, external message, submission, automation or memory write.

For p>h>=2, a B_h set has no nontrivial equal h-term sums WITH
repetitions. If S is such a set, forbidden extensions are exactly
F_h(S)=S union union_(j=1)^h j^-1(hS-(h-j)S). Cancel repeated x
and pad old summands to prove the converse; p>h permits division.
F_h(empty)=empty. Cardinality <=f_h(s)=s+sum_j s^(2h-j).
L_h(k)=f_h(k-1). Any D with |D|>L_h(k) has a greedy B_h
subset of size k. Repeated removal gives equal-size blocks and
remainder <=L_h(k)<=(h+1)k^(2h-1). k1 has no remainder.

For any |K(a,b)|<=1, two partitions with remainder sizes rA,rB
give normalized sum <=eta*(m-rA)(n-rB)/(mn)+rA/m+rB/n-rA*rB/(mn),
if eta bounds normalized complete-block sums. Therefore <=eta+L/m+L/n.
If all B_h pairs at every polynomial size exponent satisfy Paley,
choose theta=epsilon/[2(2h-1)], k=floor p^theta and apply their
bound at exponent theta/2. Remainders <=O_h(p^-epsilon/2).
This proves full all-exponents Paley equivalent to its restriction
to B_h pairs for EACH FIXED h. No restricted cancellation proved.

Moment transfer: if U bounds M_(2r)(C) on B_h sets of size k,
|S(A,B)|/(mn)<=(U/(m*k^(2r)))^(1/(2r))+L_h(k)/n.
Exact integer version with actual q blocks and actual maximum U:
(|S|-m|R|)_+^(2r)<=q^(2r)m^(2r-1)U; U0 when q0.

New sufficient SS-B: r_j unbounded fixed integers, h_j>=2 with
h_j/r_j->0, beta_j>=0->0; for each j uniformly all B_(h_j)
sets of size k_j=floor p^(1/(r_j+1)), assume
M_(2r_j)(C)<=C_j p^(1+beta_j) k_j^r_j.
Then full classical Paley follows. With theta=1/(r+1), bound is
(C_j2^r)^(1/(2r))p^((theta+beta-epsilon)/(2r))
+(h+1)p^((2h-1)theta-epsilon). Choose one j with both theta+beta
and (2h-1)theta <epsilon/2; delta=epsilon/(8r_j) eventually works.
Every j, h_j and constant is fixed before the prime grows. h_j2
or h_j=floor sqrt(r_j) are allowed. NO CONVERSE TO SS-B proved.

At fixed r,h, need epsilon>max(theta+beta,(2h-1)theta). In
particular r=h2, beta0 gives epsilon>1, so Sidon-only fourth moment
does NOT reproduce epsilon>1/3 through this extraction. The interval
example gives binom(k+h-1,h)<=h(N-1)+1 if h(N-1)<p, showing a
real size cost but not optimality or necessity of h_j=o(r_j).
Past forced-row obstruction at n=p^(1/r) for r>h does not refute
SS-B at n=p^(1/(r+1)); its normalized contribution there tends to0.
It does not follow that SS-B is true.

Also computed simplest mixed Weil quotient for p1mod4:
g(s,t)^-1g(s',t') has trace2-de with d=s'-s,e=t'-t.
Ordinary character chi(de), one-zero cases sqrt(p)chi(d or e),
both-zero p. Sum=p mn+sqrt(p)(nS(U,U)+mS(V,V))+S(U,U)S(V,V),
exactly old Hilbert-Schmidt identity; no added cancellation estimate.

Final verifier exit0:2406 exact forbidden sets;38006 candidate
extensions;33276 greedy extensions;26298 minimal-energy checks;
172 partitions;74 two-sided transfers;296 integer moment transfers
at moment orders2,4,6,8;12 arbitrary bounded kernels;91 interval
size checks;4 sufficient-sequence and203 fixed-order rational exponent
checks;112707 independently multiplied SL2 quotient cases.
The exponent ledger with r as large as39500 only checks rational
inequalities, not the unproved moment hypothesis. No floating acceptance.

Primary background HTML Shkredov arXiv2103.14670v1 archived394457bytes,
SHA256bb9624db0131540fe4feb5f414814f0631f5c05a59a8dfd32d9b5de030175f6e.
Introduction reviewed for terminology/context; no theorem imported.
No new PDF rendering/visual review. Source manifest24old+1new=25.
/tmp/paley19-before.json pins five pass18 artifacts and old manifest;
/tmp/paley19-manifest-before.json preserves old manifest bytes.

Next: prove an actual character-sensitive uniform estimate on the
SS-B slice, or another sufficient input. Additive relation exclusion
alone has not supplied it. Do not call the reduction a proof of
Paley or of SS-B. Independent review, all uniform targets, the
spectral edge/clique improvement and official prize bridge remain
open. Goal remains active and unachieved.

## Twentieth pass checkpoint: inversion strengthens SS-B

Files: research/parallel20-inversion-moments-2026-09-05.md,
research/parallel20-pass-summary-2026-09-05.md,
experiments/parallel20_inversion_moments_2026_09_05.py,
results/parallel20_inversion_moments_2026_09_05.json,
results/parallel20_pass_audit_2026_09_05.json.

Pass19 hashes checked against its audit; previous turn is progress.
Root worked locally. Most recent live worker check in pass18 found
all three terminal at usage limits, retrySeptember11. No restart,
reset, purchase, external message, submission, automation or memory write.

For p>h>=2, C sizek, at most
Q_h(k)=k+(2h-2)binom(binom(k+h-1,h),2) poles z either lie in C
or fail to make D_z={1/(c-z):c in C} a B_h set. For two distinct
h-multisets, cancel repeated common entries to get sum nu_c/(c-z),
nu nonzero integers |nu|<=h,total0,supportsize<=2h. Multiply by
the distinct denominator product. Polynomial P has degree<=2h-2
because totalnu0; P(c)=nu_c product_(d!=c)(d-c)!=0 for each
supportpoint, since p>h. Union bound over unordered multiset pairs
givesQ. Do not omit p>h or the size hypothesis. In F5,C0,1,2,
both available poles fail to give Sidon; no unconditional all-size claim.

Signs w_(1/(c-z))=chi(c-z). For x!=z,t=1/(x-z), exact
F_D,w(t)=chi(-1)chi(x-z)F_C(x); at t0 valuechi(-1)k.
Thus M_(2r)(D,w)=M_(2r)(C)-|F_C(z)|^(2r)+k^(2r)>=M_(2r)(C).
Works for p1mod4 and p3mod4, preserving the missing row exactly.

If p>f_h(2k-1), use pass19 to extend any B_h D of sizek to
D union E of size2k. LetP,N be sign classes, |P|=a,s=2a-k.
C+=P union random(k-a)-subset(E), C-=N union randoma-subset(E).
Pointwise w1_D=E1_C+ -E1_C- +(s/k)1_E. All indicated sets have
sizek and are B_h. The triangle inequality gives
M_(2r)(D,w)<=(2+|s|/k)^(2r)U0<=3^(2r)U0 when U0 bounds
unsigned B_h moments at sizek. No moment bound at size2k or at
smaller cardinalities is assumed. Signs are not probability weights.

Finite theorem: p>h, p>Q_h(k), p>f_h(2k-1) imply that a uniform
B_h sizek moment bound transfers to all sizek sets at cost3^(2r).
Restriction gives conversewithconstant1, hence equivalence up toconstant.
This transfers a hypothesized bound; it proves no independent bound.

Explicit K_h=max(2,4h(h-1),(h+1)2^(2h-1)+1). For k>=K_h,
Q_h(k)<k^(2h) and f_h(2k-1)<k^(2h). Proof uses
(h-1)/(h!)²<=1/4 and (1+(h-1)/k)^(2h)<=2 by a geometric-series
comparison; then Q<=k+k^(2h)/2<k^(2h). Completion estimate
f_h(2k-1)<=(h+1)(2k)^(2h-1)<k^(2h) is immediate fromK.

New SS-B*: on k_j=floor p^(1/(r_j+1)), assume uniform
M_(2r_j)(V)<=C_j p^(1+beta_j)k_j^r_j only for B_(h_j) sets,
where r_j unboundedfixedintegers>=3, beta_j>=0->0, and
2<=h_j<=floor((r_j+1)/2). For primep>=K_h^(r+1), finite theorem
applies. It gives unrestrictedSSwithconstant3^(2r)C_j. Original
sampling transfer proves fullclassicalPaley, choosing onej with
theta_j+beta_j<epsilon/2. delta=epsilon/(8r_j) eventually works.
The old h_j/r_j->0 restriction is superseded; boundary2h=r+1
works because Q/k^(2h)->(h-1)/(h!)²<=1/4. All orders/constants
fixed before p grows. No growing-r estimate or practical threshold claim.
Example: unproved Sidon sixth moment on p^(1/4) would imply
Paley for epsilon>1/4. Sidon fourth moment on p^(1/3) notcovered.

Final exact verifier exit0:56256 nonzero relation polynomials;
646347 rational/polynomial equivalences;1319 complete pole sets and
13168 pole candidates;5159 distinct inversion cases;112762pointwise
rows;20636 corrected moments;88 completion extensions;80 exact
minimal-energy checks;102 sign patterns;728 signed indicator
coordinates;2196 unsigned completion set moments;408 signed moment
transfers;32 completeexamples and128 original-to-B_h moment transfers.
22 shifted-polynomial threshold certificates for h2..12,4 conditional
exponent examples,240 admissible(r,h)pairs. Actual moment checks
use orders2,4,6,8; larger r in the ledger checks exponents only.
All arithmetic exact, no floating-point acceptance.

No new external theorem/source imported; manifest unchanged.
No new PDF rendering or visual review. /tmp/paley20-before.json
pins five pass19 artifacts and the whole source manifest. The script
imports the unchanged pass19 exact helper module; input hashes pin
both that helper and its proof plus the original pass5 transfer.

Next required input remains an actual uniform SS-B* upper estimate
on these highly relation-free sets, or another route. Inversion and
completion show their moment supremum is as hard up toconstants as
the unrestricted one in this range. Do not describe transfer as proof
of the moment estimate. Full classical/subgroup/spectral/prize targets
and independent mathematical review remain outstanding. Keep the full
goal active and unachieved.

## Twenty-first pass checkpoint: one-sided squarefree aggregate

Files: research/parallel21-squarefree-moments-2026-09-05.md,
research/parallel21-pass-summary-2026-09-05.md,
experiments/parallel21_squarefree_moments_2026_09_05.py,
results/parallel21_squarefree_moments_2026_09_05.json,
results/parallel21_pass_audit_2026_09_05.json.

Pass20 preserved hashes in /tmp/paley21-before.json. Previous completed
turn classified as progress. Root derived the current pass locally.
A fresh agent listing now returned only root, so root attempted three
bounded user-authorized assignments: pass21_independent_review,
squarefree_positive_upper, subgroup_prize_next_input. Their actual
outcomes must be recorded before calling them completed or running.
No reset, purchase, external message, proof submission or automation.

For signs eps with N nonzero entries and F=sum eps, H_j=j!e_j
satisfies H0=1,H1=F,H_(j+1)=F H_j-j(N-j+1)H_(j-1).
For N>=2r, H_(2r) is the characteristic polynomial of a symmetric
tridiagonal matrix, eigenvalues paired +/-lambda, |lambda|<=R,
R²<=8rN. This gives H>=-(8rN)^r,
|H|<=F^(2r)+(8rN)^r, and
F^(2r)<=2^rH+2(16rN)^r. For N<2r at actual sign sums H=0;
F^(2r)<=N^(2r)<=(2r)^r N^r handles this case including N0.
Summing actual character rows N=n or n-1 proves:
T_(2r)>=-(8r)^r p n^r/(2r)!,
|T_(2r)|<=[M_(2r)+(8r)^r p n^r]/(2r)!,
M_(2r)<=2^r(2r)!T_(2r)+2(16r)^r p n^r.
Thus one-sided T upper estimates and fixed-order moment estimates
are equivalent up to constants; negative T already has Gaussian scale.
On pass20 SS-B* slice only top squarefree positive upper bound needed.
It remains unproved for actual characters. No converse from Paley toSS.

Exact sixth ledger:
a0=15n³-30n²+16n,a2=90n²-300n+272,a4=360n-960;
M6=p a0-a2 binom(n,2)+a4T4+720T6
-15sum_(c in C)F(c)^4-15sum_(c in C)F(c)^2-n.
Do not omit these zero-row corrections. At n=floor(p^(1/4)),
termwise degree4 Weil already controls a4T4 at O(p n³).

Weighted countermodel for every prime with k=floor(p^(1/4))>=256:
P=p-k,L=floor sqrt(p),W=P-L,e=(k-1)mod2.
Extreme sum+/-k massL; zero rows massk with zero coordinate uniform,
remaining sum+/-e uniform on patterns; regular bulk massW with
v=(Pk-ke²-Lk²)/W in[k-2,k), t largest parity-k integer t²<=v,
b=t+2, lambda=(v-t²)/(b²-t²), sum+/-t,+/-b with weights1-lambda,lambda.
Every sign class uniform over patterns. Positive rational row weights.
Gram pI-J, zero mass1 percolumn, means0, odd products0,C2=-1.
C_d=L+W[(1-lambda)kappa_d(t,k)+lambda kappa_d(b,k)]
+(k-d)kappa_d(e,k-1), kappa=H_d/(N)_d.
Bounds |C4|<=2sqrt p, |C6|<=5sqrt p, M4<=3p k²,
M6/(pk³)>=k/2 ->infinity. Sidon labels exist sincef2(k-1)<k³<p.
These are rationally weighted rows, not p uniform character translates;
each zero fibre can have multiple rows. Labels do not supply character
factorization or full-field translation/projection identities. This is
NOT a counterexample to actual SS-B*, fullPaley, or an officialprize.

Exact verifier exit0:67080 recurrence/binomial comparisons,
201240 pointwise inequalities,16770 explicit H4/H6 polynomial identities,
3442 direct character product sums,22996 actual aggregate identities,
68988 actual aggregate inequalities,5749 sixth-moment ledgers;
1250 explicit weighted rows,60 columnzero/mean checks,190Gramentries,
486 distinct-subset correlations,15 explicitweightedmoments;
5 large weightedfixtures,30 individualWeil-shapedbounds,
796097 explicit Sidon unorderedpairsums;30728 factorialbounds and
3841 constantparameterchecks k256..4096. No floating acceptance.
Actual fields p3,5,7,11,13; large p4294967311,4362470461,
21743271967,68719476767,1099511627791 with k256,257,384,512,1024.
Large cases do not evaluate actual field character kernels.

Primary source McDonald-Sahay-Wyman arXiv2210.03789v2 HTML archived
339531bytes,SHA2563ae1bc2ca60ef412b271e3a03056df70fe1465fda26a89727e1929cfe9e1f622.
Section2Lemma2.1 reviewed for standard Weil comparisononly; noVCtheorem.
Manifest25old+1=26; originalbytes /tmp/paley21-manifest-before.json.
No newPDF rendering or visualreview. Fullclassical,subgroup,spectral,
prizebridge and independent mathematical review remain outstanding.
Goal active and unachieved, not blocked.

Potential next elementary model, NOT YET CHECKED or included in pass21
claims: use zero rows uniform independent signs excepttheirzero, extreme
rows allsame massL=k², W=p-k-L regularrows, and cube density
1-a e2(eps) with a=(L+1)/W relative to uniform ±1 cube. Forp>=k⁴,
positivity may follow from W>=(L+1)binom(k,2). Walshorthogonality would
give C2=-1, alloddC0, allevenC_(d>=4)=L<=sqrtp; hence all individual
orders atonce, at the cost of M4constant4 ratherthan3. On the general
slice p>=k^(r+1), extreme mass yields M_(2r)/(pk^r) growinglikek,
while every lowerfixedmoment may stay Gaussian via Rademacherpairings.
This candidate needs a separate rigorous proof and exact verification;
it remains a weightedrelaxation and supplies noactualPaleyupperbound.

Fresh worker outcome update: all three assignments ran and sent
substantive progress. pass21_independent_review completed its bounded
algebra/code review, found no mathematical corrections, and wrote
research/parallel21-independent-review-2026-09-05.md with SHA256
f399ba7212ba7b840bf0fd7c1086ee2a1180d9c7d2855c7e6ab7d36702e3d16a.
The report pins proof42bfee31ccd53abaa9081128a3dfa250c5316570c407740c439b7175759f0210,
script59094f846116bb56cd5975a2ba41986b282b3f09e6ccf44ff04d2b8fc2047269,
pass20transfer and sourceHTML. Root read the complete review. No
mathematical changes or verifier rerun were needed after that review.
This is separate-agent review, not human or formal verification; the
pinned proof's original pending annotations are superseded by this
report and current summary. Reviewworker then assigned a bounded
spectral next-input task; leave its completed reviewfile unchanged.

Current live branches after that reassignment:
- /root/pass21_independent_review now on spectral-next-input.
- /root/squarefree_positive_upper on full-translate quartic-star energy.
- /root/subgroup_prize_next_input on opposite-free sixth-energy relations.
All are user-authorized bounded tasks. Their additional notes/results
must be reviewed before treating findings as established; do not assume
reports are complete from progress messages. No usage reset or purchase.

## Twenty-second pass checkpoint: completed parallel lanes and all-order obstruction

Previous goal turn classified as progress. Root verified all five
pass21 artifact hashes against its final audit before continuing;
/tmp/paley22-before.json pins those five, the audit, and oldmanifest.
All three research workers were live on the initial authoritative
agent listing. Their completed outputs are first integrated in this
pass22 assessment, though they retain pass21 filenames.

Primary files:
research/parallel22-pass-summary-2026-09-05.md,
research/parallel22-all-orders-obstruction-2026-09-05.md,
experiments/parallel22_all_orders_obstruction_2026_09_05.py,
results/parallel22_all_orders_obstruction_2026_09_05.json,
results/parallel22_pass_audit_2026_09_05.json.
Late worker proofs:
research/parallel21-positive-upper-review-2026-09-05.md,
research/parallel21-subgroup-next-input-2026-09-05.md,
research/parallel21-spectral-next-input-2026-09-05.md.
All have corresponding exact Python verifiers and result files.
All received separate-agent cross-review, linked from the summary.
No human refereeing or formal proof certification is claimed.

Root weighted theorem: fixedr>=3,k=floor(p^(1/(r+1))),k>=2(r+1).
L=k²,W=p-k-L,alpha=(L+1)/W. Extreme vectors+/-all1 totalmassL;
cube measure totalmassW,densityrho=1-alpha e2; eachzerocolumnfibre
mass1, othercoordinates independentuniformsigns. Forp>=k⁴,k>=2,
W-(k²+1)k(k-1)>=k²(k-2)>=0 and
W-(k²+1)k>=k(k-2)(k²+k+1)>=0, so rhoin[1/2,3/2].
Gram pI-J; alloddC0,C2=-1,allevenCd=k² for d>=4. Thus every
individual distinct-product Weil-shaped bound holds through degreek.
For ALL realcoefficientvectorsb, s<=r-1,
M2s(b)<=[(3/2)(2s-1)!!+1]p||b||2^(2s).
This follows from rho<=3/2, cube pairings, and Lk^s<=p.
Forany unsignedsubsetB ofmcolumns, withU_s(m)=cube2smoment,
V_s(m)=(U_(s+1)(m)-mU_s(m))/2>=0 by covariance,
M2s(B)=Lm^(2s)+(p-L-m)U_s(m)-(L+1)V_s(m)+mU_s(m-1),
so constant(2s-1)!!+1 works. At s1 exactpm-m².
Extreme rows forceM2r/(pk^r)>=k/2->infinity; k>=2(r+1)
implies(k+1)^(r+1)<2k^(r+1). Bhlabels existfor2h<=r+1 by
f_h(k-1)<(h+1)k^(2h-1)<k^(2h)<=p. Choosemaxh forsimultaneous
lower-orderlabelproperties. LabelsdoNOTsupplyfield-differenceidentities.
Measure ispositive rationallyweighted, notpunitrowcharacters. All
lowercoefficientmoments+allindividualordersdonotforce nextmomentinthis
relaxedclass. It isNOTanactualSS-B*,Paleyorsubgroupcounterexample.

Root exactverifier exit0:508cube-densitybounds,70columnmass/mean,
203Gramentries,508allsubsetcorrelations,536unsignedformulas,
268unsignedlowerbounds,332coefficientlowerbounds,2300explicitrows;
1560binomialpairing/covariancechecks,558compressedalldegreecorrelations,
268largeWeil-shapedbounds,1922largelowerunsignedmoments,
29Bhthresholds,17genuineprime-sizefixtures r3..9,185explicitBh
multisetsums(label-onlysmallfixtures). No floating acceptance.
Independent modelreview pins proofb959300e5c9c886e2936c306a2b0846f01d47e77b795bf1c36b5325877ee926c;
no mathematicalcorrection. Its separateprogram checks37116explicitrows
in9models,501nonemptycorrelations,105coefficients,148unsignedformulas,
15444densitychecks and630Bhthresholds. Largerindependentfixturesuse
nonprimeintegerp=k^(r+1), clearlylabeled; root17largefixturesuseprimes.

Classical: defineU_C(t)=sum_(triplesQ)K(Q union{t}),
L_out=sum_(toutsideC)U_C(t)^2. Forn>=6,n⁴<=p,
L_out<=20pT6+p²n³ and20pT6<=L_out+3p²n³.
Fulltranslateidentity:sum_tU²=p sum_x e3²-T3²;
sum_x e3²=20T6+R,
R=6(n-4)T4-(n-2)(n-3)binom(n,2)+pbinom(n,3)-sum_C e2(c)²,
|R|<=pn³. OnC energy<=p²n³; degree3square<=p²n³/9.
DisjointtripleGram=pK(QunionR)-K(Q)K(R) returnstheoriginalsixroot
sum. No newupperboundorlowerarithmeticdifficultyclaimed.
Verifier216cases,27358cubicpointwisechecks,2931degree3+2762degree4+
920degree6directsums,3041starrows,74fullmatrixchecks,1022repeatedrows,
1600pairwiseGramidentities and5actualSidonslicesn6..10.
Separatecrossreviewfoundnocorrection; convenienttransferthresholdn25
andneedunboundedfixedrsequenceforfullPaleyexplicitlyrecorded.

Subgroup: for symmetricSnonzero,p>3,n=|S|,
E3=T6+(15n-60)(E2-T4)+60U-30J+R6,
T4=3n²-3n,T6=15n³-45n²+40n,
U=sum_u(r2(2u)-1),J=#{u:3u inS};R6countsopposite-freezero6tuples.
ForH, U=n(r2(2)-1),J=n1_(3inH). Positiveextensionbound
T6+R6<=E3<=T6+15n(E2-T4)+R6.
IntrinsicE2thereforeE3=T6+R6. ConditionalF4<=An²L,R6<=Dn³L
withp<=n⁴givesM<=[15+(15A+D)L]^(1/9)n^(8/9). Stillunproved.
All14order64sigmaexceptionshaveF4=1536=24n;13extensioncounts1382400,
one1397760 atp11127041(K5), allR6<2n³; thesearefinitefacts.
Symmetricnon-subgroupsbuiltfromshiftedSidonBhaveE2=T4 and
n⁴/20<=E3<=(5/2)n⁴. Q chosenfromunboundedoddprimesanddyadick;
quarticp1modnexistsfrompass4source. Multiplicativeclosureabsent.
Verifier53generalsets,44subgroups,34multiplicitychecks,5constructions,
86922admissibleweightedtriplepairs. Rootandcrossreviewcheckedalgebra;
rootfreshlycheckedprimaryThorner-Zaman3.1/3.2andexistingarchivehash.

Spectral: forp1mod4,C={chi(x)=chi(x-1)=1},m=(p-5)/4,
R_C=(J(eta,chi)^2+conj(J)^2-6p+36)/16, |J|²=p.
Quadraticrootmapz=t+A+B/tprovestheJacobi-squareidentityelementarily.
Bothanchorcorrectionsretained. UniformvectorQRtotalmass
q_C=(3p+5)/(8p)+R_C/(2m sqrtp)=3/8+O(p^-1/2),
1/4<=q_C<1/2 atp>=73. Alltwoanchoredgesbyresidueaffinedilation.
Onevectoronly; otherdirectionsandcouplinguncontrolled.
Residualvariance=(1/(16m))sum_C(L(x)-6)^2-(R_C/m)^2,
L(x)=sum_ychi(y(y-1)(y-x)). Noallvectorgap/depthimprovement.
Verifier79primeneighborhoods,1742edges,1418quadratictraceidentities,
31822quadraticmapfibres,9057rowidentities,72constantgapcertificates.

SpectralcrossreviewfoundanomittedclassicalGaussSIGNinput in the
frequencylabel, notafalseformula. DLMF20.11.1-2at(m,n)=(2,p),p1mod4,
givesG(2,p)/sqrtp=e^(-pi i/4)(1+i)/sqrt2=1. Rootsquaresconvert
this topositivequadraticGauss; henceS/√p labelsQRwithplussign.
RootverifiedofficialHTML, amendedonlynormalization/source-scopeprose,
andreranunchangedverifierexit0withsamecounts. Theoriginalreviewednote
andresultbytesarepreservedinresults/parallel22_spectral_reviewed_snapshot_2026_09_05.json.
Reviewerhashesfornote/resultresolvetothatsnapshot;otherinputscurrent.
Thefinalcombinedauditmustnotpretendreviewerreadtheamendedbytesor
thatfiniteFourier-surdchecksindependentlyprovethisGaussphase.

Sources:26oldmanifestentriespreserved+1NISTmetadata-onlyentry=27.
DirectDLMFHTMLarchive403;liveprimaryreviewworked.Noarchivehashclaimed.
/tmp/paley22-manifest-before.jsonpreservesoriginal26entrybytes.
Thorner-ZamanprimaryHTMLfreshlyreviewed;existingarchiveunchanged.
NoPDFnewlyopened/rendered/reviewed. No newexternaltheoreminrootcube.

Allthreeoriginalresearchlanesandallfourcrossreviewtaskscompleted.
Oneattemptofspectralcrossreviewfailedatmodelcapacity;anotheravailable
workercompletedthatreview. Nomodeloverride,reset,purchase,external
message,proofsubmission,automationormemorywrite. Do notcallworkers
stillrunningwithoutalivelisting. Noindependenthumanorformalreview.

Nextactualinputsremain:positiveoffCquarticenergy,classicalunbounded
momentsequence;R6(H)upperestimateusingmultiplicativeclosure;
nonconstantFourierdirectionsandcouplingorlongercomplete-spectraldepth.
Fullclassical,subgroup,spectral,cliqueandprizebridgegoalsunproved.
Goalactiveandunachieved,notblocked. Do notsubstitutethispass's
partialtheorems,finitechecksormodelobstructionsforthefullobjective.

## Pass 23 checkpoint, 2026-09-05

Continue the original full goal; do not mark it complete or blocked.
The previous completed pass was progress. Pass23 has new actual limited
estimates and a source-dependent spectral calculation, not a full proof.

Latest assessment: research/parallel23-pass-summary-2026-09-05.md.
Before central edits, 29 pass22 artifact/central hashes were verified and
pinned in results/parallel23_prior_state_2026_09_05.json, including the old
source-manifest text. All prior proof and verification artifacts are to
remain unchanged; only central documentation and the source manifest change.

Classical: for n=floor(p^1/4), the proportion of Sidon sets failing
M6<=15pn^3 is <=(216/5+o(1))n^2/p. Every good set has actual
L_out<=(2/3)p^2n^3. Proof is a complete-transform variance argument,
balanced-population moment domination and slice Poincare induction.
Variance bound V=n(sqrt(pA)+15sqrt(B))^2, mean gap G=p(30n^2-16n);
conditional failure <=min(1,(8/7)V/G^2), n>=6,n^4<=p.
Exceptional sets remain uncontrolled; no SS-B* universal transfer.
Verifier:3551 exhaustive sets,15042 insertion,1780 transform,38 slice,
19 variance,287 sampled Sidon,5 full energy,5 large analytic certificates.
Review found a report-only out-of-range conditional label. Root amended
the helper to emit applicability plus null outside theorem range, reran
successfully with unchanged counts; reviewer reconciled current hashes.
Reviewer reversed only the reporting changes in memory and recovered the
old script/result hashes exactly, proving other fields and certificates
unchanged. Current classical script hash6380ad545ff05d61ae973dfdbfa72bf8d636a27e3038f47a9d0593245f81a89f,
result23786c66eb4bfdafda63bbc8b1195b8aca39a89ad8a01c4eb97088e7e462c369.

Subgroup: B_rho=n sum_z w(z)w(z/rho), where w(z) counts
(a-1)(d-1)=z+1. B1-D3=nX exactly, D3=6n^3-9n^2+4n.
Balanced opposite-free part <=10nX << n^3(1+log n), n^2<p,
by Shkredov1504.04522v1 Theorem6, freshly checked primary HTML.
Internal opposite pairs mean fixed balanced count <=nX, not equality.
The unbalanced remainder remains open. At p1073748737,n256,X0,
R6=368640 is all unbalanced; explicit witness in results. Not a
counterexample to Gaussian-order R6. Circular bounds stay fourth power.
Exact coset r3(x)=n1_(j0)+(C^2)j0 and E3 ledger include r3(0)=nκ0.
Verifier9groups,121355 normalized weighted terms,560 ratio identities.

Spectral: actual two-anchor C, m=(p-5)/4, B=P_C, u=1/sqrt(m).
T=sum_D L²=p²-2p-3, W=sum_(s!=0,+/-1)L(s²)², Z=sum_C L².
Symmetries give 4Z=3W-2T. Grove Thm1.3 at moment2 gives W/p²→1.
Optional exact Hecke trace: W=p(p-3)-4-a8; Z=(p²-5p-6-3a8)/4.
Mean L_C=(U+6)/(4m)=O(1) from prior Jacobi calculation.
b²=||(I-uuT)Bu||²=(Z/m-meanL²)/(64p)→1/64.
Deleted cross block has norm b→1/8; uTB²u→5/32. No all-vector bound.
Verifier79fields,36465symmetry/fibre,9057rows,2369block entries.
It infers the Hecke label from W, not an independent modular computation.
Grove primary publisher HTML archived; no PDF. Normalization reviewed.

Prize: monic received degree k+a lists equal first-a coefficient/power-sum
fibres. A valid multivariate bound uses polynomial phases, not just linear
periods. If R^n<p^a, pigeonhole gives nonzero zero-constant phase degree<=a
with real sum>n(1-20/R²). At pinned p2130706433,n262144,a26215,R8,
this exceeds180224. Existence only; no production polynomial/list/MCA claim.
Exact official certificate at δ=j/n,1<=j<=131071 requires distinct
16-row remainder explanations plus 8-row nontrivial ratio-scalar maxima
<=274980728111395087. Both maxima open. Quartic-field mismatch and
trace-zero extension obstruction still apply. No current website claim.
Verifier841Newton,10monic lists,182counts,2F25centres,24general lists,
72MCA identities,11old source hashes. Prize review also checks576 cases
of larger-witness shrinking independently.

Manifest27 prior entries + Grove +4 same-commit ArkLib files=32.
Root retrieved the actual Lambda/IsMCA/mcaError definitions, source ledger
results/parallel23_prize_additional_sources_2026_09_05.json. No model
override, reset/purchase, external message, submission, automation or
memory write. One subgroup turn failed at transient model capacity and
resumed on the same model. Mathematical difficulty is not a blocker.

All four research results and all four separate-agent reviews completed.
No mathematical correction was required; classical report labels corrected.
Final audit: results/parallel23_pass_audit_2026_09_05.json, script
experiments/parallel23_final_audit_2026_09_05.py. It verifies23 result input
hashes,32 review input hashes,24 preserved prior proof/verification files,
5 Python syntax checks and local links. Manifest27 old+5 new=32.
Latest live listing has all3child agents completed, none running.
Use collaboration.list_agents for future live state; don't infer activity
from this checkpoint. The full original objective remains active/unachieved.

## Twenty-fourth parallel pass, 2026-09-05

Previous turn classified progress after verifying its authoritative
artifact hashes. New prior-state snapshot:
results/parallel24_prior_state_2026_09_05.json. Full original goal remains
active and unachieved. See research/parallel24-pass-summary-2026-09-05.md.

Four lanes ran: root unsigned inversion moments; three workers on local
exceptions, actual next spectral coefficients, and unbalanced subgroup
relations. Each has a reviewer distinct from the author. Root reviewed
the classical lane with a separate actual-row increment verifier.

New bounded result: every size-k input C has some B_h inverse D_z with
M_(2r)(D_z)<=U_r/(p−Q_h(k)) under explicit p>h,p>Q_h(k). At fixed
r=3,h=2,k=floor(p^1/4), this is (20+o(1))pk^3. Standard squared Weil
sums plus prior good-pole counting; no novelty claim. Original moment
transports to prescribed chi(d) signs, so SS-B* remains unproved.
Root verifier: 573 sets, 1,149 good-representative bounds, and 27 tuple
expansions; all signed/unsigned boundary identities and 573 integer
overflow guards pass.
Separate reviewer checks both prime congruences and a noninteger pole
fraction threshold with plain Python, without importing lane functions.

Local exceptions: exact sixth swap coefficient
1−6(p−5)/(k(p−k)); lower coefficients nonnegative on thin slice,
nonnegative iteration only k>=7. Local variance retains M10/M8 of
the actual deletion. Shell mechanism requires >k^(5/6) retained labels;
its measure is below the allowed global exception fraction. Actual sign
split loses8 at balance, termwise bound returns the old leading power.
No uniform exceptional-input control. Root separately checked 858
increment moments, 143 rows, and 3 fields at sizes outside the author's
identity coverage.

Spectral: exact S_C L=[p+6+k3(t)+k3(1−t)+k3(1−1/t)]/4.
Pinned Katz hypotheses give next diagonal a=1/2+O(p^-1/2) and next
coupling squared c²=3/64+O(p^-1/2). K2=span{u,Bu} has Rayleigh upper
(7+sqrt5)/16+O(p^-1/2), and retained-image squared norm
(15+sqrt74)/64+O(p^-1/2). Block deletion error tends sqrt3/8.
All-vectors and growing-depth input still missing. Explicit next
residual needs h^T S_C h=o(p^(7/2)), also unproved. The verifier covered
88 primes.

Subgroup: finite quartic classes n=4,8 have R6=0; n=16 has R6=480 only
at p=33713,37201,41521 and zero elsewhere. All 480 exceptional words
have pattern (4,1,1) and are unbalanced at all 10 splits. All 690
eligible prime fields and 1,770 norm determinants were verified;
E3<=15n³ in these finite classes. Exact intrinsic fibers:
I1=6n³−9n²+4n, I_nonidentity_square=18n(n−2), I_nonsquare=0.
Orbit masses equal actual r3, so there is no new general aggregate
bound. The dyadic norm bound is exponential in n and does not extend
automatically to growing n in the quartic window. The general subgroup
square-root target remains open.

No new prize lane theorem or contract update. Pass23 exact arbitrary-word
remainder-fiber and scalar-ratio maxima remain unbounded, with the same
pinned official certificate and spot-check obligation. No new worst-case
cancellation exponent, human referee claim, or formal theorem.

Audit: experiments/parallel24_final_audit_2026_09_05.py writes
results/parallel24_pass_audit_2026_09_05.json. It pins all current
proofs/results/reviews and dependencies, preserves prior proof/verification
files, and checks syntax/local links. All 32 manifest entries are unchanged.
One classical worker turn hit a selected-model capacity error and resumed
with the same model; no account-wide limit claim. No model override,
reset/purchase, external message/submission, automation, or memory write.
Use collaboration.list_agents for live activity; archived notes are not
evidence that workers remain running. Mathematical difficulty is not a
blocking condition, and the full original goal is still active.

## Twenty-fifth parallel pass, 2026-09-05

Previous pass classified progress; its authoritative artifacts were
verified before new work. Snapshot:
results/parallel25_prior_state_2026_09_05.json pins 27 files and the
32-entry source manifest. Full original goal remains active/unachieved.
Assessment: research/parallel25-pass-summary-2026-09-05.md.

Root prize projection: for any linear code over F_q, B_r=B_1 exactly
for MCA bad-scalar counts at the same radius and field. Nonzero random
row projections give survival probability >1−1/q; integer rounding
retains the B_1=q−1 endpoint. If binom(L_1+1,2)<q, L_r=L_1;
if L_1<=L<q, L_r<=floor[L(q−1)/(q−L)]. For the pinned official
q=2130706433^6 and R=274980728111395087, the integer list condition
holds, giving B8+L16<=R iff B1+L1<=R. Scalar maxima remain unbounded,
the field has not been reduced, and the spot check remains necessary.
Extension-to-base-field lists transfer by coefficient expansion; MCA
bad sets only inherit a base-field affine-line intersection bound.
No Paley-to-grand-prize equivalence or current prize-status claim.
Verifier passes 401427 exhaustive tiny-code word pairs, 25252 nonzero
list projections, extension field and affine-line checks, and the
exact official budget. The binary n=3 example is a linear repetition
code, explicitly not an RS code with three distinct F2 points.
Separate-author review of this root-authored note remains outstanding.

Worker signed inversion: exact half-moment, completed mixed-average,
and weighted fractional-linear identities. The asymmetry is negligible
on the thin slice but the common moment remains uncontrolled. The
joint good-pole upper bound retains original M6. Augmented signed M6
is preserved by repeated fractional-linear changes; no independent
sign randomization follows. Author verifier passes 666 input sets.

Worker subgroup: repeated opposite-free zero-sum six-words obey
R6rep<=15n[r4(2)−6n+8]<=15nE2, aggregated over all product ratios.
MRSS gives O(n^(69/20)(1+log n)^(1/5)) in the quartic window.
Together with prior balanced/opposite-containing estimates,
E3=T6+D6+O(n^(69/20)(1+log n)^(1/5)+n³(1+log n)), with
nonnegative error. D6 is six distinct, opposite-free, fully unbalanced
at all ten splits and remains unbounded; the error still exceeds cubic
scale. Norm/rank argument is uniform in dyadic n. Actual n=256 witness
has norm 4419283227438373820201815271409083396, valuation one at
p=1073748737 despite 128 rationally independent cyclic shifts.
Author verifier passes four quartic cases at n=64,128,256.

Worker spectral: full S3 block decomposition, orbit dimensions,
S_C²=pI/4+3A/4−3J/2, and exact off-block leakage retaining the
rank-two border. Uniform-vector polynomial iterates stay in dimension
about m/6 at any depth. Nontrivial sectors need lambda_max(A)<=
(2/3+o(1))p; elementary bound is p. Trivial-sector border persists,
and existing logarithmic depth allowance is insufficient. Author
saved proof/script before failure; root completed its 24-prime run.

All three workers ended with usage-limit errors during this pass.
The crossreviews assigned to workers did not save completed review
notes. Root continued locally and separately reviewed all three worker
deductions. Its independent verifier imports no worker functions and
checks actual sign halves and fractional-linear weights, subgroup
incidences and partition, a 128-dimensional integer Bareiss norm, and
full spectral/leakage coefficients at p=13,17,61,269. It passes.
Evidence: research/parallel25-root-review-2026-09-05.md and
results/parallel25_root_review_2026_09_05.json. Agent review is not
human refereeing or formal certification. Do not claim that the root
prize note has a distinct-author review, or that workers are running.
Use collaboration.list_agents for future live state.

Root used primary HTML only this pass: MRSS Theorem 3 energy convention
and subgroup hypothesis, Shkredov Theorem 6 with shifts −1, and the
GGR abstract for provenance only. Manifest preserves all 32 entries
and adds two HTML archives (34 total). Audit:
experiments/parallel25_final_audit_2026_09_05.py and
results/parallel25_pass_audit_2026_09_05.json. It checks hashes,
prior proof preservation, source ledgers, syntax, and local links.
No model override, reset/purchase, external message/submission,
automation, or memory write. Mathematical difficulty is not a blocker.
The full original goal remains active and no worst-case cancellation
exponent or complete proof is established.

## User-requested Prove2Me setup and local formal projection core, 2026-09-05

The user supplied a replacement API key and explicitly requested project
setup with Prove2Me, then continued research. No credential appears in
this checkpoint or project artifacts. The key was entered through an
echo-disabled prompt, authenticated against the official refresh endpoint,
and saved in the existing shared gitignored mode-600 authentication store
outside this project. Access tokens are kept only in process memory and
requests refuse redirects. Do not print or copy key/token values.

Project instructions: PROVE2ME.md. Connector scripts:
scripts/prove2me.py and scripts/prove2me_projection.py. Reused shared
workspace /Users/shawwalters/prove2me_workspace, Lean4.30.0, Mathlib
c5ea00351c28e24afc9f0f84379aa41082b1188f. Live supported environment
confirmed. Service changed0.9.6→0.9.7 during setup; version guard stopped
before publication. Root refreshed the official skill and setup/prove/
contribute/mission_solver references, inspected their differences, and
updated only those shared docs. Dirty Lean work and dependencies preserved.

Authenticated readback: four prior Paley elementary proofs ACCEPTED;
paley_two_set_conjecture still Open; proposal
b5e7121a-eb3c-48f1-a495-29fffd24e71d remains private Draft with8items.
No new mission, launch, public release, community message, vote, or rating.
Readback: results/prove2me_connection_2026_09_05.json.

New local formal result: finite incidence counting with |P|+1=qt,
q>=2,t>=1,K<q, every S-degree>=(q−1)t and P-degree<=K implies
|S|<=K, retainingK=q−1. Local proof and exact target statement compile;
proof axioms only propext,Classical.choice,Quot.sound. Private statement
submitted once under name paley_projection_survival_count. Publish job
843b8456-f9ad-4158-809d-8378028cc958 was still PENDING after15minutes.
No theorem ID or proof-submission ID existed at that point. Do not
report a hosted proof verdict. Latest authoritative queue state is
results/prove2me_projection_server_2026_09_05.json.

Resume with python3 scripts/prove2me_projection.py poll; once PUBLISHED,
run its verify action once, then poll until terminal. Scripts reuse saved
IDs and hash-pin the payload/proof. Do not duplicate the queued job.
The local private target carries an intentional sorry placeholder;
the solution is complete and imports no target or open child.

While waiting, root formally checked the next algebraic core locally:
(1) cardinalV=cardinalker(f)*cardinalF for a nonzero functional;
(2) exact positive kernel/survival counts for finiteF,V;
(3) for nonzero functionals indexed byS, if every nonzero vector
detects<=K<|F| of them then |S|<=K;
(4) any nonempty family of<=|F| nonzero functionals has one nonzero
vector on which all are nonzero. All have only the three standard
axioms, including copied dependencies. Sources/copies and results:
prove2me/Check_paley_functional_kernel_card.lean,
prove2me/Check_paley_nonzero_functional_count.lean,
results/prove2me_functional_counts_local_2026_09_05.json,
results/prove2me_projection_core_local_2026_09_05.json.
These two check files have not been submitted to the server.

research/prove2me-projection-core-2026-09-05.md explains the simpler
written MCA proof: there are at mostqbad scalars, so their quotient
witness functionals admit one common projection preserving all of them.
It also permits binom(M,2)<=q in the list-distinctness step. The earlier
strict bound remains correct and official budget unchanged. Next actual
formalization input is the quotient witness construction and survival
for the actual code definitions, followed by comparison of maxima.
Do not call those code-specific steps Lean-proved. No new cancellation
bound, scalar prize maximum bound, human review, or full proof follows.

experiments/prove2me_setup_audit_2026_09_05.py verifies setup, exact
proof/statement bytes, standard axiom records, copied local core files,
credential permission metadata without reading the secret, and local
links. results/prove2me_setup_audit_2026_09_05.json records queue status
separately from local verification. Pass25's research audit was rerun
after central-document updates. Three parallel workers remain errored
with usage limits; Prove2Me setup does not reset their limits. Root has
continued locally. Full original goal remains active and unachieved.


## 2026-09-05 pass 26: code projection compiled; attribution corrected

Previous turn classification: progress (setup and local algebraic proofs).
This turn: progress in formal verification and source accuracy; no new
cancellation estimate. Goal remains active and unachieved.

Standalone prove2me/Check_paley_mca_projection.lean compiles in pinned
Lean4.30.0/Mathlibc5ea003. Seven principal statements have only standard
axioms. It proves simultaneous bad-scalar preservation by one nonzero
row projection, input-failure equivalences, and exact equivalence of
uniform count bounds, specialized to codeword agreement on any allowed
family of coordinate sets. A source-level comparison matches official
MCA witnesses; there is no direct ArkLib import theorem or full build.

Crucial attribution: the previously archived ArkLib Errors.lean already
contains mcaError_interleaved_eq at commit e65197892890b8fd9b0dc05b8980273cf1d595cc.
Root inspected its full proof and the pinned primary file. Pass25 MCA
invariance is an independent rediscovery, not a new mathematical advance.
Historical proof bytes remain unchanged; central assessment is corrected.
The pass25 list transfer still lacks separate-author review and Lean
verification. The coefficient-pigeonhole threshold idea was also already
proved in research/list-to-winning-set.md; no new obstruction is claimed.

Private Prove2Me job843b8456-f9ad-4158-809d-8378028cc958 remains PENDING;
no theoremID, proofsubmissionID, or serververdict. Do not duplicate it.
Use scripts/prove2me_projection.py poll, then verify once it is PUBLISHED.
Three workers remain usage-limit errors. No model/reset/credits change.
Root can continue locally; goal is not blocked. No public release,
mission launch, automation, message to others, or memory write occurred.

Pass26 summary, detailed formalization, compiler log/result, source ledger,
and final audit contain the current evidence. Two pinned code-source files
are added; the 34 existing manifest entries are preserved. Remaining
work targets actual uniform scalar MCA/list bounds, signed moments,
growing-order subgroup remainders, and the full spectral operator.


## 2026-09-05 pass27: exact access to the six-distinct remainder

Previous goal turn: progress in formalization and source attribution.
Current goal turn: proved unique balanced split for distinct opposite-free
six-words, exact rational unmarking of repetitions, and a sparse D6 count.
The full goal remains active and unproved; no uniform D6 bound follows.

Two balanced splits force c=d or c=-d by subtracting c times one balance
equation and d times the other. Three statements in
prove2me/Check_paley_balanced_split_collision.lean pass Lean with standard
axioms. The arbitrary-partition relabeling and full combinatorial count
are ordinary proofs. Distinct-balanced union is exactly10B*, whereB*
is the fixed-split distinct opposite-free count. Repeated count is
15n times sum1/nu(1,1,a,b,c,d) over opposite-free four-sum=-2 candidates.
E3=n²w0²+n sumcosetW². SubtractT6,J6,repeated,10B* to getD6 exactly.
Coset label is z^n, notz^((p-1)/n); initial implementation error was
caught by independent checks and corrected before saving final results.

Candidate cost O(n²logn+r4(2)+E×(H-1)) yields O(n^(49/20)polylogn)
field/comparison operations under existing primary energy bounds, withH
supplied. StorageO(n²). This is a computational improvement only.
Every distinct opposite-free six-set has freeH-scaling orbit, so its
ordered-word count is divisible by720n. This checks bothD6 andbalanced.

Sparse results cover14cases. Six new cases atN512,1024 choose the first
eligible prime above each ofN^4/4,N^4/2,3N^4/4. All six haveE3=T6+D6.
D6_scaling_orbits atN512 are5,2,0 and atN1024 are8,3,2. No uniform
conclusion or worst-prime assertion follows. A separately implemented
C++ normalized-six-word enumeration matches every category in all14
fields, including the six new large cases. Both programs are root-authored,
so this is independent implementation rather than separate-author review.

Results and proof: research/parallel27-six-distinct-2026-09-05.md,
research/parallel27-pass-summary-2026-09-05.md, andresults/parallel27*.
All36source entries and prior proof artifacts preserved. No new hosted
submission; private incidence job843b8456-f9ad-4158-809d-8378028cc958
was polled and remainsPENDING. Do not duplicate it. No public release,
mission launch, reset, credits change, or memory write occurred.

Final pass27 coverage extends14 to17cases with dense subgroups(p,n)=(5,4),(17,8),(97,16). These test the exact formula beyond its quartic cost application. Final sparse and independent direct records cover all17.

Prove2Me update in pass27: the existing private statement job is PUBLISHED, theoremID d025929b-75ed-4c1c-b661-8e491017b153. The locally checked incidence proof was submitted once as27ff60e7-7d6d-48da-8c34-7b123c6c0a3f. Current hosted proof status: ACCEPTED, with theorem readback Proved in the pinned environment. This is terminal; no further polling or resubmission is required for this lemma.


## Pass 28: mixed moments and complete finite D6 census

Previous goal turn: progress. Full Paley conjecture and prize unproved;
goal remains active. See research/parallel28-pass-summary-2026-09-05.md.

The ordinary proof in research/parallel28-moment-recurrence-2026-09-05.md
combines MRSS E3 << n^4 log n, Shkredov invariant-set energy (published
Lemma 9, equation 23), and Konyagin (3,6) to derive
M << n^(71/72) (log n)^(1/18), uniformly in n^4/4<=p<=n^4, n>=4.
A nonnegative weighted lifting gives a log-free energy recurrence with
constant 82952 C_*. The origin is explicitly removed and then bounded.
The published recurrence (28) independently gives the same power with
log exponent 1/6. No literature novelty/current-best claim, independent
author review, Lean proof of this analytic theorem, or new hosted verdict.
The prior project power 2849/2880 improves by 1/320. The recorded
recurrence/interpolation/mixed-moment envelope supplies at most 1/72
saving; square-root cancellation and all full targets remain open.

The complete census covers 28,774 quartic-window prime/order pairs at
n=4,8,16,32,64. D6 maxima are 0,0,0,23040,184320, respectively.
All satisfy D6<=n^3. 28,753 overlapping E2/E3 records match the old sigma
census. Separate C++ direct counting matches 75 fields, covering every
observed category-count tuple and all positive-D6 maximizers.
The three workers last reported terminal usage-limit errors. Root
completed this pass without resets, model changes or credit purchases.
Prove2Me still has five accepted elementary proofs; the fifth is terminal
ACCEPTED and must not be resubmitted. No public release, mission launch,
community message or memory write occurred.


## Pass 29: exact centering and closure of the moment estimates

Previous goal turn: progress through live Lean performance diagnosis.
The single-file check finished; other proximity builds were still active
under heavy contention. No process was stopped or restarted by this pass.
No new Lean run was launched. Existing three agents were inspected live
and still report terminal usage-limit errors. Root can continue locally.

See research/parallel29-centered-recurrence-2026-09-05.md and
research/parallel29-pass-summary-2026-09-05.md. For T_s=E_s-n^(2s)/p,
the ordinary proof gives T_(2s)<=D n^(2s-1/2)T_s with D absolute and
independent of s. It combines published point--plane Theorem8 with an
elementary exact variance identity, then handles signed weights and
the origin explicitly. No unqualified centered-function theorem is used.
No independent-author review or Lean verification of this proof is claimed.

The exact-centered recurrence removes the old principal-term deficit cap.
The mixed gate still supplies best power71/72 from the available seeds.
Feeding that amplitude bound back into moments gives a closed envelope
equal to the earlier deficit through order24 and G(s)=s/36+17/6 thereafter.
It is closed under interpolation, recurrence, amplitude feedback, and
the exact coset-product inequality T_(s+t)<=(p/n)T_sT_t. All finite and
infinite mixed-gate endpoints give saving at most1/72. This is a limit
of these inequalities, not of all methods or of the actual period.

The one-process verifier checks incidence encoding, plane variance,
centered convolutions, signed weights, origin absorption, exact exponent
arithmetic, and separately labeled approximate Fourier identities.
Prime-subfield examples rule out naively replacing p by q in the invariant
set estimate. The full conjectures and scalar prize targets remain open.
The main source manifest and prior proof artifacts are preserved. The
accepted Prove2Me incidence theorem is terminal; no polling or duplicate
submission occurred. Goal remains active and unachieved.


## Pass 30: a triangle-generated remainder counterexample

Previous goal turn: progress. Current pass found a counterexample to the
literal D6<=n^3 extrapolation beyond the complete finite census. The full
goal remains active and unproved. No cancellation exponent improves.

At p=215535361, n=128, g=25525303, g^64=-1 and 1+g+g^19=0 modulo p.
Primality and the quartic window are checked exactly. The zero triangle
alone proves D6>=360n(n-20)=4976640>n^3. Sparse and separate direct
counts give D6=10967040=(5355/1024)n^3. Its triangular and primitive
parts are4976640 and5990400, or54 and65 scaling orbits of unordered
six-sets. Both parts separately exceed n^3. This refutes neither an
unspecified O(n^3) constant nor either Paley target or the prize.

For a general dyadic subgroup, kappa=|H intersect (1-H)| and
tau=kappa-3*1_(2 in H) count distinct zero triangles through tau/6
free scaling orbits. Every distinct zero-sum six-set has at most one
zero-triple partition. Pairing triangle orbits proves
10n max(0,n-28)tau² <= D6_tri <=10n²tau².
Thus uniform D6=O(n³) would require kappa=O(sqrt(n)). The primitive
remainder has no proper nonempty zero-sum subset and remains unbounded.

For one triangle orbit A, bad relative scalings are exactly
R union (-R) union (-R²), R=A/A, of size at most20. The explicit
example attains20 bad and108 good relative scalings, giving the short
lower certificate independently of full enumeration. The norm of
1+X+X^19 modulo X^64+1 is exactly p, verified by an independent integer
determinant. Its principal ideal equals the evaluation kernel by index,
so all relation polynomials in this example are integral multiples in
the quotient ring. Independent prime constraints cannot be assumed.

See research/parallel30-triangle-remainder-2026-09-05.md and the pass30
summary. The verifier records13201 exact checks and compares all
triangle-pair counts against direct six-word enumeration in five fields.
The norm search is discovery only, not a completeness claim. New proof
arguments have author review and exact finite checks; separate-author
review and Lean verification remain outstanding. No new Lean process
or hosted proof job was launched. Existing three agents were inspected
and remain in terminal usage-limit error states; root can continue.
Prior artifacts, the main source manifest, and accepted Prove2Me record
are preserved. No public release or memory write occurred.


## Pass 31: exact minimum triangle derivations

Previous goal turn: progress through an arithmetic counterexample and
structural triangle bounds. This pass classifies all 119 remainder
orbits of that fixed subgroup by their minimum triangle derivation
length. It yields no uniform upper bound or new cancellation exponent.

In R=Z[X]/(X^64+1), f=1+X+X^19 has determinant p=215535361 and
fR equals the evaluation kernel at g=25525303. There is one orbit of
distinct zero triangles. The unique integral quotient q=F_S/f has
minimum triangle term count ||q||_1, allowing all scales, signs and
repetitions. This is proved by injectivity and the exact absolute-sum
lower bound, not inferred from failed searches for shorter derivations.

The 54 triangular orbits have minimum 2. Among the 65 primitive orbits,
20 have minimum 4; 1 has minimum 38; 9 have 40; 17 have 42; and 18 have
44. Nine of the last group have one coefficient of magnitude 2, so
their squared norm is 46; the other quotients have coefficients 0,+1,-1.
The primitive exponent set {0,8,16,45,62,100} has its 38-term quotient
written explicitly in research/parallel31-short-multiples-2026-09-05.md.

Every minimum cancellation graph for a primitive six-set is connected:
a disconnected component would either have zero formal boundary,
contradicting injectivity and minimality, or give a proper zero-sum
subset. With L triangle terms it has (3L-6)/2 edges and cycle rank
L/2-2. Thus 45 actual primitive orbits require ranks 17 through 20,
while the 20 four-term quotients have tree derivations. These are not
being identified with the separate spectral necklace graphs.

The classifier recovers all 714 normalized six-sets and 119 scaling
orbits, constructs an integer adjugate fB=p, and verifies every integral
quotient and an explicit cancellation graph. A separate checker uses
trinomial row equations, checks every proper subset and product
partition, and compares with the prior independent C++ count. The two
records contain 1148 and 1026 checks. No new Lean process or hosted
proof job was launched; independent-author review is still outstanding.

The initial higher-energy source check distinguished repeated-difference
energies from sum moments. Direct coset norm splitting only gives
E_3(Q)<<|Q|^6 n^-2 log n for Q contained in F_p^*. Cauchy plus the
centered doubling recurrence gives deficit update 41/40+d_s/2, already
inside the previous envelope. No exhaustive impossibility claim is made
for other higher-energy arguments. The next input must count short
outputs of the larger cancellation networks or control the spectral
quantity by other means. All full targets remain unproved; goal active.
