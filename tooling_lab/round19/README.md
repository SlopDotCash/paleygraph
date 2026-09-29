# Round19: an obstruction to short bounds, and a general membership codec

This round resolves the immediate generalization gap left by [round18](../round18/README.md). It provides an exact membership algorithm for general short-relation encodings, including the p2013265921,N64 input that the generator-2 adapter could not handle. It also proves why simply adding more short-relation digit bounds cannot make that input's membership test exact.

The [derivation](DERIVATION.md) separates the ordinary mathematical arguments from the finite experiments. The [source audit](prior_art.md) identifies existing symbolic dynamics, separating inequalities and algebraic coset encoding. Historical novelty is unestablished; neither prize is being proved.

## A whole family of small bounds misses a one-unit error

Take the actual scalar a=(p+1)/2. Its centered vector has first coefficient-(p-1)/2. Add p to that coefficient, moving it just one unit outside the centered range while preserving every scalar congruence. Let delta be the smallest distance of any other coefficient's magnitude from the centered boundary.

The [certificate](relation_barrier.json) proves that this escaped vector passes every coefficient bound arising from every nonzero kernel relation f with ||f||_1<=delta+1. This includes arbitrarily many such relations. Its encodings remain integral. A separate symbolic dual-identity check verifies the general inequality argument; this is stronger than failure on a sampled relation list.

| p | N | All relations up to this L1 budget fail to separate | Lower bound for any separating relation's L1 | Explicit separator's L1 |
|---:|---:|---:|---:|---:|
|17|4|5|6|10|
|257|8|65|66|130|
|65537|16|16385|16386|32770|
|6700417|32|26174|26175|52348|
|2013265921|64|11607941|11607942|23215882|

The last column is constructed from the subgroup's least nonidentity odd residue and checked by exact polynomial multiplication. This brackets the minimum separator cost within a factor two; it does not determine the exact minimum. It is specific to these L1-derived digit inequalities, not every possible arithmetic constraint.

For the largest input, the current short relation has L1 norm13. A finite portfolio of16 LLL-discovered relations, with L1 norms13–20, leaves127 of128 tested single-coordinate aliases of the actual a1 vector. The universal certificate above explains a stronger failure at the chosen near-boundary scalar. The [separate review](barrier_review.json) checks all portfolio decisions, the universal slack identities, and25 exact resultants for the five actual/escaped pairs. [Separation bounds](separation_bounds.json) give the explicit binomial separators.

## A small nonlinear interface handles the general case

For E(a)=fF_a/p and u=g^-1, calibrate beta=E(1)(u) modulo p. The exact identity

```
E(a)(u)=a*beta mod p
```

follows because E(a)-aE(1) is an integral multiple of f. If beta is nonzero, a candidate D yields a=D(u)/beta modulo p. Reconstruct centered F_a and accept exactly when fF_a/p equals D.

The [codec](projection_codec.py) therefore uses a modular scalar projection followed by exact sparse re-encoding. All five saved input configurations have nonzero beta, including the general N64 relation with cofactor9985208709332560769028097. That cofactor is not needed by the fast branch. At N64, one word uses64 projection terms,64 reconstructed coefficients and576 coefficient products for the9-term relation. These are operation counts, not measured wall-clock gains.

Projection alone is insufficient: each hard escaped vector has the same projected scalar as its centered partner. Re-encoding rejects all five escaped examples. This explicitly retains the nonlinear centered section that small relation bounds discard.

When beta is zero, a fallback extracts the scalar with one adjugate coefficient row modulo k*p, then performs the same re-encoding. It preserves exact membership without a dense inverse per word or a cofactor-sized state table. Large integer weights and the cost of preparing that row remain explicit. The cases f=p, f squared and f=p squared distinguish a p-dividing cofactor from a degenerate calibration.

The [independent codec review](codec_review.json) uses full rational polynomial inverses and integer multiplication matrices. It checks:

- All484 saved positive encodings and450 probes.
- Every6561-word small cube decision, accepting exactly257 including zero.
- Five hard actual/escaped pairs.
- 204 degeneracy records, including51 intended encodings and153 mutations.
- Three malformed digit inputs and12 rejected corrupt certificates or outputs.

The [family controls](codec_family_controls.json) exhaust80 nonzero kernel relations and117200 relation-specific digit words at p17,41,97. Each result is compared with a complete direct scalar codebook. Relation coefficients range over[-2,2]^4 at p17 and41. At p97 that box contains no nonzero kernel relation, so its explicitly reported nonempty test family uses[-3,3]^4. This is additional small-parameter coverage, not a complete large-field census.

## Why raw windows were the wrong generalization

The [window hierarchy](window_hierarchy.json) exhausts the N4 and N8 ternary cubes. For N8, the surviving counts are:

| Maximum window depth | Windows only | Add odd support or zero | Add digit sum+1,-1 or zero |
|---:|---:|---:|---:|
|1|6561|3281|2033|
|2|1153|577|501|
|3|561|305|297|
|4|385|257|257|
|5|321|257|257|
|6|289|257|257|
|7|273|257|257|
|8|257|257|257|

The sharp general statement over Q=2^N+1 is: raw windows need full depth N; with odd support, depth N/2 is sufficient and necessary. Explicit witness families demonstrate sharpness through N128, and the [review](window_review.json) independently recovers every window from lifted coordinates. These large moduli may be composite. A constant sign-memory automaton from round18 retains information that growing raw windows miss; this is a specific application of established symbolic-dynamics ideas.

## What comes next

The membership problem is now inexpensive for a supplied word. General enumeration, prefix completion and aggregate norm bounds remain open tooling tasks. A finite-field projection plus re-encoding still represents p possible scalars; it does not by itself compress that search space or count its rare subsets. The next prototype should return complete families consistent with partial digit information and expose the state or branching cost before any asymptotic inference.

```sh
/opt/miniconda3/bin/python3 tooling_lab/round19/window_hierarchy.py
/opt/miniconda3/bin/python3 tooling_lab/round19/window_review.py
/opt/miniconda3/bin/python3 tooling_lab/round19/relation_barrier.py
/opt/miniconda3/bin/python3 tooling_lab/round19/barrier_review.py
/opt/miniconda3/bin/python3 tooling_lab/round19/separation_bounds.py
/opt/miniconda3/bin/python3 tooling_lab/round19/projection_codec.py
/opt/miniconda3/bin/python3 tooling_lab/round19/codec_review.py
/opt/miniconda3/bin/python3 tooling_lab/round19/codec_family_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round19/verify_round.py
```

The [manifest](manifest.json) binds sources, results and inputs and preserves the earlier rounds. Reviews are root-run separate implementations, not external peer review or Lean certification.
