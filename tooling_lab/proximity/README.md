# Proximity obstruction microscope

This lab implements two connected research instruments: a **syzygy survival
filtration along specified root motions**, and a **stack incidence ledger that
retains decoded codewords, their full agreement sets, and affine codeword tracks**.
They diagnose what a rank defect is made of and what happens when it is put into
an actual Reed–Solomon pencil. They do not attempt the prize proof.

The first iteration exposed a failure of first-order diagnostics. The second
exposed a failure of triple-only diagnostics. These are useful negative results
for tool design, accompanied by executable replacements.

## Reproduce

Python 3.10+; standard library only. From this directory:

```sh
python3 deformation_microscope.py
python3 stack_provenance.py --samples 4
python3 verify_artifacts.py
```

`--quick` runs a smaller deformation experiment; `--output FILE` chooses the JSON
output. The full deformation experiment took approximately 4 seconds and the
16-stack experiment approximately 15 seconds on this machine. Seeds are fixed.
Elapsed time is the only intentionally nondeterministic output field.

Files:

- `deformation_microscope.py`: finite-field algebra, coefficient-matrix jets,
  projected-kernel filtration, all admissible one-root replacements.
- `stack_provenance.py`: exact subset-parity census, deduplicated decoded
  codewords, affine-track hypergraph, overlap accounting.
- `verify_artifacts.py`: readback verification of syzygies, codewords,
  agreement masks, certificate multiplicities, and affine tracks.
- `iteration1_results.json`: first bounded experiment before the stack-aware
  iteration.
- `deformation_results.json`: 14 profiles up to subgroup order 96 and fields
  above 2^31, plus a characteristic sweep of the known F41 witness.
- `stack_results.json`: full certificates for 16 sampled stacks. Every stack
  is completely censused over finite affine scalars; the stacks themselves are
  sampled, not exhaustively enumerated. The projective point at infinity is
  deliberately excluded and must be added for any comparison to a projective
  budget.

## What we actually have

Read-only source audit: `/Users/shawwalters/proximityprize`, HEAD
`25de4107388ce3e371c2f56252279ab6455eb1a9`. The checkout contains untracked SW1
work dated 2026-09-05. Nothing in that checkout was modified or built here.

The starting executable witness is the `(6,6,6)` linear-cofactor syzygy on
`mu_20` over F41 in
`scripts/probes/probe_syz71_verify_f41.py`. Existing tools already compute
finite-field rank, exact bad-scalar sets by parity checks, prime sweeps, and
cyclotomic norm certificates. Reimplementing those under a new name is not an
invention.

Two newer source files materially change how to use those objects:

- `ArkLib/Data/CodingTheory/ProximityGap/Frontier/_SW1_F1_UniformSylvesterRefuted.lean`
  supplies a binomial/coset family refuting a triple-only injectivity gate. Its
  source claims a field-uniform and production-shaped construction; this lab
  independently verifies its small finite instances but did not check Lean.
- `docs/kb/deltastar-sw1-nec-2026-09-05.md` audits the distinction between a
  sufficient Fourier route and a proved necessity of that route. Its finding
  is a warning against treating the older dossier's “same BGK object” prose as
  a theorem that every possible prize approach must pass through.

Thus this lab targets a missing **joint diagnostic**, not another proposed
triple-only sufficient condition. None of its outputs determines delta-star.

## Instrument 1: projected-kernel survival

For three root sets, write `W_i(X)=prod_(a in S_i)(X-a)`. Fix the maximum
cofactor degree `q`. Construct `M` with columns `X^j W_i`, `0 <= j <= q`, in
the monomial basis. Its kernel is exactly the allowed cofactor syzygies.

Move each root on a specified straight trajectory `a(t)=a+t*v_a`, and compute

```
M(t) = M_0 + t M_1 + ...
```

using polynomial arithmetic modulo `t^(ell+1)`. No floating point or numerical
differentiation enters. Let `T_ell` be the lower triangular block matrix with
block `(i,j)=M_(i-j)` when `j<=i`. It expresses

```
M_0 c_0 = 0
M_0 c_1 + M_1 c_0 = 0
M_0 c_2 + M_1 c_1 + M_2 c_0 = 0
...
```

Define

```
sigma_ell = dim ker(T_ell) - dim ker(T_(ell-1)),   dim ker(T_-1)=0.
```

This is exactly the dimension of the possible leading syzygies `c_0` that
extend through the chosen order. **Proof:** projection from `ker(T_ell)` to
its first coefficient has a kernel consisting of tuples with `c_0=0`.
Deleting that zero identifies this kernel with `ker(T_(ell-1))`. Apply
rank-nullity. The surviving spaces are nested, so the sequence is nonincreasing.
This elementary identity is proved here for correctness, not claimed as new.

The first-order direction finder computes `L M_1(v) R`, where `R` spans
`ker M_0` and `L` spans the left kernel. Its kernel gives directions that
preserve all initial syzygies through first order.

**Crucial scope:** survival is along the specified root trajectory. Failure
at order two does **not** prove that a tangent vector cannot be integrated
along a curved trajectory with nonzero root accelerations. Also, these root
motions usually leave the multiplicative subgroup. Since `p` does not divide
`n`, the derivative of `X^n-1` is nonzero at each root: the subgroup has zero
root-motion tangent dimension. The tool reports this explicitly. Admissible
subgroup edits are separately enumerated as discrete replacements.

## Experiment and first revision

On the known F41 witness, `M` is 8 by 6 with rank 5. The ambient root
parameter space has dimension 18; the first-order preserving space has
dimension 15. But all 24 tested basis/random tangent directions had

```
sigma = [1, 1, 0, 0].
```

The translation and scaling controls both have `[1,1,1,1]`, as required by
their exact change-of-variable symmetry. A seeded random partition has no
linear syzygy. All 36 single-root replacements by an unused point of `mu_20`
destroy the F41 syzygy.

This invalidates “large tangent space means stable useful obstruction” as a
tool heuristic in this experiment. It motivated retaining the entire
projected-kernel filtration and adding explicit domain admissibility.

The F41 support *index pattern* has kernel dimensions `[1,0,0,0,0,0]` at
`p=[41,61,101,1021,100361,2147484041]`, with a deterministically chosen
primitive root in each field. This is a transport audit, not a proof of its
complete exceptional-prime set.

For the structural control, take three cosets of `mu_d` in `mu_(4d)`, so
`W_i=X^d-c_i`. With constant cofactors:

| n | d | Kernel dimension | Ambient tangent dimension | One-root replacements destroying the kernel |
|---:|---:|---:|---:|---:|
| 16 | 4 | 1 | 9 | 48 / 48 |
| 24 | 6 | 1 | 13 | 108 / 108 |
| 48 | 12 | 1 | 25 | 432 / 432 |
| 96 | 24 | 1 | 49 | Not enumerated |

Each row was checked at three characteristics, including one above 2^31.
The dimension pattern is `2d+1`; it has an elementary explanation. Simple
root motion identifies the `d` root velocities of each polynomial with all
coefficient perturbations of degree below `d`. The unique constant syzygy
requires their weighted sum to lie in `span(1,X^d)`. Since the perturbation
has degree below `d`, this imposes `d-1` independent constraints on `3d`
parameters, leaving `2d+1`.

The edit fragility also has a direct explanation. Replacing one root in
`W_i` changes it by a nonzero multiple of `W_i/(X-a)`, a nonconstant polynomial
of degree `d-1` when `d>=2`. The two unchanged binomials span
`span(1,X^d)`, so no constant syzygy survives. **Field-uniform structural
identities can be maximally fragile under one-root edits.** A fragility score
by itself cannot distinguish them from accidental finite-field witnesses.

## Instrument 2: realized stack provenance

For a stack `u_0+z*u_1`, enumerate every size-`s` agreement subset `S`. Its
Reed–Solomon parity checks give vectors

```
a_0 = H_S u_0,   a_1 = H_S u_1.
```

The equation `a_0+z*a_1=0` either supplies one candidate scalar, supplies
none, or holds for the whole field when both vectors vanish. This is the
existing exact, field-size-independent method from the repository's SYZ53
probe. The additional output here is every distinct decoded codeword, its
**maximal agreement support**, and every affine track through a pair of
scalar/codeword points. Duplicate subset certificates are compressed by
their decoded support, preserving the exact original multiplicity.

For a track `f_z=A+z*B`, let `t` be the number of coordinates at which
`u_0=A` and `u_1=B`. Every other coordinate agrees at at most one scalar.
Therefore `m` distinct track points each with agreement at least `s` satisfy

```
m (s-t) <= n-t     when t<s.
```

The tool checks this exact allocation law and reports overlaps among tracks
with at least three scalar/codeword points. Summing track capacities is
**not** justified: tracks overlap, and a scalar may admit multiple
codewords. Preserving that overlap is the purpose of the ledger.

## Second experiment: identical triple signatures, different incidence

We constrained each sampled pencil to agree on three prescribed cores at
scalars `0,1,2`, then took four seeded linear combinations of the resulting
stack kernel. Each line is exactly censused over all `C(n,s)` subsets.

| Configuration | Parameters | Exact finite bad-scalar counts of the four sampled stacks |
|---|---|---|
| F41 linear middle witness | n=20, k=10, s=14, p=41 | 5, 4, 4, 4 |
| Coset, small field | n=16, k=8, s=11, p=17 | 13, 7, 8, 12 |
| Same coset shape, medium field | n=16, k=8, s=11, p=1009 | 3, 3, 3, 3 |
| Same coset shape, large field | n=16, k=8, s=11, p=2147483713 | 3, 3, 3, 3 |

No sampled stack had a size-`s` common agreement set. At p=17, one sample
contains a five-scalar affine codeword track, and another has 20 three-scalar
tracks. At both larger fields every sampled stack has only its three forced
scalar/codeword points, with no three-point affine track.

The retained overlap is substantial: the p=17 sample with 20 three-point
tracks has only 15 distinct nodes on those tracks, representing 12 scalars.
Its 60 node memberships therefore contain 45 repeats, while one additional
bad scalar lies outside all three-point tracks. In the five-point-track
sample the track's common agreement set has size 10, just one below `s=11`;
its elementary capacity is 6 and its allocation slack is 1. These are exact
local explanations that a bad-scalar count or a syzygy rank alone does not
retain. Neither permits adding capacities across overlapping tracks.

These are observations on 16 explicit stacks, not global maxima or field
transfer theorems. The incidence census examines 38,760 subsets per F41
stack and 4,368 per coset stack. The same triple defect and tested jet/edit
signature coexist with very different realized pencil behavior. This
shows that these triple signatures do not determine the observed pencil
incidence; extra stack information is needed.

## Prior art and originality boundary

The underlying ingredients are established mathematics:

1. Ben-Sasson, Carmon, Ishai, Kopparty, Saraf, *Proximity Gaps for Reed–Solomon
   Codes* (2020), [primary paper](https://eprint.iacr.org/2020/654). The work
   analyzes algebraic decoding over function fields. Exact parity and affine
   pencil reasoning are in this methodological neighborhood, and already
   implemented locally.
2. Košir and Sethuraman, *Determinantal Varieties Over Truncated Polynomial
   Rings* (2002), [primary paper](https://arxiv.org/abs/math/0212051). Higher
   order determinantal equations over truncated rings are old; calling this
   a jet filtration does not make that theory new.
3. Teixeira, *Syzygy gap fractals—I. Some structural results and an upper
   bound* (2010), [primary paper](https://arxiv.org/abs/1008.0583). Characteristic
   dependence and determinant-based syzygy-gap analysis are also established.

Repository searches covered `tangent`, `deformation`, `jet`, `cokernel`,
`cotangent`, and `Smith normal` in the probes and knowledge base. Existing
work includes several unrelated deformation probes and an étale-domain
caution; no instance of this specific joint executable instrument was
identified in that bounded search. Web queries included combinations of
Reed–Solomon, syzygy, deformation, and affine pencil. This is a scoped search,
not proof of historical absence. **The potentially original contribution is
the combined diagnostic and certificate schema for this research workflow;
mathematical novelty is unestablished.**

## Path toward the full problem

The currently scalable part is the low-cofactor coefficient-matrix
filtration. Even there, generic large-degree dense matrices and root lists
eventually become infeasible; the production coset case must be represented
symbolically by binomials and quotient orbits. The stack census costs
`C(n,s)` and emphatically does not scale to `n=2^30`.

The next useful work is to replace the exhaustive agreement-subset generator
by a quotient-orbit or support-family representation **with a proof that no
relevant agreement support is omitted**, while retaining the codeword/track
provenance. Test that compressed generator against this exact oracle on
n=16,20, then n=24,32. A useful advance would be a certified way to account
for cross-track overlaps at scales where the naive census cannot run.
It would still require worst-case control over stacks and the precise MCA
predicate before it could support any prize statement.

The prototype does not establish such a compression theorem, a new bound,
a complete exceptional-prime classification, a production-scale algorithm,
or a proof that these tools were never used before.

## Third iteration: actual extension-field arithmetic

The production profile uses an extension field; computing modulo a large
prime does not implement it. `extension_field_probe.py` adds a separate
polynomial-basis adapter for `GF(p^e)` and leaves the two earlier prototypes
unchanged. Run:

```sh
python3 extension_field_probe.py
```

Elements are integers whose base-`p` digits encode polynomial coefficients.
Multiplication reduces modulo a supplied monic irreducible polynomial.
The adapter validates irreducibility by excluding **every** monic factor of
degree at most `e/2`; this bounded method is exact but exponential in `e`
and unsuitable for constructing large production fields. The chosen moduli
are `X^2+1` for F9 and `X^6+X^5+X^4+1` for F729, both over F3.

The resulting `extension_results.json` contains seven experiments, all
checked by comparing subset-syndrome enumeration with an independent
scalar/codeword census. For dimension-one codes, the latter uses the exact
constant-codeword agreement histogram instead of repeating equivalent
comparisons against every constant.

| Test | Exact bad-scalar count | Common agreement? |
|---|---:|---|
| Same embedded n=3,k=1,s=2 base stack in F3 | 2 | No |
| Same embedded stack in F9 | 2 | No |
| Same embedded stack in F729 | 2 | No |
| Nonbase stack on mu4 in F9, k=2,s=3 | 3 | No |
| Correlated control on mu4 in F9 | 9 | Yes |
| Correlated control on mu4 in F729, k=1,s=3 | 729 | Yes |
| Noncorrelated stack on mu4 in F729, k=1,s=3 | 1 | No |

The final example's unique bad scalar is the polynomial-basis element
`2X`, genuinely outside the prime subfield. All eight nonzero inverses in
F9 and all 728 in F729 were checked. Associativity and distributivity were
checked on every one of the 729 triples in F9, and on 729 deterministic
triples in F729; irreducibility supplies the general algebraic justification.
Every reported subset certificate is also matched against its unique
decoded codeword and full agreement support.

This implements degree-six extension arithmetic at tiny characteristic and
domain size. It does **not** run the full problem's characteristic, domain,
code dimension, or worst-case stack census, and the earlier larger-domain
experiments still use prime fields. The next integration step is to place
this tested arithmetic interface beneath the full provenance instrument,
then require identical prime-field outputs before scaling extension cases.
