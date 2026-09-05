# The prize reduction and an exact MCA calculation

The Paley conjectures and the grand Proximity Prize challenges remain open
in this workspace. This note checks the proposed connection against the
primary definition, proves a general elementary MCA bound, and gives exact
counterexamples to unrestricted domination by pairs of monomials.
These are ordinary mathematical proofs with integer computation, not new
Lean certificates. No claim of novelty is made.

## Verified source and version

The recovered local copy of Arnon–Boneh–Fenzi,
[*Open Problems in List Decoding and Correlated Agreement*](https://eprint.iacr.org/2026/680),
is dated **April 8, 2026** on its title page. The local filename
`sources/abf-2026-680-local-june.pdf` describes when it was archived, not
the paper's edition. Its SHA-256 is in `sources/manifest.json`.
The following checks concern the April copy; pages 5, 17, and 23 were
visually inspected after rendering. The July 6 revision has since been
retrieved and its relevant attack/list pages checked in
[list-to-winning-set.md](list-to-winning-set.md). That later note establishes
an unsafe suffix for the pinned winning-set target without closing the
Paley or general reduction problems.

Page 5 asks for the MCA error threshold and, separately, an interleaved
list-size threshold on smooth Reed–Solomon domains at rates
`{1/2,1/4,1/8,1/16}`, with target error such as `2^-128`.
It does not restrict the field to primes in a quartic window.
No Paley reduction was found in the inspected paper. This is a source
finding, not a theorem that no relationship exists.

Theorem 5.1, on page 23, states that list size at most `L` at radius `δ`
implies

\[
\varepsilon_{\rm mca}\!\left(C,1-\sqrt{1-\delta+\eta}\right)
\le\frac{L^2\delta n+1/\eta}{|F|}.
\]

The loss in radius matters: capacity list decoding for random Reed–Solomon
codes gives MCA up to the Johnson radius through this theorem. The pinned
note `deltastar-abf26-ld-mca-bridge-2026-06-13.md` incorrectly describes that
consequence as MCA up to capacity.

## The exact event

Let `F` be a finite field of size `q`, let `D⊆F` contain `n` distinct
evaluation points, and put `C=RS[F,D,k]`, where `1≤k<n`.
Words in `C` are evaluations of polynomials of degree less than `k`.
By Definition 4.3, a scalar `γ` is MCA-bad for words `f,g∈F^D` at radius
`δ∈(0,1)` exactly when some `S⊆D`, with `|S|≥(1-δ)n`, has both properties:

1. The folded word `f+γg` agrees with a codeword on `S`.
2. There are no two codewords agreeing respectively with `f` and `g` on
   that same set `S`.

The maximum over `f,g` of the fraction of such scalars is
`ε_mca(C,δ)`. Merely counting close folded words is insufficient unless
the no-joint-explanation condition has also been established.

## A circuit upper bound for every pair of words

For every `(k+1)`-subset `T⊆D`, define the linear functional

\[
\ell_T(w)=\sum_{x\in T}\frac{w(x)}{\prod_{y\in T\setminus\{x\}}(x-y)}.
\]

All denominators are nonzero. This is the coefficient of `X^k` in the
unique degree-at-most-`k` interpolant of `w|_T`. In particular,

\[
w|_T\text{ extends to degree }<k\iff\ell_T(w)=0,
\qquad \ell_T(X^k)=1.
\tag{1}
\]

If `γ` has an MCA witness `S`, the direction `g|_S` cannot extend to degree
less than `k`: otherwise, subtracting `γ` times that extension from the
folded codeword also extends `f|_S`, contradicting property 2. Thus `|S|>k`.
Choose any `k` points of `S` and interpolate `g` there by degree less than
`k`. Some further point of `S` violates this interpolant. Together they
form a set `T` with `|T|=k+1` and `ℓ_T(g)≠0`. The folded word still
agrees on `T`, so

\[
\gamma=-\ell_T(f)/\ell_T(g).
\tag{2}
\]

Writing `W=binom(n,k+1)`, every bad scalar therefore belongs to a union
of at most `W` singleton sets, independently of the radius. This proves

\[
\boxed{\varepsilon_{\rm mca}(C,\delta)\le\min(1,W/q).}
\tag{3}
\]

When `δ≥(n-k-1)/n`, every `(k+1)`-set is an allowed witness. Equation (2)
then describes the bad-scalar set **exactly**, using only sets with
`ℓ_T(g)≠0`. A set with `ℓ_T(f)=ℓ_T(g)=0` is jointly explainable and
contributes no bad scalar.

## Attainment over sufficiently large fields

**Theorem.** If `q>binom(W,2)` and `δ≥(n-k-1)/n`, with `δ∈(0,1)`, then

\[
\boxed{\varepsilon_{\rm mca}(C,\delta)=W/q.}
\tag{4}
\]

**Proof.** Choose the direction `g(x)=x^k`. For every `T`, equation (1)
makes its bad scalar `-ℓ_T(f)`. Distinct sets `T,U` give distinct linear
functionals on `F^D`: choose a coordinate in `T\setminus U`, where one
has a nonzero coefficient and the other has coefficient zero. Hence the
equation `ℓ_T(f)=ℓ_U(f)` cuts out a proper hyperplane of `q^(n-1)` offsets.
There are `binom(W,2)` such hyperplanes. Their union has size less than
`q^n`, so an offset exists outside all of them. Its `W` bad scalars are
distinct. Every corresponding set is allowed at the stated radius, and
`ℓ_T(g)=1` excludes joint explanation. This attains (3). ∎

This is an elementary finite-field argument. In particular, no
additive-character estimate is used. The subsequently retrieved July edition
adds Lemma 4.16, a lower bound linear in the number of allowed errors; it does
not assert the circuit-attainment formula above. No novelty claim is made.
Also, `W` grows rapidly with `n`, so this theorem does not settle the sponsor's
large-code parameter range.

## Why monomial pairs have at most 40 bad scalars at length eight

Let `D=μ_8⊆F*`, `k=4`, and `δ≥3/8`. The field has odd characteristic.
For a five-point witness `T`, write `U=D\setminus T`, a three-point set,
and denote its elementary symmetric sums by `e_j(U)`. The four non-code
monomials have circuit values

\[
(\ell_T(X^4),\ell_T(X^5),\ell_T(X^6),\ell_T(X^7))
=(1,-e_1(U),e_2(U),-e_3(U)).
\tag{5}
\]

Indeed, Lagrange interpolation gives
`ℓ_T(X^(4+j))=h_j(T)` for `j≥0`, where `h_j` is the complete homogeneous
symmetric polynomial. Its generating function is

\[
\sum_{j\ge0}h_j(T)z^j
=\frac1{\prod_{x\in T}(1-xz)}
=\frac{\prod_{u\in U}(1-uz)}{1-z^8}.
\]

Comparing the coefficients of degrees 0 through 3 proves (5).
The interpolation identity itself follows by expanding
`Σ_{x∈T} x^4/[∏_{y≠x}(x-y)(1-xz)]` in partial fractions; it equals
`1/∏_{x∈T}(1-xz)`.

There are exactly 24 triples `U` containing an opposite pair: choose one
of four pairs `{u,-u}`, then choose `w` from the six remaining points.
A triple cannot contain two opposite pairs, so this does not overcount.
For such `U={u,-u,w}`, equation (5) becomes

\[
(1,-w,-u^2,u^2w)\in(\mu_8)^4.
\]

Consequently, for every ordered pair of non-code monomials, all bad
scalars contributed by these 24 witnesses lie in `μ_8`, which has only
eight elements. Each of the remaining `56-24=32` witnesses contributes
at most one scalar by (2). Hence **every monomial pair has at most 40
MCA-bad scalars**. If the direction has degree less than 4, there are
no bad scalars; if only the offset does, the only possible bad scalar
is zero. Thus the same bound covers every monomial pair.

Combining this bound with (4) proves a general obstruction: for every
finite field of size `q>1540=binom(56,2)` containing `μ_8`, the maximum
MCA error at `δ≥3/8` is `56/q`, while the maximum over monomial pairs
is at most `40/q`. This gives arbitrarily large fields at fixed length
eight, not a growing-length asymptotic counterexample.

## Explicit smooth prime-field counterexamples to monomial pairs

Take `n=8`, `k=4`, and `δ=3/8`. For either `p=2017` or `p=65537`, use
the subgroup `D=μ_8⊆F_p*` and the pair

\[
f=X^5+2X^6+3X^7,\qquad g=X^4.
\]

The exact results are:

| Prime | Domain generator | Bad scalars for this pair | Maximum over all monomial pairs |
|---:|---:|---:|---:|
| 2017 | 438 | 56 | 40 |
| 65537 | 16 | 56 | 40 |

For `p=2017`, the domain is
`{1,229,438,548,1469,1579,1788,2016}`. This example also lies in the
literal quartic window `n^4/4≤p≤n^4`.
In each field, every one of the 56 five-point subsets gives a distinct
scalar `γ=-ℓ_T(f)`. Polynomial division by `∏_{x∈T}(X-x)` independently
produces a cubic agreeing with `f+γg` on exactly those five points.
The monomial direction `X^4` cannot agree with a cubic at five distinct
points, so the pairs are not jointly explainable on these witnesses.
The upper bound (3) proves the **global** maximum is exactly 56, not merely
a lower bound found by searching offsets.

The computation exhausts all 64 ordered pairs `(X^a,X^b)`, `0≤a,b<8`.
Exactly `(a,b)=(4,5),(5,4),(6,7),(7,6)` attain 40. Arbitrary nonzero
coefficients multiplying the two monomials only rescale the parameter `γ`;
codeword additions do not change the MCA event. Every nonnegative exponent
reduces modulo 8 on this domain. Thus this exhausts the usual monomial-pair
class, including those harmless changes of representation.

The conclusion is precise: **a maximum need not be attained by a pair of
monomials**. It does not disprove a claim that only the direction can be
chosen monomial while the offset remains arbitrary; the attaining example
already has a monomial direction. Nor does this finite example establish
anything about a proposed asymptotic restriction to other radii.

Reproduce with:

```sh
python3 experiments/mca_circuit_certificate.py
```

The script checks primality and subgroup order, all circuit functionals by
independent polynomial remainders, all monomial pairs, and every explicit
cubic witness. It writes both scalar lists and all 112 interpolation
certificates to `results/mca_circuit_certificate.json`, with a source hash.
It uses only Python's standard library and exact modular arithmetic.

## Consequence for the proposed Paley route

At the pinned repository commit
`5b00e50c3c51b3c944201a1749a1a8132e5ce167`, the introductory chain in
`DyadicLacunaryDeltaStar.lean` attributes monomial-pair extremality to
`FarLineIncidenceEquivariance.lean`. The latter's actual theorems establish
invariance under code-preserving coordinate permutations. Such invariance
allows taking a quotient by orbits; it does not force a maximizer to be
fixed by the symmetry or to consist of monomials. The examples above
explicitly rule out that unrestricted conclusion.

The same repository already records other monomial-domination failures in
`MonomialDominationKilled.lean` and `MonomialDominationBoundaryRefuted.lean`.
The current finding corroborates the gap on a larger prime field, including
a quartic-window instance. It is not a claim that this issue was previously
unknown to the project.

Accordingly, proving the thin-subgroup exponential-sum bound would still
require an independently verified theorem connecting it to **arbitrary**
MCA word pairs at the sponsor's radii and field sizes. Neither the orbit
invariance theorem nor the monomial-pair calculation supplies that step.
The two-set quadratic Paley conjecture is a further distinct statement.

The later [subset-sum note](subset-sums-and-lists.md) supplies one explicit
special-family connection: cancellation bounds the subset counts that
exactly describe the lists of `X^{k+1}` and the bad scalars of the pair
`X^{k+1},X^k` at radius `1-(k+1)/n`. It also gives a finite list lower
bound for the newer pinned official profile. These statements preserve
the distinction between a particular pair and the maximum over all pairs;
they do not repair the unsupported extremality step above.
