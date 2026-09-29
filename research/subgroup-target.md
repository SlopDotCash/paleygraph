# Primary target recovered from the prize repository

## Provenance and scope

On 2026-09-04, the old `elizaOS/proximityprize` repository URL resolved to
`SlopDotCash/proximityprize`. Its current `main` commit was
`5b00e50c3c51b3c944201a1749a1a8132e5ce167`.
Source files were downloaded at that exact commit; see
`sources/prior-work/manifest.json`. The current dossier still describes the
production result as open. This is a live verification of the repository
state, not reliance on the earlier memory snapshot.

The relevant [subgroup attack note](https://github.com/SlopDotCash/proximityprize/blob/5b00e50c3c51b3c944201a1749a1a8132e5ce167/docs/kb/deltastar-464-attack-01-paley-direct-smooth-subgroup.md)
targets additive sums over a dyadic subgroup. Its
[double-sum transfer note](https://github.com/SlopDotCash/proximityprize/blob/5b00e50c3c51b3c944201a1749a1a8132e5ce167/docs/kb/deltastar-464-paley-double-sum-singleton-gate-2026-06-25.md)
explicitly says an ordinary large-set Paley double-sum theorem does not
directly give the required individual-period bound.

Thus there are two different mathematical targets here:

| Target | Character and summation set | Current local status |
|---|---|---|
| Classical two-set Paley | Quadratic multiplicative character `χ(a-b)`, arbitrary large `A,B` | Standard bounds reconstructed; generic Gaussian moment route refuted |
| Prize campaign subgroup bound | Additive character `e_p(bh)`, one specific dyadic multiplicative subgroup `H` | Exact centered moment formulation and finite checks; required upper bound open |

The Gaussian counterexample for arbitrary `B` in the first row is **not** a
counterexample to Gaussian period moments in the second row. Both the
character and the allowed sets have changed.

The historical essay uses `p≈n⁴` and sometimes `n≈2³⁰`. These are not
complete quantifiers. A precise asymptotic target could fix constants
`0<c₋≤c₊`, require `c₋n⁴≤p≤c₊n⁴`, `n=2^μ`, `n|(p-1)`, and seek a
constant independent of `μ,p,b`. The sponsor's concrete code parameters
and required numerical constant would still have to be matched separately.
This suggested formulation has not been declared sponsor-equivalent.

## Exact moment formulation

Let `H≤F_p*` have size `n`, let `m=(p-1)/n`, and write

\[
\eta_b=\sum_{h\in H}e_p(bh),\qquad M=\max_{b\ne0}|\eta_b|.
\]

For an integer `r≥1`, define

\[
E_r(H)=\#\{(x_1,\ldots,x_r,y_1,\ldots,y_r)\in H^{2r}:
                 x_1+\cdots+x_r=y_1+\cdots+y_r\}.
\]

Additive-character orthogonality gives

\[
\sum_{b\in\mathbb F_p}|\eta_b|^{2r}=pE_r(H).
\]

Since `η_0=n` and `η_b` is constant on each multiplicative coset of `H`,

\[
\boxed{Q_r(H):=\frac{pE_r(H)-n^{2r}}n
       =\sum_{bH\in\mathbb F_p^*/H}|\eta_b|^{2r}.}
\]

In particular,

\[
M^{2r}\le Q_r(H)\le mM^{2r},\qquad
M=\lim_{r\to\infty}Q_r(H)^{1/(2r)}.
\]

These follow simply because `Q_r` is a sum of `m` nonnegative terms,
one of which is the maximum. For `r=1`, `E_1=n`, so

\[
\sqrt{\frac{n(p-n)}{p-1}}\le M\le\min(n,\sqrt{p-n}).
\]

The lower bound is asymptotic to `√n` when `n/p→0`. Replacing it by an
exact `M≥√n` without additional hypotheses is invalid; the full subgroup
`H=F_p*` has `M=1`.

## The sufficient estimate, without hiding its difficulty

If one could prove, with an absolute `K`,

\[
Q_r(H)\le m(Krn)^r
\tag{SG}
\]

at `r=ceil(log m)` for `m≥e` in the target family, then

\[
M^2\le m^{1/r}Krn\le eK n(\log m+1)
\le2eK n\log m.
\]

This would supply the desired square-root bound with a logarithmic factor.
Neither the elementary identities nor the finite calculations establish (SG).
A bound with constants depending arbitrarily on `r` is insufficient when
`r` grows with `p`; its growth must be controlled.

### A necessary lower-order consequence

The same target already requires a substantial additive-energy estimate.
Since `Σ_(b≠0)|η_b|²=n(p-n)`, for every integer `r≥1` one has

\[
pE_r-n^{2r}\le M^{2r-2}n(p-n).
\]

Thus `M²≤C²nL`, with `L=ln(p/n)`, would imply in particular

\[
E_2(H)\le\frac{n^4}{p}+C^2n^2\left(1-\frac np\right)L.
\tag{NE}
\]

In a fixed quartic window this is `E_2(H)=O(n² log n)`. The structural
Gaussian-coefficient counterexamples below do not contradict this looser
necessary bound. No uniform proof of (NE) is supplied in this workspace.
Conversely, this energy bound alone would not give the target maximum:
the fourth-moment inequality still loses a factor from the number of cosets.
It is a necessary checkpoint, not a replacement objective.

One particularly strong sufficient candidate is the real Gaussian comparison

\[
pE_r(H)-n^{2r}\le(p-1)(2r-1)!!\,n^r.
\tag{GC}
\]

Since `(2r-1)!!≤(2r)^r`, (GC) would imply (SG) with `K=2`.
The literal universal version of (GC) is **refuted below**, including a finite
example in an explicit quartic window. A version with a larger constant,
additional arithmetic hypotheses, or a sufficiently-large-size qualifier is
not refuted by those finite examples. The weaker (SG) remains open.

The later [positive-product criterion](positive-product-moments.md)
gives another sufficient route to (SG): control the positive product
of periods from the two halves of each dyadic subgroup at one even
logarithmic moment depth. Its implication has an explicit constant,
and the corresponding centered balanced-count condition is verified
through the known order-64 tower. The uniform hypothesis is unproved;
extra raw balanced relations are forced by the principal-frequency
term at larger orders.

The literal amplitude constant `C=√2` in `M≤C√(n ln(p/n))` also fails
in the quartic-window example below. An analytic Lean certificate proves
`Re η₁>43` and `ln(p/64)<12`. Exact moments through depth 12 bound `M²≤1970`
and satisfy (SG) with `K=2` at every tested depth. See
[the finite spectral obstruction](finite-spectral-obstruction.md) for the
proof, computational upper certificate, and distinction from `M²≤2n ln p`.
The later [polynomial certificate](period-polynomial-certificate.md) identifies
the exact maximizing coset in this field and encloses `M=η₁` to width `10^-18`.
It uses all signed moments through order 24, independently checked against
the carry counts. No uniform bound follows from this finite computation.

## Exact failure of the Gaussian coefficient in a quartic window

Take the prime `p=6700417`, `n=64`, and `H=⟨2⟩⊂F_p*`.
The element 2 has order 64. This prime is a factor of
`2³²+1=641·6700417`, and

\[
n^4/4=4194304\le6700417\le16777216=n^4.
\]

An exact count gives

\[
\#\{(a,b)\in H^2:1+a-b\in H\}=201.
\]

Normalize an additive quadruple `(x₁,x₂,y₁,y₂)` by its nonzero first
coordinate: `(a,b)=(x₂/x₁,y₁/x₁)` and `y₂/x₁=1+a-b`.
This is a bijection between additive quadruples and `H` times the displayed
set. Therefore `E₂(H)=64·201=12864`, and

\[
pE_2(H)-n^4=86177387072
>82334711808=3(p-1)n^2.
\]

This disproves the exact coefficient in (GC) at `r=2`. It does **not**
disprove a bound with an unspecified absolute constant or an asymptotic
statement allowing finitely many exceptions. In particular it does not
disprove the primary period conjecture or either prize challenge.

`SubgroupQuarticCounterexample.lean` checks primality, all 64 nonzero
representatives, multiplication closure, their 64th powers, the collision
count 201, the quartic window, and the displayed integer inequality.
The standard finite-subgroup and normalization interpretations are explained
here in ordinary mathematics; that bijection has not been formalized in this
local Lean file. The independent Python convolution computes the same energy.

The earlier dense calculation at `(p,n)=(65537,64)` similarly has
`E₂=19776` and violates (GC), but its prime is much smaller than `n⁴`.

### A structural reason the exact coefficient can fail

For any multiplicative subgroup `H` containing `-1` and `2`, of size `n`,
in characteristic `p>5`, the sharpened bound is

\[
E_2(H)\ge3n^2+9n.
\]

**Proof.** The usual permutation solutions to `a+b=c+d` number
`2n²-n`. The zero-sum solutions `(a,-a,c,-c)` number `n²`;
their overlap with the permutation solutions numbers `2n`. Their union
therefore has size `3n²-3n`.

Use the bijection from additive quadruples to zero-sum words obtained by
negating the last two coordinates. For each `u∈H`, the word
`(u,u,2u,-4u)` sums to zero. Its three values are distinct for `p>5`, so
it has `4!/2!=12` different permutations. The repeated value identifies `u`,
making these sets of permutations disjoint for different `u`. None of the
words contains an opposite pair, so none belongs to the earlier degenerate
classes. They supply `12n` additional solutions and prove the bound.

The earlier estimate `3n²+n` counted only four of these twelve placements;
it remains valid for `p>3`, but the stronger result applies to the present
quartic-window example. That example attains equality in the stronger bound.

Consequently, whenever additionally `9p>n³-3n`,

\[
pE_2-n^4-3(p-1)n^2\ge9pn-n^4+3n^2>0.
\]

Thus the Gaussian coefficient fails for **every** subgroup satisfying those
conditions. An unbounded family of such primes in the prescribed quartic
window has not been proved here. The finite example above supplies one.
This identifies an arithmetic resonance that a valid moment upper bound
must allow; it does not rule out (SG) with a larger constant.

The subsequent [geometric-lift analysis](geometric-lift-and-alias.md) proves
a uniform moment bound for the full geometric cycle and isolates the signed
error introduced by reduction modulo the prime. An independent carry
computation shows that doubling alone accounts for this example's fourth
and sixth moments; additional prime relations first occur at order 8.

The subsequent [coset-coherence analysis](coset-coherence.md) derives the
exact nonlinear identity `x_j²=n+Σ_s k_s x_(j+s)` and its autocorrelations.
It constructs synthetic vectors satisfying the latter, and all Gauss-sum
magnitude constraints, with almost maximal spikes. The synthetic vectors
fail the nonlinear identity and are not actual period counterexamples.

## Correction to a broad impossibility claim in the historical essay

The old phase-cancellation essay observes that
`E_r≥n^{2r}/p`, so the uncentered upper estimate
`M≤(pE_r)^(1/(2r))` has a right side at least `n`.
That specific observation is correct. It does **not** rule out all arguments
using magnitudes, energies, or higher moments.

The displayed identity for `Q_r` subtracts `n^{2r}` exactly; its upper
estimate converges to `M`. Thus the **full magnitude distribution determines
the maximum exactly**. Fixing the second moment does not fix that distribution.
For example, the magnitude vectors `(2,0,0,0)` and `(1,1,1,1)` have the same
squared norm but different maxima and fourth moments.

This correction does not supply a new upper bound. It prevents an argument
about one uncentered estimate from being treated as a universal prohibition
on the centered moment approach. Conversely, merely writing `Q_r` instead of
`E_r` is not a breakthrough; the uniform estimate (SG) remains the obstacle.

## Reproducible finite checks

`experiments/subgroup_moments.py` constructs `H` with a primitive root and
uses integer cyclic convolution to count ordered sums. It then evaluates
`E_r` and `pE_r-n^{2r}` exactly. It checks total mass, invariance under
multiplication by a subgroup generator, the second moment, and an independent
difference-count identity at fourth-moment order.

The initial dense suite covers ten `(p,n)` pairs and 82 moment levels,
including `(257,4)` and `(65537,16)` where `p=n⁴+1`, with depth up to 10.
The sparse suite reaches `p=6700417` without allocating a field-sized array.
These are small instances. They do not reach the production scale or
establish a limiting law.


The [signed quotient representation](signed-quotient-operators.md) gives a
sparse symmetric integer matrix S_H whose eigenvalues are exactly the
nonprincipal periods, once per H-coset. Its even trace is Q_r. In the working
quartic window, however, `ρ(|S_H|)≥n−2/n`; discarding its entrywise signs
therefore loses almost all the desired cancellation. The corresponding
child-product matrix and its positive part express the sufficient mixed
moment criterion exactly. No uniform signed spectral estimate is supplied.

## What remains

1. Prove or refute a uniform estimate such as (SG) for the specified family.
2. Independently audit the exact theorem chain from that analytic estimate
   to the sponsor's MCA and list-decoding statements, including all parameters
   and hypotheses. An essay describing the statements as equivalent is not
   enough evidence of that equivalence.
3. Supply both threshold directions if claiming an exact prize threshold.

One specific sufficient hypothesis is now isolated in
[`jacobi-coefficient-route.md`](jacobi-coefficient-route.md): bound the
orthogonal-polynomial recurrence coefficients by
`|α_j|≤A√(n(j+1))`, `β_j≤Bnj` through degree `ceil(ln m)`, with absolute
constants. Its implication is proved, and the original finite example passes
with `A=B=1`. The literal B=1 extension is now refuted at p=67403009, n=128.
The uniform hypothesis with other constants is not proved and is not claimed
to follow from the finite calculation or the nonlinear period identity.

The [new official profile audit](official-profile-and-trace.md) also records
a concrete mismatch between this historical quartic-window target and the
current pinned IRS benchmark. Its domain is in the KoalaBear prime field,
embedded in a sextic coefficient field. Nonzero trace-zero extension
frequencies give full-size periods. No direct transfer to that benchmark is
established; it is one official parameter point, not the entire grand prize.

The current work resolves none of these three obligations.
