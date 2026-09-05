# The current official profile and the trace obstruction

The new official benchmark is relevant context for the prize link. It
does not replace the full Paley or Proximity Prize goals. This note
audits its parameters and proves a necessary correction for any attempted
extension-field spectral argument. It does not prove that the benchmark
uses the particular character sum below, or establish a reduction to it.

## Pinned source facts

The Ethereum Foundation's
[August 20 announcement](https://blog.ethereum.org/2026/08/20/better-codes-challenge)
links the new better.codes benchmark to the broader Proximity Prize.
The public site and repository identified contract commit
`b34c0131cfa36b51111521541d7d3e35c8791082` when inspected on September 4.
The authoritative parameter source is
[`IRSProfile.lean` at that commit](https://github.com/proximity-prize/proximity-prize/blob/b34c0131cfa36b51111521541d7d3e35c8791082/ProximityPrize/Benchmark/IRSProfile.lean).

| Parameter | Value |
|---|---|
| Base field | `F_p`, `p=2³¹-2²⁴+1=2130706433` |
| Coefficient/challenge field | `F_q=F_(p⁶)` |
| Evaluation domain | `H=μ_(2¹⁸)⊂F_p*`, embedded in `F_q` |
| Domain size | `n=262144` |
| Scalar dimension | `131072` |
| Interleaving | `8` |
| Total dimension | `1048576` |
| Repetitions | `128` |
| Alphabet rate | `1/2` |

The modulus and sextic presentation come from the pinned CompPoly revision
`641694629e4557520a1539b272ec338c9f3044c7`:
[`Basic.lean`](https://github.com/zksecurity/CompPoly/blob/641694629e4557520a1539b272ec338c9f3044c7/CompPoly/Fields/KoalaBear/Basic.lean)
and [`Ext6.lean`](https://github.com/zksecurity/CompPoly/blob/641694629e4557520a1539b272ec338c9f3044c7/CompPoly/Fields/KoalaBear/Ext6.lean).
The extension is `F_p[θ]/(θ⁶+θ³+1)`.

Both `p<n⁴/4` and `q>n⁴` hold. Thus neither substituting the base field
size nor the extension field size makes this profile an instance of our
current quartic-window working target. A more general theorem or another
argument would be needed for a bridge to this profile.

The pinned [`TargetLower.lean`](https://github.com/proximity-prize/proximity-prize/blob/b34c0131cfa36b51111521541d7d3e35c8791082/ProximityPrize/Benchmark/TargetLower.lean)
requires the certified combination-round error to be at most `2^-128` at
an admissible radius, and separately scores `(1-δ)^128`.
At the pinned ArkLib revision `e65197892890b8fd9b0dc05b8980273cf1d595cc`,
[`coe_certifiedGammaError`](https://github.com/Verified-zkEVM/ArkLib/blob/e65197892890b8fd9b0dc05b8980273cf1d595cc/ArkLib/ProofSystem/ToyProblem/Impl/IRS.lean#L669)
identifies that error with the MCA term plus the two-interleaved list term
divided by q. The source does not identify it with a subgroup period
maximum. Its [`upper target`](https://github.com/proximity-prize/proximity-prize/blob/b34c0131cfa36b51111521541d7d3e35c8791082/ProximityPrize/Benchmark/TargetUpper.lean)
uses a worst-case winning-set quantity over an entire suffix of radii.

These are source inspections. The new repository's Lean targets and
dependency closure have not been rebuilt or independently kernel-audited
in this workspace. Pinned copies and checksums are retained in
[`sources/official-prize-2026-09-04/manifest.json`](../sources/official-prize-2026-09-04/manifest.json).

## Exact trace reduction for a subgroup in the prime field

Let `q=p^d`, `d>1`, and `H≤F_p*`, of size n. Write

\[
\operatorname{Tr}(b)=b+b^p+\cdots+b^{p^{d-1}},\qquad
\psi_q(b)=e_p(\operatorname{Tr}(b)).
\]

The trace is a surjective F_p-linear map to F_p. It lands in F_p by
Frobenius invariance; it is nonzero because its defining nonzero polynomial
has degree `p^(d-1)<q` and cannot vanish on all q elements. A nonzero
linear map to the one-dimensional space F_p is surjective. Its kernel
therefore has `p^(d-1)` elements.

For every `h∈H⊂F_p`, linearity gives the exact identity

\[
\eta_q(b):=\sum_{h\in H}\psi_q(bh)
=\sum_{h\in H}e_p(h\operatorname{Tr}(b))
=\eta_p(\operatorname{Tr}(b)).
\tag{11}
\]

Every nonzero trace-zero b has `η_q(b)=n`. Consequently

\[
\boxed{\max_{b\in F_q^*}|\eta_q(b)|=n,}
\]

with at least `p^(d-1)-1` maximizing frequencies. In fact these are all
the maximizers when `|H|>1`: equality in the triangle inequality forces
all `e_p(th)` to have the same phase, which for two distinct h and
nonzero t is impossible. This does not contradict the prime-field Paley
target; it identifies an extra obstruction introduced by the extension.

For an even moment, additive relations among members of H are unchanged
by the embedding. Define

\[
Q_r^{(p)}=\frac{pE_r(H)-n^{2r}}n,\qquad
Q_r^{(q)}=\frac{qE_r(H)-n^{2r}}n.
\]

These are sums over nonzero multiplicative cosets in the corresponding
fields. Then

\[
\boxed{Q_r^{(q)}=p^{d-1}Q_r^{(p)}
 +(p^{d-1}-1)n^{2r-1}.}
\tag{12}
\]

Thus removing only the zero frequency leaves many full-size periods.
The correct deletion, if the aim is to study cancellation inherited
from the prime field, is the entire trace-zero space:

\[
\sum_{\operatorname{Tr}(b)\ne0}|\eta_q(b)|^{2r}
=p^{d-1}[pE_r(H)-n^{2r}].
\tag{13}
\]

Restricting to these frequencies recovers exactly the original prime-field
period problem, with p as the harmonic-analysis field size. It supplies
no new cancellation by itself and does not automatically control the
MCA or list quantities of arbitrary extension-valued words.

## A concrete witness in the official field

For the displayed sextic, θ is nonzero and

\[
\operatorname{Tr}(1),\operatorname{Tr}(\theta),\ldots,
\operatorname{Tr}(\theta^5)=(6,0,0,-3,0,0).
\]

In particular, `b=θ` gives `η_q(θ)=262144` exactly. The trace-zero
space consists of coefficient vectors satisfying `6c₀-3c₃=0` and has
dimension five. This is an ordinary mathematical certificate, not a new
Lean proof. The existing prime-field conjecture remains unproved.

[`extension_trace_audit.py`](../experiments/extension_trace_audit.py)
independently checks the base modulus by trial division, the order of
the specified subgroup generator, and sextic irreducibility by Rabin's
criterion. Exact modular polynomial arithmetic checks every basis trace.
It also exhausts `F_25=F_5[u]/(u²-2)` with `H=F_5*`: the nonzero
frequency periods are four copies of 4 and twenty copies of -1. Direct
convolution checks (12) through moment order 16. Results are saved in
[`extension_trace_audit.json`](../results/extension_trace_audit.json).

A subsequent [subset-sum analysis](subset-sums-and-lists.md) proves a
limited prime-field character-sum connection to monomial lists and MCA
pairs. It gives a finite lower bound for the list term of the pinned
official certificate, while preserving the distinction between that
certificate and actual winning-set soundness. It does not supply the
general spectral reduction left open in this note.


The subsequent [list-to-winning-set proof](list-to-winning-set.md) establishes
an ordinary mathematical upper-track suffix statement for this exact pinned
profile: winning density exceeds 2^-128 from radius 122641/262144 up to the
minimum relative distance. It uses a base-field coefficient-pigeonhole list
and a linear projection, not a Paley estimate. The associated exact score
inequality has centibits 11649. No Lean certificate or sharp threshold is
claimed, and the general spectral-to-code bridge remains open.
