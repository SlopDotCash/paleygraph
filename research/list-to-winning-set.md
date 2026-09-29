# From a scalar list to the prize winning set

**Status: the Paley conjectures and exact grand-prize thresholds remain
unproved.** This note proves an ordinary mathematical lower bound on the
winning-set density of the pinned toy protocol. For its interleaved
Reed–Solomon profile, the density exceeds 2^-128 throughout

\[
 \boxed{\frac{122641}{262144}\le\delta<\frac{131073}{262144}.}
 \tag{1}
\]

The left endpoint is about 0.4678382874. Its corresponding upper-track
score inequality holds with **11649 centibits**. This is not a Lean
certificate, an accepted submission, a claim of the best known bound, or
an end-to-end implementation of an attacking prover. It improves what
this workspace has established about the actual winning-set quantity;
it does not identify the optimal threshold.

## Source and target checked

The July 6, 2026 edition of
[Arnon–Boneh–Fenzi, *Open Problems in List Decoding and Correlated Agreement*](https://eprint.iacr.org/2026/680)
was retrieved successfully after earlier attempts had failed. The
55-page PDF is archived as `sources/abf-2026-680-july.pdf`; its title,
metadata, and relevant rendered pages were checked. The April archive is
retained separately.

ABF Definition 6.11 on page 31 defines the winning set for Construction
6.9. Lemma 6.12 bounds it using the two-word list, under an assumption
that the maximum such list is smaller than the field. The argument below
uses a single-word list and fixes the second received word to zero.
It has no assumption that the maximum list is smaller than the field.
No novelty is claimed for this elementary list/projection method.

The pinned
[upper target](https://github.com/proximity-prize/proximity-prize/blob/b34c0131cfa36b51111521541d7d3e35c8791082/ProximityPrize/Benchmark/TargetUpper.lean)
asks for `epsilonStar < winningSetDensity` over the entire admissible
suffix, and separately a score bound. The pinned ArkLib
[definition](https://github.com/Verified-zkEVM/ArkLib/blob/e65197892890b8fd9b0dc05b8980273cf1d595cc/ArkLib/ProofSystem/ToyProblem/SoundnessBounds.lean#L1146)
is the supremum of winning fractions over violating instances. Its
comments distinguish this from a formally proved converse to optimal
protocol game error. We establish that precise combinatorial quantity
in ordinary mathematics, without claiming the missing formal converse.
The encoding map remains fixed throughout; allowing a new encoding map
with the same code range would change the inner-product constraints.

## A list-to-winning-set lemma without a maximum-list restriction

Let F have q elements, and let
`Enc:F^K→(F^s)^n` be an **injective** linear encoder with minimum
relative distance d_min. Suppose a received word w has L distinct
messages a in a set A satisfying

\[
 \Delta(\operatorname{Enc}(a),w)\le\delta_0<d_{\min}.
\]

Then there exists a single choice of v∈F^K such that, for every
`δ₀≤δ<d_min`, the instance

\[
 (f_1,f_2)=(w,0),\qquad (\mu_1,\mu_2)=(0,1)
 \tag{2}
\]

is violating and has winning-set density at least

\[
 \boxed{\frac{L}{q+L-1}.}
 \tag{3}
\]

Moreover, if `binom(L,2)<q`, there exists such a v for which the
stronger bound is L/q.

**Projection proof.** Choose v uniformly from F^K and map
`a↦⟨a,v⟩`. For distinct a,b, the nonzero linear functional
`v↦⟨a−b,v⟩` vanishes with probability exactly 1/q. Let C(v) count
ordered pairs in A with equal images, including the diagonal. Then

\[
 \mathbb E_v C(v)=L+L(L-1)/q.
\]

For some v, C(v) is at most its mean. Cauchy–Schwarz on the fibers gives

\[
 |\{\langle a,v\rangle:a\in A\}|
 \ge\frac{L^2}{C(v)}\ge\frac{qL}{q+L-1}.
 \tag{4}
\]

For the stronger assertion, the union bound over unordered pairs is
`binom(L,2)/q<1`; some projection is injective on A.

**Violation and winning proof.** Any message b satisfying
`⟨b,v⟩=μ₂=1` is nonzero. Injectivity of Enc gives a nonzero codeword,
whose distance from the second received word zero is at least d_min.
Therefore no jointly constrained pair of codewords can agree with (w,0)
on enough coordinates at any δ<d_min. The input is violating throughout
the suffix.

For each `γ=⟨a,v⟩` with a∈A, the combined received word is still w,
and the combined claim is `μ₁+γμ₂=γ`. The codeword Enc(a) satisfies
that claim and is within δ₀, hence within δ. Thus every image in (4) is
a winning challenge. The same v and w work throughout the suffix;
no general monotonicity assumption on worst-case winning density is
needed. Dividing (4) by q proves (3). ∎

This is a witness argument for ABF's winning set. It does not assert
that the MCA error equals this density, or that every large two-word
list can be substituted for the single-word list in this lemma. In fact,
the pair (w,0) has no MCA-bad scalar: any agreement witness for w has the
joint explanation (u,0). The winning instance exploits the additional
inner-product claim. Thus this argument does not lower-bound the grand
MCA challenge's error by the same quantity.

## A coefficient-pigeonhole list on any base-field domain

Let `D⊂F_p`, `|D|=n`, `F_p⊂F_q`, and let `k≤r<n`.
For the scalar code `RS[F_q,D,k]`, there exists a received word with
at least

\[
 \boxed{L_r=\left\lceil\frac{\binom nr}{p^{r-k}}\right\rceil}
 \tag{5}
\]

nearby codewords at radius `δ_r=1−r/n`.

**Proof.** For each r-subset T⊂D, form

\[
 Q_T(X)=\prod_{x\in T}(X-x).
\]

Partition these subsets according to the coefficients of degrees
`k,k+1,…,r−1`. There are at most p^(r−k) labels, so a fiber has at
least L_r subsets. All their polynomials have the same part f of
degrees at least k, with leading term X^r. Set

\[
 u_T=f-Q_T.
\]

Then deg u_T<k, and f−u_T=Q_T has exactly r roots in D. Distinct
subsets yield distinct polynomials and codewords. This proves (5).
Indeed, the entire list at f is exactly that coefficient fiber: any
nearby codeword makes f−u a monic degree-r polynomial with at least
r roots in D, hence a unique Q_T. ∎

Embedding the word and codewords into one interleaved row, with all
other rows zero, preserves column distance. Therefore (5) applies to
the pinned eight-row code. Smoothness of D and a Paley bound are not
needed for this construction.

This uses the same common-leading-coefficient principle as ABF
Appendix C. We use the consistent sign `u_T=f−Q_T`. In the July PDF,
the displayed definition in the proof of Lemma C.5 on page 49 has the
opposite sign while the following displayed difference is incompatible
with it; replacing that definition by `u_S=f−V_S` gives the stated
vanishing identity. This is a local sign correction, not a refutation
of the list-size statement.

## Exact specialization to the pinned profile

The [profile audit](official-profile-and-trace.md) fixes

\[
 p=2130706433,\quad q=p^6,
 \quad n=262144,\quad k=131072,
 \quad s=8,\quad d_{\min}=131073/262144.
\]

Choose r=139503. The number of fixed coefficients is r−k=8431,
and exact integer arithmetic gives

\[
 L_r=
 \left\lceil\frac{\binom{262144}{139503}}{2130706433^{8431}}\right\rceil
 =85677801616821870413970774.
 \tag{6}
\]

The challenge field has size

\[
 q=93571093019388561295270373781649880353786165192103559169.
\]

The following comparisons hold with integers:

\[
 \binom{L_r}{2}<q,
 \qquad 2^{128}L_r>q.
 \tag{7}
\]

Select exactly L_r messages from a fiber supplied by (5). The
injective-projection part of the lemma proves a winning fraction at
least L_r/q>2^-128, throughout (1). This gives an existential word f
and vector v. They are not enumerated at production size, and no
claim of an efficient explicit search for them is made.

For comparison, increasing the required agreement count to r+1 gives

\[
 \left\lceil\frac{\binom{262144}{139504}}{2130706433^{8432}}\right\rceil
 =35350350170772326
 <\left\lfloor q/2^{128}\right\rfloor+1
 =274980728111395088.
 \tag{8}
\]

For r≥n/2, the unrounded ratio in (5) decreases, since its successive
ratio is `(n−r)/(p(r+1))<1`. Thus the selected radius is optimal **only
within this coefficient-pigeonhole lower bound at this fixed profile**.
Failure of that lower bound at the next grid point does not establish
safety there or exclude a stronger list construction.

The upper target's score condition at this radius is verified exactly by

\[
 2^{11649}\,139503^{12800}\ge262144^{12800},
 \qquad
 2^{11648}\,139503^{12800}<262144^{12800}.
 \tag{9}
\]

Taking 100th roots proves
`2^(−11649/100)≤(139503/262144)^128`. Thus 11649 is the least integer
centibit value satisfying that score inequality for this radius. This
calculation does not prove a lower-track security guarantee of 116.49 bits.

## The earlier huge list now bounds the winning set

The previously verified [root-lift construction](subset-sums-and-lists.md)
gives L≥2^8154 for a single-word list at `δ₀=4095/8192`. Since
q<2^186, (3) yields

\[
 \boxed{\operatorname{winningSetDensity}(\delta)
       >1-2^{-7968}\quad
       (4095/8192\le\delta<d_{\min}).}
 \tag{10}
\]

This is a new consequence in the workspace of the earlier list result;
it follows from (2)–(4), rather than from the upper certificate Γ.
The earlier observation Γ>1 alone did not establish any such conclusion.

A smaller explicit received-word example also follows from the July
antipodal construction: in G=μ_128, take the singleton 1 and 32
opposite pairs excluding {1,−1}. All `binom(63,32)=916312070471295267`
subsets have size 65 and sum 1. Lifting by a=2048 gives the center
`X^133120−X^131072` and that many codewords at radius 63/128.
The improved radius (1) uses (5), not an assertion that these sums are
zero. The earlier fixed-order obstruction to *odd zero-sum subsets*
therefore remains correct within its stated scope.

## Verification and remaining work

[The script](../experiments/list_to_winning_set.py) independently enumerates
coefficient fibers, codewords, projections, and winning challenges in four
small fields. It checks all 86086 projection vectors across the cases and
their exact collision average. Two cases have L≥q; they attain all q
winning challenges, directly testing the regime excluded by the maximum-list
assumption in ABF Lemma 6.12.

| p | n | k | r | Exact list | Winning challenges |
|---:|---:|---:|---:|---:|---:|
| 5 | 4 | 2 | 2 | 6 | 5 |
| 7 | 6 | 3 | 3 | 20 | 7 |
| 13 | 6 | 3 | 4 | 3 | 3 |
| 17 | 8 | 4 | 5 | 4 | 4 |

At production size the script checks (6)–(9), minimum-distance/grid
conditions, the original list parameters for (10), and source hashes.
It computes successive binomials by exact divisions and checks the
selected and next values independently with `math.comb`. No numerical
logarithm is an acceptance condition. The JSON retains input parameters,
compact hashes of the large comparison integers, and full small witnesses.
It is [saved here](../results/list_to_winning_set.json).

These are ordinary proofs and exact calculations. The missing work
includes a Lean proof against the pinned upper target, any claim of a
sharp prize threshold, the general spectral-to-code reduction, and the
Paley conjectures themselves. None of those is established by this result.
