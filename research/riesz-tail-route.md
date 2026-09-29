# Large-spectrum entropy and the remaining Riesz-product estimate

The uniform prime-field subgroup bound remains open. This note audits a
direct large-spectrum argument: a standard dissociation estimate is too weak
for the requested scale, and a whole-subgroup Riesz product requires a
nontrivial normalization estimate. The simplest proposed normalization is
false in an exact quartic-window example. A corrected sufficient condition
is stated and proved equivalent to the target up to constants, so it is
not presented as an easier theorem or an unconditional solution.

The large-spectrum and Riesz-product methods are established. A primary
reference is [James R. Lee, *Covering the large spectrum and generalized
Riesz products*, v2](https://arxiv.org/html/1508.07109v2), particularly
Sections 2 and 3.2. The arguments below are self-contained, and no novelty
is claimed. We do not assume any independence of the subgroup elements.

## 1. A direct entropy bound for a dissociated subset

A set `Λ⊂F_p` is dissociated when the only relation
`Σ_{λ∈Λ} ε_λ λ=0`, with every `ε_λ∈{-1,0,1}`, is the zero relation.
Let `A⊂F_p`, `|A|=n`, and suppose

\[
\left|\sum_{x\in A}e_p(\lambda x)\right|\ge\delta n
\quad\text{for every }\lambda\in\Lambda,
\qquad 0\le\delta<1.
\]

Write `D=|Λ|`. Then

\[
\boxed{D I(\delta)\le\log(p/n),\qquad
I(\delta)=\frac{1+\delta}{2}\log(1+\delta)
          +\frac{1-\delta}{2}\log(1-\delta).}
\tag{E}
\]

Choose unit complex phases `z_λ` so that
`E_{x∈A} Re(z_λ e_p(λx))≥δ`, and set

\[
R_t(x)=\prod_{\lambda\in\Lambda}
 [1+t\operatorname{Re}(z_\lambda e_p(\lambda x))],\qquad 0<t<1.
\]

Every factor is positive. Expanding the product and averaging over F_p,
dissociation kills every nonconstant term, so `E_{F_p}R_t=1`. Consequently
`E_A R_t≤p/n` and Jensen gives `E_A log R_t≤log(p/n)`.
For `y∈[-1,1]`, concavity of the logarithm places its graph above the
chord joining the endpoints:

\[
\log(1+ty)\ge y\operatorname{artanh}(t)+\tfrac12\log(1-t^2).
\]

Summing over Λ and taking the average on A gives
`D[δ artanh(t)+(1/2)log(1-t²)]≤log(p/n)`. Taking `t=δ` proves
(E); the case δ=0 follows directly. Since `I'(δ)=artanh δ≥δ`, we also get

\[
D\delta^2\le2\log(p/n).
\tag{E2}
\]

This is a direct Chang-type estimate, not an assumption about random
phases. The endpoint δ=1 can be included by taking a limit.

## 2. Why one multiplicative orbit does not make this estimate sufficient

Let H be a multiplicative subgroup of size n and write
`η_b=Σ_{h∈H}e_p(bh)`. If `|η_b|=M`, all frequencies in bH have the same
magnitude M. Apply (E) to a dissociated subset `Λ⊂bH`, with `A=H` and
`δ=M/n`.

The sums of the `2^D` subsets of Λ must be distinct in F_p; otherwise their
difference is a forbidden signed relation. Therefore

\[
D\le\log_2 p.
\]

Even if Λ attained this maximal possible size, (E) would remain compatible
with `M=n/2` in the quartic window. Indeed,
`I(1/2)=(3/4)log3-log2<(1/4)log2`, because `27<32`. Thus
`D I(1/2)<(1/4)log p`, whereas `p≥n²` implies
`log(p/n)≥(1/2)log p`. The quartic window `n⁴/4≤p≤n⁴` has `p≥n²`
for `n≥2`.

Hence this use of (E) alone cannot yield even `M=o(n)`, much less
`M=O(sqrt(n log(p/n)))`. This is a limitation of the displayed argument,
not a claim that an actual subgroup with `M=n/2` exists asymptotically.
It does not rule out arguments combining large spectra with additional
sum-product or arithmetic information.

## 3. Using the whole subgroup: an exact sufficient condition

Now let H have even order n. Choose `H_+` containing one member of each
opposite pair `{h,-h}`, and set `d=n/2`, `m=(p-1)/n>1`. For `σ∈{1,-1}` and
`0<t<1`, define

\[
R_\sigma(b;t)=\prod_{h\in H_+}[1+\sigma t\cos(2\pi bh/p)],\qquad
Z_\sigma(t)=\frac1{p-1}\sum_{b\ne0}R_\sigma(b;t).
\]

The product is independent of the choice of H_+ and is constant on each
nonzero multiplicative coset bH. Every factor is positive. Thus
`R_σ(b;t)≤m Z_σ(t)`. Applying the same chord inequality as before, now
without any dissociation assumption, gives

\[
\boxed{\sigma\eta_b\le
\frac{2\log(mZ_\sigma(t))-(n/2)\log(1-t^2)}
 {\operatorname{artanh}(t)}.}
\tag{R}
\]

Here `Σ_{h∈H_+} cos(2πbh/p)=η_b/2`, which accounts for the constants.

In particular, suppose `n≥4log m`, and put `t=sqrt(log m/n)≤1/2`.
If a uniform constant `C≥0` satisfies

\[
Z_+(t), Z_-(t)\le\exp(Cnt^2),
\tag{RP}
\]

then `artanh t≥t` and `-log(1-t²)≤(4/3)t²` in (R) imply

\[
\boxed{M\le(2C+8/3)\sqrt{n\log m}.}
\]

For `n<4log m`, the trivial bound already gives `M≤2sqrt(n log m)`.
This establishes the sufficient implication with explicit constants.

There is also a converse at this scale. If `M≤A sqrt(n log m)`, then
`log(1+y)≤y` gives

\[
\log R_\sigma(b;t)\le\sigma t\eta_b/2
\le tM/2\le(A/2)nt^2.
\]

So (RP) follows with `C=A/2`. The single-scale condition is equivalent
to the original maximum bound up to constants. It has not removed the
unproved mathematical step.

## 4. The normalization and a strict counterexample to Z ≤ 1

In the group algebra `R[X]/(X^p-1)`, let

\[
a_0(t)=[X^0]\prod_{h\in H_+}
 \left[1+\frac t2(X^h+X^{-h})\right].
\]

Fourier orthogonality gives the exact identity

\[
\boxed{Z_+(t)=\frac{p a_0(t)-(1+t)^{n/2}}{p-1}.}
\tag{N}
\]

The coefficient of `(t/2)^j` in a_0 is the number of signed zero sums
using exactly j distinct opposite pairs, with at most one sign chosen
from each pair. All these coefficients are nonnegative. This differs
from the ordered relation counts in ordinary moments: repeated choices
and opposite-pair cancellations have been removed. Nevertheless it is
not legitimate to set `Z_+(t)≤1` automatically.

Use the previously verified quartic-window subgroup
`p=67403009`, `n=128`, `H=<64701253>`. All 128 multiples in H of

\[
(1,1074964,1550267,64777777)
\]

are distinct unordered zero quadruples with four distinct opposite pairs.
Direct enumeration of signed pairs verifies that these are exactly the
128 squarefree signed zero sums of size four, and that there are none
of size one, two, or three. Consequently

\[
a_0(t)\ge1+128(t/2)^4=1+8t^4\qquad(t>0).
\]

At `t=1/4`, (N) yields the rational lower bound

\[
Z_+(1/4)\ge
\frac{(33/32)p-(5/4)^{64}}{p-1}
>\frac{201}{200}>1.
\]

The first strict comparison is checked using integers. Even the weaker
comparison to 1 is exactly
`(p+32)4^64>32·5^64`. No evaluation of an exponential or trigonometric
sum is needed. The higher omitted coefficients can only increase the
left side at positive t.

Thus even after removing the principal frequency, this strong Riesz
normalization fails inside the target window. It does not refute (RP)
with an unspecified C, the target bound, or any prize claim. The value
`t=1/4` is the stated finite test value, not asserted to equal
`sqrt(log m/n)` for this field.

## Verification and the remaining arithmetic problem

[`riesz_tail_audit.py`](../experiments/riesz_tail_audit.py) checks prime
and subgroup orders and independently enumerates the signed zero triples
and quadruples. It matches the quadruples to the multiplicative orbit and
proves the displayed rational inequality. It also computes exact group-ring
products at both signs for `(p,n)=(97,32),(2017,8),(17393,16)`, checking
their total mass, multiplicative invariance, and principal deletion.
Greedy maximal dissociated subsets have sizes 5, 4, and 8; exhaustive
signed-relation checks and constant coefficients verify their normalization.
Results, including the 128 signed witnesses and source hash, are in
[`riesz_tail_audit.json`](../results/riesz_tail_audit.json).

These are ordinary mathematical proofs and exact finite checks, not Lean
certificates. There is still no uniform upper bound for the excess in (N)
at the critical scale. Replacing it by (RP) only restates the desired
cancellation in another form. A useful next argument must bound actual
arithmetic relations, or exploit them to control the tail, without assuming
independence or silently discarding the normalization.
