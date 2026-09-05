# A finite obstruction to the literal square-root-of-two constant

The main conjecture is still open. The result here refutes the stronger
universal inequality

\[
M\le\sqrt{2n\ln(p/n)}
\tag{L2}
\]

in the concrete window `n⁴/4≤p≤n⁴`. All logarithms are natural.
It does not refute an unspecified constant, a sufficiently-large-`n` version,
or the different inequality `M²≤2n ln p`.

## Explicit witness and analytic certificate

Let `p=6700417`, `n=64`, and `H=⟨2⟩⊂F_p*`. The subgroup arithmetic and
quartic window are checked in `SubgroupQuarticCounterexample.lean`.
At the nonzero frequency `b=1`, set

\[
\eta_1=\sum_{h\in H}\exp(2\pi i h/p).
\]

The new certificate `SharpConstantCounterexample.lean` proves

\[
\Re\eta_1>43,\qquad \ln(p/64)<12.
\]

Consequently

\[
M^2\ge |\eta_1|^2>43^2=1849>1536
>128\ln(p/64),
\]

which disproves (L2). The Lean file defines the actual complex exponential
sum, connects its real part to the cosine sum, and proves the norm inequality.
It also checks primality, cardinality, nonzero representatives, multiplication
closure, membership of 1 and 2, and the quartic window. The interpretation of
a finite multiplicatively closed nonempty subset of `F_p*` as a subgroup is
standard mathematics; that abstract identification is not encoded locally.
The verification log is `results/sharp-constant-lean-verification.txt`.
All six printed axiom censuses contain only the standard Lean axioms, and
the log records the exact proof-source hash and pinned mathlib revision.

The analytic proof uses only rational bounds. For `0≤u≤4`, define

\[
c=1-u^2/32,\quad d=2c^2-1,\quad
L(u)=\begin{cases}2d^2-1&d\ge0,\\-1&d<0.\end{cases}
\]

If `|x|≤u`, the inequality `cos(x/4)≥1-x²/32`, followed by two double-angle
identities, gives `cos x≥L(u)`. Squaring is used only when both compared
quantities are nonnegative; the other branch uses `cos x≥-1`.

Fold each angle using `t=min(h,p-h)` and take `u=63t/(10p)`.
Since `π<3.15`, `|2πt/p|≤u≤3.15<4`. Exact rational summation yields

\[
\sum_{h\in H}L\!\left(\frac{63\min(h,p-h)}{10p}\right)>43.
\]

For the logarithm, `p/64<2¹⁷` and `ln 2<0.694` imply
`ln(p/64)<17·0.694<12`.

For orientation only, floating evaluation gives
`Re η₁≈43.8024827976263` and
`Re η₁/√(64 ln(p/64))≈1.6104691799182`, whereas the right side of
(L2) is approximately `38.4646077121223`. These approximations are not used
in the certificate. `experiments/sharp_constant_witness.py` independently
recomputes the rational lower sum and clearly separates diagnostic decimals.

## Exact moments through logarithmic depth

Let `C_r(x)` count ordered `r`-term sums from `H` equal to `x`, with
`C₀(0)=1`. Multiplication by `H` preserves each `C_r`, so it suffices to
store one entry per multiplicative coset, plus zero. If `λ(x)` is the index
of that orbit and `a_i` is its representative, the exact recurrence is

\[
C_{r+1}(a_i)=\sum_{h\in H}C_r(a_i-h).
\]

This is additive convolution, compressed using a proved symmetry. It gives

\[
E_r=C_r(0)^2+n\sum_{i=1}^{m}C_r(a_i)^2,\qquad
Q_r=(pE_r-n^{2r})/n.
\]

`experiments/quotient_moments.py` checks every orbit label and cardinality,
total mass `Σ_x C_r(x)=n^r`, and an overflow bound before every integer-array
sum. Energies and the large zero-frequency subtraction use arbitrary-precision
integers. Independent dense convolutions agree on two small fields; the target's
second and third energies agree with the earlier sparse computation.

Here `m=104694` and `ceil(ln m)=12`. Integer root comparisons give:

| Moment depth r | Certified computational upper bound on M² |
|---|---:|
| 2 | 36695 |
| 4 | 4148 |
| 6 | 2520 |
| 8 | 2148 |
| 10 | 2022 |
| 12 | 1970 |

At the final depth,

\[
E_{12}=3360801151538389751318410002884550208,
\]
\[
Q_{12}=3406624544637182825373883182020746016505
\le1970^{12}.
\]

Together with `M²ʳ≤Q_r`, this establishes the finite interval
`43<M≤√1970<44.385`. The upper certificate is an exact integer computation
with a mathematical derivation; it has not been replayed in the Lean kernel.

This finite enclosure has since been sharpened using all signed moments
through order 24. The [polynomial certificate](period-polynomial-certificate.md)
proves that the unique maximizing coset is `H` and
`43.802482797626304198≤M=η₁≤43.802482797626304199`.
Every other coset has absolute period below 42. The new bound uses exact
rational arithmetic and a mathematical proof, not a new Lean certificate.

The Gaussian comparison fails at every tested order `r=2,…,12` here.
The weaker proposed estimate `Q_r≤m(2rn)^r` passes **all** orders `r=1,…,12`.
This demonstrates that the Gaussian coefficient and the weaker moment target
have different behavior even at the logarithmic depth. It supplies no uniform
bound for other fields or arbitrarily large subgroups.

## Matching the prior repository's conventions

The pinned file
[`GeneralizedPaleyRamanujan.lean`](https://github.com/SlopDotCash/proximityprize/blob/5b00e50c3c51b3c944201a1749a1a8132e5ce167/ArkLib/Data/CodingTheory/ProximityGap/GeneralizedPaleyRamanujan.lean)
uses `Real.log`, and defines a bound on the **squared** period with an
unspecified coefficient `C`. Thus its `C` equals the square of the amplitude
constant used in the local target. Its `DepthLogSubGaussian` bridge instead
has the envelope `M²≤2n ln p`; changing `ln p` to `ln(p/n)` while keeping
the coefficient 2 is not an equivalent bound.

Indeed the present example satisfies that different envelope. The exact
upper bound gives `M²≤1970`. Since `p²>2⁴⁵` and `ln 2>0.69`,

\[
2n\ln p>128\cdot(45/2)\cdot0.69=1987.2>1970.
\]

The prior essay's observed constants below `√2` came from its stated sampled
prime windows. This witness does not dispute those reported samples; it shows
that extrapolating them to every prime in the broader quartic window fails.
No sponsor-level implication or counterexample is claimed.

## Remaining mathematical obligation

A uniform estimate such as `Q_r≤m(Krn)^r`, with one absolute `K` and
`r≈ln m` across the precisely quantified target family, remains unproved.
The finite resonance must be accommodated, but it does not eliminate this
route. The next step must estimate the centered convolution growth uniformly;
further finite examples alone cannot establish that statement.
