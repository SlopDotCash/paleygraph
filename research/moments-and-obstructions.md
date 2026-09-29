# Character moments: proved identities and failed proof routes

Throughout, `p` is an odd prime, `χ(0)=0`, and `χ` is the quadratic character.
Write `m=|A|`, `n=|B|`, and

\[
F_B(x)=\sum_{b\in B}\chi(x-b),\qquad
S(A,B)=\sum_{a\in A}F_B(a),\qquad
M_{2r}(B)=\sum_{x\in\mathbb F_p}F_B(x)^{2r}.
\]

All proofs below are ordinary mathematical proofs, except the explicit finite
obstruction also checked in the accompanying Lean file. They do not prove the
main conjecture. Standard estimates are identified as such, and no novelty
claim is made for the identities or obstruction arguments.

## 1. Exact second moment and the standard size barrier

For distinct `u,v`,

\[
\sum_x\chi(x-u)\chi(x-v)=-1.
\]

Indeed, for `x ≠ v`, the map `t=(x-u)/(x-v)` is a bijection onto
`F_p \ {1}`. The product has character `χ(t)` because `(x-v)^2` is a
nonzero square. The missing term `x=v` is zero, and
`Σ_{t≠1} χ(t)=-1`. For `u=v`, the sum is `p-1`.
Consequently,

\[
\boxed{\sum_x F_B(x)=0,\qquad M_2(B)=pn-n^2.}
\]

Apply Cauchy–Schwarz to the centered indicator `1_A-m/p` and `F_B`:

\[
|S(A,B)|^2\le
\frac{m(p-m)n(p-n)}p\le pmn.
\]

Thus `|S|/(mn) ≤ √(p/(mn))`. If `m,n > p^ε` and `ε>1/2`, this proves
the target estimate with `δ=ε-1/2`. It gives no power saving at or below
`ε=1/2`. The exact centered inequality is verified in squared integer form
in the experiments, including empty and full sets.

This is the classical spectral / second-moment estimate, not new progress
on the small-set conjecture.

## 2. What termwise Weil estimates give at higher moments

Expanding `M_{2r}` gives a sum over ordered `2r`-tuples from `B`. Tuples in
which every element occurs an even number of times contribute at most `p`
each. Their number is at most `(2r-1)!! n^r`: pair the positions, then
assign an element of `B` to each pair. This overcounts, which is harmless.

For any remaining tuple, its polynomial
`∏_{j=1}^{2r}(X-b_j)` is not a square, so the standard multiplicative
character version of Weil's bound gives an absolute complete sum at most
`(2r-1)√p`. Hence

\[
M_{2r}(B)\le (2r-1)!!\,p n^r+(2r-1)\sqrt p\,n^{2r}.
\]

Hölder, followed by enlargement of the outer sum from `A` to `F_p`, yields

\[
\frac{|S(A,B)|}{mn}
\le \left(\frac{(2r-1)!!p}{mn^r}
             +\frac{(2r-1)\sqrt p}{m}\right)^{1/(2r)}.
\]

Taking larger fixed `r` can make the first term small. The second still
requires `m>√p` by a power of `p`. Interchanging `A,B` does not solve the
balanced small-set case. Letting `r` grow does not turn a base ≥1 into a
quantity smaller than 1.

This recovers the mechanism behind the classical Karatsuba range:
`m>p^(1/2+η)` and `n>p^η`. See Appendix A of
[Fouvry–Shparlinski–Xi (2025)](https://doi.org/10.4171/RMI/1530).
It is a limitation of this particular estimate, not a theorem excluding
all possible uses of higher moments.

## 3. An exact fourth-moment decomposition

Define

\[
E_B=\sum_{a\in B}F_B(a)^2,\quad
K(Q)=\sum_x\prod_{b\in Q}\chi(x-b),\quad
T_4(B)=\sum_{Q\in\binom B4}K(Q).
\]

Then

\[
\boxed{M_4(B)=p(3n^2-2n)-6n^3+14n^2-9n-6E_B+24T_4(B).}
\]

**Proof.** Partition the ordered quadruples by multiplicities:

| Multiplicity pattern | Total contribution |
|---|---|
| `4` | `n(p-1)` |
| `2+2` | `3n(n-1)(p-2)` |
| `3+1` | `-4n(n-1)` |
| `2+1+1` | `-6n(n-1)(n-2)-6E_B+6n(n-1)` |
| `1+1+1+1` | `24T_4(B)` |

For the fourth row, let `a` be the repeated element and let `{b,c}` be
the unordered singleton pair. Its complete sum is
`-1-χ(a-b)χ(a-c)`, since the square factor deletes the term `x=a`
from the two-point correlation. The ordering multiplicity is 12. Also,
`Σ_{b<c; b,c≠a} χ(a-b)χ(a-c) = (F_B(a)^2-(n-1))/2`.
Summing the five rows proves the formula.

**Remaining issue.** This identity isolates the signed aggregate `T_4(B)`;
it does not bound it. Applying absolute Weil bounds separately to each
quartic returns the same square-root barrier. Neither positivity of `M_4`
nor the identity permits discarding `T_4`.

## 4. Why a uniform Gaussian moment bound is false

Consider the tempting claim, for a fixed integer `r≥2`,

\[
M_{2r}(B)\le C_r p n^r\quad\text{for every }B.
\tag{G}
\]

Choose any `n` elements `q` that are nonzero squares and put `B={-q}`.
Then `F_B(0)=n`, so `M_{2r}(B)≥n^{2r}`. Therefore (G) would require
`n^r≤C_r p`. Taking `n=⌊p^α⌋` with `α>1/r` contradicts this for
arbitrarily large primes. For `r≥3`, one may take `1/r<α<1/2`, so
the failure occurs even below the square-root set-size threshold.

**Concrete Lean certificate.** At `p=1009`, take
`B={-j² mod 1009: 1≤j≤30}`. These representatives are distinct,
`n=30`, and `n²=900<1009`. The zero row alone gives

\[
M_8(B)\ge30^8=656100000000
>105\cdot1009\cdot30^4=85815450000.
\]

Here `105=7!!` is the Gaussian pairing coefficient. The actual Legendre
symbol, primality of 1009, cardinality, row value, and moment inequality are
handled in `MomentObstruction.lean`. The moment lower bound uses a single
nonnegative summand; it does not rely on a numerical evaluation of all rows.

This finite certificate refutes the displayed candidate estimate. The
preceding asymptotic argument refutes every possible constant `C_r`.
Neither argument is a counterexample to Paley: the forced row set has size 1.

## 5. Many spikes: an elementary counting obstruction

For a row set `U`, put

\[
N(U)=\{b\in\mathbb F_p:\chi(a-b)=1\text{ for all }a\in U\}.
\]

Let `d=(p-1)/2`. Every column has exactly `d` positive entries. Double
counting pairs `(U,b)` with `|U|=k` and `b∈N(U)` gives the exact identity

\[
\sum_{U\in\binom{\mathbb F_p}k}|N(U)|=p\binom dk.
\]

Thus some `U` has

\[
|N(U)|\ge p\frac{\binom dk}{\binom pk}.
\tag{1}
\]

For any integer `n` at most the right side, choose an `n`-element subset
`B⊆N(U)`. All `k` rows indexed by `U` equal `n`, and therefore

\[
\boxed{M_{2r}(B)\ge k n^{2r}.}
\tag{2}
\]

This uses only the number of positive entries per column, not Weil's bound.
It also applies to `p≡3 (mod 4)`, with the indicated orientation.

### A constant allowance for spikes still fails

The proposed repair

\[
M_{2r}(B)\le C_r(pn^r+n^{2r})
\tag{G'}
\]

is also false for every fixed `C_r` and `r≥2`. Choose a fixed integer
`k>C_r`. For fixed `k`, the right side of (1) is `(2^{-k}+o(1))p`.
Choose `n=⌊c_k p⌋` for any fixed `0<c_k<2^{-k}`. Equations (1)–(2) give

\[
\frac{M_{2r}(B)}{pn^r+n^{2r}}
\ge\frac{k}{1+p/n^r}\longrightarrow k>C_r.
\]

For `r≥3`, the same contradiction works with `n=⌊p^α⌋` and
`1/r<α<1/2`. Removing any fixed number of rows also fails to fix (G):
choose more forced rows than the number removed.

### A logarithmic allowance is necessary in a relevant range

For `1≤k≤d`,

\[
\frac{\binom dk}{\binom pk}
=2^{-k}\prod_{j=0}^{k-1}\left(1-\frac{j+1}{p-j}\right)
\ge2^{-k}\left(1-\frac{k^2}{p-k+1}\right).
\]

The last step uses `∏(1-u_j)≥1-Σu_j` for `0≤u_j≤1`.
If `k²≤(p-k+1)/2`, (1) is at least `p/2^{k+1}`.
Consequently one may choose

\[
k=\lfloor\log_2(p/n)\rfloor-1
\]

whenever `k≥1`, `k≤d`, and the displayed size condition holds.
For fixed `0<α<1` and `n=⌊p^α⌋`, these conditions hold for large `p`,
and `k=(1-α)log₂p+O(1)`. Hence there exist such sets with

\[
M_{2r}(B)\ge((1-\alpha)\log_2p+O(1))n^{2r}.
\]

If `αr>1`, the `pn^r` term is negligible relative to `n^{2r}`.
Any uniform bound of the form `C_r pn^r+D_r(p,n)n^{2r}` then needs
`D_r(p,n)` to grow at least logarithmically along these parameters.
This is a necessary lower bound, not a matching upper bound.

## 6. Two scope checks

**Prime powers cannot be silently substituted.** In `F_{r²}`, every
nonzero element of the subfield `F_r` is a square, since
`(r²-1)/2=(r-1)(r+1)/2` is divisible by `r-1`. For `A=B=F_r`,
`S(A,B)=r(r-1)`. This rules out the unrestricted prime-power analogue
for `ε<1/2`. It says nothing against the prime-field conjecture.

**A full proof would imply very small quadratic nonresidues.** If all
integers `1,...,2H-1` were residues modulo `p`, the sets
`A={1,...,H}` and `B={0,-1,...,-(H-1)}` would have `S(A,B)=H²`
(take `2H<p`). Choosing `H>p^ε` contradicts Paley's conclusion for
large `p`. Thus the conjecture implies a `p^{o(1)}` bound for the least
positive quadratic nonresidue. The elementary estimates above do not
establish that consequence either.
