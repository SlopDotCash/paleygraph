# A bound on extra relations, summed over primes

This is an unconditional elementary bound for a fixed dyadic subgroup
order, summed over the prime fields in which that subgroup exists. It
does not give the required worst-case spectral estimate. No novelty or
Lean-formalization claim is made.

## Statement

Fix a power of two `n≥2` and an integer `k≥2`. Let `T_k(n)` be the
number of ordered k-tuples of complex n-th roots of unity whose sum is
zero. For a prime `p≡1 mod n`, let `Z_k(p,n)` be the corresponding
zero-sum count for `μ_n⊂F_p*`. Then `Z_k(p,n)≥T_k(n)` and

\[
\boxed{\sum_{\substack{p\ \mathrm{prime}\\p\equiv1\pmod n}}
 [Z_k(p,n)-T_k(n)]\log p
 \le [n^k-T_k(n)]\log k.}
\tag{N}
\]

For fixed n and k only finitely many summands are nonzero. No assertion
about the number of primes in a progression is needed for this theorem.

## Proof

Put `N=n/2` and let ζ be a primitive complex n-th root. Its minimal
polynomial is `Φ_n(X)=X^N+1`: irreducibility follows from Eisenstein at
2 after `X↦X+1`, since N is a power of two. For each exponent tuple
`a=(a₂,…,a_k)∈{0,…,n-1}^(k-1)`, define

\[
f_a(X)=1+X^{a_2}+\cdots+X^{a_k}.
\]

Call a intrinsic if `f_a(ζ)=0`. This is equivalent to its reduction
modulo `X^N+1` being zero as an integer polynomial. In particular it
then vanishes at every primitive n-th root in every eligible prime
field. Multiplying every entry of a zero tuple by the inverse of its
first entry shows that the number of intrinsic normalized tuples is
`T_k(n)/n`.

For a nonintrinsic tuple, let `D_a` be the determinant of multiplication
by `f_a` on the free integer module `Z[X]/(X^N+1)`. Over C this map
diagonalizes by evaluation at the primitive roots, so

\[
D_a=\prod_{u\in(\mathbb Z/n\mathbb Z)^*}f_a(\zeta^u)
\in\mathbb Z\setminus\{0\},\qquad |D_a|\le k^N.
\]

For `p≡1 mod n`, `X^N+1` splits into N distinct linear factors modulo
p. Therefore the multiplication matrix has nullity

\[
t_{a,p}=\#\{g\in\mathbb F_p^*: \operatorname{ord}(g)=n,
                                      \ f_a(g)=0\}.
\]

An integer matrix with nullity t modulo p has determinant divisible by
`p^t`. For example, lift row elimination modulo p to integer operations
whose determinants are prime to p; the last t rows become divisible by
p, which gives the assertion on determinant valuation. Thus

\[
\sum_{p\equiv1\ (n)}t_{a,p}\log p\le\log|D_a|\le N\log k.
\]

For each fixed primitive root g in `F_p`, exponent tuples enumerate
every ordered tuple with first entry 1 exactly once. The intrinsic
subset is the same under every primitive choice. Consequently

\[
\sum_{a\ \mathrm{nonintrinsic}}t_{a,p}
=\frac{N}{n}[Z_k(p,n)-T_k(n)].
\]

There are `(n^k-T_k(n))/n` nonintrinsic normalized tuples. Sum the
determinant estimate over them and cancel `N/n` to obtain (N).
The finite product of the nonzero determinants also proves that only
finitely many primes can contribute. ∎

## Explicit intrinsic counts and an energy consequence

The Q-basis `1,ζ,…,ζ^(N-1)` shows that a zero sum must use each root
and its opposite equally often. Thus odd intrinsic counts vanish, and

\[
T_{2r}(n)=(2r)![t^{2r}]
\left(\sum_{j\ge0}\frac{t^{2j}}{(j!)^2}\right)^{n/2}
\le(2r-1)!!\,n^r.
\]

The inequality counts perfect matchings of the 2r positions, assigning
opposite roots to each matched pair. Every intrinsic tuple has at least
one such matching; repeated coverage is allowed in an upper bound.
In particular `T₂=n` and `T₄=3n²-3n`.

Because `-1∈μ_n`, `Z₄(p,n)=E₂(μ_n)`. For any real `L>1` and `t>0`,
(N) implies

\[
\#\{p\ge L:p\equiv1\pmod n,\quad
E_2(\mu_n)>3n^2-3n+tn^2\}
\le \frac{(n^4-3n^2+3n)\log4}{tn^2\log L}.
\tag{9}
\]

Taking `L=n⁴/4` bounds the number of exceptional primes in the quartic
window by `O(n²/(t log n))`. This is an absolute count of possible
exceptions. Without a separate lower bound for the number of primes
in that window and progression, it is not a proved density statement.
It also leaves individual exceptional primes entirely possible.

At order `2r`, the numerator in (N) is of size `n^(2r) log(2r)`.
This growth is too large to control the desired centered moments of
size `(Krn)^r` at `r≈log((p-1)/n)` by this argument. Moreover, (N)
counts extra relations above the complex intrinsic count; it does not
subtract the prime-field principal-frequency term `n^(2r)/p`.
These are the explicit limitations of this route.

## A sharper norm average, and its pointwise limitation

The product estimate can be sharpened using the arithmetic-geometric mean
inequality. Put `R=(n^k-T_k(n))/n`, the number of nonintrinsic normalized
tuples. For any fixed primitive complex root ζ, orthogonality gives

\[
\sum_a |1+\zeta^{a_2}+\cdots+\zeta^{a_k}|^2
=k n^{k-1}.
\]

Indeed the k diagonal terms each sum to `n^(k-1)`, while all cross
terms sum to zero. Intrinsic tuples contribute zero. Applying the
arithmetic-geometric mean inequality to the R remaining positive
squared absolute values, then multiplying over the N embeddings, gives

\[
\prod_{a\ \mathrm{nonintrinsic}}|D_a|
\le\left(\frac{k n^{k-1}}R\right)^{NR/2}.
\]

The same divisibility argument therefore proves the stronger estimate

\[
\sum_{p\equiv1\ (n)}[Z_k(p,n)-T_k(n)]\log p
\le\frac{n^k-T_k(n)}2
\log\frac{k n^k}{n^k-T_k(n)}.
\tag{N'}
\]

For k=4, this still cannot improve the elementary pointwise energy bound
`E₂≤n³` in the quartic window once `n≥16`. To see this without decimal
approximations, write `n=2^s`, `s≥4`, and `T=T₄(n)`. Keeping only one
prime in (N') gives the upper-bound expression

\[
U=T+\frac{n^4-T}{2\log p}\log\frac{4n^4}{n^4-T}.
\]

Since `p≤n⁴` and the logarithm in the numerator is greater than `log4`,

\[
U>T+\frac{n^4-T}{4s}\ge\frac{n^4}{4s}\ge n^3.
\]

The last inequality follows from `2^s≥4s` for `s≥4`, by induction.
Thus this explicit pointwise consequence is weaker than the elementary
bound. This is a limitation of (N') as used here, not an impossibility
theorem for every method involving cyclotomic norms.

## At higher moments, extra relations occur at every eligible prime

The sparse exceptional-prime picture at order four cannot be extended
unchanged to the logarithmic moment depth. Cauchy–Schwarz applied to the
r-fold additive convolution gives

\[
Z_{2r}(p,n)=E_r(H)\ge\frac{n^{2r}}p.
\]

For `p≤n⁴`, comparison with the intrinsic pairing upper bound yields

\[
Z_{2r}(p,n)-T_{2r}(n)
\ge n^{2r-4}-(2r-1)!!n^r.
\tag{10}
\]

In particular, for every dyadic `n≥1024` and every prime `p≡1 mod n`
with `p≤n⁴`,

\[
Z_{10}(p,n)-T_{10}(n)\ge n^5(n-945)>0.
\]

This is a universal obstruction to treating all extra high-order relations
as rare prime exceptions. It does not contradict square-root cancellation:
the lower bound is supplied by the principal-frequency contribution that
the desired centered estimate must remove. The next argument must control
fluctuations above that contribution, rather than merely exclude primes
with additional relations.

## Exact finite audit

[`cyclotomic_norm_audit.py`](../experiments/cyclotomic_norm_audit.py)
enumerates all normalized quadruples for `n=4,8,16`, computes the integer
multiplication determinants, factors their product, and independently
checks every relevant matrix nullity against primitive-root evaluations.
Pair-sum enumeration independently gives each prime-field energy.
The required prime product divides the norm product, and every
determinant satisfies the stated integer upper bound.

The complete lists of primes `p≡1 mod n` with extra zero quadruples are:

| n | Primes with `E₂>3n²-3n` |
|---|---|
| 4 | 5 |
| 8 | 17, 41 |
| 16 | 17, 97, 113, 193, 257, 337 |

Completeness follows from factoring every nonzero determinant, not a
search over a bounded interval of primes. Results and source checksum
are saved in [`cyclotomic_norm_audit.json`](../results/cyclotomic_norm_audit.json).

## Relation to the recently inspected polynomial results

[Yip–Yoo, arXiv:2608.02568v1, Theorem 1.1 and Proposition 2.7](https://arxiv.org/html/2608.02568v1)
concern exact decompositions `H=A+B`. Their polynomial argument first
proves `|A||B|≤|H|`, then uses the opposite inequality from the equality
of sets to force exact factorization and unique representations.
Our collision counts do not supply such a decomposition. Applying that
factorization here without an additional theorem would be unjustified.

[Kim–Yip–Yoo, arXiv:2309.09124v4, Theorem 1.1](https://arxiv.org/html/2309.09124v4)
bounds a rectangle satisfying `AB+λ⊆S_d∪{0}`. Counting many pairs
`h₁+h₂=h₃+h₄` does not verify that rectangle hypothesis. No bridge
from these statements to the required energy or recurrence bounds has
been established here. Their Conjecture 2.12 also quantifies over every
nontrivial multiplicative character, broader than the specifically
quadratic formulation separately tracked in this workspace.

The later [kernel-discriminant calculation](kernel-discriminant.md)
independently recovers the complete order-4,8,16 classifications above
and extends them to order 32 using one polynomial. Its valuation formula
and the [orbit factorization](quadruple-orbits-and-cube.md) still do not
turn the elementary height bound into adequate pointwise control.
