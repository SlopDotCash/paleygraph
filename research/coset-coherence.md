# What multiplicative symmetry does and does not control

The main subgroup bound remains unproved. This note derives exact
constraints on its period vector, then shows why the linear correlation
constraints alone permit large spikes. A further nonlinear identity excludes
the constructed examples. Neither the identity nor the counterexamples
resolve the original conjecture. No novelty claim is made for these results.

## An unconditional bound already available

Theorem 1.1 in
[Kurlberg's 2007 exposition](https://people.kth.se/~kurlberg/eprints/short_expsum.pdf)
of the Bourgain–Konyagin / Bourgain–Glibichuk–Konyagin argument states:
for every fixed `α>0`, there is `β(α)>0` such that

\[
|H|>p^\alpha\quad\Longrightarrow\quad
\max_{b\ne0}|\eta_b(H)|\ll |H|p^{-\beta(\alpha)}.
\]

In the proposed quartic window, choose, for example, `α=1/5`.
For sufficiently large parameters this gives a **fixed power saving**
`M≪n^(1-4β)` with an implied constant depending on the fixed window.
The pinned prior attack note's description of BGK as only `n^(1-o(1))`
is therefore inaccurate for a fixed exponent regime. This known result
still does not give the desired `M≪sqrt(n log(p/n))`. The statement is
included as a verified baseline, not a claim about the best current bound.

The source PDF is stored with its hash in `sources/manifest.json`; the
displayed theorem on page 2 was visually checked.

## Exact autocorrelations of the period vector

Let `p` be prime, `H≤F_p*` have even order `n`, and `m=(p-1)/n`.
Choose a primitive root `g` and write

\[
x_j=\eta_{g^j}=\sum_{h\in H}e_p(g^j h),\qquad j\in\mathbb Z/m\mathbb Z.
\]

These numbers are real because `-1∈H`. They satisfy

\[
\sum_j x_j=-1,\qquad |x_j|\le n,
\qquad \boxed{\sum_j x_jx_{j+s}=p\,1_{s=0}-n.}
\tag{1}
\]

For the last identity, sum `η_b overline(η_(g^s b))` over all `b∈F_p`.
Orthogonality makes this `p|H∩g^sH|`. Remove the principal contribution
`n²`, then divide by the multiplicity `n` of each coset. Since the
intersection has size `n` for `s=0` and zero otherwise, (1) follows.

Equivalently, the discrete Fourier transform of `(x_j)` over `Z/mZ`
has value `-1` at zero and absolute value `sqrt(p)` at every nonzero
frequency. This transform is the appropriate list of multiplicative
Gauss sums. Thus (1) retains all their magnitudes, not only one total
second moment. It does not retain their arithmetic phases.

## Exact vectors with these correlations and a large spike

Suppose additionally that `m` is a prime congruent to 3 modulo 4 and
`m≥n≥8`. Put

\[
a=\frac{n+1}{m},\qquad D=n+\frac{1-n^2}{m}>0,
\qquad
v_j=-a+n\,1_{j=0}+\sqrt D\,\chi_m(j),
\tag{2}
\]

where `χ_m(0)=0` is the quadratic character of `F_m`.
Here `p=nm+1`; the construction works whether or not this number is prime.
In particular it works at every prime instance meeting these hypotheses.

**Proposition.** The vector `(v_j)` has sum `-1`, satisfies exactly the
autocorrelations (1), and has `|v_j|≤n` for all `j`. Nevertheless,
`v_0=n-(n+1)/m≥n-2`: it permits a spike close to the trivial upper bound.

**Proof.** The standard quadratic-character identities give

\[
\sum_j\chi_m(j)=0,\quad \chi_m(-j)=-\chi_m(j),\quad
\sum_j\chi_m(j)\chi_m(j+s)=m1_{s=0}-1.
\]

For completeness, at nonzero `s` the last sum reduces by a change of
variables to `Σ_y χ_m(y²-1)`. The equation `z²=y²-1` has `m-1` solutions:
choose any nonzero value of `y-z`, then set `y+z` to its inverse.
Counting the same solutions as `Σ_y(1+χ_m(y²-1))` gives the value `-1`.

The cross terms between the delta at zero and the character cancel
because the character is odd. Thus the off-diagonal autocorrelation in
(2) is

\[
ma^2-2na-D=-n,
\]

and the diagonal value is

\[
ma^2-2na+n^2+(m-1)D=p-n.
\]

Also `Σv_j=-ma+n=-1`. Finally, `a<2`, `0<D≤n`, and
`sqrt(n)+2≤n` for `n≥8`. These give the claimed coordinate bounds. ∎

These are **synthetic vectors, not actual Gaussian periods**. They show
that the precise linear correlation identities allow almost maximal spikes;
any argument using these identities must bring in further arithmetic
information to exclude them. The construction is not a disproof of a
Paley conjecture. No unbounded family of simultaneous primes in the
specified dyadic quartic window is being claimed from this construction.

Exact examples with both `m` and `p` prime are:

| n | m | p=nm+1 | Synthetic spike |
|---:|---:|---:|---|
| 8 | 131 | 1049 | `8-9/131` |
| 16 | 1087 | 17393 | `16-17/1087` |
| 32 | 8219 | 263009 | `32-33/8219` |
| 64 | 65647 | 4201409 | `64-65/65647` |

All four lie in `n^4/4≤p≤n^4`. The arithmetic is exact, including the
correlations expressed over `Q(sqrt(D))`. The finite examples do not
refute an unspecified asymptotic constant.

## A further identity retaining arithmetic interaction

For each coset define the nonnegative integer

\[
k_s=\#\{h\in H\setminus\{-1\}:1+h\in g^sH\}.
\]

Then `Σ_s k_s=n-1`, and every actual period vector satisfies

\[
\boxed{x_j^2=n+\sum_s k_s x_{j+s}.}
\tag{3}
\]

Indeed, write an ordered pair in `H²` as `(u,hu)`. Its sum is
`u(1+h)`. The case `h=-1` contributes `n`; every other `h` contributes
the period at the coset of `g^j(1+h)`. This proves (3) for all `j`.
The argument is an exact group-ring identity, before any Fourier evaluation.

The inverse map `h↦h^-1` preserves the coset of `1+h` because
`1+h^-1=(1+h)/h`. Among the allowed `h`, its only fixed point is 1.
Consequently every `k_s` is even except the coefficient at the coset
`2H`, which is odd. In particular,

\[
\sum_s k_s^2\ge 2n-3.
\tag{4}
\]

Combining (3) with (1) gives two exact scalar consequences:

\[
\sum_jx_j^3=pk_0-n^2,\qquad
\sum_jx_j^4=p\left(n+\sum_s k_s^2\right)-n^3.
\tag{5}
\]

For the cubic identity, multiply (3) by `x_j` and sum using (1).
For the quartic identity, square (3), sum, and use the same correlations.
Equivalently the additive energy is

\[
E_2(H)=n^2+n\sum_s k_s^2.
\]

Thus (4) recovers `E_2≥3n²-3n`. This is consistent with the earlier
moment analysis, not a new improvement on that bound.

For the original resonant example `p=6700417`, `H=⟨2⟩`, `n=64`, the
kernel has 28 entries equal to 2, one equal to 3, and one equal to 4.
In particular `k_0=3`, `Σk_s²=137`, and `E_2=12864`.
This identifies the small set of cosets carrying the fourth-moment excess.

## The constructed spikes fail the arithmetic identity

For each of the four synthetic examples, the program also constructs
the **actual** subgroup and its kernel in `F_p`. At index zero, the
residual of (3) for the synthetic vector is

\[
v_0^2-n-\sum_s k_s v_s=A+B\sqrt D,
\]

where

\[
A=(n-a)^2-n+(n-1)a-nk_0,\qquad
B=-\sum_s k_s\chi_m(s).
\]

The output saves these rational/integer coefficients and verifies
`A²-B²D≠0`, so the residual is nonzero without a numerical square root.
Hence these vectors cannot be used as actual counterexamples to the
period conjecture. No integrality or Galois-conjugacy properties were
imposed on the synthetic vectors either.

The exact group-ring identity is checked independently by enumerating all
`n²` ordered sums in each of the five prime fields. Primality, coset order,
kernel parity, energy, and the synthetic correlation coefficients are also
checked. For the two smallest auxiliary primes, every quadratic-character
autocorrelation is directly enumerated in addition to the general proof.

```sh
python3 experiments/coset_coherence.py
```

Results and source hash are in `results/coset_coherence.json`.
The script needs only Python's standard library. The proofs in this note
have not been formalized in Lean.

## Remaining question

The period equation (3) is stronger than the linear correlation data. The
next useful question is whether its particular arithmetic kernel supports
a uniform bound on the maximum, rather than only reproducing low moments.
The present argument has **not** obtained such a bound, and has not proved
that (1)–(3) alone suffice. The logarithmic-depth estimate (SG) and the
connection to arbitrary MCA word pairs remain unresolved.

The subsequent [Riesz-product analysis](riesz-tail-route.md) audits a
direct large-spectrum use of the same multiplicative invariance. Its
dissociation entropy bound alone is too weak at the quartic scale. Using
the whole subgroup instead requires a normalization estimate that is
equivalent to the target up to constants; setting that mean to at most
one is refuted by an exact arithmetic witness. This does not show that
the full nonlinear period identities are insufficient.

The subsequent [mixed-period analysis](mixed-periods-and-shifted-energy.md)
constructs those full identities. Their normalized solutions are exactly
the shifted actual period vectors. It also proves that the total ordered
triple collisions among the cyclotomic matrix cells equal the nontrivial
multiplicative energy of `(H−1)\{0}`. This gives an exact arithmetic
quantity to investigate; its known upper bound does not close the needed
energy estimate, and the full moment target remains unproved.
