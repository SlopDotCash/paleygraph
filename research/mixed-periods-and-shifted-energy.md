# Mixed period identities and shifted multiplicative energy

**Status: no uniform Paley or prize bound is proved.** The mixed identities
retain more arithmetic information than the single squared-period identity.
Here we derive their matrix representation and an exact collision identity,
then certify two larger finite examples. These are tools for the ongoing
investigation, not a completed reduction or a novelty claim.

The multiplication-matrix framework is classical; see
[Hoshi–Kanai](https://arxiv.org/html/2105.14872), introduction and §2.3.
The related symmetric matrix and circularity framework appear in
[Garcia–Lorenz–Todd](https://arxiv.org/html/2112.13886), §§2–3. We give the
arguments below in the workspace's notation to make the assumptions and
normalizations explicit.

## The full mixed identities

Let p be an odd prime, H≤F_p* have even order n, and m=(p−1)/n. Choose
a primitive root g, and put H_t=g^tH, with indices modulo m. Define

\[
 C_{ts}=\#\{z\in H_t:1+z\in H_s\},\qquad
 x_j=\sum_{h\in H}e_p(g^j h).
\]

All x_j are real. Writing a pair (u,v)∈H×H_t as (u,zu) proves the
group-ring identity, and hence

\[
 x_jx_{j+t}=n\,1_{t=0}+\sum_s C_{ts}x_{j+s}.
 \tag{1}
\]

The exceptional z=−1 occurs precisely for t=0. At t=0 the row C₀s
is the kernel κ_s already studied in [coset-coherence.md](coset-coherence.md).
The changes z↦−1−z and z↦1/z give

\[
 C_{ts}=C_{st}=C_{-t,s-t},\qquad
 \sum_s C_{ts}=n-1_{t=0},\qquad
 \sum_t C_{tt}=n-1.
 \tag{2}
\]

The last equality follows from `C_tt=C_(-t),0` and the row sum.
In particular C is a sparse symmetric nonnegative integer matrix.

Put e₀=(1,0,…,0)ᵀ and let 1 denote the vector of m ones. The symmetric
matrix

\[
 T=\begin{pmatrix}0&\sqrt n\,e_0^T\\\sqrt n\,e_0&C\end{pmatrix}
 \tag{3}
\]

has eigenvector `(1,√n·1)` of eigenvalue n. For each j it has eigenvector
`(√n,x_j,x_(j+1),…,x_(j+m−1))` of eigenvalue x_j, by (1). The earlier
autocorrelation identity `Σ_t x_(j+t)x_(k+t)=p·1_(j=k)−n` shows that
these m+1 vectors are orthogonal, each of squared norm p. Thus the
complete spectrum of T is n and the m nonprincipal periods, each counted
once. Removing the rank-one projection on the first eigenvector removes
the principal frequency exactly. This is consistent with the existing
quotient-convolution computations, not a new bound on the spectrum.

There is also an integer representation

\[
 L=C-ne_0\mathbf1^T.
 \tag{4}
\]

Its eigenvectors are the shifted period vectors, since their coordinate
sum is −1. Their Gram matrix is `pI−nJ`, whose eigenvalues are p and 1,
so they form a basis. Consequently L's eigenvalues are exactly the x_j.

This establishes a useful completeness statement. Suppose a complex
vector y satisfies `Σy_t=−1` and all m equations

\[
 y_0y_t=n\,1_{t=0}+\sum_s C_{ts}y_s.
 \tag{5}
\]

Then `Ly=y₀y`, so y₀ is an actual period. The periods are distinct:
equality of two distinct coset sums would give a nonzero polynomial of
degree at most p−1, constant coefficient zero, vanishing at ζ_p. Its
minimal polynomial `1+X+⋯+X^(p−1)` makes that impossible. The minimal
polynomial follows from Eisenstein after X↦X+1. Simplicity of the
eigenvalue and the normalization `Σy=−1` now force y to be the
corresponding shifted period vector.

Thus the full mixed system has no extra normalized solutions. Proving
a bound on its solutions would prove the period bound, but solving or
restating this exact system does not itself supply that bound.

## Global double and triple collisions

Write `(c)_r=c(c−1)⋯(c−r+1)`. Partition F_p\{0,−1} into the cells
`S_ts={z:z∈H_t,1+z∈H_s}`. Then

\[
 \sum_{t,s}C_{ts}=p-2,\qquad
 \sum_{t,s}(C_{ts})_2=(n-1)(n-2).
 \tag{6}
\]

For the second identity, an ordered pair of distinct points x,y in a
cell gives `a=y/x∈H\{1}`, `b=(1+y)/(1+x)∈H\{1}`, with a≠b. Conversely
each such (a,b) gives uniquely `x=(b−1)/(a−b)`, `y=ax`; neither point is
0 or −1. There are (n−1)(n−2) choices. This also gives
`ΣC_ts²=p+n²−3n`, agreeing with the spectrum of (3).

For triples, remove zero from the shifted subgroup:

\[
 R=(H-1)\setminus\{0\},\qquad
 \mathcal X(H)=E^\times(R)-\bigl(2(n-1)^2-(n-1)\bigr).
\]

Here `E×(R)=#{(r₁,r₂,r₃,r₄)∈R⁴:r₁r₂=r₃r₄}`. The subtracted term
counts the two trivial pair matchings, subtracting their overlap. The
exact identity is

\[
 \boxed{\mathcal X(H)=\sum_{t,s}(C_{ts})_3.}
 \tag{7}
\]

**Proof.** For distinct x,y,z in one cell, set

\[
 a=y/x,\quad b=(1+y)/(1+x),\quad
 c=z/x,\quad d=(1+z)/(1+x).
\]

All four belong to H\{1}; moreover

\[
 (a-1)(d-1)=(b-1)(c-1).
 \tag{8}
\]

Neither trivial matching holds, since x,y,z are distinct. Conversely,
take a,b,c,d∈H\{1} satisfying (8), excluding `a=b,d=c` and `a=c,d=b`.
The nonzero shifted factors imply a≠b, a≠c, c≠d, and b≠d. Define

\[
 x=(b-1)/(a-b),\qquad y=ax,\qquad z=cx.
\]

Equation (8) gives `(1+z)/(1+x)=d`; the other three ratio equations
follow directly. These points are distinct and avoid 0,−1, so form an
ordered triple in one cell. The constructions are inverse. ∎

Consequently 𝒳≥0 and is divisible by 6. More strongly, 𝒳=0 if and only
if every C_ts≤2. This is equivalent to circularity: any two distinct
sets aH+b,cH+d with nonzero a,c have intersection of size at most two.
If their centers agree, distinct cosets are disjoint; if the centers
differ, translating and scaling their intersection reduces to a cell
S_ts. This gives a certificate of the property without enumerating all
affine copies over the ambient field.

## What (7) gives for the necessary fourth-energy estimate

Recall `E₂(H)=n²+nΣ_s κ_s²`, and that every κ_s is even except κ at
the coset 2H, which is odd. Let q be that odd entry and set

\[
 D=\sum_s\kappa_s^2-(2n-3)
  =(q-1)^2+\sum_{s\ne2H}\kappa_s(\kappa_s-2)\ge0.
\]

For an even integer c≥0, `c(c−2)≤(c)₃/3`. For odd q≥1,
`(q−1)²≤(q)₃/3+2`, as one checks at q=1,3 and then from
`3(q−1)≤q(q−2)` for q≥5. Therefore

\[
 0\le D\le\frac13\sum_s(\kappa_s)_3+2
       \le\frac{\mathcal X(H)}3+2,
\]

and hence the unconditional comparison

\[
 \boxed{E_2(H)\le 3n^2-n+\frac n3\mathcal X(H).}
 \tag{9}
\]

When 𝒳=0, parity gives the sharper exact value `E₂=3n²−3n`.
A uniform bound `𝒳=O(n log n)` would give the necessary
`E₂=O(n² log n)`. It would still not control logarithmic-depth moments.
It is not shown here to be a necessary consequence of the desired
spectral bound, and no such uniform upper bound on 𝒳 is proved.

[Shkredov's shifted-energy result](https://arxiv.org/html/1504.04522)
gives `E×(H−1)≪n²log n` for n<√p. Its convention includes zero; our
R omits it. This known bound also bounds 𝒳 but is quantitatively too
weak in (9). The task is to control the nontrivial excess, or its
concentration in row zero, more sharply. Subtracting the trivial
solutions from a coarse upper bound does not create a power saving.

## Exact finite certificates

The script independently counts products in R and sums in H, reconstructs
every nontrivial collision as an ordered triple, checks the inverse map,
and recovers every matrix cell of size at least three. In five small
fields it independently constructs the entire matrix from residues and
checks all mixed group-ring identities. In two of those it also checks
traces of L through degree eight against direct additive convolution.

| p | n | 𝒳(H) | E₂(H) | Circular |
|---:|---:|---:|---:|:---:|
| 1049 | 8 | 0 | 168 | yes |
| 2017 | 8 | 0 | 168 | yes |
| 17393 | 16 | 0 | 720 | yes |
| 6700417 | 64 | 114 | 12864 | no |
| 67403009 | 128 | 720 | 51840 | no |
| 1073748737 | 256 | 0 | 195840 | yes |
| 17179869697 | 512 | 0 | 784896 | yes |

Every row in this table lies in `n⁴/4≤p≤n⁴`. Two additional small
checks use (p,n)=(17,8),(97,8), outside that window. The final two
subgroups have generators 1064280392 and 13395504394, respectively.
The script verifies primality by trial division and exact dyadic orders
by modular powers. No numerical eigenvalue calculation is used.

The finite circular examples do not prove circularity throughout the
window: the two resonant examples already refute that assertion. Even
the minimum fourth energy alone does not give the target maximum bound.
The full subgroup conjecture, classical Paley conjecture, and prize
reduction remain open.

Run `python3 experiments/mixed_period_collisions.py`; exact results and
all cells of size at least three are in `results/mixed_period_collisions.json`.
The three inspected primary HTML sources are archived with checksums in
`sources/mixed-periods-2026-09-04/manifest.json`. No Lean formalization or
best-current claim is made.

Subsequent work in [kernel-discriminant.md](kernel-discriminant.md)
classifies fourth-energy exceptions through order 32 by an exact integer
discriminant. [quadruple-orbits-and-cube.md](quadruple-orbits-and-cube.md)
then separates repeated-entry and four-distinct zero-sum orbits. Neither
step bounds the global shifted-energy excess or proves (SG).
