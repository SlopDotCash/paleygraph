# Separate-agent review of the second Krylov coefficient and leakage

**Verdict: no mathematical correction or unresolved source-dependent
hypothesis found in the stated fixed-subspace result.** The asymptotics
for \(a,c\), the two-dimensional Rayleigh bound, and the norm bound
with leakage retained follow from the cited published sheaf results
and the supplied exact identities. This is not an all-vector spectral
bound or a proof of any Paley or prize target.

I inspected the pinned primary texts, reopened both primary PDFs,
checked the specialized hypotheses, and independently recalculated
the normalization, local obstructions, cohomology dimensions, anchor
terms, and constants. I also performed a separate exact finite check
in four actual fields, including \(p=1201\), outside the author's
recorded range. This is agent review, not human refereeing or formal
verification. Only this review file was written.

## Reviewed bytes

| Input | SHA-256 |
| --- | --- |
| [Proof](parallel24-spectral-operator-2026-09-05.md) | `7abbf64aeaa1804d4cbad5e9daa1f1a4dd90945919d05d6d3f5e4fd8fbe1b857` |
| [Verifier](../experiments/parallel24_spectral_operator_2026_09_05.py) | `42f1d718cd12915e4aa62b11c4e344d91fc6b826b0cf3ac2b6da04579e6be0dd` |
| [Results](../results/parallel24_spectral_operator_2026_09_05.json) | `02ca9d465bb403a7e17d264d57b017922c0ee33d6322098211343763902944e5` |
| [Pass-23 coupling](parallel23-spectral-coupling-2026-09-05.md) | `6d8cd29185b0c86b33961a340d11f5b2313697f2950172807678fbf69c568ed2` |
| [Earlier Jacobi calculation](parallel21-spectral-next-input-2026-09-05.md) | `2bc5a06bebfd1cdb2cd7f0d61d9f95c8448fe7bd64ace017f37edbaa7093df98` |
| [Earlier seeded result](parallel6-seeded-kernels-2026-09-04.md) | `33333657fe78eae92a70b5a2ea43a2cdc39c45c63f29f434214efa9b7c3b50c6` |
| `sources/katz-g2-hypergeometric.pdf` | `0bf485bcc9dde2afebc8268af679566bcff65aa6d6ed8c9f4b0387d7e9e954c1` |
| `sources/katz-g2-hypergeometric.txt` | `b627860e16214478e5565effcf4affd4cebb0e1f58b6f51f34005ccbd2123a40` |
| `sources/katz-gauss-kloosterman-monodromy.pdf` | `8711424f8edb14f38c0e61606ef49d6efc67d5a93ac9f5bd9fa3e13b7c7d9ffe` |
| `sources/katz-gauss-kloosterman-monodromy.txt` | `8f78cbd117172142014520943769a319d5180c610e0f0bd24861efe6fa6716bc` |

Every result-file input hash matches the current bytes. The result
explicitly records its post-run metadata refresh; it does not claim
that the original long numerical run was repeated after that refresh.
The current verifier parses as Python. My finite check below is a
separate calculation, not a rerun of that script.

## 1. Rank and weight normalization

[Katz, *G₂ and hypergeometric sheaves*, §2, printed pp. 3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
allows repeated characters within each list, provided the upstairs
and downstairs lists are disjoint. For type \((r,r)\), it supplies
geometric irreducibility, rank \(r\), weight \(2r-1\), and tame
local descriptions. The trace sign is \((-1)^{2r-1}=-1\).
Thus the specialization to trivial upstairs and quadratic downstairs
characters is valid for both \(r=2,3\) at every prime in this task.

Here is the normalization directly in the problem's convention.
In the defining raw sum, put \(x_i=u_i y_i\); then \(\prod_i u_i=t\).
For each \(u\), summation over \(y\ne0\) gives

\[
 \sum_{y\ne0}\psi((u-1)y)\chi(y)
 =G_p\chi(u-1)=G_p\chi(1-u),
\]

including zero when \(u=1\). The last equality uses
\(p\equiv1\pmod4\). Therefore the raw sum is \(G_p^r\) times
the \(r\)-fold multiplicative convolution defining \(k_r\).
Twisting once over \(\mathbb F_p\) by the constant Frobenius scalar
\(G_p^{-r}\) gives base-field trace \(-k_r\) and weight \(r-1\),
since \(|G_p|=\sqrt p\). It changes neither geometric irreducibility
nor inertia.

The phrase “trace \(-k_r\)” should be read here on
\(\mathbb F_p\)-points. Over a degree-\(e\) extension, the fixed
constant twist contributes \(G_p^{-re}\); it is not a newly selected
factor \(G_{p^e}^{-r}\). The exact extension trace consequently
contains the ratio \(G_{p^e}^r/G_p^{re}\). Its absolute value is one,
so purity has the stated normalization. The argument uses base-field
trace sums and does not assume that this extension ratio is always
one. The proof's explicit fixed-twist qualification is sufficient.

## 2. The invariant and coinvariant tests are complete

The local data specialize to \(J_r(1)\) at zero and \(J_r(\chi)\)
at infinity. For rank three, the pseudoreflection at one has
eigenvalues \((-1,1,1)\): its determinant is \(\chi^3=\chi\).
Its fixed space has dimension two, so this is a semisimple reflection,
not a length-three unipotent block. These are the local forms supplied
by Katz's endpoint and pseudoreflection descriptions cited above.

For the twist \(\chi(t)^e\chi(t-1)^f\), local scalar characters at
zero, one, and infinity have exponents \(e,f,e+f\) modulo two.
This checks all four cases without assuming that a nonconstant trace
function automatically defines a geometrically nontrivial twist.

An invariant in a tensor of two irreducible factors gives a nonzero
map from the dual of one factor to the other; a coinvariant gives the
dual Hom condition. Either nonzero map is an isomorphism. Consequently:

- Ranks two and three cannot pair into a constant factor, with any of
  the four twists.
- For a same-rank square and a nontrivial twist, if \(e=1\), the
  zero-endpoint eigencharacter changes from trivial to quadratic.
  If \(e=0,f=1\), the infinity eigencharacter changes from quadratic
  to trivial. Neither can match the dual factor. This also covers
  \(e=f=1\) through the first test.
- A single twisted rank-three sheaf remains geometrically irreducible
  of rank three and has no constant subobject or quotient.

For the three parameter maps \(t,1-t,1-1/t\), the preimage of the
reflection point is respectively \(1,0,\infty\). Each map is a
degree-one automorphism of the same three-punctured curve. At the
reflection location of one pullback, either distinct pullback has
a single length-three Jordan block. Dualizing and applying a scalar
Kummer twist preserve Jordan block lengths. The block partition
\((1,1,1)\) therefore cannot match \((3)\), even where a twist is
ramified. This excludes every twisted cross-pullback isomorphism.

No assertion that the entire normalized sheaf is arithmetically
self-dual is needed for these tests. The dual's rank and local block
data suffice; inverse trivial and quadratic eigencharacters equal
the originals. No independence or large-monodromy assumption is
being substituted for the explicit local tests.

## 3. Compact-support cohomology gives the claimed constants

[Katz, *Gauss sums, Kloosterman sums, and monodromy groups*,
§§2.3.1–2.3.3 and the weight step in §3.6](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf)
provides the Euler-characteristic and trace formulas and the
compact-support weight bound used here. On
\(U=\mathbb P^1\setminus\{0,1,\infty\}\), \(\chi_c(U)=-1\).
The reviewed sheaves are lisse and tame there. Compactly supported
degree-zero sections vanish on this nonproper connected curve;
the absence of geometric coinvariants gives \(H_c^2=0\).
Thus a rank-\(R\) tensor has \(\dim H_c^1=R\). A sheaf of weight
\(w\) contributes at most \(R p^{(w+1)/2}\) in absolute value.

Applied to the factors just checked, this gives the following full
\(U(\mathbb F_p)\)-sum constants:

| Trace product, with the specified allowed twists | Rank | Weight | Bound |
| --- | ---: | ---: | --- |
| \(k_3\) | 3 | 2 | \(3p^{3/2}\) |
| \(Lk_3\) | 6 | 3 | \(6p^2\) |
| \(L^2\), nontrivial twist | 4 | 2 | \(4p^{3/2}\) |
| \(k_3^2\), nontrivial twist | 9 | 4 | \(9p^{5/2}\) |
| Two distinct rank-three pullbacks | 9 | 4 | \(9p^{5/2}\) |

The curve has exactly the three omitted geometric points throughout;
the degree-one pullbacks create no additional ramification. The mask
on \(U(\mathbb F_p)=\mathbb F_p\setminus\{0,1\}\) is an average
of its four displayed twists, so averaging does not enlarge those
constants. There is no use of a boundary-stalk purity estimate to
silently restore zero or one. Untwisted diagonal sums are handled
separately by exact arithmetic, not by the zero-coinvariant argument.

These steps use published cohomological theorems as imported inputs.
I found no further unproved hypothesis that requires labeling the
fixed-subspace conclusions conditional on a new conjecture.

## 4. Missing anchors and explicit error constants

I independently obtained \(SL=pf_0+1\). In the inverse substitution
for \(S(D_0D_1L)(t)\), the missing parameter zero contributes
\(+1\), giving \(\chi(t)k(1-1/t)+1\). The two anchor corrections
in the exact neighborhood mask contribute another \(+1\) on \(C\),
because \(L(0)=L(1)=-1\). Hence

\[
 S_CL=(p+6+K)/4,
\]

with the constant six exactly as written. The symmetry substitutions
permute \(C\), and invariance of \(g\) gives the factors three in
\(\sum K=3s_0\) and \(\sum gK=3s\). This is valid on short
orbits as well; no orbit-size division is involved.

The complete quadratic correlation gives
\(\sum_U L^2=p^2-2p-3\). For the multiplicative convolution,
the Jacobi Mellin eigenvalue has absolute value one at its two
exceptional characters and \(\sqrt p\) at the other \(p-3\).
Thus Parseval gives

\[
 \sum_{t\ne0}k(t)^2=\frac{(p-3)p^3+2}{p-1}
 =p^3-2(p^2+p+1).
\]

The restriction to \(U\) subtracts \(k(1)^2\). Here
\(k(1)=\sum f_0L=U_J\), and the self-contained earlier Jacobi
calculation gives \(|U_J|\le2p\). Zero is already absent from the
Mellin sum. These facts validate both diagonal main terms in (11).
The six cross terms in \(\|K\|^2\) are retained and bounded by
the cross-pullback estimate, giving main term \(3p^3/4\).

The mean bound \(|\ell|\le4\) follows from
\((2p+6)/(p-5)\le4\) for \(p\ge13\). Consequently

\[
 |d-p^2/4|\le3p^{3/2}+(18p+3)/4.
\]

Dividing by \(p^2\), the right side is decreasing for
\(p\ge1024\), and at 1024 is less than \(1/8\). This proves the
stated positive denominator \(d\ge p^2/8\). Also
\(|s|\le6p^2+12p^{3/2}\le(51/8)p^2\), so
\(|3s/d-\ell|\le157<160\), proving
\(|a-1/2|\le20/\sqrt p\). Positivity of \(d\) at this threshold
is not a claim that \(c\) is already positive at every such size.

## 5. Actual leakage and independent finite verification

The row identity gives \(b=\sqrt d/(8\sqrt{mp})>0\), so the sign
in the limiting two-by-two matrix is fixed by the chosen direction
\(v=g/\sqrt d\). Projecting the recurrence onto the complement of
\(1,g\) gives exactly \((I-\Pi)Bv=h/(8\sqrt{pd})\). The two
subtracted projections in \(\|h\|^2\) have order \(p^2\), below
the \(3p^3/4\) main term. This verifies \(c^2=3/64+O(p^{-1/2})\).

Since \(Bu=qu+bv\) and \(Bv=bu+av+w\), with \(w\perp u,v\)
and \(\|w\|=c\), the Gram matrix of the actual images is

\[
 \begin{pmatrix}q^2+b^2&b(q+a)\\b(q+a)&b^2+a^2+c^2\end{pmatrix}
 \longrightarrow\frac1{64}\begin{pmatrix}10&7\\7&20\end{pmatrix}.
\]

Its top eigenvalue is \((15+\sqrt{74})/64\). The compressed
Rayleigh matrix has top limit \((7+\sqrt5)/16\). These are different
quantities; the leakage is correctly present in the first one.
The deleted cross-block has precisely the two possible nonzero
eigenvalues \(\pm c\), establishing its norm without discarding the
orthogonal complement.

For an independent finite check I defined \(k\) by multiplicative
convolution with \(\chi(1-t)\), rather than by the author's full
\(SD_0\) calculation. I then formed the actual compressed matrix
images \(B1\) and \(Bg\) directly, using exact pairs
\(a+b\sqrt p\), and projected those images directly to check both
Gram diagonal entries, the off-diagonal entry, and the leakage.
The calculation passed at \(p=13,37,61,1201\). Useful exact values
are:

| p | d | \(c^2\) |
| ---: | --- | --- |
| 13 | 0 | next direction undefined; not divided by d |
| 37 | 384 | 0 |
| 61 | \(3072/7\) | \(7/244\) |
| 1201 | \(8890368/23\) | \(7713562969/130762559668\) |

These independently checked finite identities do not certify the
imported sheaf theorems or derive the asymptotic limits from data.
The small vanishing values of \(c\) are consistent with eventual
nonvanishing, which is the asymptotic assertion.

## 6. Scope relative to the earlier seeded result

The earlier seeded note already supplied hypergeometric kernel
correlation tools and quantitative bounds on specified spans. It also
showed an actual missing inversion-odd spectral sector, with an
extremal vector outside the old span at \(p=13\). The present work
does not overturn that obstruction or establish a new general sheaf
independence theorem.

The additional evaluated quantities are the next actual recurrence
coefficient and its residual norm, using three particular Möbius
pullbacks. The endpoint/reflection mismatch is what justifies those
new correlations; the earlier disjoint-seed/anchor hypotheses should
not be applied to them without this check. The resulting upper bound
is for the actual two-dimensional Krylov domain and includes its
nonzero outgoing coupling.

All these directions remain in the \(S_3\)-invariant sector. The next
\(h^TS_Ch\) estimate, growing-depth control, and the operator on
other directions remain unproved. Thus the reviewed result is a
valid fixed-subspace advance, with no conclusion about the full
spectral edge or a worst-case Paley exponent.
