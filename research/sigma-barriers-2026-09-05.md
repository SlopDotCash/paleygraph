# Barrier map for the two-set Paley conjecture

**Status: a synthesis of PROVED elementary facts and CITED published
obstructions. Nothing here proves or refutes the conjecture. Its purpose is
to state, with proofs where elementary, exactly which properties any proof
of the conjecture must use and which proof shapes are already excluded.**

Throughout, `p` is an odd prime, `χ` the Legendre symbol with `χ(0)=0`,
and for `A,B ⊆ F_p`

    S(A,B) = Σ_{a∈A} Σ_{b∈B} χ(a+b).

The conjecture (Satake 2020, Conjecture 7) asserts: for every `ε∈(0,1)` there
are `δ>0`, `p₀` with `|S(A,B)| ≤ p^{-δ}|A||B|` whenever `p>p₀` and
`|A|,|B|>p^ε`. The `a−b` form is the same statement with `B` replaced by `−B`.

## 1. Four analytic barriers, each pinned to a single lossy step

### B1. The Fourier / Cauchy–Schwarz barrier (PROVED)

For `c≠0`, `Σ_x χ(x)χ(x+c) = −1` (substitute `x ↦ x(x+c)`; the map
`x ↦ 1+c/x` is a bijection from `F_p^*` onto `F_p∖{1}`). Hence for any `B`,

    Σ_x ( Σ_{b∈B} χ(x+b) )² = |B|(p−1) − |B|(|B|−1) = |B|(p−|B|),

and by Cauchy–Schwarz over `x∈A`,

    |S(A,B)| ≤ √( |A| |B| (p−|B|) ).

This is below the trivial bound `|A||B|` only when `|A||B| > p−|B|`. The lossy
step is the extension of the sum over `x∈A` to all of `F_p`: it discards the
information that `A` is small. No argument that only uses the second moment of
`F_B(x)=Σ_b χ(x+b)` over all `x` can pass `|A||B| ≈ p`.

### B2. The Hölder–Weil barrier (PROVED derivation, standard)

Write `S(A,B) = Σ_{a∈A} F_B(a)` and apply Hölder with exponent `2k` and then
extend to all `x`:

    |S(A,B)|^{2k} ≤ |A|^{2k−1} Σ_{x∈F_p} F_B(x)^{2k}
                 = |A|^{2k−1} Σ_{b_1,…,b_{2k}∈B} Σ_x χ(Π_i (x+b_i)).

By the Weil bound the inner sum is at most `(2k−1)√p` unless the multiset
`{b_i}` has every element of even multiplicity, in which case it is at most
`p`. The number of such ordered tuples is at most `(2k−1)!!·|B|^k` (choose a
perfect matching of the `2k` positions and a value for each pair; for `k=2`
the exact count is `3|B|²−2|B|`, so the constant `k!` sometimes quoted is
wrong, as the verifier shows). Thus

    |S(A,B)| ≤ |A|^{1−1/(2k)} ( p (2k−1)!! |B|^k + (2k−1)√p |B|^{2k} )^{1/(2k)}.

The second term is `≤ p^{−δ}|A||B|` only if `|A| > (2k−1)^{…} p^{1/2+δ'}`.
This is Karatsuba's regime `|A| > p^{1/2+ε}`, `|B| > p^{ε}`. The lossy step
is again the extension of `Σ_{x∈A}` to `Σ_{x∈F_p}`; the Weil input is sharp
and cannot be improved for a single polynomial. Any proof for
`|A| ≤ p^{1/2}` must retain the restriction to `A` inside the amplification,
which is exactly the "local-to-average transfer" problem examined in
[the crux note](sigma-crux-2026-09-05.md).

### B3. The polynomial-method barrier (CITED)

Hanson–Petridis (Refined estimates concerning sumsets contained in the roots
of unity, Proc. LMS 2021) prove via Stepanov's method that a clique in the
Paley graph on `p` vertices has size at most `√(p/2)+1`; the Randomstrasse101
problem list (Problem 25) records `(1+o(1))√(p/2)` as the best known bound.
The method constructs an auxiliary polynomial of degree about `|A|²` that
vanishes to high order on the clique; the degree must stay below `p` for the
derivative argument, which pins the exponent at `1/2`. No variant of the
method that uses only vanishing multiplicities has gone below this.

### B4. The subfield obstruction (PROVED, elementary)

Let `q=p²` and let `χ_q` be the quadratic character of `F_q`. The subgroup
`F_p^*` has index `p+1` in the cyclic group `F_q^*`, and `p+1` is even, so
`F_p^* ⊆ (F_q^*)²`, i.e. `χ_q(u)=1` for every `u∈F_p^*`. Therefore with
`A=B=F_p ⊂ F_q` (so `|A|=|B|=q^{1/2}`),

    Σ_{a,b∈F_p} χ_q(a−b) = p(p−1) = |A||B|(1−q^{−1/2}).

The verbatim two-set statement over `F_q` is false for every `ε<1/2`.
Hence a proof over prime fields must use a property that distinguishes `F_p`
from `F_{p²}`: the absence of intermediate subfields, which enters analytic
arguments only through sum-product / incidence theorems (Bourgain–Katz–Tao and
their descendants). Every step of a would-be proof that is "field-agnostic"
(Weil bounds, orthogonality, Fourier analysis, moment identities, spectral
identities of the Paley matrix, the polynomial method as used in B3) cannot
by itself pass `ε=1/2`, because all of those steps hold verbatim over `F_{p²}`.
This is the single most useful diagnostic for a proposed proof: locate the
step that fails over `F_{p²}`; if there is none, the proof is wrong.

The verifier constructs `F_{p²}` explicitly and confirms `χ_q ≡ 1` on `F_p^*`
for every odd prime `p ≤ 60`.

### B5. The short-interval (Burgess) special case (PROVED reduction, see crux note)

Take `A=B={1,…,N}` with `N=⌊p^ε⌋`. If every integer in `[2,2N]` were a
quadratic residue mod `p` then `S(A,B)=N²`, so the conjecture forces the
least quadratic non-residue `n_p ≤ 2N ≤ 2p^ε` for large `p`: Vinogradov's
conjecture, open since 1919 (Burgess 1957 gives `n_p ≪ p^{1/(4√e)+o(1)}`).
More generally `S(A,B) = Σ_n w(n)χ(n)` with the triangle weight
`w(n)=#{(a,b)∈[1,N]²: a+b=n}`, so the conjecture implies cancellation in
smoothly weighted character sums over intervals of length `p^ε`, beyond the
Burgess range `N > p^{1/4+ε}`. Any general proof must therefore contain a proof
of a beyond-Burgess estimate; conversely, no method that is known to be
stuck at the Burgess exponent for intervals can prove the conjecture.

## 2. Two excluded proof shapes (CITED)

### B6. Low-degree sum-of-squares certificates

Kunisky–Yu (A degree 4 sum-of-squares lower bound for the clique number of
the Paley graph, arXiv:2211.02713, CCC 2023) prove that the degree-4 SOS
relaxation of the clique number of the Paley graph on `p` vertices has value
`Ω(p^{1/3})`. Since the two-set conjecture implies a clique bound
`p^{o(1)}`, no proof of the conjecture can be encoded as a degree-4 SOS
certificate for the clique number; Randomstrasse101 Problem 28 records that
degree 4 can at best improve the exponent from `1/2` to `1/3`. This does not
exclude certificates of degree growing with `log p`, which is the regime the
workspace's spectral line (traces of depth `j` with `j/log p → ∞`) targets.

### B7. The restricted-isometry "square-root bottleneck"

Bandeira–Mixon–Moreira (A conditional construction of restricted isometries,
arXiv:1410.6457) show that the Paley equiangular tight frame satisfies the
restricted isometry property with sparsity `Ω(M^{1/2+ε})` conditionally on a
number-theoretic conjecture about quadratic residues (their exact
formulation is not quoted here; UNVERIFIED CITATION for its form).
Satake (On the restricted isometry property of the Paley matrix,
arXiv:2011.02907, local copy `sources/satake-2011.02907.txt`, Theorem 10)
proves the sharper statement: if Conjecture 7 holds then for every
`0<α<1/2` there is `β₀=β₀(α)>0` such that the Paley matrix `Φ_p` has the
`(p^{τ+β₀}, p^{τ−1/2+o(1)})`-RIP for `max(α+β₀, 1/2−β₀) < τ < 1/2`, so
`K=Ω(M^γ)` with some `γ>1/2`. Randomstrasse101 Problem 29 states the
unconditional version as open. The two-set conjecture therefore sits
upstream of the RIP problem; all deterministic RIP constructions being stuck
at the square-root bottleneck is the same phenomenon as B1–B3 in
signal-processing language, not independent evidence about the conjecture.

## 3. What follows for the programme in this workspace

1. The spectral/necklace line (passes 2–11) is field-agnostic in every
   step that has been written down: the identities `S²=pI−J`, the projection
   `P`, the necklace traces and the Katz-type bounds all hold over `F_{p²}`.
   By B4 the line cannot reach `ε<1/2` unless a step is added that fails over
   `F_{p²}`. The growing-depth criterion `j/log p → ∞` is exactly where such
   a step would have to enter; fixed-depth trace bounds are field-agnostic and
   therefore cannot suffice on their own.
2. The subgroup line (Gauss sums over `μ_n`, `p≈n⁴`) is a genuinely
   prime-field statement (Bourgain–Glibichuk–Konyagin use sum-product), which
   is why it admits `p^{−δ}` savings while the arbitrary-set problem does not.
3. The only known mechanism that passes B1–B4 for arbitrary sets is the
   sum-product route of Chang (Duke 2008) with exponent `4/9`; see
   [the sum-product note](sigma-sumproduct-2026-09-05.md) for the current
   bottleneck.
4. A candidate proof should be tested first against the interval case (B5)
   and against `F_{p²}` (B4). Both tests are cheap and decisive.

## Verification

`experiments/sigma_barriers_2026_09_05.py` checks, with exact integer
arithmetic: the shift orthogonality and second-moment identity of B1 for all
`B` over `p ≤ 13` and random `B` over `p ≤ 199`; the Chung bound; the Weil-side
tuple count in B2 for `k=2,3` at small `p` (the exact number of even-multiplicity
tuples and the exact inner sums); the subfield identity `χ_q ≡ 1` on `F_p^*`
for all odd `p ≤ 60`; and the interval implication in B5 by exhibiting, for
each `p ≤ 500`, the equality `S([1,N],[1,N]) = N²` for `N = ⌊(n_p−1)/2⌋`.
Results go to `results/sigma_barriers_2026_09_05.json`.
