# All symmetry sectors and an exact elliptic-operator reduction

**Status: an exact reduction covering the full actual operator at every
depth. No uniform character-sensitive norm saving was obtained.**
The uniform-seed Krylov method is confined to a subspace of asymptotic
relative dimension 1/6, even with unlimited depth. Three explicit
invariant blocks cover the missing sectors. On the nontrivial blocks,
the desired edge becomes an upper bound 2p/3 for an averaged elliptic
kernel; the elementary estimate still gives p. On the trivial block,
the exact nonvanishing rank-two correction is retained.

All statements concern primes p≡1 mod4, p≥13, and the actual set

    C={x:χ(x)=χ(x−1)=1}, m=(p−5)/4,
    S_C=(χ(x−y))_(x,y∈C), B=(I−J/p+S_C/√p)/2.

They transfer to every anchor edge by residue affine dilation. This
is not a Paley proof, a growing-depth moment estimate, or a prize result.

## 1. A complete block decomposition

Let R and I₀ be the permutation operators on C induced by
x↦1−x and x↦1/x. They are orthogonal involutions. Put U=RI₀, so
U³=I and RUR=U^(−1). Multiplicativity proves that S_C commutes with
both R and I₀: simultaneous inversion multiplies χ(x−y) by χ(xy)=1.
The all-ones matrix commutes with every permutation, so B does too.

Define the following rational orthogonal projections:

    P_cyc=(I+U+U²)/3,
    P_+=(I+R)P_cyc/2, P_-=(I−R)P_cyc/2,
    P_st=I−P_cyc,
    P_st,±=(I±R)P_st/2.                            (1)

Their four ranges are pairwise orthogonal and sum to R^C. Every
range is invariant under B and S_C. The operator

    V=(U−U²)/√3

commutes with B, satisfies VᵀV=P_st and VR=−RV, and maps the two
standard ranges isometrically onto one another. Thus the full operator
is the direct sum

    B ≅ B_+ ⊕ B_- ⊕ B_st ⊕ B_st,                 (2)

where B_st is its restriction to the minus standard range. In
particular the full norm, spectrum and traces at **every** power are
recovered from these three blocks. Equation (2) is not a truncation.

Here are their exact dimensions. Let ε₂=1 if χ(2)=1 and zero
otherwise, and let ε₃=1 if p≡1 mod3 and zero otherwise. The exceptional
S₃ orbits in C are

    {−1,2,1/2}, size 3, present iff ε₂=1;
    {x:x²−x+1=0}, size 2, present iff ε₃=1.

All other orbits have size six. For the second assertion the roots
have order six; when p≡1 mod4 and 3 divides p−1, they are squares,
and x−1=x² is a square as well. The exceptional orbits are disjoint
for p>3. Writing g=(m−3ε₂−2ε₃)/6 gives

    d_+=g+ε₂+ε₃, d_-=g+ε₃, d_st=2g+ε₂,
    m=d_++d_-+2d_st.                              (3)

Equivalently, trace R=ε₂ and trace U=trace U²=2ε₃ give the same
dimensions directly from (1). The largest block has dimension
m/3+O(1); the trivial block has dimension m/6+O(1).

Since the uniform vector u is in range P_+, every polynomial in B
applied to u is in that range. Thus the entire uniform-seed Krylov
space has dimension at most d_+, at any depth. Its orthogonal omitted
space has dimension at least 5m/6+O(1). This is more specific than
the earlier [single-kernel inversion-even obstruction](parallel6-seeded-kernels-2026-09-04.md):
it concerns Krylov iteration from u, not every unsymmetrized kernel
in the old family. The [pass24 directions](parallel24-spectral-operator-2026-09-05.md)
remain in this trivial sector.

Both other sector types can contain the maximum in actual fields.
At p=13, C={4,10} and the top B eigenvalue is in the sign block.
At p=17, C={2,9,16}, S_C=I−J, and the top eigenvalue
(1+1/√17)/2 is in the two-dimensional standard sector. These finite
examples reject a universal trivial-sector extremizer assumption;
they make no asymptotic assertion about the location of the extremum.

## 2. A character-sensitive identity on the whole space

Let S be the full prime-field sign matrix, D₀=diag(χ(t)),
D₁=diag(χ(t−1)), and define actual restricted kernels

    K₀=(S D₀ S)_(C,C), K₁=(S D₁ S)_(C,C),
    K₀₁=(S D₀D₁ S)_(C,C).

Writing L(t)=Σ_yχ[y(y−1)(y−t)], homogeneity gives the concrete
elliptic kernel

    (K₀)_(x,y)=L(y/x), x,y∈C.                     (4)

The diagonal is included: L(1)=−1. Reflection and inversion give

    K₁=R K₀ R,
    K₀₁=I₀ R K₀ R I₀−J.                         (5)

For the second identity substitute t=1/z in its defining complete
sum. The missing z=0 term of K₁(1/x,1/y) is one, because x,y∈C.
It must be subtracted, producing −J. Also K₀ commutes with I₀,
either by the same substitution or L(1/t)=χ(t)L(t).

Consequently

    A=[K₀+R K₀ R+I₀ R K₀ R I₀]/3                 (6)

is the full S₃ conjugacy average of K₀. It commutes with all four
projections in (1). The exact common-neighbor mask is

    E_C=(I+D₀+D₁+D₀D₁)/4−(E_{0}+E_{1})/2.

Using S²=pI−J and (S(E_0+E_1)S)_(C,C)=2J, then (5), proves

    S_C²=pI/4+3A/4−3J/2.                          (7)

The coefficient −3/2 retains both the inversion boundary and the
two anchor columns. Its diagonal checks exactly:
m−1=p/4−3/4−3/2=(p−9)/4.

The elementary ambient norm estimate is

    ||K₀||≤||S||²||D₀||≤p, hence ||A||≤p.         (8)

For every nontrivial block J=0. Equation (7) then gives the complete
operator identity

    S_σ²/p=I/4+3A_σ/(4p), σ=−,st.                (9)

Thus A_σ≥−pI/3, but the available upper bound is still A_σ≤pI.
The desired two-anchor edge on either block is exactly the stronger
requirement

    λ_max(A_σ)≤(2/3+o(1))p.                       (10)

Indeed (10) is equivalent to ||S_σ||≤(√3/2+o(1))√p and hence
||B_σ−I/2||≤√3/4+o(1). This identifies a required one-third saving
against (8), on actual elliptic kernels. Equation (10) is not proved.
It is a reduction of the all-vector arithmetic estimate, not an
independently stronger known theorem.

On the sign block the compression of each conjugate in (6) is the
same, so A_-=P_-K₀P_- on its range. The analogous statement holds
on the trivial block. On the standard multiplicity block (6) retains
the average of the two equivalent standard components. Replacing
that average by its largest summand only recovers (8); no uniform
sector saving was obtained by this route.

## 3. Keep the trivial-sector border explicitly

The normalized full-field matrix T=S/√p−J_all/p is an orthogonal
involution. Its C compression is 2B−I. Therefore the actual leakage
Gram matrix is

    E=T_(C^c,C)ᵀT_(C^c,C)=I−(2B−I)²=4B(I−B).

The desired full edge is equivalent to E≥(1/4−o(1))I. Substituting
(7), with products in the C coordinate space, gives

    E=3(I−A/p)/4+(3/(2p)−m/p²)J
                      +(S_CJ+JS_C)/(p√p).        (11)

This includes every vector and every symmetry sector. On the
nontrivial sectors it reduces to (9)–(10). On the trivial sector the
last two terms cannot be silently dropped.

More explicitly, use the actual pass24 u,v,b, put ζ=m/p and
a₀=uᵀS_Cu=(ℓ−6)/4, and let Q=uuᵀ. Since
S_Cu=a₀u+2√p b v, the exact correction in (11) is

    κQ+2ζb(uvᵀ+vuᵀ),
    κ=ζ(3/2−ζ)+2ζa₀/√p.                          (12)

By the proved pass24 estimates, κ→5/16 and 2ζb→1/16.
Thus this is a genuine rank-two correction, with a nonvanishing
off-diagonal coefficient. At small primes where the first residual
vanishes, (11) still applies and (12) is read without the v term.
The decomposition does not make span{u,v} invariant under A; that
would discard the further coupling already computed in pass24.

## 4. Why this has not overcome the depth limitation

For Z=4(B−I/2)/√3, the exact block decomposition gives, for every j,

    tr T_(2j)(Z)=tr T_(2j)(Z_+)+tr T_(2j)(Z_-)
                                      +2tr T_(2j)(Z_st),          (13)

where T_k is the Chebyshev polynomial. If a block of dimension d>0
has tr T_(2j)(Z_σ)≤M with M≥0, then

    ||Z_σ||≤cosh(log(2(M+d))/(2j)).                (14)

To prove this, an eigenvalue outside [−1,1] contributes
cosh(2j arcosh|λ|), while every other even-Chebyshev term is at
least −1. If none is outside, (14) is automatic. Thus logarithmic
size of the trace allowance still requires j/log p→∞ for an edge
1+o(1), since d_+,d_-,d_st are all Θ(p).

The recorded [pass13 principal-rank allowance](parallel13-principal-budget-2026-09-05.md)
has, for two anchors,

    Ψ₂=(19+5√13)/6>1, Θ₂=25/3,
    j≤min{(1/2−ε)log p/log Ψ₂,
           (1−ε)log p/log Θ₂}.

This is a comparison with that existing source-dependent allowance,
not a new proof of its convolution machinery. It remains O(log p).
Its principal correlation error budget grows like Ψ₂^j/√p.
The present block decomposition changes dimensions by fixed factors;
it supplies no cancellation against that exponential rank allowance.

One can also transfer a whole-space Chebyshev bound to a sector:
every other sector contributes at least minus its dimension, so
isolating one term of (13) costs only O(m). Thus a full O(p) bound
stays O(p) per sector, not o(1) or an exponentially smaller quantity.
Conversely, proving all three sector bounds at a sufficiently growing
depth would control the entire operator and would avoid the missing
sectors of uniform-seed iteration. That analytic step is still open.

Powers of (6) expand into three-conjugate elliptic words. No improved
rank/conductor or signed aggregate bound for those words is established
here. Applying (8) at every step gives no saving. One more fixed
Lanczos coefficient would not resolve this precise remaining issue.

## 5. Verification and source scope

The new [integer verifier](../experiments/parallel25_spectral_full_operator_2026_09_05.py)
and [results](../results/parallel25_spectral_full_operator_2026_09_05.json)
check the actual orbit types, all rational projectors after clearing
denominators, equivalence of the standard blocks, elliptic kernels,
the boundary coefficient in (7), the full leakage identity, and
sector traces and uniform-seed containment in selected fields.
No floating-point spectral fit supplies a theorem.

The new identities and decomposition are elementary and proved above;
no new cohomological or Hecke theorem is imported. The limiting border
coefficients use the separately reviewed pass24 theorem. The pass6
and pass13 records are used to locate the earlier missing sectors and
depth limitation, respectively. No new determinant certificate is
claimed as a constant edge, and no prime-field counterexample to the
desired edge is constructed. Human and formal review remain open.
