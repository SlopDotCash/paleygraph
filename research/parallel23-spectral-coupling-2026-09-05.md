# The uniform direction has nonvanishing spectral coupling

For actual prime-field two-anchor neighborhoods, the direction controlled
in pass 22 is **not** an approximate eigenvector. Its coupling to its
orthogonal complement tends to 1/8 in the compressed Fourier projection.
This evaluates a previously unestimated quantity; it supplies no upper
bound on the full spectrum and no Paley or prize proof. The asymptotic
input is an existing elliptic-family moment theorem, not a new theorem
proved here about that family.

## Statement

Let p≡1 mod4 be prime, p≥13, let χ be the quadratic character with
χ(0)=0, and put

    C={t:χ(t)=χ(t−1)=1},  m=(p−5)/4,  u=1_C/√m,
    L(t)=Σ_y χ[y(y−1)(y−t)],  S_C=(χ(x−y))_(x,y∈C).

Use the already normalized projection from the
[previous spectral note](parallel21-spectral-next-input-2026-09-05.md),
compressed to C:

    B=(I_C−J_C/p+S_C/√p)/2,  q=u^TBu,  Q=uu^T,
    b²=||(I_C−Q)Bu||².

Here J_C is the all-ones matrix. Along these primes,

    b²→1/64,               b→1/8,                 (1)
    u^T B²u→5/32.                                  (2)

Consequently deleting just the cross terms between the uniform vector
and its orthogonal complement has a nonvanishing operator error:

    ||B−[QBQ+(I_C−Q)B(I_C−Q)]|| = b → 1/8.         (3)

The statement is uniform over every actual anchor edge, because residue
affine dilation carries its neighborhood to C and preserves S.

## Source input and exact restricted second moment

Write D=F_p\{0,1}, and define

    T=Σ_(t∈D)L(t)²,
    W=Σ_(s∈F_p\{0,1,−1})L(s²)²,
    Z=Σ_(t∈C)L(t)².

The needed published input is W/p²→1. Brian Grove's
[Hypergeometric moments and Hecke trace formulas](https://link.springer.com/article/10.1007/s11139-026-01372-y)
(published 6 April 2026), Theorem 1.3 with moment 2, supplies it:
the introductory point-count formula identifies H_p(t)=−L(t) on D
when p≡1 mod4. The omitted s=0,±1 each have squared value 1.
Theorem 2.3(2), with k=2 and G_4(z,p)=z²−p, also records the
Ahlgren–Ono exact trace formula

    W=p(p−3)−4−a_8(p),
    a_8(p)=Tr(T_p | S_4(Γ_0(8))).                  (4)

Only this fixed second moment and optional exact identity are imported.
No claim of a self-contained modular-form proof or a uniform growing
moment order is made. The [publisher HTML](../sources/grove-hypergeometric-moments-2026.html)
was retrieved and its theorem statements and normalization reviewed;
no PDF was used. Archive metadata is in the
[source record](../results/parallel23_spectral_source_2026_09_05.json).

Here is an elementary transfer from that source to C. The changes of
variable y↦1−y and y↦y/t prove, for t∈D,

    L(1−t)=L(t),      L(1/t)=χ(t)L(t).             (5)

Thus L(t)² is invariant under t↦1−t and t↦1/t. The three twisted
sums with weights χ(t), χ(t−1), and χ(t(t−1)) are equal. Call their
common value A. The first change interchanges the first two weights;
the inverse change sends the second to the third. Both changes are
bijections of D, including the exceptional short permutation orbits.

The exact indicator on D is

    1_C(t)=[1+χ(t)+χ(t−1)+χ(t(t−1))]/4.

Each t∈D has 1+χ(t) square roots in F_p\{0,1,−1}, so

    4Z=T+3A,   W=T+A,   4Z=3W−2T.                (6)

To evaluate T, let f(y)=χ(y(y−1)). Then Σf=−1 and Σf²=p−2.
The exact translate identity Σ_t χ(y−t)χ(z−t)=pδ_(y,z)−1 gives

    Σ_(t∈F_p)L(t)²=p(p−2)−1.

Directly L(0)=L(1)=−1, hence

    T=p²−2p−3.                                    (7)

Equations (6), (7) and W/p²→1 yield Z/p²→1/4. In particular,

    Z/(mp)→1.                                     (8)

Keeping (4) gives the exact refinement

    A=−a_8(p)−p−1,
    Z=[p²−5p−6−3a_8(p)]/4.                       (9)

The source theorem implies a_8(p)=o(p²). No explicit numerical error
constant or prime threshold is extracted here. In the finite verifier,
a_8(p) is defined from W via (4); this does **not** independently
compute a Hecke trace or verify the imported trace formula.

## Variance, compression, and the block error

Set ℓ=(1/m)Σ_(t∈C)L(t). The earlier exact anchor/Jacobi calculation
gives, with U=J(η,χ)²+conj(J(η,χ))²,

    Σ_(t∈C)L(t)=(U+6)/4,     |U|≤2p.

Therefore ℓ=O(1). Also (S_C1_C)(t)=(L(t)−6)/4. Subtracting the
mean removes the −6 exactly, and gives

    ||(I_C−Q)S_Cu||² = [Z/m−ℓ²]/16.              (10)

Both I_Cu and J_Cu are parallel to u. Consequently

    b² = [Z/m−ℓ²]/(64p) → 1/64,                  (11)

by (8). An optional exact expression from (9) is

    b²=1/64−(6+3a_8(p))/(64p(p−5))
             −(U+6)²/(64p(p−5)²).                (12)

The previous result gives q→3/8. Since Bu=qu+(I_C−Q)Bu is an
orthogonal decomposition, u^TB²u=q²+b²→9/64+1/64=5/32.

For (3), write w=(I_C−Q)Bu. Symmetry gives the deleted matrix
uw^T+wu^T. If w≠0, its restriction to span(u,w) is the matrix
with zero diagonal and both off-diagonal entries ||w||; it vanishes
on the remaining orthogonal space. Its norm is therefore ||w||=b.
The formula also holds if w=0, which can occur at small primes.

## Interpretation and verification limits

The first Rayleigh quotient remains controlled. This calculation shows
that replacing its direction by an isolated scalar block cannot have
o(1) error in operator norm. It does not exclude a method retaining
the coupling, controlling a larger subspace, or bounding the remaining
operator directly. It also gives no adversarial-vector estimate and
no control of the spectrum from finitely many moments.

The [exact finite verifier](../experiments/parallel23_spectral_coupling_2026_09_05.py)
checks the full-field, square-map, symmetry, restricted-moment, row-sum
and variance identities on actual prime fields. It independently
forms neighborhood row sums and exact rational projection entries.
Its [results](../results/parallel23_spectral_coupling_2026_09_05.json)
distinguish these identities from the imported asymptotic theorem.
No floating-point fit proves (1), and no modular coefficient is
computed independently. Separate-agent review is recorded separately;
human refereeing and formal verification are not claimed.
