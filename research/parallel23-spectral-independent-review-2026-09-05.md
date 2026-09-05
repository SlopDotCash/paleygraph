# Separate-agent review of the two-anchor spectral coupling

Verdict: **no mathematical correction is required** in the pinned coupling
note. Its stated limits hold for the displayed compression on actual
two-anchor neighborhoods as primes p≡1 mod4 tend to infinity:

    b²→1/64,  b→1/8,  uᵀB²u→5/32,
    ||B−[QBQ+(I−Q)B(I−Q)]||=b.

The proof uses an existing published second-moment theorem. It neither
proves that theorem anew nor obtains an upper bound on all vectors, a
spectral-edge estimate, a growing-anchor result, or the Paley or prize
conjecture. This is a proof review by a separate agent, not independent
human refereeing or formal verification.

## 1. Inputs and source scope

The proof, verifier, recorded results, and prior Jacobi/anchor derivation
were read. The exact archived primary publisher HTML of Brian Grove,
[*Hypergeometric moments and Hecke trace formulas*](https://link.springer.com/article/10.1007/s11139-026-01372-y),
was inspected at the introductory point-count normalization, Theorem 1.3,
the definition of G_k, Theorem 2.3(2), and Lemma 3.4. The archive and
source record agree on the HTML hash. No PDF was used. This review does
not claim a fresh publisher retrieval, current prize-status verification,
or a review of every theorem or reference in Grove's paper.

SHA-256 values checked against the actual bytes during this review:

| Input | SHA-256 |
| --- | --- |
| `research/parallel23-spectral-coupling-2026-09-05.md` | `6d8cd29185b0c86b33961a340d11f5b2313697f2950172807678fbf69c568ed2` |
| `experiments/parallel23_spectral_coupling_2026_09_05.py` | `dc2bf169ce8ae0b9669ac9858375e66a0c958876b6e27a0bec2103ff67f68de3` |
| `results/parallel23_spectral_coupling_2026_09_05.json` | `cf0d31dd2d2620c17cf2c36092b8b0c62997c7385300e4710cb0cc755b65c55f` |
| `sources/grove-hypergeometric-moments-2026.html` | `1174074d2d3104414f03446b060a245d935c0fe2ca80beb2c14a01bbe10ec440` |
| `results/parallel23_spectral_source_2026_09_05.json` | `c83900cd4933c33257f92561f3d36348266ef4f4b120e0f35a42dad104000099` |
| `research/parallel21-spectral-next-input-2026-09-05.md` | `2bc5a06bebfd1cdb2cd7f0d61d9f95c8448fe7bd64ace017f37edbaa7093df98` |
| `results/parallel21_spectral_flat_direction_2026_09_05.json` | `c88cac808f6352902b61b943298ecf2ade43858604439b29066f98607fd2f29a` |

All six input hashes recorded by the coupling verifier match these files.
No pre-existing artifact was edited during this review.

## 2. Published normalization and the optional exact trace formula

For nonsingular Legendre parameters t∉{0,1}, the affine point count and
the single point at infinity give

    #E_t(F_p)=p+1+L(t),  a_t^Leg(p)=−L(t).

Grove's introduction states H_p(t)=χ(−1)a_t^Leg(p). Thus H_p(t)=−L(t)
on this domain when p≡1 mod4. This is the unscaled trace normalization;
there is no extra factor p or √p.

Theorem 1.3 at its fixed moment m=2 says

    p^(−2) Σ_(s∈F_p) H_p(s²)² → C(1)=1.

The introduction defines H_p(0)=1. Lemma 3.4, equation (3.4), with d=2
gives H_p(1)=χ(−1)=1 for these primes. The excluded parameters s=0,1,−1
therefore contribute exactly 3, and W=Σ_s H_p(s²)²−3. Restricting the
published limit over primes to p≡1 mod4 is valid. Consequently W/p²→1.

For the optional exact refinement, the definition of G gives
G_4(z,p)=z²−p. Theorem 2.3(2) at k=2 reads

    −Tr(T_p | S_4(Γ_0(8)))
      =4+Σ_(s=2)^(p−2) [a_(s²)^Leg(p)²−p].

There are p−3 summands, not p−2. The sign and endpoints yield precisely
W=p(p−3)−4−a_8(p). This formula is not required for the limit once
Theorem 1.3 is imported. The deduction a_8(p)=o(p²) follows from this
formula and that limit; the note gives no numerical rate or threshold.

## 3. The transfer to the actual neighborhood

Let D=F_p\{0,1}. Substituting y=1−z in L(1−t) gives χ(−1)L(t)=L(t).
Substituting y=z/t in L(1/t) gives χ(t^(−3))L(t)=χ(t)L(t).
Thus L² is invariant under both bijections t↦1−t and t↦1/t on D.

For weights w_0(t)=χ(t), w_1(t)=χ(t−1), w_2(t)=χ(t(t−1)),
the first bijection interchanges w_0 and w_1, while the second sends
w_1 to w_2. Therefore all three weighted sums of L² equal A.
This uses bijections rather than division by an orbit size, so the
short S_3 orbits introduce no exceptions.

For each t∈D, z²=t has exactly 1+χ(t) solutions, including the
nonsquare case of zero solutions. None can be 0 or ±1. The mask on D
and these root counts give

    4Z=T+3A,  W=T+A,  hence 4Z=3W−2T.

The full-field correlation identity has value p−1 on equal arguments
and −1 otherwise. Applied to f(y)=χ(y(y−1)), it gives

    Σ_t L(t)²=p Σ_y f(y)²−(Σ_y f(y))²
             =p(p−2)−1.

Here Σf=−1 and Σf²=p−2. Directly L(0)=L(1)=−1, so
T=p²−2p−3. Combining these exact identities with W/p²→1 gives
Z/p²→1/4. Since m=(p−5)/4, this implies Z/(mp)→1.

The optional refinement also has the correct constants:

    A=−a_8−p−1,  4Z=p²−5p−6−3a_8.

## 4. Anchors, mean, variance, and the operator norm

The full-field mask is

    1_C=(1+χ(t)+χ(t−1)+f(t))/4−(δ_0+δ_1)/2.

The two anchor corrections imply m=(p−5)/4. They can also be used to
check the mean independently of the prior formula for R_C. Quadratic
correlation gives

    Σ_(t∈F_p) L(t)=0,
    Σ_t χ(t)L(t)=Σ_t χ(t−1)L(t)=1,
    Σ_t f(t)L(t)=U.

Since L(0)+L(1)=−2, the mask therefore gives
Σ_C L=(U+2)/4+1=(U+6)/4.
The earlier elementary Jacobi evaluation U=J²+conj(J)², with |J|²=p,
supplies |U|≤2p. Thus ℓ=(U+6)/(p−5)=O(1), uniformly over these primes.

For t∈C, convolution of this mask with χ(t−y) gives a zero constant
term, two correlation terms −1, a term L(t), and anchor subtraction
−1. Hence the row sum is (L(t)−2)/4−1=(L(t)−6)/4.
The factor and the constant in the note are correct. Centering the rows
removes −6 and yields

    ||(I−Q)S_Cu||²=(Z/m−ℓ²)/16.

Since the I and J parts of Bu are parallel to u, its perpendicular
part is (I−Q)S_Cu/(2√p). Thus

    b²=(Z/m−ℓ²)/(64p)→1/64.

The nonnegative square root gives b→1/8. Substituting the optional
expression for Z/m gives exactly equation (12), including both
denominators and the squared term (U+6)². Also the average row sum is
(ℓ−6)/4=O(1), so directly

    q=(1−m/p)/2+(ℓ−6)/(8√p)→3/8.

Symmetry of B and orthogonality then give
uᵀB²u=||Bu||²=q²+b²→5/32.

Finally w=(I−Q)Bu is perpendicular to the unit vector u. The deleted
matrix is uwᵀ+wuᵀ. When w≠0 its matrix in the orthonormal basis
(u,w/||w||) is [[0,b],[b,0]], with eigenvalues ±b, and it is zero on
the orthogonal complement. Its operator norm is exactly b. The zero
case is included and is allowed at small primes. This is an exact
operator-norm proof, not an inference from the Frobenius norm alone.

For every anchor edge {a,b}, χ(b−a)=1 and the affine bijection
t↦a+(b−a)t identifies its common-neighbor sign matrix with S_C.
It also preserves I, the all-ones matrix, and the uniform vector.
The assertions are therefore identical for every edge at a fixed p;
no averaging over anchors is present. All statements require p≡1 mod4
and p≥13, so m>0. No claim is made for all prime powers or more anchors.

## 5. Verifier and coverage limitations

The verifier was read completely and parsed successfully. Its direct
polynomial character sums and direct neighborhood rows are distinct
computations. The ranges for D, the square preimages, and C are correct.
The variance uses exact fractions. The representation of Bu and q as
rational-plus-√p pairs has the correct normalization and no floating
point. Prior Gaussian-integer Jacobi values are used for finite
cross-comparisons; they are archived inputs, not freshly recomputed here.

The recorded run covers 79 primes p≡1 mod4 from 13 through 997,
36,465 parameter records checking the two symmetries and square fibres,
9,057 actual neighborhood rows, and 79 exact projection variances.
Its 2,369 deleted-block entry checks are restricted to p≤101; they
check the Frobenius identity and action on the constant vector. The
general operator norm follows from the rank-two proof above. The
recorded 79 prior Jacobi/variance comparisons are correctly distinguished
from the direct row computation.

This review checked the recorded input hashes, source statements,
algebra, and implementation; it did not repeat the full finite run.
The finite verifier does not prove the published asymptotic, compute
Hecke traces, or independently verify Theorem 2.3. In particular,
its variable a8 is defined by rearranging the imported formula, and
later checks involving that variable verify algebraic consistency only.
The code, result labels, and note state this limitation explicitly.

There is no hidden inference from a finite numerical trend to a limit.
The established conclusion excludes an o(1)-operator-error replacement
by the two diagonal blocks for this particular uniform/complement
decomposition. It leaves methods retaining the coupling, larger
subspaces, and the remaining spectral upper estimate open.
