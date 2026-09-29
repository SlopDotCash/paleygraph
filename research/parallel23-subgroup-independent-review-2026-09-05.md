# Separate-agent review of the subgroup partial upper bound

Verdict: **no mathematical correction is required** in the pinned note.
The product-ratio bijection, zero-product ledger, balanced-subset bound,
coset convolution, and explicit circular-subgroup bounds are correct.
The imported shifted-energy theorem has the required hypotheses and
energy convention. Focused independent computations also confirm the
actual circular witness and the strict fixed-partition inequality.

The proved third-power bound applies only to the product-balanced part
of the opposite-free sixth energy. The unbalanced remainder remains
uncontrolled at that scale. No total-energy exponent, cancellation
exponent, Paley conjecture, or official-prize result is established.
This is a review by a separate agent, not independent human refereeing
or formal verification.

## 1. Inputs and source scope

Read the complete proof, verifier, result records, and lane audit, as well
as the relevant collision bijection and earlier sixth-energy decomposition.
The archived primary [Shkredov HTML](https://arxiv.org/html/1504.04522v1)
was inspected directly at §2's multiplicative-energy definition and §4,
Theorem 6. No PDF was used, and no fresh network retrieval or review of
the entire published proof is claimed here.

SHA-256 values checked against the actual files:

| Input | SHA-256 |
| --- | --- |
| `research/parallel23-subgroup-upper-2026-09-05.md` | `1b97e491bc1c0b2cfe391388b0340d89d1f2a7168822db116513b63f4148bd41` |
| `experiments/parallel23_subgroup_upper_2026_09_05.py` | `ea2023a2cf0b78434efebb9913249c6dc74aec9247da2ad5b8d8a2000251414f` |
| `results/parallel23_subgroup_upper_2026_09_05.json` | `d675ecc47be55cde21ba176fbd5d3bfdf9c0da35af79c093079d31e47cb82770` |
| `results/parallel23_subgroup_upper_audit_2026_09_05.json` | `87df5f5ec779256a1fde49faf6ce4c740f08ed6d6b9cc8fcd15d75675bcfe3c4` |
| `research/parallel21-subgroup-next-input-2026-09-05.md` | `84005053f3954ee2d36eccf85b2786f1fe00d6bd956c76b23b1b40e82b991c22` |
| `research/parallel22-subgroup-independent-review-2026-09-05.md` | `3915e01c73041265eac37be84b9f9756ffc57fa72ba257bf546ab2e6b23fe693` |
| `research/mixed-periods-and-shifted-energy.md` | `de7d7202e55253955b62f73c973d2bdee42b48123eb344c875d5f0527a21f4aa` |
| `sources/mixed-periods-2026-09-04/shkredov-1504.04522.html` | `be13a230cc3cfdb780fd8dca0606d35a6f561221cb3ba67a860be6f759a1f459` |

The recorded input hashes and the lane audit's three artifact hashes
match the current bytes. The lane audit is accurately described as
syntax/hash/link checking, not another mathematical verification run.

## 2. Product-ratio fibers and the exact ledger

For an ordered pair of triples counted by B_ρ, all denominators in
λ=x₁/(y₂y₃), a=λy₂, d=λy₃, b=λx₂/ρ, c=λx₃/ρ are nonzero and belong
to H. Direct substitution gives λx₁=ad, and the product equation gives
λy₁=ρbc. Thus the scaled triples are exactly

    λx=(ad,ρb,ρc),  λy=(ρbc,a,d).

Their equal sums are equivalent to ad−a−d=ρ(bc−b−c). Conversely these
displayed triples, divided by any λ∈H, have product ratio ρ and recover
x₁/(y₂y₃)=λ. The inverse is unique, including repeated coordinates.
The free scaling coordinate contributes exactly n, not n².

Consequently B_ρ=nΣ_z w(z)w(z/ρ). This includes z=0; division by ρ
is valid even at that summation index. Every equal-sum triple pair has
a unique product ratio in H, so E₃=nΣ_(ρ∈H)N_ρ. Dilation preserves the
counting measure on all F_p, which proves N_ρ≤N₁ by Cauchy–Schwarz.
At ρ=1, translation z↦z+1 gives N₁=E×(H−1), with zero included.

The zero-product bucket for H−1 has size 2n−1. Its squared contribution
is (2n−1)². Writing r=n−1, the two trivial product matchings in the
nonzero shifted set have union size 2r²−r. Hence

    E×(H−1)=(2n−1)²+2(n−1)²−(n−1)+𝒳
            =6n²−9n+4+𝒳.

The ordered permutation-pair count is

    D₃=36 binom(n,3)+9n(n−1)+n=6n³−9n²+4n.

The middle term chooses the repeated value, the different value, and
three placements in each triple. Thus B₁−D₃=n𝒳 is exact. The previous
collision bijection identifies 𝒳 with Σ_(i,j)(C_ij)₃, so it is
nonnegative, and 𝒳=0 is equivalent to all actual cells being at most two.

## 3. Opposite-free restriction and the analytic input

Two equal-sum, equal-product triples sharing one value become two pairs
with the same sum and product after cancellation. Their monic quadratic
polynomials agree, so their unordered pairs agree. This proves that a
nonpermutation triple pair has disjoint cross-support. It does not
exclude opposite entries within an individual triple.

Fix an unordered 3+3 partition of the six positions and orient it, for
example by requiring the first half to contain position zero. Keeping
that half and negating the other converts a balanced zero-sum word to
an equal-sum, equal-product pair. An opposite-free word cannot produce
a shared entry across this pair, and hence cannot produce a permutation
pair. The conversion is injective for the fixed partition. Therefore

    R_(6,bal,I)≤B₁−D₃=n𝒳.

There are exactly binom(6,3)/2=10 unordered partitions. The verifier's
choice of two other positions to join position zero implements these
ten without duplication. The union bound gives R_(6,bal)≤10n𝒳;
it is not an equality when a word is balanced at several partitions.

Theorem 6 applies to prime fields, multiplicative subgroups Γ,Π with
|Γ||Π|<p, and arbitrary nonzero shifts. Its definition in §2 counts
all product equalities in subsets of F_p, including zero. Setting
Γ=Π=H and both shifts to −1 therefore gives

    E×(H−1) ≪ n² log n+2n² ≪ n²(1+log n).

The requested quartic window with n≥4 implies n²<p. No restriction on
the shift lying outside −H is required. Thus the displayed balanced
bound follows with an absolute implied constant for this family. The
two finite out-of-window cases are not used to prove that analytic claim.

Summing the bound N_ρ≤N₁ over n ratios yields only
E₃≪n⁴(1+log n). Obtaining third-power total energy through this identity
requires an additional bound on the sum of the correlations. Neither
the individual shifted-energy estimate nor the balanced bound supplies
that missing aggregate estimate.

## 4. Coset convolution and the circular bounds

Use cyclic coset indices and −1∈H. The substitutions z↦−1−z and
z↦1/z give C_ij=C_ji=C_(−i,j−i), with row sum n−1_(i=0).
For x∈H_j, substituting z=−a/x in a+b=x gives

    r₂(x)=C_(−j,−j)=C_(j,0)=κ_j,  r₂(0)=n.

For convolution by one more H element, x−h=0 contributes n when j=0.
For x−h∈H_s, the same substitution counts
C_(−j,s−j)=C_(j,s). Thus

    r₃(x)=n1_(j=0)+(Cκ)_j=n1_(j=0)+(C²)_(j,0),
    r₃(0)=nκ₀.

Squaring the zero term and summing n copies of each nonzero-coset value
gives

    E₃=n³+2n²||κ||²+n²κ₀²+n(C⁴)_(0,0).

Here symmetry is used twice: (C²)₀₀=||κ||² and
||(C²)e₀||²=(C⁴)₀₀. These are exact actual-coset identities; they do
not require an abstract matrix model or a spectral approximation.

Under circularity, every κ_j≤2. Swapping the two summands in r₂(x)
leaves one fixed ordered pair exactly when x∈2H. Therefore the entry
at that coset is one and the remaining entries are zero or two. Since
Σκ=n−1, one obtains Σκ²=2n−3.

The dyadic proof of κ₀=0 is valid: S={x∈H:1−x∈H} has at most two
elements and is invariant under x↦1−x and x↦1/x. Their potential
fixed points force the three distinct elements −1,2,1/2 unless S is
empty. Otherwise both involutions on a nonempty S must swap the same
two elements, forcing x²−x+1=0 and order six. This is impossible in a
dyadic subgroup for p>3.

For v=Cκ the exact row sums and symmetry yield

    max v_j≤2(n−1),  Σv_j=n(n−1)−κ₀=n(n−1).

Hence ||v||²≤2n(n−1)², and substitution yields exactly
E₃≤2n⁴+n³−4n². Also E₂=n²+n(2n−3)=3n²−3n=T₄.
The prior positive opposite-pair decomposition then gives E₃=T₆+R₆,
and subtraction yields R₆≤2n⁴−14n³+41n²−40n.
There is no omitted U/J correction in this specialization: equivalently,
circularity gives r₂(2)=1 and excludes 3∈H, since otherwise the pairs
(−1,3),(3,−1),(1,1) would all sum to 2.

## 5. Focused independent checks and actual witnesses

The full verifier was read and parsed; its recorded run was not repeated.
Its normalization fixes the first coordinate to one using the free H
scaling action. The pair and triple multiset weights restore exactly
the orders of the remaining five positions. Both opposite-freeness and
the property of having some balanced partition are invariant under
these within-block permutations and scaling. The fixed split uses
ρ=−ab/(h₄h₅h₆), with the correct minus sign from negating three entries.

Two small focused computations were made separately from that verifier,
using the recorded subgroup elements after independently checking
primality, cardinality, and that all elements are roots of X^n−1.
Having n distinct roots identifies the actual order-n subgroup.

At p=1073748737,n=256, independently counting nonzero shifted products
reconfirmed 𝒳=0. The displayed tuple

    (1,914267366,972187974,9468345,10646993,240926795)

has all six entries in that subgroup and integer sum 2p. All 15 pairs
were checked not to sum to zero, and all ten product partitions were
checked unbalanced. The quartic window was checked by integer
inequalities. This independently confirms an actual unbalanced circular
relation. This focused check did not independently recount the full
recorded R₆=368640 or E₃=249088000.

At p=7204033,n=64, an independent enumeration of unordered triple
multisets with their permutation weights produced

    B₁=1580032, D₃=1536256, B₁−D₃=43776,
    fixed-partition opposite-free count=25344.

It also found the equal-sum, equal-product nonpermutation pair

    x=(1,265781,3254010),
    y=(2405335,3519792,4798698).

The second triple contains the opposite pair
2405335+4798698=p. This is a concrete reason that equality in the
opposite-free subset bound would be false, within the requested window.

The original recorded checks cover nine groups and all 560 product
ratios; four full coset matrices give 397 nonzero-coset triple-count
checks and four rooted-walk checks. Direct triple sum/product buckets
cover five small groups. The two cases (17,8) and (97,8) are correctly
marked outside the quartic window. Large circular certificates use the
exact shifted collision identity rather than a constructed full coset
matrix. These finite checks establish neither circularity throughout
the family nor an asymptotic bound for the unbalanced remainder.

Only this review file was written. The original proof, verifier,
results, source archive, and central documents were left unchanged.
