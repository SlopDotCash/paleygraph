# A limited character-sum connection to Reed–Solomon lists and MCA

The Paley and prize goals remain unproved. This note supplies an explicit
connection for one family of words: a bound on prime-field subgroup periods
controls subset sums, which count nearby Reed–Solomon codewords and the bad
scalars of a particular monomial pair. Root lifting also gives a finite list
lower bound for the pinned official profile. None of these statements controls
the worst-case MCA error at a fixed positive gap below capacity.

Subset-sum counting with character sums and its decoding connections are
established subjects. Relevant primary sources are
[Li–Wan, *On the subset sum problem over finite fields*](https://arxiv.org/abs/0708.2456)
and [Zhu–Wan, *An Asymptotic Formula For Counting Subset Sums Over Subgroups Of Finite Fields*](https://arxiv.org/html/1101.0289).
The latter gives stronger specialized estimates using a distinct-coordinate
sieve. We give the elementary bound needed here, with no novelty claim.

## 1. From additive characters to fixed-cardinality subset counts

Let `D⊂F_p`, where p is prime and `|D|=n<p`. Put

\[
M(D)=\max_{b\ne0}\left|\sum_{x\in D}e_p(bx)\right|,\qquad
N_s(c)=\#\{S\subseteq D:|S|=s,\ \sum_{x\in S}x=c\}.
\]

For an integer `U≥max(1,M(D))`, `0≤s≤n`, and `t=min(s,n-s)`, we have

\[
\boxed{\left|N_s(c)-\frac1p\binom ns\right|
\le\frac{p-1}{p}\binom{U+t-1}{t}.}
\tag{S}
\]

**Proof.** Let
`F_b(z)=∏_{x∈D}(1+z e_p(bx))`. Fourier inversion gives

\[
N_s(c)=\frac1p\sum_{b\in F_p}e_p(-bc)[z^s]F_b(z).
\]

The principal term is `binom(n,s)/p`. For `b≠0`, the formal logarithm
through degree s is

\[
\log F_b(z)=\sum_{j=1}^s
 \frac{(-1)^{j-1}}j\left(\sum_{x\in D}e_p(jbx)\right)z^j
 \pmod {z^{s+1}}.
\]

Since `s<p`, every `jb` in this expression is nonzero. Taking absolute
values term by term in the exponential series bounds the coefficient by

\[
|[z^s]F_b(z)|\le[z^s]\exp\left(U\sum_{j\ge1}\frac{z^j}j\right)
=[z^s](1-z)^{-U}=\binom{U+s-1}s.
\]

The triangle inequality proves (S) with s in place of t. Complementation
gives `N_s(c)=N_{n-s}(ΣD-c)`, establishing the stated version. In particular,

\[
N_s(c)\ge L:=\max\left(0,\left\lceil
\frac{\binom ns-(p-1)\binom{U+t-1}t}{p}\right\rceil\right).
\tag{S+}
\]

If `L>0`, every field element occurs as a sum of s distinct elements of D.
This is a prime-field statement. Replacing p by an extension size while
keeping D in the base field would be false; the sums still lie in F_p.

For context, the usual Gauss-sum bound gives `M(H)≤√p` for a multiplicative
subgroup H: expand its indicator in the m multiplicative characters trivial
on H. The principal Gauss sum is -1 and the others have magnitude √p, so
`M(H)≤(1+(m-1)√p)/m≤√p`. For `H=F_p*`, the exact maximum is 1.

## 2. Exact list and MCA identities at one step below capacity

Let `F_p⊂F_q`, `D⊂F_p`, `|D|=n`, and `1≤k<n`. Let
`C=RS[F_q,D,k]` denote evaluations of polynomials of degree strictly less
than k. Set `s=k+1` and `δ₀=1-s/n`.

For the received word `v=X^{k+1}`, the list at relative distance at most
δ₀ has cardinality exactly

\[
\boxed{|\operatorname{List}_C(v,\delta_0)|=N_{k+1}(0).}
\tag{L}
\]

Indeed, if u is a nearby codeword, the monic polynomial `X^{k+1}-u` has
at least k+1 distinct roots in D. Its degree forces it to equal
`Q_S(X)=∏_{x∈S}(X-x)` for a unique size-(k+1) subset S. Its X^k coefficient
is zero, so `ΣS=0`. Conversely this condition cancels both the degree k+1
and degree k coefficients of `X^{k+1}-Q_S`, leaving degree less than k.
The root set is exactly S, and distinct S give distinct codewords. The
argument works over F_q and forces these particular codewords to have
base-field coefficients.

For the pair `f=X^{k+1}`, `g=X^k`, the MCA-bad scalars at the same radius are

\[
\boxed{\{\gamma\in F_q: \gamma\text{ is bad for }(f,g)\}
=\{-\sum S:S\subseteq D,\ |S|=k+1\}\subseteq F_p.}
\tag{M}
\]

The same factorization now has X^k coefficient γ, which must be `-ΣS`.
The no-joint-explanation requirement in the MCA definition is automatic:
`X^k` cannot agree with a degree-<k polynomial on k+1 distinct points.
Any larger witness set contains such a subset. Conversely the displayed
factorization provides a size-(k+1) witness for each listed scalar.

Thus `L>0` in (S+) implies that this pair has exactly p bad scalars and
MCA probability `p/q`. In the prime field this is 1. It is **the probability
for this pair**, not an upper bound on the supremum over arbitrary pairs.

For `n→∞`, `s/n→ρ∈(0,1)`, `log p=o(n)`, and `M(D)=o(n)`, (S) gives
uniformly in c

\[
N_s(c)=\frac1p\binom ns(1+o(1)).
\]

To see this, take `U=max(1,ceil M)`. The error binomial has logarithm
`o(n)` since `U=o(n)`, while `log binom(n,s)=n h(ρ)+o(n)` for the binary
entropy in natural units. The relative error in (S) is therefore
`exp(-n h(ρ)+o(n))`. In the quartic window, the conjectured Paley bound
implies these hypotheses, but so would any fixed power saving
`M≤np^{-ε}`. Square-root cancellation is not necessary for this consequence.
The radius here approaches capacity with gap `1/n`; it is not a
counterexample at any fixed positive capacity gap.

## 3. Root lifting gives a list lower bound farther below capacity

Let `H=μ_n⊂F_p*`, let `a|n` and `a|k`, with `k+a≤n`, and put
`G=H^a=μ_N`, `N=n/a`, `s'=k/a+1`. For any size-s' subset `T⊂G` of sum zero,
define

\[
Q_T(X)=\prod_{z\in T}(X^a-z),\qquad
v(X)=X^{k+a},\qquad u_T=v-Q_T.
\]

The degree-k coefficient of Q_T is `-ΣT=0`; all exponents are multiples
of a. Hence `deg u_T≤k-a<k`. Each z has exactly a preimages in H, so
Q_T has exactly `a s'=k+a` distinct roots there. Distinct T give distinct
u_T, proving

\[
\boxed{|\operatorname{List}_{RS[F_q,H,k]}(X^{k+a},\delta_a)|
\ge N_{s'}^G(0),\qquad \delta_a=1-(k+a)/n.}
\tag{Lift}
\]

For `a>1` we assert an injection, not a classification of every nearby
codeword. The same lower bound holds for any number of interleaved rows:
place u_T in one row and zero in every other row, using v in the same
row of the center. Column agreement is unchanged.

There is a useful exact obstruction to extrapolating this construction.
Fix dyadic `N≥4`. An odd-sized subset of the complex Nth roots of unity
cannot sum to zero. Its 0/1 polynomial F of degree below N would otherwise
be divisible by `Φ_N=X^{N/2}+1`, forcing equal coefficients in each
opposite pair, and thus an even number of terms. For an odd subset of
size t, its nonzero integer cyclotomic norm consequently satisfies

\[
1\le|\operatorname{Norm}_{\mathbb Q(\zeta_N)/\mathbb Q}F(\zeta_N)|
\le t^{N/2}.
\]

If that subset becomes a zero sum after reduction to `μ_N⊂F_p`, the
evaluation homomorphism at a primitive root modulo p shows that p divides
this norm. Therefore no such zero sum exists when `p>t^{N/2}`.

At rate 1/2, the complement of a size-`N/2+1` zero-sum subset has odd
size `t=N/2-1`. Thus (Lift) produces no such subsets once
`p>(N/2-1)^{N/2}` for fixed N. A fixed gap `a/n=1/N` keeps N fixed.
This is a limitation of this construction, not a list-decoding theorem
for arbitrary received words.

## 4. Exact finite consequence for the pinned official profile

Use the [audited profile](official-profile-and-trace.md) at contract commit
`b34c0131cfa36b51111521541d7d3e35c8791082`: `p=2130706433`, `q=p^6`,
`n=262144`, scalar dimension `k=131072`, eight interleaved rows.

Take `a=32`, so `N=8192`, `s'=4097`, `t=4095`, and
`δ_a=4095/8192=1/2-1/8192`. An exact fourth-moment computation gives

\[
E_2(G)=203464704,\qquad
Q_2(G)=\frac{pE_2(G)-N^4}{N}=52370599862533.
\]

The latter is the sum of fourth powers of absolute nonprincipal periods
over multiplicative cosets; hence `M(G)^4≤Q₂(G)<2691^4`.
The energy computation uses

\[
\kappa_c=\#\{x\in G\setminus\{-1\}:1+x\in cG\},\qquad
E_2(G)=N^2+N\sum_c\kappa_c^2.
\]

The nonzero kernel entries consist of one 1, 4029 entries of 2, and 33
entries of 4. Independently, direct membership tests give
`#{(x,y)∈G²:1+x-y∈G}=24837`, whose product with N is the same energy.
Cosets are identified exactly by the label `(1+x)^N mod p`.

Substituting `U=2691` into (S+) yields the integer bound

\[
L=\left\lceil\frac{\binom{8192}{4097}
 -(p-1)\binom{6785}{4095}}p\right\rceil
\ge 2^{8154}>q.
\]

By (Lift), the scalar list at `X^{131104}` at radius `4095/8192` has
at least this many codewords. So does the twice-interleaved official code
(sixteen scalar rows), by the injection above. At the pinned ArkLib source,
[`coe_certifiedGammaError`](https://github.com/Verified-zkEVM/ArkLib/blob/e65197892890b8fd9b0dc05b8980273cf1d595cc/ArkLib/ProofSystem/ToyProblem/Impl/IRS.lean#L669)
identifies the lower-track certificate as
`Γ(δ)=MCA(C,δ)+Λ(C^{⋈2},δ)/q`. Nonnegativity of MCA and monotonicity
of list size therefore imply

\[
\boxed{\Gamma(\delta)>1>2^{-128}
\quad\text{for every }\delta\ge4095/8192.}
\]

Consequently the pinned
[`ProtocolClaim` lower certificate](https://github.com/proximity-prize/proximity-prize/blob/b34c0131cfa36b51111521541d7d3e35c8791082/ProximityPrize/Benchmark/TargetLower.lean)
cannot use such a radius. This calculation **alone** does not lower-bound
winning-set soundness: Γ is an upper certificate and can exceed 1. A later
[list-to-winning-set argument](list-to-winning-set.md) uses the same
single-word list to prove a winning density greater than `1−2^(-7968)`
throughout this radius suffix below minimum distance. That is a separate
ordinary proof, not an inference from Γ>1. The official dependency closure
has not been rebuilt or kernel-audited here.

For comparison, applying (S+) directly to H gives at least `2^{262103}`
subsets of size `k+1` for every prescribed base-field sum. By (M), the pair
`X^{k+1},X^k` at radius `1/2-1/n` has probability exactly `p/q=p^{-5}`,
which is below `2^{-128}` because `p^5>2^{128}`. Its small MCA probability
coexists with a huge list. No worst-case MCA upper bound follows.

## Reproduction and remaining obligation

Run `python3 experiments/subset_sum_list_certificate.py` from the project.
The [script](../experiments/subset_sum_list_certificate.py) verifies the
eleven pinned source checksums, the base prime and generator orders, exact
kernel energies for nine lifting factors, and the independent pair count
at a=32. It saves [JSON results](../results/subset_sum_list_certificate.json)
with integer bit lengths and hashes, avoiding enormous decimal output.

Small-field checks cover all 3490 size/residue combinations at `(p,n)=(17,16)`
and `(97,32)`. They compare the first case to exhaustive enumeration of
11440 subsets, reconstruct 114 polynomial MCA witnesses, and verify the
exact zero-sum list sizes 672 and 5832096. The root-lifting test exhausts
the image subsets at a=2,4,8 for the length-32 code, verifying all 64
codewords supplied by a=2 and their complement norms divisible by 97.
These are ordinary proofs with exact computational checks, not new Lean
certificates.

The missing prize bridge still concerns arbitrary received words and
arbitrary MCA pairs at the required radius. For a general gap a without
the root-lifting restriction, cancellation of the top a coefficients of
a vanishing polynomial imposes multiple elementary-symmetric conditions
on its root subset. A linear additive-character estimate only addresses
the first such condition. Neither (S) nor the finite list example provides
the higher-condition estimate or the uniform centered period bound (SG).
