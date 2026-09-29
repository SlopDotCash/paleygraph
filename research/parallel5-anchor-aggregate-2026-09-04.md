# Direct anchor weights on the hypergeometric kernel span

**Status: the translated-anchor loss identified in the previous kernel
note is removed for this kernel span. The full necklace aggregate and
Paley conjecture remain unproved.** A direct Kummer twist gives a uniform
bound for any product of distinct anchor weights, without taking the
Mellin L1 norm of the weight. A consequence is uniform quadratic
equidistribution of the kernel span over prescribed adjacency patterns
with nearly half-logarithmically many anchors.

Let p≡1 mod 4 be prime, G=F_p*, χ(0)=0, and

\[
 K_r=S(DS)^{r-1},\quad k_r(t)=K_r(1,t),\quad
 S_{xy}=\chi(x-y),\quad D=\operatorname{diag}(\chi(x)),
 \qquad h_r(t)=p^{-(r-1)/2}k_r(t).
\]

The [rank-kernel construction](parallel3-necklace-2026-09-04.md)
and [pairwise audit](parallel4-kernel-aggregate-2026-09-04.md) give
geometrically irreducible rank-r sheaves F_r with trace −k_r at F_p
points. Their generic weight is r−1. At zero their tame inertia has
only the trivial eigencharacter, and at infinity only the quadratic
eigencharacter. At one the actual middle-extension stalk has dimension
r−1. Their only singularities are 0,1,∞, all tame. These facts follow
from [Katz, Section 2, printed pp.3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf).
The infinity formula on printed page 5 was visually checked again.
Euler–Poincare, the trace formula, and the upper weight bound are the
same imported [Katz-GKM inputs](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf)
used in the preceding proofs. No new source theorem is assumed.

## 1. A product of distinct anchor weights

For a nonempty set T⊂G, put m=|T| and

\[
 w_T(t)=\prod_{a\in T}\chi(t-a),\qquad
 d_T(r,s)=
 \begin{cases}
 mrs,&1\in T,\\
 mrs+r+s-1,&1\notin T.
 \end{cases}
\]

For all positive r,s and every multiplicative character ρ,

\[
 \boxed{\left|\sum_{t\in G}\rho(t)w_T(t)k_r(t)k_s(t)\right|
       \le d_T(r,s)p^{(r+s-1)/2}.}                  \tag{1}
\]

In particular the previously problematic weight χ(1−t) has the bound
rs p^((r+s−1)/2), because χ(−1)=1 and T={1}. No Mellin inversion of
w_T is used, and there is no restriction on the ranks relative to p.

**Proof.** On G_m form the tensor of the actual middle extensions

\[
 \mathcal T=F_r\otimes F_s\otimes\mathcal L_\rho
                 \otimes\mathcal L_{\chi(\prod_{a\in T}(t-a))}.
\]

The last factor has zero stalk at every a∈T, since the roots are
distinct and its local monodromy there is the nontrivial quadratic
character. The displayed tensor therefore has precisely the trace in
(1), including all its zero values. On
U=P¹\(\{0,1,∞\}∪T) its rank is rs, its weight is r+s−2, and its
local monodromies are tame. It has no punctual subsheaf: the tensor of
the actual invariant stalks injects into the invariant stalk of the
tensor, just as in the preceding pass.

We must exclude compactly supported H². If r≠s, rank mismatch between
the irreducible factors does so. If r=s and T contains a≠1, F_r and
F_s are unramified at a but the last Kummer factor has scalar inertia
−1 there. No invariant or coinvariant is possible. This also shows
that restriction to the smaller U has not created a homomorphism.
The remaining case is T={1}, r=s. A possible isomorphism
F_r≅F_r^∨⊗L_(ρ^(-1))⊗L_(χ(t−1)) first forces ρ=1 by inertia at
zero. At infinity the right side then has only the trivial
eigencharacter, whereas the left side has only the quadratic one.
This is impossible. Thus H_c²(U,mathcal T)=0 in every case.

Euler–Poincare on U gives dimension mrs when 1∈T, and (m+1)rs
otherwise. At points in T the actual tensor stalk is zero. When
1∉T its stalk at one has dimension (r−1)(s−1), so the open–closed
sequence subtracts exactly that dimension. Consequently

\[
 H_c^0(G_m,\mathcal T)=H_c^2(G_m,\mathcal T)=0,
 \qquad \dim H_c^1(G_m,\mathcal T)=d_T(r,s).
\]

Its upper weight is r+s−1. The trace formula proves (1). A larger
stalk obtained by middle-extending the tensor on U has not been
substituted for the actual tensor stalk.

## 2. A signed quadratic bound for every coefficient vector

Let R≥1, let U_R=Σ_(r=1)^R r²=R(R+1)(2R+1)/6, and define

\[
 G^{\rho,T}_{rs}=\frac1p\sum_{t\in G}\rho(t)w_T(t)h_r(t)h_s(t).
\]

Since r+s−1≤rs for positive r,s, equation (1) implies the entry bound

\[
 |G^{\rho,T}_{rs}|\le
 \bigl(m+1-1_{1\in T}\bigr)rs/\sqrt p.
\]

For any complex vectors x,y, bound their bilinear form by the last
entrywise majorant, then use
Σr|x_r|≤sqrt(U_R)||x||₂. This proves

\[
 \boxed{\|G^{\rho,T}\|\le
       \frac{(m+1-1_{1\in T})U_R}{\sqrt p}.}         \tag{2}
\]

This does not assume that ρ, or the matrix, is real. As h_r is real,
for every complex coefficient vector z one obtains

\[
 \boxed{\left|\frac1p\sum_t\rho(t)w_T(t)
             \left|\sum_{r=1}^Rz_rh_r(t)\right|^2\right|
       \le\frac{(m+1-1_{1\in T})U_R}{\sqrt p}\|z\|_2^2.}
                                                               \tag{3}
\]

For T={1} this is U_R/√p, compared with the O(R²) bound from the
previous Mellin L1 argument. It tends to zero for R=o(p^(1/6));
all logarithmic ranks are included. With several anchors it suffices
that (m+1)R³=o(√p).

## 3. Exact adjacency cells, with anchors removed

Fix A⊂G, a=|A|≥1, and signs ε_b∈{−1,1}. Define the actual cell

\[
 C(A,\epsilon)=\{t\in G\setminus A:
                   \chi(t-b)=\epsilon_b\text{ for every }b\in A\}.
\]

Let

\[
 c_R=(3R^2-R-2)/2,\quad e=1_{1\in A},\quad
 \Gamma_A=a/2+1-2^{-a}-e/2,
\]
\[
 \mathcal E_{A,R}=
 \frac{2^{-a}c_R+\Gamma_A U_R}{\sqrt p}
 +\frac{2^{1-a}+aU_R/2}{p}.
\]

Then, uniformly over A, its signs, and all complex z,

\[
 \boxed{\left|\frac1p\sum_{t\in C(A,\epsilon)}
                 \left|\sum_{r=1}^Rz_rh_r(t)\right|^2
           -2^{-a}\|z\|_2^2\right|
       \le\mathcal E_{A,R}\|z\|_2^2.}               \tag{4}
\]

To retain the boundary correctly, first write

\[
 P_A(t)=2^{-a}\prod_{b\in A}(1+\epsilon_b\chi(t-b))
       =2^{-a}\sum_{T\subseteq A}
                         \left(\prod_{b\in T}\epsilon_b\right)w_T(t).
\]

Off A this is the exact cell indicator. At a point b∈A it can instead
equal 1/2, and it must be removed. The empty-subset term has error at
most 2^(-a)(c_R/√p+2/p)||z||² by the preceding kernel aggregate.
Summing (3) over nonempty T costs

\[
 2^{-a}\sum_{\varnothing\ne T\subseteq A}
        (|T|+1-1_{1\in T})=\Gamma_A.
\]

Finally |h_r(t)|≤r for t∈G. Away from one this follows from purity
and rank; at one it follows from the mixed upper weight of the
middle extension and its stalk rank r−1. Thus
|Σz_r h_r(t)|²≤U_R||z||². Since 0≤P_A(b)≤1/2, deleting all the
anchor points costs at most aU_R/(2p)||z||². This proves (4).

For any fixed ε>0, R=O(log p), and

\[
 a\le(1/2-\epsilon)\log_2p,
\]

the relative error 2^a mathcal E_(A,R) tends to zero, uniformly in
the anchors and their signs. In particular (4) then gives

\[
 \sum_{t\in C(A,\epsilon)}|\textstyle\sum_r z_rh_r(t)|^2
 =(1+o(1))p2^{-a}\|z\|_2^2.
\]

These are simultaneous quadratic bounds, not just pointwise estimates
for each kernel. They also imply injectivity of evaluation of this
R-dimensional kernel span on each cell once the relative error is
less than one.

## 4. The remaining spectral gap

Equation (4) controls only the span of the R specified kernels. It is
not an operator estimate on all functions supported on the cell, and
does not bound the largest eigenvalue of the Paley adjacency matrix
restricted to that cell. In particular, no argument currently puts its
extremal eigenvectors in this kernel span or proves that this span
approximates them. The full word aggregate in the
[spectral transfer](parallel2-spectral-transfer-2026-09-04.md) still
requires control. Taking R of order p in (2) makes the error useless;
the dimension gap cannot be removed by that substitution.

The new result closes the particular next-anchor estimate left open
in the fourth kernel note. It supplies a reusable analytic input,
but does not prove the desired clique, two-set, or thin-subgroup bound.

## Verification

The [exact verifier](../experiments/parallel5_anchor_aggregate_2026_09_04.py)
checks the full kernels, character-weighted pair sums, and every Mellin
mode through exact convolution-operator positive-semidefinite
certificates. Equality cases are allowed. It separately checks signed
quadratic combinations and the full adjacency-cell expansion, including
the half-valued anchor corrections, with integer-cleared bounds.
Ranks at least p and complex coefficient vectors are included.
[Results](../results/parallel5_anchor_aggregate_2026_09_04.json) record
finite counts and input hashes. The uniform theorem rests on the proof
and its stated primary inputs, not on finite tests alone.
