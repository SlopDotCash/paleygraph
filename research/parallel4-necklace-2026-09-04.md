# Two arbitrary blocks followed by the third necklace label

**Status: a growing family is proved, with an exact boundary correction.**
Let p≡1 mod 4 be prime, χ(0)=0, A={0}, B={1}, and C={0,1}.
For r,h≥1 put k=r+h+1 and s=min(r,h). Then

\[
 \boxed{|N(A^rB^hC)|\le
 (s+1)p^{(k+1)/2}+|t_s-(-1)^s|,\qquad
 t_j=\frac{\operatorname{tr}(C_A^j)}{p-1}.}
 \tag{1}
\]

In particular, the remainder in (1) is at most p^(s/2)+1. There is no
rank<p restriction. This covers a family with arbitrarily many B labels,
and gives |N(AABBC)|≤3p³+2. It resolves 30 of the 90 length-five words
with multiplicities (2,2,1), after the previously proved label transfers.
The other 60 such words remain unresolved here.

The result does not bound all words or produce cancellation in the weighted
aggregate needed for the spectral edge. The last section gives a precise,
zero-retaining reduction for the remaining three-marked-vertex family and
explains why the present one-variable identification does not cover it.

## 1. Exact weighted average and its zero coordinate

Write S_xy=χ(x−y), D_b=diag(χ(x−b)), C_b=D_bS, D=D_0, and
C*=C_0 restricted to F_p*. For h≥1 define the full p by p matrix

\[
 Q_h=\sum_{b\in F_p}\chi(b)C_b^hD C_b,
 \qquad Q_h^*=Q_h|_{F_p^*\times F_p^*}.
 \tag{2}
\]

For a nonzero second anchor b, scaling every necklace variable by b
multiplies N(A^rB_b^hC_b) by χ(b): the original summand contains
2(r+h+1)+1 character factors. Therefore

\[
 (p-1)N(A^rB^hC)=\operatorname{tr}(C_0^r Q_h).
 \tag{3}
\]

The term b=0 in (2) is zero. Simultaneous scaling shows that Q_h(cx,cy)
=Q_h(x,y) for c≠0. Indeed each C_b is invariant under joint scaling
of b and its two coordinates, and χ(b) and D contribute two cancelling
quadratic factors. Thus Q_h* is multiplicative convolution. It is a
normal real matrix, though symmetry is not assumed.

Each C_b kills the constant vector, so Q_h does too. Its row indexed
by zero is constant off the diagonal. The exact entries are

\[
 (Q_h)_{00}=(p-1)((-1)^h-t_h),\qquad
 (Q_h)_{0y}=t_h-(-1)^h\quad(y\ne0).
 \tag{4}
\]

For the first equality, expand (2) at (0,0), scale z=bv for b≠0,
and use χ(z)²=1 except at z=0. The result is

\[
 (Q_h)_{00}=(p-1)\sum_{v\ne0}(C_1^h)_{0v}\chi(v-1).
\]

For q_v=χ(v−1), S q=p e_1−1 and C_1 q=−q. Also
(C_1^h)_{00}=t_h, because every diagonal entry of C*^h equals t_h.
Removing v=0 gives (4); its second equality follows from the row sum.

Since C*1=−1, for r≥1 the blocks of C_0^r, with zero first, are

\[
 C_0^r=\begin{pmatrix}0&0\\(-1)^{r-1}\mathbf1&(C^*)^r\end{pmatrix}.
\]

Consequently the exact necklace identity, retaining the zero coordinate, is

\[
 \boxed{N(A^rB^hC)=
 \frac{\operatorname{tr}((C^*)^rQ_h^*)}{p-1}
       +(-1)^{r-1}[t_h-(-1)^h].}
 \tag{5}
\]

## 2. The Mellin factorization

All multiplicative characters below, including the trivial character, are
extended by zero at zero. Put

\[
 K_h=S(DS)^{h-1},\qquad k_h(u)=(K_h)_{1u},\qquad
 J(\psi,\chi)=\sum_{v\ne0}\psi(v)\chi(1-v).
\]

The [preceding pass](parallel3-necklace-2026-09-04.md) proves that
k_h is the h-fold multiplicative convolution of χ(1−u).
For a multiplicative character ψ define

\[
 I_\psi(u)=\sum_b\chi(b)\chi(1-b)\psi(u+(1-u)b).
 \tag{6}
\]

The eigenvalue of Q_h* on the vector (ψ(x))_(x≠0) is exactly

\[
 \boxed{\lambda_\psi(Q_h^*)=
 J(\psi,\chi)\sum_{u\ne0}\chi(u)k_h(u)I_\psi(u).}
 \tag{7}
\]

To prove it, expand the row with x=1 in (2), then sum the last coordinate
y against ψ(y). For z≠0,
Σ_(y≠0)ψ(y)χ(z−y)=ψ(z)χ(z)J(ψ,χ); z=0 is killed by D.
The row of C_b^h at 1 vanishes when b=1. For b≠1 set
u=(z−b)/(1−b). If z=b its C_b factor vanishes. Otherwise
(C_b^h)_(1z)=k_h(u). The remaining weight is
χ(u)χ(b)χ(1−b)ψ(u+(1−u)b), giving (7).
In particular the exclusion z=0 is retained by the zero extension of ψ.

## 3. A rank-two hypergeometric family, including u=1

The imported input is [Katz, *G₂ and hypergeometric sheaves*, Section 2,
printed pp. 3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf):
raw sums, disjoint-list irreducibility, rank/weight, the actual stalk trace
at 1, and tame local monodromy. These are the same inputs as in the
[previous proof](parallel3-necklace-2026-09-04.md), together with the
curve trace formula, Euler–Poincaré, and Deligne's bound on H_c¹ weights.
The [archived PDF](../sources/katz-g2-hypergeometric.pdf) has SHA256
`0bf485bcc9dde2afebc8268af679566bcff65aa6d6ed8c9f4b0387d7e9e954c1`.

Suppose ψ≠1. Take the rank-two raw hypergeometric sum with upper
characters (χψ,1) and lower characters (χ,ψ), using the convention
that lower characters are conjugated in the raw summand. Then

\[
 \boxed{\operatorname{RawHyp}[(\chi\psi,1);(\chi,\psi)](u)
       =p I_\psi(u)\quad(u\ne0).}
 \tag{8}
\]

Here is a direct normalization check. Substituting t=(b−1)/b in (6)
and then v=1/t gives

\[
 I_\psi(u)=\psi(-1)\sum_{v\ne0}
 (\chi\psi)(v)\bar\psi(1-v)\psi(1-u/v).
 \tag{9}
\]

Every deleted value has a zero character factor. In particular (9)
holds also at u=1. In the raw sum set x_i=v_i y_i. Its two rank-one
factors give G(ψ) and G(ψ̄), with their factors ψ(−1) and ψ̄(−1).
Their product leaves G(ψ)G(ψ̄) times the convolution on the right of
(9). Since G(ψ)G(ψ̄)=ψ(−1)p, (8) follows with a positive factor p.

The same calculation over every extension E/F_p gives #E times I_(ψ_E),
with norm-compatible characters. Thus a single Tate twist (1) of this
hypergeometric sheaf, denoted G_ψ, has trace −I_ψ over every extension.
It is geometrically irreducible of rank two and weight one on
U=P¹\{0,1,∞}. The lists are disjoint exactly when ψ≠1. Its monodromy
at 1 is a tame unipotent pseudoreflection, and its actual stalk there has
dimension one and trace 1, agreeing with I_ψ(1)=−1. Its semisimple
0-monodromy characters are 1 and χψ.

Let F_h be the normalized rank-h sheaf from the preceding pass. Its
F_p trace is −k_h, its weight is h−1, its 0-monodromy has only the
trivial eigencharacter, and its actual stalk at 1 has dimension h−1.
Form, on G_m, the tensor of actual middle extensions

\[
 T_{h,\psi}=F_h\otimes G_\psi\otimes L_\chi.
\]

Its trace sum is the sum in (7). Its generic rank is 2h, weight h,
and all three missing-point monodromies are tame. H_c²(U,T|_U)=0:
for h≠2 a nonzero global coinvariant would identify irreducibles of
different ranks. For h=2 it would instead require
G_ψ≅F_2^∨⊗L_χ geometrically. At zero the latter has eigencharacters
(χ,χ), whereas G_ψ has (1,χψ). These multisets cannot be equal.
**Thus the equal-rank case has no omitted invariant contribution.**
This also handles ψ=χ: then G_χ is the normalized quadratic rank-two
family itself, but the extra Kummer twist χ prevents an invariant.

T has no punctual subsheaf: the natural map to the middle extension of
its restriction is injective, because the tensor of invariant subspaces
injects into the invariant subspace of the tensor. Euler–Poincaré on U
gives dim H_c¹(U,T)=2h. The open–closed exact sequence removes the
actual stalk at 1, of dimension h−1, giving

\[
 \dim H_c^1(G_m,T)=h+1.
\]

This uses the tensor's actual stalk, not the generally larger invariant
stalk obtained by middle-extending the tensor on U. Deligne's weight
bound and the trace formula consequently give

\[
 \left|\sum_{u\ne0}\chi(u)k_h(u)I_\psi(u)\right|
 \le(h+1)p^{(h+1)/2}\qquad(\psi\ne1).
 \tag{10}
\]

The rank-one boundary h=1 works in the same calculation: the actual
stalk is zero and H_c¹ has dimension two. All local data permit arbitrary
h; in particular h≥p does not require another argument.

## 4. The trivial character and the final bound

The trivial character is not covered by the disjoint-list theorem and
is evaluated directly. For u≠0,

\[
 I_1(u)=-1-\chi(u)+1_{u=1}.
\]

The correction at u=1 is necessary. Since
Σk_h=Σχk_h=(−1)^h and k_h(1)=t_h, (7) gives

\[
 \boxed{\lambda_1(Q_h^*)=2(-1)^h-t_h.}
 \tag{11}
\]

Also I_χ=k_2, so λ_χ(Q_h*)=−t_(h+2), an independently checkable
exact value. The latter follows from
k_2(u^(-1))=χ(u)k_2(u) and convolution at 1.

The elementary Jacobi bound |J(ψ,χ)|≤√p, including its exceptional
characters, combines with (10). For (11) use |t_h|≤p^(h/2), obtained
from ||C*||≤√p. Hence every character satisfies

\[
 \boxed{\|Q_h^*\|\le(h+1)p^{(h+2)/2}.}
 \tag{12}
\]

Because both restricted operators are normal multiplicative convolutions,
the dimension-(p−1) trace inequality and (5) prove (1) first with s=h.
Reflection swaps A,B exactly, and reversal then interchanges the block
lengths, proving N(A^rB^hC)=N(A^hB^rC). Choose the shorter block.

The earlier projective label transfer adds at most kp^(k/2) for any
permutation of A,B,C. The normalized bound for all these permutations is

\[
 \frac{|N|}{p^{k/2+1}}
 \le\frac{s+1}{\sqrt p}+\frac{k}{p}
       +p^{-(k-s+2)/2}+p^{-k/2-1}.
 \tag{13}
\]

It tends to zero uniformly at logarithmic length, and for k=o(√p).
This is a family-wise bound; no cancellation between different words
has been proved.

At length five, t_2=−1 makes (5)'s correction +2 for AABBC.
Its dihedral/reflection orbit has ten words. All label permutations
give 30 distinct words, out of the 90 with multiplicities (2,2,1).
At most four positions change under the required inversion transfer,
so all 30 obey |N|≤3p³+4p^(5/2)+2≤5p³ for p≥5.

## 5. The exact remaining three-gap reduction and its obstruction

For a,b,c≥1 consider w=B A^(a−1) B A^(b−1) C A^(c−1), of
length k=a+b+c. Define

\[
 H(x,y,z)=\sum_v\chi(v)\chi(v-x)\chi(v-y)\chi(v-z).
\]

The same weighted anchor average, followed by scaling its C coordinate
to 1, gives the exact two-variable sum

\[
 N(w)=\sum_{x,y\in F_p}K_a(x,y)k_b(y)k_c(x)H(x,y,1).
 \tag{14}
\]

No factor of p−1 remains: the nonzero C coordinate supplies precisely
the p−1 equivalent scalings. Its zero value was killed by its A factor.
For nonzero x,y, inversion of v and translation by 1 give

\[
 H(x,y,1)=\chi(xy)K_2(x^{-1}-1,y^{-1}-1)-1.
 \tag{15}
\]

The subtracted 1 is the missing inverse-coordinate zero term. This is
a cross-ratio kernel, not a function of y/x. For example, at p=13,
H(2,4,1)=−3 and H(3,6,1)=1, although both ratios equal 2 and all four
nonzero parameters avoid 1. Thus replacing H in (14) by a one-variable
k_j(y/x), or assuming independent hypergeometric factors, is false.

Every zero-coordinate contribution to (14) can also be given explicitly.
Use K_a(0,y)=(−1)^(a−1)χ(y), K_a(0,0)=0, and
H(0,y,1)=p·1_(y=1)−1−χ(y). Restricting (14) to x,y≠0 therefore
requires the correction

\[
 \boxed{(-1)^{a+c-2}[p t_b-2(-1)^b]
       +(-1)^{a+b-2}[p t_c-2(-1)^c].}
 \tag{16}
\]

The other two (2,2,1) bracelet representatives are AABCB and ABABC.
The natural weighted-anchor graph (cycle edges, label spokes, and edge
between the two anchors) is nonplanar for each, so the existing planar
Fourier duality cannot be applied directly. Here are explicit K_(3,3)
minor certificates, numbering cycle vertices 0,…,4 in the displayed
word and anchors 5,6. Contract the edge 0–1. The branch partitions are

- AABCB: left {0,1},{3},{6}; right {2},{4},{5}.
- ABABC: left {0,1},{3},{5}; right {2},{4},{6}.

All nine cross adjacencies are present; surplus edges may be deleted.
This obstructs that direct plane-graph route only. It does not rule out
other identities or prove that these word bounds are impossible.

## Verification

The [verifier](../experiments/parallel4_necklace_2026_09_04.py) constructs
Q_h by its full literal anchor average. It checks (4), (5), (7), (11),
homogeneity, normality, and the norm bound by exact Bareiss positivity.
It evaluates (8) from the raw four-variable hypergeometric definition in
cyclotomic integers, including u=1 and every nontrivial character in its
small test fields. It independently checks (14)–(16), counts the word
orbits, verifies the graph minors, and directly sums original small-field
necklace tuples. Ranks at least p are included. No floating-point value
is used as an acceptance condition.

[Results](../results/parallel4_necklace_2026_09_04.json) record exact
values, counts, and hashes. These finite checks supplement the proof;
the imported cohomological theorems and the still-missing spectral
aggregate are not established by the tests.
