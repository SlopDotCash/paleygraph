# Arbitrary gaps between two exceptional necklace labels

**Status: the general two-gap family is proved, using Katz's hypergeometric theorem and Deligne's weight bound.** Let p≡1 mod 4 be prime, let A={0}, B={1}, C={0,1}, and let

\[
 w_{j,m}=B A^{j-1} C A^{m-1},\qquad k=j+m,\qquad s=\min(j,m),
 \quad j,m\ge1.
\]

Then

\[
 \boxed{|N(w_{j,m})|\le(s+1)p^{(k+1)/2}+1.}
 \tag{1}
\]

Thus the two exceptional labels need not be adjacent or separated by only one vertex. The constant is explicit and uniform in both gaps, including ranks/gaps at least p. This controls this word family at logarithmic depth; it does not control all words or their weighted aggregate, and it does not prove a Paley conjecture.

## Primary inputs and source boundary

The hypergeometric input is [Katz, *G₂ and hypergeometric sheaves*, Section 2, printed pp. 3–5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf). Page 3 defines the raw sum; page 4 states disjoint-list irreducibility, rank, weight, the trace formula including the stalk at 1, and tame pseudoreflection monodromy at 1. Page 5 specifies the endpoint monodromies. These statements review *Exponential sums and differential equations*, Chapter 8. The source permits arbitrary list lengths and repeated characters, with no rank<p hypothesis. We apply the disjoint lists (1,…,1) and (χ,…,χ) of equal length j; χ≠1 because p is odd.

The other input is the curve cohomology/weight machinery, as used in [Katz, *Gauss sums, Kloosterman sums, and monodromy groups*, §2.3.1–2.3.2, printed p. 32; §3.6, printed pp. 38–40](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf): Euler–Poincaré, the trace formula, and Deligne's bound on weights of H_c¹. The underlying theorems are imported. The finite-sum identification, exact Betti count for the tensor used here, and matrix consequences are derived below.

The PDF locations are zero-based pages 2–4 for the first source; the second PDF contains two printed pages per PDF page, with §2.3 on zero-based page 20 and the weight argument of §3.6 on page 24. The first source is [archived locally](../sources/katz-g2-hypergeometric.pdf). Local source files are managed centrally; this lane writes only its note, script, and result. The result records the URLs and independently fetched byte hashes.

## Kernels and a hypergeometric identification without exceptional quotients

Set

\[
 S_{xy}=\chi(x-y),\quad D=\operatorname{diag}(\chi(x)),
 C_0=DS,\quad K_j=S(DS)^{j-1},\quad C^*=C_0|_{F_p^*\times F_p^*}.
\]

The function f(t)=χ(1−t) on F_p* is the multiplicative convolution kernel of C*. Its j-fold convolution is

\[
 k_j(t)=(K_j)_{1t}
       =\sum_{u_1\cdots u_j=t}\prod_{i=1}^j\chi(1-u_i).
 \tag{2}
\]

Indeed C*_(xy)=f(y/x), and the restriction of K_j is D* (C*)^j. This proves (2) directly from matrix multiplication. It also gives the Mellin transform J(ψ,χ)^j, with the usual zero-at-zero convention for all multiplicative characters. The argument below does not divide by a Gauss sum associated to ψχ and therefore does not discard the exceptional characters.

Fix a nontrivial additive character e_p and let G_p=Σ_(y≠0)χ(y)e_p(y). Let RawHyp_j(t) be Katz's raw hypergeometric sum with j trivial upstairs characters and j quadratic downstairs characters. In its definition put x_i=u_i y_i. Then

\[
 \begin{aligned}
 \operatorname{RawHyp}_j(t)
 &=\sum_{u_1\cdots u_j=t}
   \prod_i\left(\sum_{y_i\ne0}\chi(y_i)e_p((u_i-1)y_i)\right)\\
 &=G_p^j k_j(t).
 \end{aligned}
 \tag{3}
\]

The inner sum is zero when u_i=1 and is G_p χ(u_i−1) otherwise. Here χ(−1)=1. Formula (3) includes t=1 exactly, not just smooth fibers.

Let H_j be Katz's geometrically irreducible rank-j middle-extension sheaf of type (j,j). Its trace at every t∈F_p*, including t=1, is −RawHyp_j(t). Define F_j by the geometrically constant twist whose Frobenius eigenvalue is G_p^(−j). It follows that

\[
 \operatorname{tr}(\operatorname{Frob}_t\mid F_j)=-k_j(t),
 \qquad t\in F_p^*.
 \tag{4}
\]

F_j has weight j−1 on U=P¹\{0,1,∞}, since the hypergeometric weight is 2j−1 and |G_p|=√p. It remains geometrically irreducible. It is tame at all three missing points, and its invariant subspace at 1 has dimension j−1.

The constant twist is defined over F_p, not chosen afresh at each extension. More explicitly, for E/F_p of degree d, with G_E the norm/trace-compatible quadratic Gauss sum,

\[
 \operatorname{tr}(\operatorname{Frob}_{t,E}\mid F_j)
 =-\left(\frac{G_E}{G_p^d}\right)^j k_{j,E}(t).
 \tag{5}
\]

Thus (4) is a statement at F_p points; no assertion that its minus sign persists unchanged over every extension is needed. The ratio in (5) has absolute value one. Purity and geometric monodromy come from the sheaf and its single specified constant twist, not from guessing signs in convolution powers.

## The gap operator and its actual boundary stalk

On F_p* define the real symmetric matrix

\[
 (V_j)_{xy}=\chi(xy)(K_2)_{xy}(K_j)_{xy}.
 \tag{6}
\]

By homogeneity it is multiplicative convolution with kernel χ(t)k₂(t)k_j(t). Therefore its eigenvalues are

\[
 \lambda_\psi(V_j)=\sum_{t\ne0}\psi(t)\chi(t)k_2(t)k_j(t).
 \tag{7}
\]

For j≠2 we prove

\[
 \boxed{\|V_j\|\le(j+1)p^{(j+1)/2}.}
 \tag{8}
\]

Let G_m=P¹\{0,∞}, and form on G_m the tensor product

\[
 T_{j,\psi}=F_2\otimes F_j\otimes L_{\psi\chi}.
\]

It is the tensor product of the actual middle extensions, rather than the middle extension of the tensor product on U. By (4), its F_p trace sum is exactly (7). Its generic rank is 2j, its weight on U is j, and it is tame at 0,1,∞. Its actual stalk at 1 has dimension

\[
 \dim (T_{j,\psi})_1
 =\dim(F_2)_1\dim(F_j)_1=1\cdot(j-1)=j-1.
 \tag{9}
\]

In particular, tensoring and middle extension have not been interchanged. That interchange can create an additional invariant at 1 when j is even.

The rank mismatch gives H_c²(U,T_(j,ψ)|_U)=0 for j≠2. By duality a nonzero such group would give a nonzero homomorphism between the irreducible rank-two factor and a Kummer twist of the dual rank-j factor. Nonzero homomorphisms between irreducibles are isomorphisms, impossible when their ranks differ. This works for every ψ, including the trivial and quadratic characters. Also H_c⁰(U,T|_U)=0. Since χ_c(U)=−1 and every Swan conductor is zero, Euler–Poincaré gives

\[
 \dim H_c^1(U,T|_U)=2j.
\]

The natural map T→ι_*(T|_U), for ι:U→G_m, is injective. At 1 this is the elementary injection

\[
 (F_2)^{I_1}\otimes(F_j)^{I_1}
 \hookrightarrow(F_2\otimes F_j)^{I_1}.
\]

Thus T has no nonzero section supported at finitely many points, and H_c⁰(G_m,T)=0. The open–closed exact sequence now supplies

\[
 0\longrightarrow T_1\longrightarrow H_c^1(U,T|_U)
 \longrightarrow H_c^1(G_m,T)\longrightarrow0.
\]

It follows, retaining the actual boundary stalk, that

\[
 \boxed{\dim H_c^1(G_m,T)=2j-(j-1)=j+1.}
 \tag{10}
\]

Deligne's bound makes H_c¹(U,T|_U) mixed of weights at most j+1. Its quotient in (10) has the same upper weight bound. The trace formula therefore bounds (7) by (j+1)p^((j+1)/2), proving (8). The only rank-based exception is j=2; it is treated explicitly below.

For j=1 the same argument is valid without an omitted endpoint: F₁ has rank one and zero stalk at 1, T has rank two, and H_c¹ has dimension two. It gives ||V₁||≤2p. No condition j<p was used anywhere; the local data and dimension calculation are explicit for every j.

## Exact necklace trace and the zero coordinate

The weighted second-anchor average from [the preceding pass](parallel2-necklace-2026-09-04.md) gives

\[
 (p-1)N(w_{j,m})
 =\sum_{x,y\in F_p}\chi(y)(K_2)_{xy}(K_j)_{xy}(K_m)_{xy}.
 \tag{11}
\]

For x,y≠0 use (K_m)_(xy)=χ(x)(C*^m)_(xy). The y=0 terms vanish. For x=0,y≠0, use (K_a)_(0y)=(-1)^(a−1)χ(y); their total is (p−1)(−1)^(k−1). Therefore

\[
 \boxed{N(w_{j,m})
 =\frac{\operatorname{tr}((C^*)^mV_j)}{p-1}+(-1)^{k-1}.}
 \tag{12}
\]

Since ||C*||≤√p, (8) and (12) give

\[
 |N(w_{j,m})|\le(j+1)p^{(k+1)/2}+1\qquad(j\ne2).
 \tag{13}
\]

Reversal exchanges the two arcs and preserves the necklace. Thus one can choose j=s=min(j,m).

## The equal-rank exception j=2

Write v=(χ(x))_(x≠0). The exact operator correction proved in [the second pass](parallel2-necklace-2026-09-04.md) is

\[
 H=V_2+pI-pvv^T,\qquad \|H\|\le2p^{3/2},
 \qquad V_2v=(p^2-2p-2)v.
\]

It uses the actual Sym²(Legendre) stalk at 1, with trace correction a(t)=k₂(t)²−p+p·1_(t=1). The large quadratic-character eigenvalue is retained and subtracted, not put under the smaller norm bound. Since C*v=−v, the exact consequence of (12) is

\[
 N(w_{2,m})
 =\frac{\operatorname{tr}((C^*)^mH)}{p-1}
       -p t_m+(p-1)(-1)^m.
\]

Thus, with k=m+2,

\[
 |N(w_{2,m})|\le2p^{(k+1)/2}+p^{k/2}+p-1.
\]

If the shorter arc is two then k≥4. For p≥5 this is at most 3p^((k+1)/2)+1, because p^(k/2)(√p−1)≥p−2. This supplies (1) in the sole exceptional case and completes the proof.

## Scope, label permutations, and verification

The result covers every arrangement with exactly one B, exactly one C, and all remaining labels A. For k≥3 these are genuine three-label words. The [previous label transfer](parallel-necklace-2026-09-04.md) also covers every permutation of the three labels, at the cost of at most kp^(k/2). Hence the normalized permuted-family bound is

\[
 \frac{|N|}{p^{k/2+1}}
 \le\frac{s+1}{\sqrt p}+\frac{k}{p}+p^{-k/2-1}.
\]

It tends to zero uniformly when k=o(√p), including logarithmic length. The two singleton exceptional positions may now be separated by arbitrary arcs. Words with more occurrences of both exceptional labels lead to higher anchor correlations; neither those individual estimates nor the weighted sum needed for spectral edges is proved here.

[The script](../experiments/parallel3_necklace_2026_09_04.py) verifies (2), (3), (6), (11), and (12) using exact integer matrices and literal sums. Hypergeometric sums are checked as cyclotomic integers using coefficient vectors modulo the p-th cyclotomic polynomial, including t=1. The operator bound is checked by exact Bareiss positivity certificates for (j+1)²p^(j+1)I−V_j² when j≠2; the known corrected operator is used for j=2. Checks include ranks j≥p. [Results](../results/parallel3_necklace_2026_09_04.json) retain check counts, values, and local/source hashes. These tests supplement the uniform proof; they do not certify the imported cohomological theorems or prove anything about the unresolved word families.
