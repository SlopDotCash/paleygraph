# A uniform three-label family from a twisted Legendre operator

**Status: a proved growing family, using an explicitly identified literature input.** For every prime p≡1 mod 4 and every integer r≥1,

\[
 \boxed{|N(A^rBC)|\le 2p^{(r+3)/2}+1,}
 \tag{1}
\]

where A={0}, B={1}, C={0,1}. The length is k=r+2, so the normalized estimate is at most 2/√p+p^(−k/2−1), uniformly in r. This handles three distinct labels at every length, rather than just another finite-depth check. It does not handle arbitrary runs or imply control of the spectral edge: that requires cancellation in a weighted aggregate of necklaces.

The proof imports the Mellin bound for a particular twisted Legendre family. A second bound below treats A^rBAC at every length, with its large exceptional character removed explicitly. All matrix reductions, exceptional characters, and zero-coordinate corrections are proved below. No novelty or Lean-formalization claim is made.

## Definitions and the exact reduction

Retain S_xy=χ(x−y), D=diag(χ(x)), C_a=D_a S, and K₂=S D S. For clarity, the letter C in a word means the doubleton label, whereas the matrix C_a always has a singleton anchor a. Define the real symmetric matrix

\[
 W=S\circ K_2,
 \tag{2}
\]

where the circle denotes entrywise multiplication. Write C*=C_0|_(F_p*×F_p*) and W*=W|_(F_p*×F_p*). Then the exact identity is

\[
 \boxed{N(A^rBC)=\frac{\operatorname{tr}((C^*)^{r+1}W^*)}{p-1}-(-1)^r.}
 \tag{3}
\]

To see this, let the second anchor vary over b∈F_p. The corresponding word sum is

\[
 N_b=\operatorname{tr}(C_0^r C_b D C_b).
\]

For b≠0, simultaneous scaling of every summation variable gives N_b=χ(b)N_1: there are r+2 cycle edges and r+3 anchor factors. The weight χ(b) is necessary. Hence

\[
 (p-1)N_1=\sum_b\chi(b)N_b.
 \tag{4}
\]

The elementary identity

\[
 \sum_b\chi(b)\chi(x-b)\chi(y-b)=(K_2)_{xy}
\]

shows, by computing matrix entries, that

\[
 \sum_b\chi(b)C_b D C_b=W D S=W C_0.
\]

Thus the right side of (4) is tr(C_0^(r+1)W). There is still a zero-coordinate term in this trace. With coordinate zero first, let **1** denote the vector of p−1 ones. The blocks are

\[
 C_0=\begin{pmatrix}0&0\\\mathbf1&C^*\end{pmatrix},
 \qquad
 W=\begin{pmatrix}0&-\mathbf1^T\\-\mathbf1&W^*\end{pmatrix},
 \qquad C^*\mathbf1=-\mathbf1.
\]

Indeed (K₂)_(0x)=−χ(x) for x≠0. Therefore

\[
 \operatorname{tr}(C_0^{r+1}W)
 =\operatorname{tr}((C^*)^{r+1}W^*)-(p-1)(-1)^r,
\]

which proves (3), including its sign and every zero-coordinate contribution.

## The needed operator bound

The following bound is uniform in p:

\[
 \boxed{\|W^*\|\le2p.}
 \tag{5}
\]

Let

\[
 E(t)=\sum_z\chi(z(z-1)(z-t)),\qquad
 h(t)=\chi(1-t)E(t).
\]

Homogeneity gives W_(x,y)=h(y/x) for x,y≠0. Multiplicative characters form an orthogonal eigenbasis, with eigenvalues

\[
 \lambda_\psi=\sum_{t\ne0}\psi(t)h(t).
 \tag{6}
\]

**Imported input and its precise application.** [Katz, *Sato-Tate theorems for finite-field Mellin transforms*, Theorem 15.1, pp. 52–53, and Corollary 4.2, p. 22](https://web.math.princeton.edu/~nmk/mellin186.pdf) treats the twist TwLeg by χ(1−t), with trace −h(t). The object TwLeg(1)[1] has Mellin dimension two and weight zero. Its unipotent local monodromy at 0 and ∞ makes every nontrivial ψ good; Corollary 4.2 then gives |λ_ψ|≤2p. The zero stalk at 1 agrees with h(1)=0. The shift and Tate twist give trace h(t)/p. This is an imported cohomological theorem.

The trivial character, excluded from that application, is handled exactly. Since S²=pI−J and S**χ**=p e_0−**1** for the vector **χ**=(χ(x)),

\[
 \sum_y W_{xy}=(S K_2)_{xx}
 =((pI-J)D S)_{xx}=1-p\,1_{x=0}.
\]

For x≠0, W_(x,0)=−1. Thus W* **1**=2**1**, and λ_1=2. All eigenvalues in (6) are therefore bounded by 2p, proving (5). No unbounded or omitted exceptional-character contribution remains.

Finally, ||C*||≤√p follows either by compression of D S or from the exact normality identities in [the earlier kernel note](localized-necklace-identities.md). The dimension p−1 trace inequality gives

\[
 |\operatorname{tr}((C^*)^{r+1}W^*)|
 \le (p-1)\|C^*\|^{r+1}\|W^*\|
 \le2(p-1)p^{(r+3)/2}.
\]

Together with (3), this proves (1).

## Label symmetry and growing-depth scope

Rotation and reversal preserve the necklace exactly. Reflection exchanges A,B exactly, and the [previous inversion transfer](parallel-necklace-2026-09-04.md) exchanges B,C with error at most k p^(k/2). Every permutation of three labels uses at most one inversion transposition together with reflections. Therefore, for every ordering of three distinct labels L,M,Q and r≥1,

\[
 |N(L^r M Q)|\le2p^{(k+1)/2}+1+kp^{k/2},\qquad k=r+2.
 \tag{7}
\]

Its normalized bound is at most 2/√p+k/p+p^(−k/2−1). In particular it tends to zero at logarithmic length, and even for k=o(p). This does not count the other words needed at that length, and it does not estimate their weighted sum.

The method uses the adjacency of the two exceptional labels B,C. A weighted second-anchor average then produces precisely S∘K₂, whose Mellin norm is controlled by Katz's fixed twisted Legendre family. If those labels are separated, the corresponding exact average is a product of three kernels rather than this one operator. For positive arc lengths ℓ,m, ℓ+m=k, putting B at one end and C at the other gives

\[
 \boxed{(p-1)N
 =\sum_{x,y}\chi(y)(K_2)_{xy}(K_\ell)_{xy}(K_m)_{xy},
 \qquad K_j=S(D S)^{j-1}.}
 \tag{8}
\]

This follows by the same χ(b)-weighted average and is checked independently. Formula (3) is the case min(ℓ,m)=1. The next section treats min(ℓ,m)=2. For arbitrary larger ℓ,m, the TwLeg norm theorem is not a bound for the new entrywise kernel products. A bound for all gaps must be proved separately; replacing them by W* would be an unjustified extension. More exceptional labels create higher anchor correlations.

## One separating vertex: remove the large quadratic mode exactly

For r≥1 and k=r+3 there is a second uniform estimate:

\[
 \boxed{|N(A^rBAC)|\le2p^{(k+1)/2}+p^{k/2}+p-1.}
 \tag{9}
\]

On F_p* let v=(χ(x)) and define the symmetric matrices

\[
 V_{xy}=\chi(xy)(K_2)_{xy}^2,\qquad
 H=V+pI-pvv^T.
 \tag{10}
\]

Their exact properties are

\[
 \|H\|\le2p^{3/2},\qquad
 N(A^rBAC)=\frac{\operatorname{tr}((C^*)^mH)}{p-1}
             -p t_m+(p-1)(-1)^m,\quad m=r+1.
 \tag{11}
\]

First apply (8) with ℓ=2 and m=r+1. At x,y≠0, use (K_m)_(xy)=χ(x)(C*^m)_(xy). The terms with y=0 vanish, while x=0 contributes (p−1)(−1)^(m−1), since (K_j)_(0y)=(-1)^(j−1)χ(y). Symmetry of V therefore gives

\[
 N(A^rBAC)=\frac{\operatorname{tr}((C^*)^mV)}{p-1}+(-1)^{m-1}.
\]

As C*v=−v and v^Tv=p−1, substitution of (10) proves the trace identity in (11).

Here is the uniform norm argument. V is multiplicative convolution with kernel χ(t)E(t)^2. For a character ψ, put ρ=ψχ. Define a(t)=E(t)^2−p for t≠0,1, and a(1)=1. Then the corresponding eigenvalue is

\[
 \lambda_\psi(V)=\sum_{t\ne0}\rho(t)a(t)
                   +p\bigl((p-1)1_{\rho=1}-1\bigr).
 \tag{12}
\]

Equivalently, for every t≠0, a(t)=E(t)²−p+p·1_(t=1). The positive p correction at t=1 must be retained; using E(t)²−p there would change the operator and the exceptional mode.

The function a is the trace of Sym²(Leg). [Katz, Theorem 15.3, p. 53](https://web.math.princeton.edu/~nmk/mellin186.pdf), gives Sym²(Leg)(3/2)[1] Mellin dimension two and weight zero. Its local monodromy at 0 and ∞ is unipotent, so nontrivial ρ are good. Corollary 4.2 yields |Σρ(t)a(t)|≤2p^(3/2). At t=1 its invariant trace is 1, agreeing with the definition of a.

The remaining character ψ=χ is not estimated by that theorem. Instead, the elementary two-point correlation gives

\[
 \sum_t E(t)^2=p(p-2)-1,\qquad E(0)=-1,
 \qquad \lambda_\chi(V)=p^2-2p-2.
\]

For the first equality, write E(t)=Σ_z f(z)χ(z−t), with f(z)=χ(z(z−1)), and use Σf²=p−2 and Σf=−1. Equation (10) changes this exceptional eigenvalue to −2. For every other ψ it changes (12) to Σρ(t)a(t). Hence every eigenvalue of H is bounded by 2p^(3/2), proving the first part of (11). In particular, applying the smaller bound directly to V would have been false at large p: its exceptional eigenvalue has order p².

The trace inequality and |t_m|≤||C*||^m≤p^(m/2) now give (9), since m=k−2. Its normalized bound is at most 2/√p+1/p+(p−1)p^(−k/2−1), uniformly in r. Rotations, reversal, and permutations of the labels give further instances; permutations cost at most k p^(k/2) by the earlier transfer. Thus both one- and two-edge separating arcs are controlled at logarithmic length. The general case min(ℓ,m)≥3 and weighted aggregates of all words remain outside the proof.

## Exact verification

[The script](../experiments/parallel2_necklace_2026_09_04.py) uses Python integers and NumPy object matrices. It checks the weighted matrix average literally over every anchor, the block matrices, row sums, both growing families through r=16, (8) for arc lengths through six, and label-permuted bounds. It verifies finite-field norm inequalities without floating point: Bareiss elimination certifies positive leading principal minors of 2pI±W* and 4p³I−H². In small fields it also directly sums original necklace tuples. [The result file](../results/parallel2_necklace_2026_09_04.json) records the counts, values, and hashes of local inputs.

All finite tests supplement the uniform proof. The full degree-two necklace conjecture, the weighted aggregate needed for spectral edges, both Paley targets, and the prize reduction remain open.
