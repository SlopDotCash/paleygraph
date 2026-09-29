# A smaller principal-rank budget and separated boundary correlations

**Status: a source-dependent improvement of the complete aggregate's
logarithmic depth range, checked locally.** This does not reach the
long-depth spectral criterion or improve a clique bound. Independent
mathematical review of the convolution arguments remains outstanding.

This pass audits the distinctions used in the preceding proof and then
uses them more efficiently. Principal constituents have smaller total
rank than the entire raw convolution. The error constituents also have
trivial inertia at the moving endpoint, unlike the principal constituent.
These two facts improve the complete energy bound without discarding
the error, exceptional fibers, constant direction, or anchor correction.

## 1. What the source review checks

The relevant complete pages of Katz,
[*Rigid Local Systems*](../sources/katz-rigid-local-systems.pdf), were
rendered and visually rechecked: PDF pages 49,50,63,64,67,107,109,139.
They support arithmetic convolution duality; exact compact convolution
with the Kummer property-P object; the constant infinity kernel;
geometric inversion on the applicable category; the tame local quotient
and rank rules; and the pure middle image. This is a targeted review
of those inputs, not independent validation of every preceding proof.
No new BBD page or new weight theorem is claimed here.

In particular, compact convolution and middle convolution remain distinct.
The former retains a constant kernel from infinity. Zero extension
after a twist can also create punctual kernels at the selected anchors.
The [ninth-pass weight gap](parallel9-all-degrees-2026-09-05.md) and
[tenth-pass raw accounting](parallel10-word-aggregate-2026-09-05.md)
retain both. Their arithmetic self-duality and word-recovery arguments
are used below with their original hypotheses and review limitation.

## 2. The error has no singularity at the moving endpoint

Fix anchors A of size a and y∉A. For a word w of length m, write

\[
 0\longrightarrow E_w\longrightarrow K_w
          \longrightarrow F_w[1]\longrightarrow0,
 \qquad e_w=\operatorname{rank}H^{-1}(E_w)=r_w-d_w.
                                                               \tag{1}
\]

The finite conductor c_y is additive on this exact sequence. The raw
recurrence gives c_y(K_w)=1, and the moving-point pseudoreflection
gives c_y(F_w[1])=1. Thus c_y(E_w)=0. Each geometric simple
constituent has nonnegative conductor there, so each has conductor
zero: it is lisse at y and is not punctual there. Extensions of these
local systems are lisse at y. Consequently H⁻¹(E_w) extends lisse
across y. This statement concerns the x-variable; its arithmetic
Frobenius data may still depend on the chosen y.

For m≥1, F_w is geometrically irreducible with nontrivial inertia at y.
It therefore cannot be a subquotient of an error constituent. For any
two length-m words u,v, the tensor F_u⊗H⁻¹(E_v) has neither
geometric invariants nor coinvariants on
U=P¹\(A∪{y,∞}). Its rank is d_u e_v, its usual weights are at
most 2m−1, and it is tame at the a+2 punctures. The curve trace
formula and weight bound give

\[
 \left|\frac1p\sum_{x\in U(\mathbb F_p)}
        h_u(x)\,\eta_v(x)\right|
       \le\frac{a d_u e_v}{p},                        \tag{2}
\]

where h_u=p^(−m/2)Tr(F_u), and
η_v=q_v−h_v=p^(−m/2)Tr(H⁻¹(E_v)) with the earlier raw-path sign
normalization. These functions are real. The H_c¹ dimension is
a d_u e_v and its weights are at most 2m, explaining the factor
1/p rather than 1/√p. The assertion also covers e_v=0.

## 3. Sum the principal local monodromies instead of the raw ranks

Let h=2^(a−1), N=2h−1, and let

\[
 D_m=\sum_{|w|=m}d_w,\qquad W_m=N^m.
\]

At a chosen finite anchor, let A_(m,k) and B_(m,k) denote the
numbers of trivial and quadratic Jordan blocks of size k, summed
over all length-m words. Anchor symmetry makes these sums the same
at every anchor. At infinity their roles are reversed: the quadratic
counts are A_(m,k) and the trivial counts B_(m,k). This holds initially
and is preserved because exactly h nonempty labels contain a specified
anchor and exactly h nonempty labels have odd cardinality; the two
complementary counts are both h−1. Thus

\[
 D_m=\sum_{k\ge1}k(A_{m,k}+B_{m,k}),\qquad
 A_{0,1}=1,\quad B_{0,k}=0.
\]

Before convolution, the aggregate twisted block counts are

\[
 G^+_{m,k}=(h-1)A_{m,k}+hB_{m,k},\qquad
 G^-_{m,k}=hA_{m,k}+(h-1)B_{m,k}.
\]

The rank rule and finite quotient rule give the exact recursion

\[
 D_{m+1}=aND_m+NW_m-(a+1)\sum_kG^+_{m,k},             \tag{3}
\]

\[
 A_{m+1,k+1}=G^-_{m,k}\ (k\ge1),\qquad
 B_{m+1,k-1}=G^+_{m,k}\ (k\ge2),                     \tag{4}
\]

with A_(m+1,1) determined by the total rank. In particular,

\[
 A_{m+1,1}=(a-1)ND_m+NW_m
       -[(a+1)h-a]\sum_k A_{m,k}
       -[(a+1)h-1]\sum_k B_{m,k}.                    \tag{5}
\]

These are sums of actual local types, not a classification of sheaves
by their local types. Individual words remain distinct by the earlier
inverse-convolution argument. Equations (3)–(5) compute D_m at any
depth without enumerating N^m words.

## 4. The exact exponential rate of the principal rank sum

Assume a≥2, and let λ_a be the larger root of

\[
 \lambda^2-h(a+a^{-1})\lambda+N=0.                   \tag{6}
\]

Then

\[
 N<\lambda_a<h(a+1)-1=:\lambda_{\rm raw}.             \tag{7}
\]

For the first inequality evaluate (6) at N, obtaining
−hN(a−1)²/a<0. For the second, the root is smaller than the sum
of the two positive roots h(a+a⁻¹), which is at most h(a+1)−1
for a≥2.

Put

\[
 u=\frac{(a-1)N}{\lambda_a-N},\quad
 A_* =\frac{a(u-1)}{a-1},\quad B_* =\frac{u+a}{a-1},
\]

\[
 f_+(k)=uk-A_*(1-a^{-k}),\qquad
 f_-(k)=uk-B_*(1-a^{-k}).                             \tag{8}
\]

Both functions vanish at zero, and f_+(1)=1. We have u>1, and
β:=f_−(1)=((a−1)u−a)/a>0. For a≥3, (7) gives
u>N/h≥7/4>a/(a−1); for a=2, u=(1+√13)/2>2.
Successive differences in (8) then show

\[
 k\le f_+(k)\le uk,\qquad
 \beta k\le f_-(k)\le uk.                            \tag{9}
\]

Define the positive weighted count

\[
 Q_m=\sum_{k\ge1}\bigl(f_+(k)A_{m,k}+f_-(k)B_{m,k}\bigr).
\]

Substitution into (3)–(5), using (6), gives the exact scalar identity

\[
 Q_{m+1}=\lambda_a Q_m+N^{m+1},\qquad Q_0=1,
 \quad Q_m=\frac{\lambda_a^{m+1}-N^{m+1}}{\lambda_a-N}. \tag{10}
\]

For an explicit coefficient check of this substitution, set
c_A=(a+1)h−a and c_B=(a+1)h−1. The coefficient of A_(m,k)
on the left of (10), after using (5), is

\[
 (a-1)Nk-c_A+h f_+(k+1)+(h-1)f_-(k-1)
       =\lambda_a f_+(k).
\]

For B_(m,k), it is

\[
 (a-1)Nk-c_B+(h-1)f_+(k+1)+h f_-(k-1)
       =\lambda_a f_-(k).
\]

These identities include k=1 because f_−(0)=0. They follow by
comparing the coefficients of k, 1 and a^(−k) in (8) and reducing
λ_a² with (6). The remaining term is NW_m f_+(1)=N^(m+1).

Let c=min(1,β)>0. Equations (9)–(10) prove, with explicit constants,

\[
 \boxed{\frac{Q_m}{u}\le D_m\le\frac{Q_m}{c},
             \qquad D_m\asymp_a\lambda_a^m.}         \tag{11}
\]

The comparison constants depend only on a. Write

\[
 \Psi_a=\lambda_a^2/N,\qquad
 \Theta_a=\lambda_{\rm raw}^2/N,\qquad N<\Psi_a<\Theta_a.
\]

For a=2,

\[
 \lambda_2=(5+\sqrt{13})/2,\quad
 \Psi_2=(19+5\sqrt{13})/6\approx6.17129,
 \qquad\Theta_2=25/3\approx8.33333.
\]

The first principal-rank sums are 1,6,33,156,729,3258,14505.
The raw-rank sums at the same depths are 1,7,41,223,1169,6007,30521.

## 5. The improved full energy bound

Use the actual normalized transfer matrix and all notation of the
[complete-energy proof](parallel11-full-energy-2026-09-05.md):
T=(S/√p−J/p)(2^a E−I)/√N, H_j=||T^j||_F². Let

\[
 L_j=\sum_{|w|=j}r_w
     =\frac{a\lambda_{\rm raw}^j-N^j}{a-1},\quad
 E_j=L_j-D_j,\quad F_j=L_j^2/N^j.
\]

The principal Gram estimate, applied to the actual equal coefficients,
costs (aD_j²/N^j+1)/√p. Equation (2) bounds the principal/error
cross term, while the pointwise lower-weight estimate bounds the
error square. The generic-column bound is therefore

\[
 V_j=1+\frac{aD_j^2/N^j+1}{\sqrt p}
               +\frac{2aD_jE_j+E_j^2}{pN^j}.          \tag{12}
\]

No arithmetic value of a boundary constituent has been guessed.
Restore the omitted finite rows and exceptional columns exactly as
before. Keeping the constant direction lost by multiplication by S,
and then the full-power J and anchor correction, gives

\[
 U_j=pV_j+4(a+1)F_j+(a+1)N^j,
\]

\[
 d_0=\frac{2^a}{2\sqrt N}\sqrt{a(1-1/p)}
              +\sqrt{\frac{|C|N^2+p-|C|}{pN}},
\]

\[
 \boxed{|\operatorname{tr}T^{2j}|\le H_j
     \le\left(\sqrt{U_j}+j d_0N^{(j-1)/2}\right)^2.}  \tag{13}
\]

Only the generic-column budget changed; all other terms from pass
eleven remain present. Either this bound or the preceding one may
be used at a finite parameter; no claim of a strict improvement at
every such parameter is needed.

For fixed a≥2 and 0<ε<1/2, (11)–(13) give H_j≤(1+o(1))p
uniformly over anchor cliques when

\[
 \boxed{j\le\min\left\{
       \frac{(1/2-\epsilon)\log p}{\log\Psi_a},
       \frac{(1-\epsilon)\log p}{\log\Theta_a}
                              \right\}.}             \tag{14}
\]

Indeed, D_j²/N^j=O_a(Ψ_a^j), and the error and finite boundary
budgets are O_a(Θ_a^j). The two conditions make their respective
ratios to √p and p tend to zero. The perturbation is o(√p)
because N<Θ_a and j=O_a(log p). Both upper limits in (14) exceed
the earlier (1/2−ε)log p/log Θ_a. For a=2 the first limit
is controlling: Ψ₂²>Θ₂, and the leading depth constant improves
from 1/(2log(25/3)) to 1/(2log((19+5√13)/6)).

This is still a bounded multiple of log p. It does not meet the
required j/log p→∞, remove the [bootstrap obstruction](parallel12-bootstrap-obstruction-2026-09-05.md),
or imply a new clique bound. The arbitrary-two-set, uniform subgroup,
exceptional-prime, and exact official Reed–Solomon prize obligations
remain separate and unproved.

Equation (11) also shows precisely why further improvements in counting
these ranks cannot finish this route. The principal/principal error
budget aD_j²/(N^j√p) is comparable to Ψ_a^j/√p. Since Ψ_a>1,
it diverges when j/log p→∞ at fixed a. This is a lower bound on
the size of the current error allowance, not on the true correlations
or spectral trace. Closing the goal requires cancellation beyond this
rank-based allowance, or a different argument for the extreme eigenvalues.

## 6. Reproduction and limits of the checks

The [verifier](../experiments/parallel13_principal_budget_2026_09_05.py)
compares (3)–(5) with full word enumeration, checks the mirrored
infinity counts, and checks that every retained raw-error constituent
has zero conductor at y. It verifies the weighted coefficient identities
and (10) in exact quadratic fields through depth 80 at a=2,…,10.
It compares (13) with the exact full energies retained from pass eleven,
and includes a nonvacuous p=10009 example. Results are in
[the exact output](../results/parallel13_principal_budget_2026_09_05.json).
The algebraic derivation proves the arbitrary-depth rank identity;
finite checks alone do not prove weights, geometric irreducibility,
or the asymptotic claim. No independent or formal verification is claimed.
