# Independent review of the quartic-star reformulation

**Verdict: no mathematical defect found in the stated proposition, overlap
ledger, zero-row corrections, or equivalence directions.** The inequalities
hold under the stated conditions n>=6 and n^4<=p, using the ordinary
squarefree multiplicative-character Weil bound. The note correctly describes
an equivalent upper-bound problem; it does not establish the required
upper bound. Its elliptic summands do not remove the degree-six difficulty,
which returns exactly in their pairwise correlations.

This is a separate-agent mathematical and code review. It is not human
refereeing, a Lean certificate, or independent certification of the Weil
bound. The reviewed verifier was read but **not executed again**: no
specific concern justified rerunning it. Its recorded input hashes match
the current files.

## 1. Scope and pinned inputs

The principal reviewed document is
[the quartic-star note](parallel21-positive-upper-review-2026-09-05.md).
The reviewer also read its
[verifier](../experiments/parallel21_positive_upper_review_2026_09_05.py),
[recorded result](../results/parallel21_positive_upper_review_2026_09_05.json),
and the relevant transfer statements in passes 20 and 21.

| Input | SHA-256 |
|---|---|
| research/parallel21-positive-upper-review-2026-09-05.md | 832805eb897c8a8235d7d04909adcb4c99c263a8d5b86f7261266ae5c198e9ea |
| experiments/parallel21_positive_upper_review_2026_09_05.py | 3fc374343b96f008b8f403dd80d16a78e480148fbc1b6c7832cbc267f8023598 |
| results/parallel21_positive_upper_review_2026_09_05.json | 6047ae01dc4df8710f5b716032548eb7882e4e6c0bf4da690ec2fd9e1e957fc0 |
| research/parallel21-squarefree-moments-2026-09-05.md | 42bfee31ccd53abaa9081128a3dfa250c5316570c407740c439b7175759f0210 |
| research/parallel20-inversion-moments-2026-09-05.md | 31ec520dcb38decced35eb111128d65383f4a4f6fb2eca545d8bc5406b55922f |

The Weil statement was freshly checked in the primary HTML of
[McDonald--Sahay--Wyman, Lemma 2.1](https://arxiv.org/html/2210.03789v2#S2).
The applications here use quadratic characters and polynomials with three
or four distinct simple roots. Such polynomials cannot be a constant
multiple of a square; the usual hypotheses are satisfied. No theorem
about VC dimension is used. No prior artifact or central assessment was
edited by this review.

## 2. Full translate transform and overlap coefficients

Retain the reviewed notation e_j, T_j, U_C(t), L_all, L_in and L_out.
The full character kernel identity is

\[
 \sum_t\chi(x-t)\chi(y-t)=p1_{x=y}-1.
\]

For x=y there are p-1 nonzero terms. For x!=y, deleting t=y and
substituting u=(x-t)/(y-t) gives every u except 1 exactly once; the
character sum is therefore -1. There is no assumption that chi(-1)=1.
Consequently

\[
 L_{all}=p\sum_x e_3(x)^2-T_3^2.
\]

The orientation U_C(t)=sum_x chi(x-t)e_3(x) is consistent with this
identity in both odd prime congruence classes. The factor p, the negative
rank-one term, and the absence of an extra chi(-1) all check out.

For N signs taking values +/-1, each ordered pair of three-element
subsets contributes according to its intersection size:

| Intersection size | Surviving elementary coefficient | Multiplicity |
|---:|---|---:|
| 0 | e_6 | 20 |
| 1 | e_4 | 6(N-4) |
| 2 | e_2 | (N-2)(N-3) |
| 3 | 1 | binom(N,3) |

In the third row, the two surviving indices can be ordered in two ways,
and the common pair can be chosen in binom(N-2,2) ways. Thus the
coefficient is 2 binom(N-2,2), not binom(N-2,2). In the second row, the
common index is outside the four surviving indices and the four indices
split into an ordered pair of two-element sets in six ways. This verifies
all coefficients in the reviewed equation (7).

## 3. Both exceptional-row corrections are correct

At x=c in C there is exactly one zero character value. Averaging over
replacing it by +/-1 adds e_2(c)^2 to the square of e_3(c), while the
averages of e_2, e_4 and e_6 do not change. The required correction to
the n-sign ledger is therefore **minus** e_2(c)^2.

An independent algebraic check avoids the random-sign argument. The
difference between the n-sign and (n-1)-sign versions of the right side
of the cubic overlap identity is

\[
 6e_4+2(n-3)e_2+\binom{n-1}{2}=e_2^2,
\]

where the last equality is the quadratic overlap identity on the n-1
nonzero signs. This confirms the same sign and magnitude. After summing
rows and using T_2=-binom(n,2), the exact remainder is

\[
 R=6(n-4)T_4-(n-2)(n-3)\binom n2+p\binom n3
       -\sum_{c\in C}e_2(c)^2.
\]

For an omitted parameter t in C and a triple Q={t,a,b}, the repeated
root contributes

\[
 \sum_x\chi(x-t)^2\chi(x-a)\chi(x-b)
   =-1-\chi((t-a)(t-b)).
\]

The second term removes the row x=t from the complete pair correlation.
Summing those triples gives exactly the correction
-binomial(n-1,2)-e_2(t) in equation (12). The on-set identity (14) follows
because each four-element subset has four choices of distinguished t.
No row or repeated-root term has been silently discarded.

## 4. Uniform constants on n>=6, n^4<=p

Every inequality needed for the proposition checks algebraically:

| Quantity | Verified upper bound | Size condition used |
|---|---|---|
| abs(6(n-4)T_4) | (3/4)pn^3 | abs(K_4)<=3 sqrt(p), n^2<=sqrt(p) |
| (n-2)(n-3)binom(n,2)+sum_c e_2(c)^2 | (1/4)pn^3 | n^2+2n<=p |
| p binom(n,3) | (1/6)pn^3 | binomial bound |
| T_3^2 | (1/9)p^2 n^3 | abs(K_3)<=2 sqrt(p), n^3<=p |
| abs(U_C(t)), t in C | sqrt(p)n^3 | 2<=n sqrt(p) |
| L_in | p^2 n^3 | n^4<=p |

For the second row, the preceding bound is n^4/2+n^5/4. It is at most
pn^3/4 precisely when n^2+2n<=p, which follows from n>=6 and n^4<=p.
The first row uses
18(n-4)binom(n,4)<=3n^5/4. The omitted-parameter bound follows from

\[
 |U_C(t)|\le3\sqrt p\binom{n-1}{3}
                 +2\binom{n-1}{2}
 \le\tfrac12\sqrt p\,n^3+n^2\le\sqrt p\,n^3.
\]

The signs in R give

\[
 -pn^3\le R\le(11/12)pn^3\le pn^3.
\]

These estimates hold throughout the claimed range, including n=6; no
unmentioned eventual-size assumption is needed for this proposition.
Combining the exact identities gives

\[
 L_{out}+L_{in}=20pT_6+pR-T_3^2.
\]

For the first direction, subtracting the nonnegative L_in and T_3^2 and
using R<=pn^3 gives

\[
 L_{out}\le20pT_6+p^2n^3.
\]

For the other direction, solving for 20pT_6 gives

\[
 20pT_6=L_{out}+L_{in}-pR+T_3^2
 \le L_{out}+(1+1+1/9)p^2n^3
 \le L_{out}+3p^2n^3.
\]

Thus the sharper 19/9 intermediate constant and the rounded constant 3
are both valid. The lower bound T_6>=-pn^3/20 also follows correctly
from nonnegativity of L_all and T_3^2. It provides no positive upper
estimate.

## 5. Equivalence directions and the exact remaining difficulty

For any family of sets in this size range, the two implications are

\[
 T_6\le A pn^3\quad\Longrightarrow\quad
 L_{out}\le(20A+1)p^2n^3,
\]

\[
 L_{out}\le D p^2n^3\quad\Longrightarrow\quad
 T_6\le(D+3)pn^3/20.
\]

These are equivalences of the existence of uniform constants, with the
displayed changes of constants. They are not an identity of suprema and
do not assert that either constant is already available. They apply to
Sidon sets simply by restricting the family; Sidon is not used in the
proof.

The genus-one description of each off-C summand is valid: its four roots
are distinct and the characteristic is odd. But pairwise correlations
of these summands return the original six-root polynomial. Indeed, for
triples Q,R and I=Q intersection R, D=Q symmetric-difference R,

\[
 f_Q(x)f_R(x)=f_D(x)1_{x\notin I},
\]

so the full transform gives

\[
 \sum_t U_Q(t)U_R(t)
 =p\left[K(D)-\sum_{a\in I}f_D(a)\right]-K(Q)K(R).
\]

For Q=R the bracket is p-3. For disjoint triples it is exactly K(Q union
R), a squarefree degree-six character sum. Thus an unproved assertion of
near orthogonality for the elliptic family would be assuming the missing
signed cancellation in a different form. Individual quartic Weil bounds
do not justify it. Termwise degree-six Weil bounds give scale
sqrt(p)n^6 for T_6, a factor of order p^(1/4) too large at
n=floor(p^(1/4)).

The reviewed document explicitly acknowledges this return of the
six-root problem. **No circular inference appears in its actual proof.**
Calling the interface a proved upper bound or a reduction in the
arithmetic difficulty would, however, exceed its content.

There are two further scope boundaries:

- The proposition itself starts at n=6. The cited pass-20 Sidon-to-all-set
  transfer has its own finite hypotheses. Its convenient threshold for
  h=2 is K_2=25; on the sixth-moment slice it applies eventually when
  n>=25. The five smaller numerical fixtures do not certify that transfer.
- A single sixth-moment estimate would not prove the full classical
  conjecture. The pass-21 sufficient criterion for the full conjecture
  requires an unbounded sequence of fixed moment orders. Neither the
  subgroup square-root target nor the official-prize bridge follows from
  this interface.

These qualifications are compatible with the reviewed note's stated
limitations; they are not defects requiring a correction there.

## 6. Verifier review and evidence limits

The code constructs the actual prime-field quadratic character, computes
e_j by a descending elementary-symmetric recurrence, and independently
sums distinct-subset character products for degrees 3, 4 and 6. Its full
translate transform uses chi(x-t), matching the mathematical orientation.
It checks every overlap type of the cubic-feature Gram identity in the
specified small examples, plus the complete-energy and exceptional-row
identities. The large cases are actual Sidon sets on the stated size
slice.

The recorded result contains 216 cases, 216 overlap ledgers, 1,022
repeated-root row checks and five slice-constant checks. Its three pinned
input hashes match the present verifier and the two cited transfer notes.
The Python source parsed successfully. These are checks of artifact
consistency and code logic; this review does not relabel the historical
run as a new independent execution.

The integer arithmetic is safe for the recorded fixtures. They have
p<=10007 and n<=10, so the initial quadratic-residue table and pointwise
cubic polynomial are well within int64. The symmetric coefficient sums
are bounded by p binom(n,d), and the code checks that bound before their
reductions. The matrix transform bounds the absolute sum in each row;
the dot-product helper bounds the sum of absolute products before each
integer dot product. The remaining moment sum uses Python integers.
No floating-point tolerance decides an identity or inequality.

No mathematical or implementation finding required a verifier rerun or
an edit to the reviewed artifacts. The output remains a valid exact
reformulation with an unproved uniform upper bound.
