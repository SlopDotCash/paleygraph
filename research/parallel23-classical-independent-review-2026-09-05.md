# Separate-agent review of the almost-all Sidon sixth-moment upper bound

**Verdict: no mathematical correction found in the stated restricted
theorem.** The proof establishes the finite conditional probability
bound, the asymptotic constant \(216/5\), and the quartic-star payoff
\(2/3\) for the sets meeting the stated sixth-moment threshold. It does
not establish the uniform input required by the Paley transfer.
The initially identified reporting issue is corrected in the current
verifier and results pinned below.

This is a separate agent's derivation and source-code review, not human
refereeing or formal certification. I read the proof, verifier, and
recorded results; independently recalculated the main identities and
constants; checked Python syntax and all recorded input hashes; and
checked the relevant prior quartic-star definitions and transfer
threshold. I did not rerun the long finite checks. No author, subgroup,
or central artifact was changed by this review. After root reran the
corrected verifier, I inspected the fresh results and verified that the
change affects only reporting labels and their associated script hash.

## Reviewed bytes

| Input | SHA-256 |
| --- | --- |
| [Proof](parallel23-classical-upper-2026-09-05.md) | `08a79fd9109396f1c31d0efdd70eba3bcb77968f9b9b02daf40ab2b6b674cc90` |
| [Verifier](../experiments/parallel23_classical_upper_2026_09_05.py) | `6380ad545ff05d61ae973dfdbfa72bf8d636a27e3038f47a9d0593245f81a89f` |
| [Recorded results](../results/parallel23_classical_upper_2026_09_05.json) | `23786c66eb4bfdafda63bbc8b1195b8aca39a89ad8a01c4eb97088e7e462c369` |
| [Quartic-star interface](parallel21-positive-upper-review-2026-09-05.md) | `832805eb897c8a8235d7d04909adcb4c99c263a8d5b86f7261266ae5c198e9ea` |
| [Pass-20 transfer](parallel20-inversion-moments-2026-09-05.md) | `31ec520dcb38decced35eb111128d65383f4a4f6fb2eca545d8bc5406b55922f` |
| [Earlier fourth-moment result](parallel6-classical-2026-09-04.md) | `b9c1cea861c8df6b99c0966f376c0f709bc974819650891ec903d0a43f28a9c6` |

The recorded hash of the prior quartic-star independent review also
matches the current file:
`970e079b585d649b2233d5bf20b0cb836e773a1b39b75ae7e4234eef5323f5fb`.

## 1. Balanced sign sampling really includes the zero

Fix a row \(x\). Its population over the possible set elements is
exactly one zero and \((p-1)/2\) copies of each sign. A uniform
\(n\)-subset therefore samples the stated population without
replacement.

After choosing the position of zero, let \(m=(p-1)/2\). Every fixed
balanced assignment on the remaining positions is compatible with
exactly \(m!\) perfect matchings: each positive position is matched
bijectively to a negative position. Each compatible matching has
exactly one choice of its independent pair orientations that produces
that assignment. Consequently the random matching construction gives
the uniform balanced assignment, with no weighting depending on the
chosen row or sample.

Conditional on the matching and zero position, selected pairs with
both endpoints cancel. Every selected singleton endpoint contributes
an independent fair sign. The resulting sum uses \(\ell\le n\)
independent signs. Jensen's inequality applied when adding each further
mean-zero sign shows \(U_s(\ell)\le U_s(n)\). Thus the even-moment
domination holds with the exceptional zero present; no independent-row
assumption is made.

The sixth moment is exactly

\[
 U_3(n)=15n^3-30n^2+16n.
\]

The general pairing bound \(U_s(j)\le(2s-1)!!j^s\) has the correct
direction: each surviving even-multiplicity word is compatible with
at least one pairing. Words with \(s\) distinct indices, each used
twice, give leading coefficient \((2s-1)!!\); every other surviving
word contributes lower polynomial degree. Both the finite mean bound
and the later asymptotic expansion are justified.

## 2. The slice Poincare induction has the correct constants

Use the author's convention that \(D_n(f)\) is half the expected
squared difference across an oriented Johnson-graph edge. Put
\(q=p-n\). Given a uniform \((n-1)\)-set \(A\), its complement has
\(q+1\) elements. Two independent insertions are distinct with
probability \(q/(q+1)\). Conditional on being distinct, they produce
a uniform oriented edge of the \(n\)-set slice. This proves

\[
 \mathbb E_A\operatorname{Var}_{b\notin A}f(A\cup\{b\})
 =\frac q{q+1}D_n(f).                                \tag{1}
\]

For adjacent \((n-1)\)-sets write
\(A=S\cup\{a\}\), \(A'=S\cup\{b\}\), where \(|S|=n-2\).
The term \(f(S\cup\{a,b\})\) cancels in their insertion averages.
There are exactly \(q\) remaining paired terms, indexed by the elements
outside \(S\cup\{a,b\}\), with denominator \(q+1\). Cauchy–Schwarz
gives

\[
 D_{n-1}(g)\le\frac{q^2}{(q+1)^2}D_n(f).
\]

The induced \(n\)-set edge is uniform: each edge with common
\((n-1)\)-set arises once for each choice of its distinguished common
element. This count is independent of the edge. The induction base is
the exact identity \(\operatorname{Var}_1f=(p-1)D_1(f)/p\).

The variance decomposition applies because sampling \(A\) and then
one insertion makes the resulting \(n\)-set uniform. The induction
coefficient simplifies as

\[
 \frac q{q+1}+rac{(n-1)(q+1)}p\frac{q^2}{(q+1)^2}
 =\frac{nq}p,
\]

using \(p=n+q\). Combining this with (1) gives precisely
\(n(p-n+1)/p\) times the expected insertion variance. The argument
works for all \(1\le n<p\), including the induction endpoints. It
does not substitute an independent-coordinate variance inequality.

## 3. Character orientation and insertion corrections

With \((Qh)(b)=\sum_x\chi(x-b)h(x)\), the relevant matrix identity is

\[
 Q^TQ=pI-J,
\]

not an assertion about \(Q^2\). For \(x\ne y\), the substitution
\(u=(x-b)/(y-b)\), with the zero term \(b=y\) omitted, gives the
off-diagonal correlation \(-1\); the diagonal is \(p-1\). Thus

\[
 \sum_b(Qh)(b)^2=p\sum_xh(x)^2-\Bigl(\sum_xh(x)\Bigr)^2
\]

for both prime congruence classes. In particular, the fact that the
kernel is skew-symmetric when \(p\equiv3\pmod4\) introduces no sign
change into this norm identity. No quadratic Gauss-sum sign is needed.

Writing \(z=\chi(x-b)\), the even powers satisfy
\(z^{2k}=1-\mathbf1_{x=b}\) and the odd powers equal \(z\).
The binomial expansion of \((F+z)^6\) therefore gives exactly

\[
 M_6(A\cup\{b\})=M_6(A)+15M_4(A)+15M_2(A)+p-1
 +Q(6F^5+20F^3+6F)(b)-15F(b)^4-15F(b)^2.
\]

Both negative exceptional-row terms and the constant \(p-1\) are
required and correctly retained. Passing from variance on the
\(p-n+1\) admissible insertions to the full sum of squared insertion
terms is an upper bound. Its denominator cancels the matching factor
in the slice inequality, leaving exactly \(n/p\).

The polynomial squares have coefficients
\((36,240,472,240,36)\) and \((1,2,1)\) as stated. Moment domination
and the full transform give
\(\mathbb E\|QH\|_2^2\le p^2A\) and
\(\mathbb E\|Z\|_2^2\le pB\). The triangle inequality in the joint
space of subsets and field coordinates then yields

\[
 \operatorname{Var}M_6(C)
 \le n(\sqrt{pA}+15\sqrt B)^2=V(p,n).
\]

The argument uses no unavailable quartic-family cancellation estimate
and no assumed uniform sixth-moment bound. There is no circular upper
bound concealed in this variance calculation.

## 4. The exact threshold, constants, and Sidon conditioning

The threshold is the finite value \(15pn^3\). The proved mean bound
places it at distance at least

\[
 G(p,n)=p(30n^2-16n)>0
\]

above the actual mean. Consequently Chebyshev gives
\(\mathbb P(M_6>15pn^3)\le V/G^2\) without replacing the finite
threshold or its gap by an asymptotic expression.

Since \(j=n-1\), the leading terms are
\(A=34020n^5(1+O(1/n))\) and
\(B=105n^4(1+O(1/n))\). The ratio
\(\sqrt B/\sqrt{pA}=O((pn)^{-1/2})\) tends uniformly to zero
under \(n^4\le p\). Therefore

\[
 \frac V{G^2}=
 \left(\frac{34020}{900}+o(1)\right)\frac{n^2}p
 =\left(\frac{189}5+o(1)\right)\frac{n^2}p.
\]

There are \(p(p-1)/2\) three-label relations \(2a=b+c\), with the
pair \(b,c\) unordered. There are \(p(p-1)(p-3)/8\) relations
between two disjoint unordered pairs. The factor eight accounts for
the two orientations inside each pair and exchange of the two pairs.
After canceling the inclusion-probability denominators, their union
bound is exactly

\[
 \frac{(n)_3}{2(p-2)}+\frac{(n)_4}{8(p-2)}
 =\frac{n(n-1)(n-2)(n+1)}{8(p-2)}.
\]

The final \(1/8\) bound is valid: the numerator is
\(n^4-2n^3-n^2+2n\le n^4-2\le p-2\) for \(n\ge3\).
This proves a positive Sidon probability, at least \(7/8\), throughout
the theorem's finite range. Dividing by that lower bound gives

\[
 \mathbb P(M_6>15pn^3\mid C\text{ Sidon})
 \le\min\{1,(8/7)V/G^2\}.
\]

The constant is \((8/7)(189/5)=216/5\). On the exact slice
\(n=\lfloor p^{1/4}\rfloor\),
\(n^2/p=p^{-1/2}(1+o(1))\). The proof is uniform in the prime and
does not average over primes or require a congruence restriction.

## 5. The quartic-star payoff is correctly normalized

For the actual character row there are \(N(x)=n-\mathbf1_{x\in C}\)
nonzero signs. Newton's identity gives

\[
 6e_3(x)=F_C(x)^3-(3N(x)-2)F_C(x).
\]

Since \(n\ge2\), \(3N(x)-2\ge0\). Squaring, dropping its
nonpositive fourth-power term, and bounding
\((3N(x)-2)^2\le9n^2\) gives

\[
 \sum_xe_3(x)^2\le\frac{M_6(C)}{36}+rac{n^2M_2(C)}4.
\]

The exact second moment is \(M_2(C)=pn-n^2\). Thus a set satisfying
the finite sixth-moment threshold obeys the slightly stronger bound

\[
 \sum_xe_3(x)^2\le\frac23pn^3-\frac14n^4
 \le\frac23pn^3.
\]

The quartic-star transform in the prior interface is exactly
\(U_C=Qe_3\), with the same orientation. Its full norm is
\(p\sum_xe_3(x)^2-T_3(C)^2\). Removing the nonnegative on-set
energy can only lower it, proving
\(L_{\mathrm{out}}(C)\le(2/3)p^2n^3\). This last implication
actually needs no Sidon hypothesis once the moment threshold holds.
All zero and on-set rows remain accounted for.

## 6. Verifier scope and remaining limitation

The source implements the proved Rademacher recursion, both insertion
corrections, the correctly oriented full transform, and the slice
factor. Its radical certificate uses
\(n[pA+225B+30\lceil\sqrt{pAB}\rceil]\), which is an upper enclosure
of \(V\). Its comparisons with the unrounded radical first separate
the nonpositive case before squaring, preserving the inequality's
direction.

The recorded exact checks cover 3,551 exhaustive small slice sets,
15,042 insertion identities, 1,780 transform identities, 38 slice
variance checks, and 19 sixth-variance checks. The larger actual
fields contain 287 sampled Sidon sets, with five full quartic-star
energy fixtures. The five still larger prime certificates evaluate
the analytic probability formula only. These are distinct scopes.
For the five NumPy fixtures, \(n\le10\) and \(p\le10007\), so the
row powers are safely within int64; the checked bound
\(p\max|e_3|<2^{63}\) also controls every partial dot-product sum.
Subsequent energies use Python integers.

The initial reporting issue is corrected. The helper now sets
`sidon_certificate_applicable` to the exact predicate
`n >= 6 and n**4 <= p`, and emits `sidon_bad_fraction_upper: null`
outside that range. I checked all 19 exhaustive small-field slices:
every one has applicability false and a null conditional certificate.
The five actual in-range slices have applicability true and retain
their previous mathematical certificates; the five large analytic
certificates are also unchanged.

For a strict comparison without another field computation, I reversed
only the reporting guard in memory and recovered the previously
reviewed script hash
`f1c3ebaa086d728dab88e76ad16642cda8d40b738e60bd2e614400a56e8b3916`.
Reversing only the applicability/null-label changes and the corresponding
script-hash entry in the new result recovered the previously reviewed
result hash
`9c1710a9ea3ba6dc109f3438d3562d9385fb898d108251569d6bbcc30a702f72`.
Thus every other result field, including all verification counts and
all in-range coefficients, is unchanged. The proof note's hash is
unchanged as well. These comparisons wrote no author files and did
not rerun the verifier.

The pass-20 convenient threshold specializes to
\(K_2=\max\{2,8,25\}=25\), as the proof states. That transfer needs
the input for every Sidon set at the prescribed size. An exceptional
proportion \(O(p^{-1/2})\) does not supply this universal quantifier,
and the sampled successes do not remove the gap. No uniform bound on
the exceptional family or on the subgroup target follows from this
review.
