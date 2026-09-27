# Pass 38: independent audit of the incidence step

Status: a specific missing hypothesis is verified. This note does not prove a counterexample to the unrestricted numerical inequality, does not import the old 22/9 bound, and does not establish the Paley conjecture.

## Sources inspected

- Shkredov, *Some new inequalities in additive combinatorics*, arXiv:1208.2344v3, Lemma 7 on PDF page 6 and Theorem 34 on PDF pages 21-23. The actual pages were rendered and inspected; source: `sources/parallel37-energy-source/shkredov-1208.2344v3.pdf`.
- Its published version, recovered independently by the source-audit lane: `sources/parallel38-energy-source/shkredov-2013-published.pdf`. PDF page 10, printed page 197, Lemma 2, and PDF page 29, printed page 216, Theorem 8, were visually inspected in this lane.
- Shkredov-Vyugin, *On additive shifts of multiplicative subgroups*, [arXiv:1102.1172v1](https://arxiv.org/pdf/1102.1172v1), Corollary 5.1(3), formula (49), on page 15, and its stated antecedent Lemma 4.1 on pages 11-13. This earlier corollary states the unrestricted three-set estimate; it is presented as a simple consequence, without a separate proof. The associated note points to the union of normalized Cartesian products. The PDF itself is stamped `arXiv:1102.1172v1 [math.NT] 6 Feb 2011`; an exact-version download was byte-identical to the initially retrieved PDF. Raw bytes are archived at `sources/parallel38-energy-source/shkredov-vyugin-1102.1172v1.pdf` (194491 bytes, SHA256 `8a8368e912b16cd7ef8693014879c03dff2b5160e2af53c61fcd5aced9329080`). These page and equation numbers refer specifically to v1, not unversioned HTML.

## What changed in the published incidence lemma

Write `t = |H|`. Besides the two size assumptions

\[
|Q||Q_1||Q_2|\ll t^5,
\qquad |Q||Q_1||Q_2|t\ll p^3,
\]

published Lemma 2 requires

\[
\left|(Q_1\times Q_2)\cdot\Delta^{-1}(Q)\right|
=\frac{|Q_1||Q_2||Q|}{t}. \tag{J}
\]

Here the set on the left consists of the pairs `(q1/q, q2/q)` with `qi` in `Qi` and `q` in `Q`. This equality is absent from arXiv v3 Lemma 7. The published bound has the same right-hand side as the old incidence estimate, but its scope is narrower.

## Exact multiplicity calculation

Let `G = F_p^*/H` be the quotient group. Let `A`, `B`, and `C` be the coset sets underlying `Q1`, `Q2`, and `Q`, respectively, so `|Q1|=t|A|`, and similarly for the other sets. For `u,v` in `G`, define

\[
m(u,v)=|\{c\in C:cu\in A,\ cv\in B\}|,
\quad U=\{(u,v):m(u,v)>0\}.
\]

The normalized field-pair image in (J) is the disjoint union of the `H x H` cosets indexed by `U`. Consequently

\[
\left|(Q_1\times Q_2)\cdot\Delta^{-1}(Q)\right|=t^2|U|,
\qquad
\frac{|Q_1||Q_2||Q|}{t}=t^2|A||B||C|
=t^2\sum_{u,v}m(u,v).
\]

Thus (J) holds **if and only if** `m(u,v)=1` for every pair in `U`. Equivalently, the map

\[
A\times B\times C\longrightarrow G\times G,
\qquad(a,b,c)\longmapsto(a/c,b/c)
\]

is injective.

This multiplicity affects the exact incidence count. Put

\[
\rho(u,v)=|\{(x,y)\in uH\times vH:y-x=1\}|.
\]

For any fixed coset triple `(a,b,c)`, normalizing an equation `y-x=z` by `z` gives `t rho(a/c,b/c)` solutions: there are exactly `t` choices of `z` in its coset. Therefore

\[
\sum_{z\in Q}(Q_1\circ Q_2)(z)
=t\sum_{u,v}m(u,v)\rho(u,v). \tag{M}
\]

Replacing a list of normalized pairs by its support `U` loses the factor `m`. Distinct-coset-pair estimates alone therefore do not justify that replacement. This is an identity and an identified missing inference, not a claimed counterexample to an asymptotic bound.

For example, if `A=B=C=D` is a subgroup of `G` of size `M`, then `U=D x D` and `m(u,v)=M` throughout `U`. More generally, for **any** common coset set `A=B=C=D`, one has `m(1,1)=|D|`. In particular, when `Q1=Q2=Q=S` is a union of two or more `H` cosets, condition (J) necessarily fails.

The root script `experiments/parallel38_coset_multiplicity.py` includes a nested-subgroup example at `p=97`, `H=<5^12>` of order `8`, and `S=<5^4>` of order `24`. Its calculation was independently checked: the direct incidence count is `48`, the weighted count is `48`, and discarding the quotient multiplicities gives `16`. There are nine normalized coset pairs, all of multiplicity three. Both `H` and `S` are invariant under additive negation, and the direct counts for `x+y in S` and `y-x in S` both equal `48`; the script's addition convention therefore agrees with the correlation convention used here.

A useful conditional repair keeps the weights. Suppose a distinct-pair estimate

\[
\sum_{(u,v)\in V}\rho(u,v)\le C t^{2/3}|V|^{2/3}
\]

is valid for every level set `Uj={(u,v):m(u,v)>=j}` under consideration. Writing `M=max m` and `L=sum m`, layer decomposition and Holder give

\[
\sum m\rho
=\sum_{j=1}^{M}\sum_{U_j}\rho
\le C t^{2/3}\sum_{j=1}^{M}|U_j|^{2/3}
\le C t^{2/3}M^{1/3}L^{2/3}.
\]

Together with (M), this yields

\[
\sum_{z\in Q}(Q_1\circ Q_2)(z)
\le C t^{-1/3}M^{1/3}(|Q_1||Q_2||Q|)^{2/3}.
\]

This implication is exact. It retains a potential `M^(1/3)` loss, with `M <= min(|Q|,|Q1|,|Q2|)/t`. It is conditional on the stated distinct-pair bound and all its hypotheses for the level sets; this note does not independently import a stronger unrestricted theorem.

## Consequence for the old Theorem 34 proof

The algebra preceding the three-set incidence application can be checked independently. Let

\[
\psi(x)=(H\circ H)(x),\quad E=\sum_x\psi(x)^2,
\quad F_3=\sum_x\psi(x)^3.
\]

The matrix `M_{h,h'}=psi(h-h')` on `H` is positive semidefinite with constant row sum `E/t`. Hence `tr(M^3) >= E^3/t^3`, giving the old formula (75). The contribution with a zero edge is at most `3t F3`. For any fixed choice of an edge, removing terms where that edge has weight below `d` loses at most `d E^2/t`, since each row sum is `E/t`. Taking `d` a sufficiently small absolute multiple of `E/t^2` and using Cauchy-Schwarz, with `sum C3(H)^2=F3`, yields the old lower bound

\[
\sum_{\alpha,\beta\ \mathrm{retained}}
\psi(\alpha)^2\psi(\beta)^2\psi(\alpha-\beta)^2
\gg \frac{E^6}{t^6 F_3},
\]

unless the removed zero-edge contribution already proves a stronger energy estimate.

The dyadic sets `Si` have size `|Si| << t^3/(2^(3i)d^3)` by the valid special case with the other two sets equal to `H`. That special case satisfies (J): its map is `c -> (1/c,1/c)` on the quotient, which is injective. Under the contradictory assumption `E >> t^(22/9) log^(2/3)t`, the two old size hypotheses for every triple of dyadic sets also hold in the quartic regime. In the notation `E=t^3/K`, they follow from `|Si| << Kt`, `K << t^(5/9) log^(-2/3)t`, and

\[
|S_i||S_j||S_k|\ll t^{14/3}\log^{-2}t,
\quad t|S_i||S_j||S_k|\ll t^{17/3}\log^{-2}t\ll p^3.
\]

The failure occurs at the next step: the old proof applies the three-set bound to **all** triples `(Si,Sj,Sk)` in formula (80), checking only those two size conditions. Its hypotheses do not establish (J). In fact, for each diagonal triple `(Si,Si,Si)` with `Si` containing more than one `H` coset, (J) fails by the calculation above. The sets are symmetric under additive negation, so changing the displayed convolution to the correlation in the lemma does not remove this obstruction.

The root lane supplied a concrete dyadic example, independently reproduced here by exact finite-field enumeration: take `p=1153`, `H=<75>` of order `8`. Then `E=168`, `d=E/(16t^2)=21/128`, and the old level `S4={x != 0:8d < psi(x) <= 16d}` has `24` elements, comprising three `H` cosets. Its normalized field-pair image has size `1600`, whereas (J) requires `24^3/8=1728`. Of its field-pair fibers, `1536` have size `8` and `64` have size `24`. This explicitly disproves automatic satisfaction of the missing hypothesis for the old proof's own level sets; it is not a numerical counterexample to the energy estimate.

If the unrestricted incidence bound were available, the remaining calculation would correctly give `E^6 << t^(44/3) log^4 t`, hence the advertised `E << t^(22/9) log^(2/3)t`. The identified defect is the incidence hypothesis, not a miscalculated final exponent or a quartic size restriction.

## Published replacement and project implication

The visually inspected published Theorem 8 contains the energy bound

\[
E(H)\ll t^{32/13}\log^{41/65}t
\qquad(t\ll p^{2/3}),
\]

as one arm of its displayed minimum. Since `32/13 > 49/20 > 22/9`, this published replacement removes the historical conflict with the later `49/20` result. The older arXiv statement must not be used to improve the project's accepted exponent. This lane has not audited the full published proof or changed any central bounds.
