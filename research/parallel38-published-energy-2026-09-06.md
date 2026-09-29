# Published subgroup energy statement: source comparison

Date: 2026-09-06. Bounded primary-source audit; no proof or build claims.

## Result

The publisher's version of I. D. Shkredov, *Some new inequalities in additive combinatorics*, gives exponent **32/13**, not **22/9**, for the small-subgroup additive-energy term. Its invariant-incidence lemma also contains an additional hypothesis absent from arXiv:1208.2344v3. These are directly observed version differences, not a conclusion inferred from later publication dates.

The old publisher site works over HTTP even though its HTTPS endpoint fails. The retrieval chain is:

- [Publisher volume 3, issue 3-4 page](http://mjcnt.phystech.edu/en/published_volumes.php?i=3&v=3).
- [Publisher article 66](http://mjcnt.phystech.edu/en/article.php?id=66).
- [Publisher full-text PDF](http://mjcnt.phystech.edu/en/download.php?id=66).

Exact downloaded bytes, the two decisive rendered pages, the publisher HTML pages, and SHA-256 hashes are archived in `sources/parallel38-energy-source/manifest.json`. The publisher PDF SHA-256 is `79dbfd8841ae220880ad8c692684bcd0889e417b700a96e2ed64f3ac7983c206` (8,617,044 bytes; 52 PDF pages including a cover).

## Bibliographic resolution

The published title page, PDF page 2, identifies **Moscow Journal of Combinatorics and Number Theory, 2013, vol. 3, issues 3-4, pp. 189-239 [pp. 425-475]**. The publisher issue page independently identifies volume 3, issue 3-4 (2013), with a printed-date field `10.05.2014`. The title page's received/revised date fields are internally out of chronological order; they should not be used to infer the revision history.

The live [author homepage](https://www.math.purdue.edu/~ishkredo/) instead lists 3:2 (2013), 237-288. That metadata does not match the publisher's actual reprint. For this audit, the reprint's title-page data controls.

## Exact published energy statement

On printed page **216 [452]**, PDF page **29**, section 7, **Theorem 8, equation (75)** assumes that \(p\) is prime, \(\Gamma\subseteq\mathbb F_p^*\) is a multiplicative subgroup, and \(|\Gamma|\ll p^{2/3}\). Setting \(n=|\Gamma|\), its conclusion is

\[
 E(\Gamma)\ll
 \min\left\{
 n^{32/13}\log^{41/65}n,
 n^3p^{-1/3}\log n+p^{1/26}n^{31/13}\log^{8/13}n
 \right\}.
\]

This statement was read directly from the complete page image, archived as `shkredov-2013-published-page29.png`. The proof continues through printed page 220 [456], PDF page 33. Printed page 219 independently repeats the term \(E\ll t^{32/13}\log^{41/65}t\). The audited energy normalization, printed page 193, is the ordinary number of ordered additive quadruples; equation (9) also gives \(E(A)=\sum_x(A\circ A)(x)^2\). Logarithms are base 2, as specified on printed page 194.

For comparison, [arXiv:1208.2344v3](https://arxiv.org/abs/1208.2344v3), dated 2012-11-06, has **Theorem 34, equations (72)-(74)**, on its page 21. It states a \(22/9\) exponent, including

\[
E(\Gamma)\ll n^{22/9}\log^{2/3}n
\quad\text{when}\quad n\ll p^{3/5}\log^{-6/5}p.
\]

The current arXiv submission-history page still ends at v3. Therefore its latest arXiv PDF is not interchangeable with this published theorem statement.

Arithmetic comparison: \(32/13-49/20=3/260>0\), whereas \(49/20-22/9=1/180>0\). Thus the published 2013 exponent is consistent with a later improvement to \(49/20\); the apparent exponent conflict arose from using the different arXiv statement. This is a numerical comparison, not independent verification of either proof.

## Exact additional incidence hypothesis

The preprint's **Lemma 7, equation (18)**, page 6, has its published counterpart at printed page **197 [433]**, PDF page **10**, **Lemma 2, equation (18)**. The complete page is archived as `shkredov-2013-published-page10.png`.

For \(\Gamma\)-invariant \(Q,Q_1,Q_2\subseteq\mathbb F_p^*\), the published lemma lists **three** conditions:

\[
 |Q||Q_1||Q_2|\ll|\Gamma|^5,
 \qquad |Q||Q_1||Q_2||\Gamma|\ll p^3,
\]
\[
 \left|(Q_1\times Q_2)\cdot\Delta^{-1}(Q)\right|
 =|Q_1||Q_2||Q||\Gamma|^{-1}.
\]

Its conclusion is

\[
 \sum_{x\in Q}(Q_1\circ Q_2)(x)
 \ll |\Gamma|^{-1/3}(|Q||Q_1||Q_2|)^{2/3}.
\]

The displayed cardinality equality is an explicit extra assumption in the publisher's text. The arXiv v3 Lemma 7 lists only the preceding two size conditions. The published paper defines \(\Delta_k(A)=\{(a,\ldots,a):a\in A\}\) on printed page 194. In this multiplicative expression, the set on the left consists of pairs \((q_1/q,q_2/q)\).

The lemma's lead-in cites Stepanov's method [25] and reference [24]. Printed page 238 identifies [24] as I. D. Shkredov and I. V. V'yugin, *On additive shifts of multiplicative subgroups*, Mat. Sbornik 203:6 (2012), 81-100. This audit did not duplicate the separate antecedent-paper audit.

## Limits and consequence for proof reuse

No explicit erratum or author statement explaining the change was found in this bounded lookup. The supported conclusion is narrower and concrete: the published energy theorem differs, and the published incidence lemma adds a condition. This audit does not assign a cause to the arXiv proof discrepancy.

Any argument that invokes the preprint's general invariant-incidence lemma must now establish the publisher's additional cardinality condition or supply another valid incidence result. The published source by itself does not justify treating the unrestricted arXiv lemma or its \(22/9\) energy consequence as a published theorem.

No repository build, external message, submission, or proof-status promotion was performed. Files written by this lane are confined to this note and `sources/parallel38-energy-source/`.
