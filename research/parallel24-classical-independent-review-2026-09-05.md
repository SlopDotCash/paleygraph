# Independent root review of the exceptional-set lane

**Verdict: the stated local identities and quantified losses check out.
They do not supply an upper bound for every exceptional Sidon set.**
The reviewer is the root agent, distinct from the lane author. This is
mathematical and computational agent review, not human refereeing or
formal verification.

I read the complete final [note](parallel24-classical-exceptions-2026-09-05.md),
its verifier, and results. I independently expanded the six moments of
the row increment Δ=v−u from the actual inside/outside character
populations. Their sums are F and −F and their nonzero counts are
N and p−1−N. In particular, the even increment numerators are
S+2F², S+8F²+6N(p−1−N), and S+32F²+30N(p−1−N).
The odd ones are −pF, −(4p−3)F, and −(16p−15)F. Substitution into
the binomial expansion independently reproduces all four coefficients
of the sixth-moment drift, including the row-zero correction.

The two lower-degree positivity arguments are valid for n≥6,p≥n⁴.
The iteration must distinguish n=6: the leading coefficient there is
negative, and the note now explicitly restricts useful nonnegative
iteration to n≥7. This is a lower bound on a conditional mean, with
no upper estimate implicit in the coefficient.

The localized variance uses the exact full-field Gram matrix and a
triangle inequality before restriction to b∉A. Its numerator retains
the actual tenth and eighth moments of A. A mean bound over random A
does not apply to every deletion of an exceptional input, as the note
correctly states.

The shell sampling identity follows from the two inclusion
probabilities, t/n and (n−t)/(p−n). For t≥1 on the specified slice,
0≤λ_t≤t/n. The separate negative value λ_0 is harmless: its absolute
sixth power times the uniform Weil ratio is below the target. The
necessary retained-label threshold exceeds n^(5/6). The factorial
moment union bound retains its p dependence as (2n²/p)^t before the
coarser (2/n²)^t estimate is used on n⁴≤p<(n+1)⁴. Thus no comparison
at an arbitrarily large p with fixed n has been silently made. The
conclusion excludes this evaluated propagation argument, and does not
exclude stronger arithmetic propagation.

For the sign split, χ(d)=χ(c−z) holds in both odd-prime congruence
classes. The minority size is (n−|F_C(z)|)/2 because the sum of weights
equals χ(−1)F_C(z). The triangle inequality gives the displayed factor
64, hence the normalized balanced-split loss of 8. Inserting the
individual Weil estimate yields a leading sixth-power term 5n;
there is no better exponent hidden in this bound. Finally, the second
moment only gives an upper count of strongly unbalanced poles. A dense
set of unsigned-good poles need not intersect that sparse set on the
basis of those cardinalities alone.

The existing Weil input is the monic quadratic-character specialization
of [McDonald–Sahay–Wyman, Lemma 2.1](https://arxiv.org/html/2210.03789v2),
read live by root in this pass. Its repeated-root application is valid
by removing even multiplicities and retaining the missing-root errors,
as written out in the separate inversion note. No stronger character
estimate is imported by this review.

## Independent computations

The [review verifier](../experiments/parallel24_classical_review_2026_09_05.py)
imports no functions from the lane. It explicitly enumerates the
inside/outside row increments at p=17,29,97, with input sizes 6,7,6.
It verifies 858 increment moments, 143 complete sixth-row expansions,
and three complete swap averages. These sizes extend beyond the
author's small-field identity coverage; they remain outside the quartic
slice and are used for exact identities only. All computations are
Python integers. The [results](../results/parallel24_classical_review_2026_09_05.json)
record the actual sets and sums. Their agreement checks the algebra;
the general positivity and asymptotic volume statements are assessed
from the proof, not inferred from these examples.

The parent suggested explicitly addressing λ_0 and keeping the
p-dependent tail before publication of the lane. The final note
includes both clarifications. No asserted theorem required weakening.
The full classical estimate and goal remain unproved.

## Reviewed artifacts

| Input | SHA-256 |
| --- | --- |
| `research/parallel24-classical-exceptions-2026-09-05.md` | `1844eb2a77942c3e719acb84b34688c8cc8250316a9d8859cf435eab5ef07d86` |
| `experiments/parallel24_classical_exceptions_2026_09_05.py` | `1c1468664e4486258fa7e172b5bf8dbbbe153fba748bde7cf5865ee24899b39b` |
| `results/parallel24_classical_exceptions_2026_09_05.json` | `8bbaa2adc847779220c5a7c5d168a109f984a1cc241c92c512814fb522501b47` |
| `experiments/parallel24_classical_review_2026_09_05.py` | `b9cfea3271a95a6593dcae6f2e5a5adf9a822101fce75e91f5bc1179b45b05e3` |
| `results/parallel24_classical_review_2026_09_05.json` | `ad9b6c289989970fb49da91ae7963187b14c0c88f4a63c6889af7328f292dd96` |
