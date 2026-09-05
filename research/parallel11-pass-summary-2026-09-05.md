# Eleventh pass: the complete spectral aggregate at logarithmic depth

**Progress: the bound now applies to the full spectral trace in a stated
logarithmic depth range, with every correction retained.** The
[written proof](parallel11-full-energy-2026-09-05.md) is checked locally
and depends on the preceding source-dependent necklace arguments.
Independent mathematical review remains outstanding. The full Paley
and official prize goals remain unproved.

For an a-anchor clique, put b=2^a, N=b−1, let E select its common
neighborhood, and set T=(S/√p−J/p)(bE−I)/√N. The positive energy
H_j=||T^j||_F² bounds the exact even trace: |tr(T^(2j))|≤H_j.
It also obeys two exact identities:

\[
 H_{j+1}=H_j/N+(N-N^{-1})\|ET^j\|_F^2,
 \qquad
 H_{j+1}-2H_j+H_{j-1}=\frac{(N-1)^2}{N}\operatorname{tr}T^{2j}.
\]

The actual word coefficients are all one in the soft-mask expansion.
Their raw-rank budget has explicit growth rate
Θ_a=[2^(a−1)(a+1)−1]²/(2^a−1). For a=2 this is 25/3.
Combining that budget with the preceding correlation and weight-gap
bounds gives a nonasymptotic estimate for H_j. The proof handles all
finite singular rows, exceptional anchor columns, the constant direction
missing from S, the rank-one J term, and the exact anchor correction
inside the powers.

Consequently, for every fixed a≥2 and ε∈(0,1/2), uniformly over the
anchor cliques,

\[
 j\le\frac{(1/2-\epsilon)\log p}{\log\Theta_a}
 \quad\Longrightarrow\quad
 H_j\le(1+o(1))p,\quad |\operatorname{tr}T^{2j}|\le(1+o(1))p.
\]

This is the complete normalized aggregate, extending the previous
open-curve word-span result. It does not provide an improved clique
bound: the required edge criterion needs depths with j/log p→∞,
which are outside the proved range.

The next sufficient estimate can now be written in terms of actual
row-energy fractions. With ρ_j=b||ET^j||_F²/H_j−1, the logarithm
of H_j/p is exactly Σ_(l<j)log(1+(N−1)ρ_l/N). An o(j) upper
bound for that sum at the required long depths would supply the
spectral criterion. This is still unproved. A converse comparison
between H_j and |tr(T^(2j))| shows that the positive formulation
retains the original difficulty, up to explicit factors. No random
row distribution, positive bias, or monotonic energy is assumed.

The [verifier](../experiments/parallel11_full_energy_2026_09_05.py)
and [results](../results/parallel11_full_energy_2026_09_05.json)
record 28 clique cases, 168 energy recurrences, 140 curvature identities,
168 trace/energy comparisons in each direction, 168 full upper bounds,
168 constant-direction checks and 672 boundary checks. Thirty exact
all-word matrix sums and 117,120 pointwise bounds test the raw-path
normalization and singular points. An abstract Jordan-block example
checks the exact quadratic energy at the spectral edge; a separate
p=10009 case has a certified full upper bound below 2p.
The [artifact audit](../results/parallel11_pass_audit_2026_09_05.json)
pins final files and sources and confirms preservation of pass ten.

The subgroup exponent remains 8/9 on the existing density-one quartic
prime class. Its uniform target and exceptional primes, the classical
arbitrary-two-set conjecture, and the exact official prize bridge remain
open. This pass was completed locally; no external submission or formal
verification is claimed.
