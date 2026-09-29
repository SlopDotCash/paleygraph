# Separate-agent review of the opposite-pair subgroup lane

Result: no mathematical correction found in the exact decomposition,
its conditional Konyagin payoff, or the symmetric Sidon construction
with its stated parameter conditions. The construction drops
multiplicative closure and does not refute any subgroup conjecture.
The lane proves no new subgroup cancellation exponent.

This is a separate agent's proof and source-code review, not
independent human refereeing or formal verification. I recalculated
the combinatorial factors and parameter inequalities, read the verifier
and its recorded results, checked Python syntax and recorded input
hashes, and checked the algebra in the cited Konyagin proof and
prime-count specialization. I did not rerun the expensive fixture
enumerations or freshly inspect the external analytic paper. Root
reports a fresh primary-source check of Thorner–Zaman Corollary 3.1
and equation (3.2); the external theorem remains an imported input.

## Reviewed bytes

| Input | SHA-256 |
| --- | --- |
| [Subgroup proof](parallel21-subgroup-next-input-2026-09-05.md) | `84005053f3954ee2d36eccf85b2786f1fe00d6bd956c76b23b1b40e82b991c22` |
| [Verifier](../experiments/parallel21_subgroup_next_input_2026_09_05.py) | `5404b71d3189f8454635ddb25baaf0f8d3baa49feb4d83748efcdb947c10d359` |
| [Recorded results](../results/parallel21_subgroup_next_input_2026_09_05.json) | `29b79704857e65bba3e45911d6991b1f6ce087a865c4f92fd05a9d78938243aa` |
| [Konyagin derivation](sigma-subgroup-2026-09-05.md) | `8a5ef969af5bc1c058c53a2e21f19cc440181928de2aaf3533339de9ef98119f` |
| [Pass-4 prime-count specialization](parallel4-subgroup-2026-09-04.md) | `e3f76ed5f91e08470244f69ab115af2933754c731648f41254f88eda17dc4fe0` |

## 1. Deletion really is unique at the multiset level

A nonintrinsic six-term zero-sum word cannot have two disjoint
opposite pairs. Their deletion would leave two entries summing to
zero, so all six entries would be intrinsically paired. Opposite
pairs belonging to different value pairs are necessarily disjoint.
Thus all opposite pairs in a nonintrinsic word use the same value
pair {u,−u}, and the smaller of its two multiplicities is exactly one.
Deleting any occurrence of that opposite pair leaves the same
four-element multiset.

The residual four-multiset is nonintrinsic; otherwise adding the
deleted pair would make the original word intrinsic. A four-term
zero-sum word with an opposite pair is intrinsic because its remaining
two entries must also be opposite. Consequently the residual is
opposite-free. These facts justify counting extensions of multisets
without either a positional overcount or a choice of a second base.

For n=2m, the intrinsic counts are correctly

    T_4=24 binom(m,2)+6m=3n²−3n,
    T_6=720 binom(m,3)+180m(m−1)+20m
       =15n³−45n²+40n.

The middle sixth-order term chooses which of two occupied value
pairs has multiplicities (2,2); its ordered weight is 6!/(2!2!)=180.
No further factor of two is needed.

## 2. Extension multiplicities and U/J corrections

For an opposite-free four-multiset with multiplicities a_v, its
ordered weight is 4!/∏a_v!. Adding a previously absent opposite
value pair multiplies this by 30. If one of the two values is already
present a times, the multiplier is 30/(a+1). The base cannot contain
both signs.

The three possible patterns give the following sums over all n/2
opposite value pairs:

    1111: 30(n/2−4)+4·15 =15n−60,
    211:  30(n/2−3)+10+2·15 =15n−50,
    31:   30(n/2−2)+30/4+15 =15n−75/2.

Pattern 22 would force opposite values and pattern 4 would force
4u=0, excluded by p>3 and u≠0. At small n, a formal coefficient
can be negative only for a pattern which cannot occur: for example
1111 needs at least four occupied opposite value pairs. The proof of
the upper and lower bounds uses actual extension counts, so it is
valid also for n=2 and n=4. It does not rely on treating every formal
coefficient as positive independently of its multiplicity class.

A pattern-31 word has unique tripled value u and remaining value
−3u. Thus its count is A_31=4J. For a pattern-211 word, fix u
as the repeated value. Among the ordered representations of −2u,
the inadmissible equal pair is exactly (−u,−u). When 3u∈S,
the only further forbidden representations are (u,−3u) and
(−3u,u). These three deletions are distinct under p>3. Choosing
the repeated positions gives

    A_211=6Σ_u[r_2(2u)−1−2·1_(3u∈S)]=6U−12J.

The coefficient differences from the 1111 class are 10 and 45/2.
Hence

    10A_211+(45/2)A_31=60U−30J,

which confirms the sign and every coefficient in equation (1).
For a subgroup, multiplication by u maps representations of 2
bijectively to representations of 2u, and 3u∈H iff 3∈H.
Therefore U=n(r_2(2)−1) and J=n·1_(3∈H), proving the stated
compressed formula.

Each extension factor is at most 30, and there are n/2 value pairs.
This directly gives the positive-count bounds in equation (4). In
particular intrinsic fourth energy forces the entire opposite-containing
nonintrinsic sixth count to vanish, so E_3=T_6+R_6 exactly.

## 3. Conditional analytic payoff and finite exception interpretation

From the two P6 assumptions,

    E_3≤[15+(15A+D)L]n³=:B n³,

because T_6≤15n³. The imported inequality at orders (3,3) is

    M≤(p E_3²)^(1/18)n^(1/3).

Using p≤n⁴ gives precisely B^(1/9)n^(8/9). The lower end of
the quartic interval is not needed for this particular substitution.
I checked the cited proof's weighted Hölder, Fourier inversion and
bilinear orthogonality steps; they do not require an even order or
positivity of the inverse Fourier transform of |η|³.

The payoff remains conditional on uniform A,D and the specified
growth of L. With a fixed logarithmic power it preserves the exponent
8/9 up to logarithmic factors; it does not yield the square-root target.
The recorded finite computations cannot
promote either P6 assumption to a uniform theorem.

The finite exception arithmetic is consistent. At n=64, the baseline
coefficient is 15n−60=900, so an excess E_2−T_4=1536 contributes
1,382,400 before the K correction. At p=11127041, K=5 contributes
60·64·4=15,360 more, producing 1,397,760. The resonant case with
excess 768 and K=3 gives 900·768+60·64·2=698,880. The fourteen
listed R_6 values are all below 2·64³. These observations explain
the listed fixtures only and do not classify all prime exceptions.

## 4. Symmetric Sidon construction and asymptotic quantifiers

The integer labels satisfy 0≤c_t<2q². Equality of two pair sums
recovers the ordinary index sum modulo 2q with no carry, then the
square sum modulo q. Since q is odd, the sum and square sum recover
the product. The unordered residue pair is therefore determined, and
the indices lie in [0,q), making the unordered integer pair unique.
Translating by 8q² preserves this Sidon property.

The parameter chain is valid, including its endpoint:

    6max B<60q²<240k²<4k⁴=n⁴/4≤p,

since k≥8 implies k²>60. Thus all signed sums of at most six
positive representatives have absolute value less than p and cannot
vanish solely from modular wraparound. Also 8q²≤b_t<10q².
This makes B and −B disjoint nonzero sets, excludes 1, and forces
two positive and two negative entries in a signed four-term relation.
The Sidon property then forces intrinsic pairing, giving E_2=T_4.

For six-term relations, 4min B>2max B excludes four positive and
two negative entries, and the more imbalanced cases are excluded as
well. Every relation has exactly three positive and three negative
entries. There are exactly binom(6,3)=20 choices for the positive
positions, giving E_3(S)=20E_3(B). The ordered triple sums of B
lie in an integer interval with at most 6q²+1≤25k² values. Therefore

    E_3(S)≥20k⁶/(25k²)=n⁴/20.

The claimed fourth-power scale also has an elementary upper bound:
Sidon implies r_2^B(x)≤2, hence r_3^B(x)≤2k and
E_3(B)≤2k·Σ_x r_3^B(x)=2k⁴. Thus E_3(S)≤(5/2)n⁴.
This additional observation concerns the constructed sets only.

This is compatible with intrinsic fourth energy but makes E_3/n³
at least n/20. It exceeds every fixed power of log n along an
unbounded parameter sequence. The lower estimate for R_6 subtracts
T_6=O(n³); it becomes substantive for sufficiently large n, and is
not asserted positive from its displayed expression at every finite n.

The construction theorem is conditional on its stated q and p.
To obtain an unbounded dyadic sequence without asserting a new prime
gap estimate, choose arbitrary odd primes q tending to infinity and
take k as the largest power of two at most q. Then k≤q<2k and
n=2k tends to infinity through a dyadic subsequence. The imported
pass-4 specialization supplies p≡1 mod n in the quartic interval
for every sufficiently large such n. Its parameter inequality is
h/φ(n)=(3/2)n³≥(n⁴)^(7/12+1/12)=n^(8/3).
No numerical threshold, a prime-prefix scan, or theorem about primes
between consecutive fourth powers is being substituted for that input.

Absence of 1 is a decisive lost subgroup property. The examples
retain symmetry, prime-field order, dyadic cardinality, the quartic
window and intrinsic fourth energy, but they do not retain
multiplicative closure. Thus they establish only the insufficiency
of those relaxed properties; no subgroup counterexample follows.

## 5. Verifier audit and coverage limits

The energy calculation convolves ordered pair counts and squares
ordered triple counts. The R_6 calculation is independent: it uses
unordered triples with factorial ordering weights, rejects opposite
pairs within each triple, and rejects opposite pairs across the two
halves by a bit-mask intersection. The two cross-intersection tests
are equivalent by negation, so only one is needed.

For sums s≠−s, restricting to one of the two buckets and multiplying
by two restores the two orders of the halves. For s=0, nested loops
already retain both half orders, so the factor one is correct.
Sequential integer division in the factorial weight is exact because
the full product of multiplicity factorials divides the tuple factorial.

The recorded 53 general decompositions, 44 subgroup specializations,
34 multiplicity checks, five Sidon constructions and 86,922 admissible
weighted triple-pair terms match the verifier's parameter lists and
result file. Only n≤32 Sidon constructions receive the additional
direct R_6 enumeration; for larger constructions R_6 is obtained from
the proved identity. The lower bound on E_3 is checked directly in
all five fixtures, but the finite program is not the proof of its
all-parameter interval-count argument.

Python syntax and all input hashes already recorded by the result
file passed this review's checks. One provenance limitation is that
the result file does not itself pin the accompanying proof note or
the pass-4 prime-count dependency. Both are pinned in this review
table. This is an artifact-scope observation, not a mathematical gap
or a request to alter historical files.

No correction is required within the stated claims. Separate human
review of the analytic input and formal verification remain distinct
from this completed agent review. The uniform subgroup square-root
target, its exceptional primes, the classical Paley conjecture and
the exact official-prize bridge remain unproved.
