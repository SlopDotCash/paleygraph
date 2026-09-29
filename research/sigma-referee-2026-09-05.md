# Sigma pass, direction `referee`: refereeing parallel9 (all-degree necklace bound) and parallel10 (word aggregate)

**Status:** Neither refereed note is REFUTED: 1,417,148 exact checks (full run,
309 s; Section 8) found no integer witness against parallel9's inequality (T)
(dense traces to p = 3389, a = 1 to p = 2*10^6), against its fixed-variable
inner-sum bound (to p = 1,000,033), or against parallel10's Gram bound (10)
(four non-vacuous cases at p = 100049 and 1,000,033). Both notes are classified **PLAUSIBLE BUT UNREVIEWED
GAPS**: every step is either re-derived here (PROVED), or a verbatim-quoted Katz
theorem (CITED), or a standard Deligne/BBD weight statement whose local source
(the BBD scan) has no text layer and which this pass could not re-verify
(UNVERIFIED CITATION). No mathematical error was found; the gaps are
citation-level (arithmetic vs. geometric statements, Weil II 1.8.4 / 3.3.1, BBD
5.1.14 / 5.3.1). PROVED here independently: the exact indicator expansion of
Task 3; the rank-growth lemma (R)/(R') and the pseudoreflection structure as a
consequence of Katz 3.3.6-3.3.7; parallel10's word recovery (1)-(2); the raw-rank
recurrence (6) by an Euler-characteristic computation; the dimension obstruction
(13). OPEN: Kunisky's Conjecture 1.14 at degree a>=2 remains proved only
conditionally on the cited inputs; nothing here touches the two-set Paley
conjecture.

Definitions and notation follow the two refereed notes; `p` is always a prime
with p = 1 mod 4, `chi` the Legendre symbol with chi(0)=0, `S_xy = chi(x-y)`,
`d_Z(x) = prod_{z in Z} chi(x-z)`, `D_Z = diag(d_Z)`.

## 1. What the two notes claim, and how it compares with Kunisky

### 1.1 Claim (a): parallel9-all-degrees-2026-09-05.md

Let A be a set of a >= 1 distinct anchors (a < p) and Z_1,...,Z_k nonempty subsets
of A. Put N = tr(D_{Z_1} S D_{Z_2} S ... D_{Z_k} S). The note claims, for k >= 3,

    (T)   |N| <= 3a (2a+2)^{k-2} p^{(k+1)/2},

together with N = 0 for k = 1 and |N| <= a^2 p for k = 2, and states that (T)
proves Kunisky's Conjecture 1.14 at every fixed degree a. Its proof: fix y not in
A; build the "principal" sheaves F_0 = L_chi(x-y), G_i = ME(F_{i-1} tensor
L_{d_{T_i}}), F_i = MC_chi(G_i) (quadratic middle convolution), show rank F_i =
d_i >= i+1 with nondecreasing increments (R'), track the difference between the raw
compact-convolution objects K_i (whose trace functions are exactly the matrix
chains R_i = S D_{T_i} R_{i-1}, with all zero masks) and F_i[1] through an exact
sequence 0 -> E_i -> K_i -> F_i[1] -> 0 with E_i of perverse weight <= i (W), bound
the total constituent mass by c_i <= (2a+2)^i (C), and close the necklace by
writing N = sum_y d_{Z_k}(y) I(y), I(y) = sum_x d_{Z_1}(x) chi(x-y) R_{k-2}(x,y),
bounding |I(y)| <= (2a-1) c_{k-2} p^{(k-1)/2} for y not in A and summing over y
trivially.

### 1.2 Claim (b): parallel10-word-aggregate-2026-09-05.md

With the same objects, for a word w = (T_1,...,T_m) put U = F_p \ (A u {y}),
h_w(x) = p^{-|w|/2} Tr(Frob_x | F_w), q_w(x) = (-1)^{|w|} p^{-|w|/2} R_w(x,y),
d_w = rank F_w, r_w = generic rank of H^{-1}(K_w). The note claims:
(1)-(2) the word is recovered from the local monodromy of F_w, hence F_u and F_v
are geometrically non-isomorphic for u != v; (3) F_w^vee = F_w(|w|); (4) for any
finite set W of words and any complex vector z,
|p^{-1} sum_U |sum_w z_w h_w|^2 - ||z||^2| <= eps_W ||z||^2 with
eps_W = (a D_W + 1)/sqrt p, D_W = sum d_w^2; (5) a cell version with s fresh
adjacency conditions; (6) the raw-rank recurrence c_y = 1, r = 1 + sum_{v in A} c_v,
c'_v = r for v in T and c_v otherwise; (7)-(8) closed recurrences for the word sums
W_m, R_m, Q_m = sum r_w^2, Z_m with Perron growth Lambda_a (Lambda_2 =
(9+sqrt 65)/2); (9) |q_w - h_w| <= (r_w - d_w)/sqrt p on U; (10) the same
quadratic bound as (4) for the raw q_w with error eps_W + 2 gamma_W sqrt(1+eps_W)
+ gamma_W^2, gamma_W = sqrt(L_W/p), L_W = sum (r_w-d_w)^2; (11)-(12) cell
versions; (13) if (2^a-1)^m > p-a-1 the evaluation map has a kernel, so no
all-coefficient isometry with error < 1 can hold.

### 1.3 Kunisky verbatim (sources/kunisky-2303.16475v1.html, Section 1.5-1.6)

Definition 1.13: "Let k >= 1 and Z_1,...,Z_k subseteq F_p. The associated necklace
character sum is Sigma(Z_1,...,Z_k) := sum_{x_1,...,x_k in F_p}
chi(x_2-x_1)...chi(x_k-x_{k-1})chi(x_1-x_k) prod_{i=1}^k prod_{z in Z_i}
chi(x_i-z)", and (15): "Sigma(Z_1,...,Z_k) = Tr(D_1 S_{G_p} ... D_k S_{G_p})".

Conjecture 1.14 (Necklace character sum estimates, degree a): "For all k >= 1,
lim_{p -> infinity} p^{-(k/2+1)} max_{Z_1,...,Z_k subseteq F_p, Z_i != emptyset,
|Z_1 u ... u Z_k| <= a} |Sigma(Z_1,...,Z_k)| = 0."

Trivial bound (17): "|Sigma(Z_1,...,Z_k)| <= p^{k/2+1}."

Remark 1.16: "when Z_1 = ... = Z_k = {z} ... there is a 'spurious degree of
freedom' in this sum, and it can in fact be rewritten as ... (p-1) sum_{x in
F_p^{k-1}} chi(q~(x)) ... Accordingly, we will only be able to show that this sum
(for fixed a and k) is of order O(p^{(k+1)/2}). We will see in Appendix B.1 that
the same happens when Z_1 = ... = Z_k = {z_1,z_2} ... It is reasonable to
conjecture that there are a few special cases for which |Sigma| = O(p^{(k+1)/2}),
and that, outside of those cases, |Sigma| = O(p^{k/2})."

Theorem 1.17: "If Conjecture 1.14 holds at degree a, then Conjecture 1.8 holds at
degree a." Theorem 1.18: "Conjectures 1.14, 1.8, and 1.9 all hold at degree
a = 1." Kunisky's proof of 1.18 (Section 5, (81)/(83)) gives
|Sigma({z},...,{z})| <= k p^{(k+1)/2} via Katz's Gauss-sum equidistribution
(his Theorem 3.9 = [Kat88, 9.6]); Theorem B.1 gives
|Sigma({z,z'},...,{z,z'})| <= k p^{(k+1)/2} + 2^k p^{k/2}.

### 1.4 Differences between the notes' definitions and Kunisky's

1. Same sum. Since S is symmetric for p = 1 mod 4, N = Sigma(Z_1,...,Z_k)
   exactly (checked literally at p = 13 for three words, Section 8).
2. Same range. Kunisky maximizes over nonempty Z_i with |Z_1 u ... u Z_k| <= a;
   the note takes nonempty Z_i inside a fixed a-set A. These coincide (enlarge the
   union to size a; (T) is monotone in a).
3. Exponent. (T) gives p^{-(k/2+1)} |N| <= 3a(2a+2)^{k-2} p^{-1/2} -> 0 at fixed
   (a,k), so (T) does imply Conjecture 1.14 at degree a. Note that (T) is NOT the
   square-root bound p^{k/2} that Remark 1.16 expects for generic words; it is the
   weaker uniform exponent (k+1)/2 that Kunisky himself proves at a = 1 and for
   the constant two-anchor word.
4. Constants. At a = 1 Kunisky's own Theorem 1.18 gives the constant k, versus
   3 * 4^{k-2} in (T); Kunisky's result is stronger at a = 1 for all k >= 3. The
   note does not say so; it should.
5. Kunisky restricts nothing about p beyond p = 1 mod 4 (needed for the Paley
   graph); the note adds a < p, harmless.
6. Claim (b) has no counterpart in Kunisky; its h_w, q_w, U and the cell C are the
   note's own objects. The trace functions q_w are literally the columns of the
   matrix chains, so (10) is a statement about the Paley matrices.

## 2. Numerical audit of claim (a)

### 2.1 Data

Previous worker (this direction, killed mid-way; artifacts reused, not
recomputed): explore1 (dense exact traces, 8 patterns, 27 primes 101..1493, k <= 8);
explore2 (a = 1 constant word, exact via (p-1) x FFT inner sum, 45 primes up to
1,999,957, k <= 8 for p < 3*10^5, else k <= 7); explore3 (fixed-y inner sums for
7 patterns with a = 2,3, 70 primes 1009..493,657, k <= 8, three random y each);
explore4 (dense exact traces, 9 patterns incl. random words and anchors, primes
1601..3413, k <= 8).

This pass (experiments/sigma_referee_2026_09_05.py): Part A, dense exact traces
for 11 patterns (a = 1,2,3; structured and random words; random anchor positions
per prime), 18 primes 101..3389, k = 1..8 (the 2^53 exactness guard allows
k = 8 up to p = 3800); Part B, exact fixed-variable inner sums by integer FFT
convolution at p = 10009, 100049, 1000033 for the same patterns, with generic
and anchor values of the fixed variable. All arithmetic is integral; the
float64 matrix products are exact under the asserted guard p^{(j+1)/2} < 2^53
and were cross-checked against int64 and Python-integer products.

### 2.2 Results for (T)

- No violation. Every tested (a,k,p,word) satisfies N^2 <= (3a(2a+2)^{k-2})^2
  p^{k+1} exactly. Combined with explore1/explore4 this covers all listed
  patterns at every prime p = 1 mod 4 sampled below 3413, and a = 1 up to
  2*10^6.
- Non-vacuity. (T) beats the trivial bound p^{k/2+1} iff p > p_0(a,k) :=
  (3a(2a+2)^{k-2})^2: p_0(1,k) = 144, 2304, 36864, 589824, 9.4*10^6, 1.5*10^8
  for k = 3..8; p_0(2,3) = 1296, p_0(2,4) = 46656, p_0(2,5) = 1.7*10^6;
  p_0(3,3) = 5184, p_0(3,4) = 3.3*10^5. Hence among dense computations only
  (a,k) in {(1,3),(1,4),(2,3)} are non-vacuous below p = 3500, and (3,3) only
  beyond 5184 (not reached by dense matrices here). The a = 1 data of explore2
  is non-vacuous for k = 3,4,5,6 (p up to 2*10^6); the largest observed ratio
  |N| / bound is 0.167 at (a,k) = (1,3) (i.e. |N| ~ 2 p^2 against 12 p^2).
- Exponents. Log-log slopes of |N| against p (nonzero traces):
  a = 1 constant word, explore2 (p up to 2*10^6): 2.03, 2.49, 2.88, 3.55, 4.00,
  4.55 for k = 3..8 against the claimed (k+1)/2 = 2, 2.5, 3, 3.5, 4, 4.5; the
  normalized values |N|/p^{(k+1)/2} stay in [0.01, 3.0] with many values near 2
  for k = 3, i.e. (T)'s exponent is sharp at a = 1 (as Kunisky's (20) shows
  exactly: N = (p-1) x a (k-1)-variable sum of size p^{(k-1)/2}).
  a = 2 constant word {z_1,z_2}: traces agree with the a = 1 word to relative
  accuracy ~10^{-3} at every prime (e.g. p = 1601, k = 3: -1.996 vs -1.993 in
  units of p^2), consistent with Kunisky's Theorem B.1: the Moebius map
  x -> x/(x-1) conjugates D_{01} S D_{01} to D_0 S D_0 up to one row/column.
  Other words (Part A, 18 primes to 3389, medians over p >= 1000 of
  |N|/p^{(k+1)/2}): the alternating word ({0},{1},{0},{1},...) is degenerate at
  even length (medians 0.79 at k = 6 and 0.73 at k = 8, maxima 1.64 and 1.17)
  but far below at odd length (medians 6*10^{-4}, 6*10^{-4}, 1.2*10^{-3} at
  k = 3, 5, 7; at k = 3 the values are ~1.3 p, i.e. O(p)); a random two-anchor
  word has medians <= 0.03 for all k; the three-anchor words have medians in
  [8*10^{-4}, 0.26] (a3_full k = 8: 0.26, max 0.39; random three-anchor words
  8*10^{-4} to 0.025) with ratios to p^{k/2} between 0.03 and 13, i.e. sizes
  between p^{k/2} and p^{(k+1)/2}. For the constant
  three-anchor word the odd-length traces vanish identically for 14 of the 18
  (p, anchor) configurations (the cyclic Moebius permutation of the three
  anchors conjugates D S to -D S when its multiplier is a non-square; not
  pursued). Full table: results/sigma_referee_2026_09_05.json,
  A_dense_traces.summary.
- Conclusion of Section 2.2: (T) is numerically consistent; the exponent (k+1)/2
  is attained by the a = 1 word and its Moebius conjugates and is an upper
  envelope for the others. No refutation is possible at accessible p because the
  constant 3a(2a+2)^{k-2} exceeds every observed normalized value by a factor
  >= 6.

### 2.3 Results for the inner sums (the actual content of the proof)

The proof bounds, for each fixed value y of one necklace variable, the
(k-1)-variable sum I_k(y) (here realised as the diagonal entry
[S D_{Z_2} S ... D_{Z_k} S]_{yy}, a cyclic rotation of the note's I(y)) by
(2a-1)(2a+2)^{k-2} p^{(k-1)/2} for y not in A, and by the trivial p^{k/2} for
y in A. Part B and explore3 test exactly this.
- No violation at p = 10009, 100049, 1000033 (Part B: 528 generic-y checks,
  270 anchor-y checks) nor in explore3.
- The bound (2a-1)(2a+2)^{k-2} p^{(k-1)/2} is non-vacuous against p^{k/2} iff
  p > ((2a-1)(2a+2)^{k-2})^2: a = 2: 324, 11664, 419904, 1.5*10^7 for k = 3..6;
  a = 3: 1600, 102400, 6.6*10^6 for k = 3..5. So Part B tests (2,3),(2,4),(2,5)
  and (3,3),(3,4) non-vacuously, explore3 the same set.
- Exponents: explore3 slopes of max_y |I| against p (70 primes) are within 0.1 of
  (k-1)/2 for every pattern and every k = 3..8 (a2_full: 1.07, 1.48, 1.96, 2.59,
  2.95, 3.36; a3_cyc: 0.95, 1.47, 2.03, 2.54, 3.05, 3.56), with
  max_y|I|/p^{(k-1)/2} having medians 1.0-2.9 and maxima <= 6.3. Part B at
  p = 10^4, 10^5, 10^6 gives generic-y ratios in [0.07, 3.7] for every pattern
  and k (a3_full k = 4: 2.96, 3.66, 1.48; a2_alt k = 7: 1.19, 2.87, 1.43), with
  three-point slopes within 0.3 of (k-1)/2 except where a small value at
  p = 10009 distorts the fit (see B_inner_sums.summary).
- Anchor fibres are genuinely degenerate, and the note's trivial bound for them
  is sharp. Part B evaluates the fixed-variable sum also at anchors. Killed
  anchor (y in Z_1, contributes nothing to N): for the word ({0},{0,1},{0},{0,1})
  and y = 0 the exact value is (u_4)_0 = p(p-2) + O(p) at p = 10009 and 100049
  (ratio to p^{k/2} = p^2 equal to 1.000). Surviving anchor (y in A \ Z_1, which
  the note bounds by the operator norm p^{k/2}): for the three-anchor word
  ({0,1},{1,2},{0,2},{0,1},{1,2},{0,2}), k = 6, y = anchor 2, the value is
  p^3 (1 + O(1/p)) at p = 10009 and 100049 (ratios 100.005 and 316.293 to
  p^{5/2}, i.e. 1.0000 to p^{k/2}); in every other tested (word, k, anchor)
  case the ratio to p^{(k-1)/2} is at most 3.3. Mechanism: with y an anchor,
  L_chi(x-y) is itself an anchor Kummer sheaf; the first twist containing y
  removes the singularity at y, and u_2(x) = sum_t chi(x-t) chi(t-z') 1_{t != y}
  has the spike u_2(z') = p-2 at another anchor z', which then propagates. So
  the exponent (k-1)/2 is NOT uniform over y in F_p: the note's three-way split
  (generic y by cohomology, surviving anchors by the trivial p^{k/2} bound,
  killed anchors by the mask) is exactly what is needed, the trivial bound is
  attained, and (T) is unaffected because (a-1) p^{k/2} is p^{-1/2} below the
  main term. The pass-9 verifier never evaluated anchor fibres; this is new.

Conclusion: the one genuinely new analytic input of (a), "every fixed-variable
sum is O_a,k(p^{(k-1)/2}) uniformly", is numerically solid up to p ~ 10^6 for all
patterns tested, including the words where the full necklace is much smaller.

## 3. Elementary reformulation (Task 3)

### 3.1 The exact identity (PROVED; verified exactly, Part C)

For Z subset F_p let E_Z(x) = 1_{x not in Z} prod_{z in Z} (1 + chi(x-z))/2 be
the common-neighbourhood indicator. Then, pointwise on F_p,

    E_Z = 2^{-|Z|} sum_{W subseteq Z} D_W  -  B_Z,
    B_Z(z) = (1/2) 1[ z adjacent to every z' in Z \ {z} ] for z in Z, B_Z = 0 off Z.

Proof: off Z expand the product; at z in Z the product sum_{W} d_W(z) equals
2^{|Z|-1} prod_{z' != z}(1 + chi(z-z')) since chi(0) = 0 kills every W containing
z, while E_Z(z) = 0. When Z is a clique, B_Z = (1/2) 1_Z, which is parallel2's (8).
Multilinearity then gives, for any diagonal chain,

    tr(E_{Z_1} S ... E_{Z_k} S) = 2^{-sum |Z_i|} sum_{W_i subseteq Z_i}
        Sigma(W_1,...,W_k)  +  (terms with at least one B_{Z_i}),

with empty W_i allowed. Each boundary term contains a factor B supported on at
most a points, so it is at most (a/2) ||S(...)S||_op <= (a/2) p^{k/2}: boundary
terms are O_a(p^{k/2}), one power of sqrt p below (T) and of the same order as
the generic necklace size. Empty labels reduce by S^2 = pI - J to shorter
words plus rank-one J terms (verified). Part C verifies the pointwise identity
(88 instances), the clique case (69), and the trace expansion (30 random
(p, anchors, word) cases with p <= 53, a <= 3, k <= 4) by exact integer
enumeration of all 2^{sum|Z_i|} subset words.

### 3.2 Does p^{(k+1)/2} follow from standard multivariable bounds?

Sigma(W_1,...,W_k) = sum_{x in F_p^k} chi(q(x)) with q a product of linear
forms: the character sum of the rank-one local system L_chi(q) on the
complement M of a hyperplane arrangement in A^k. The standard sufficient
condition for square-root cancellation |Sigma| <= C p^{k/2} is that H^i_c(M,
L_chi(q)) vanish for i != k, which for arrangements is the non-resonance
condition of Esnault-Schechtman-Viehweg / Schechtman-Terao-Varchenko: for every
dense edge (intersection of hyperplanes) the sum of the exponents of the
hyperplanes through it must not be an integer (UNVERIFIED CITATION: not fetched
in this pass; stated for orientation only). Here every exponent is 1/2, so every
edge through an even number of hyperplanes is resonant. Resonant edges exist for
every word: the point x_1 = ... = x_k = z (z in every W_i) lies on the k edge
hyperplanes and on all sum |W_i| anchor hyperplanes; for a = 1 the whole
arrangement is central at that point with 2k hyperplanes, and the resonance is
exactly the G_m-symmetry that produces the factor (p-1) in Kunisky's (20). So the
standard non-resonant bounds do not apply directly to any necklace; the
genuinely degenerate words are those with a resonant central point of the
whole arrangement (a = 1, and its PSL_2-conjugates such as the constant
two-anchor word), where H^{k-1}_c is nonzero of weight k-1 and the sum is
p^{(k+1)/2} up to constants. For the other words the data of Section 2 show
sizes between p^{k/2} and p^{(k+1)/2}, and for some (alternating word at k = 3)
even O(p); their cohomology is not "generic" either, but the top weight k+1
never appears beyond the a = 1 degeneration in the data.

### 3.3 Does the cyclic structure force the loss of exactly one sqrt p?

Not for all words; the loss is a feature of the proof method, not of the sums.
Inductive one-variable convolution (the Katz middle-convolution chain of the
note) proves that every (k-1)-variable fixed-y sum I(y) is O(p^{(k-1)/2}), which
is the correct square-root bound for those sums (they are sums over the
complement of a non-central arrangement whose resonance is resolved one
variable at a time). The last variable y is then summed with absolute values,
which costs a factor p instead of sqrt p: this is where the single sqrt p is
lost, and it is unavoidable for the method because the map y -> I(y) is the
trace function of a sheaf of rank ~ (2a+2)^{k-2} on the y-line whose
cancellation over y is exactly the statement that the k-variable sum is
O(p^{k/2}); nothing in (a) controls it, and for the a = 1 word there is no
cancellation to be had: scaling x_i -> y x_i gives I(y) = chi(y)^{2k-1} I(1) =
chi(y) I(1) on F_p^*, so sum_y chi(y) I(y) = (p-1) I(1) exactly (Kunisky's (20);
checked exactly against dense traces in explore2). So: the loss of one
sqrt p is forced for the degenerate words of 3.2, and is a method artefact for
the others. Passing from (T) to Kunisky's expected p^{k/2} for the generic words
would require the H^{k}_c weight analysis of the full k-variable arrangement,
which neither note attempts.

## 4. Line-by-line classification of the argument of (a)

Katz citations are quoted from sources/katz-rigid-local-systems.pdf via the
text dump katz-rls.txt; section numbers are Katz's.

A1. Trace identity (F): Tr(Frob_x | K_i) = (-1)^{i+1} R_i(x,y) at every x,
including the finite singular set. PROVED (standard): Grothendieck-Lefschetz
for R pi_! plus proper base change; the sign follows from Tr(L_chi[1]) = -chi and
the [2] shift of the compact convolution. The pass-9 verifier checked the sign
convention on integer matrices; my Part B re-derives N = sum_y d_{Z_1}(y)(u_k)_y
against dense traces (22 words, p = 101, 401).

A2. Exactness of Phi_T and preservation of the stratification A u {y} and of
tameness. CITED: Katz 2.3.3 ("for a perverse object K on Y, both j_!K and Rj_*K
are perverse on X, and as functors ... both j_! and Rj_* are exact"), 2.6.1-2.6.2
("property P_! ..., the functor L -> L*_!K is an exact functor from Perv(G) to
itself"), 2.9.1 ("j_*L_chi[1] on A^1 has property P"), 3.3.6 ("suppose in
addition that K is tamely ramified everywhere. Then so is K*mid j_*L_chi[1]").
The stratification claim for the raw objects needs the case analysis of 3.3.3
(constant -> 0; Kummer translate -> punctual; punctual -> Kummer translate), which
the note performs; PROVED given 3.3.3.

A3. Irreducibility and "type 2d" of every G_i, hence of F_i. CITED: Katz 3.3.3,
which is stated for "K perverse irreducible on A^1" over an algebraically closed
field with "2d) If K is a perverse irreducible in P which is not of type 2a, 2b,
or 2c, then K*mid j_*L_chi[1] is a perverse irreducible in P which is not of type
2a, 2b, or 2c", and 2.9.7 (MC_chi and MC_chibar "are inverses of each other").
G_1 is rank one with |T_1|+1 >= 2 finite singular points, so not a Kummer
translate (2c/2b) nor punctual (2a); for i >= 2, rank G_i = d_{i-1} >= 2. PROVED
given the citations and (R').

A4. Local monodromy rules and rank formula. CITED verbatim: Katz 3.3.6 "(1) at
s in A^1 - U: (F(s)/F(s)^{I(s)}) tensor L_chi(x-s) = G(s)/G(s)^{I(s)}. (2) at
infinity: there exists a tame I(infinity)-representation M with F(inf) =
M/M^I, G(inf) = M tensor L_chi/(M tensor L_chi)^I ... rank M = sum_s rank
(F(s)/F(s)^I)"; 3.3.7 "rank G = sum_s rank(F(s)/F(s)^I) - rank((F(inf) tensor
L_chi)^I)". The note's block-count formulas i_v^+(F) = delta + i_v^+(G),
i_v^-(F) = i_v^+(G) - j_v^+(G) and the infinity analogues are consequences (a
tame I-representation with eigenvalues +-1 is determined by its Jordan blocks;
V/V^I shrinks unipotent blocks by one and V is recovered from V/V^I and dim V).
PROVED: re-derived in this note's Part D implementation and checked on every
transition to depth 14/7/5/4 for a = 1/2/3/4 (77,140 transitions in the full run).

A5. Pseudoreflection at y: F_i(y) = chi + 1^{d_i-1} (i even), J_2 + 1^{d_i-2}
(i odd), invariant codimension one. PROVED (induction from F_0(y) = chi through
A4; verified on every transition).

A6. Rank growth (R): d_{i+1} - d_i = -delta_i + sum_{v in F(T)} t_v(F_i) >=
(|F(T)|-1) delta_i >= delta_i, and (R'). PROVED: the untwisted rank formula uses
the involution 2.9.7 (rank MC(F_i) = rank G_i = d_{i-1}); twisting exchanges the
two block counts at v in T and at infinity iff |T| odd; t_v(F_i) >= delta_i by
A4. I re-derived the identity and the verifier confirms (R) on every transition
and d_i >= i+1, delta nondecreasing, d_{i+1} <= a d_i + 1.

A7. Weight gap (W): 0 -> E_i -> K_i -> F_i[1] -> 0 with E_i of weights <= i,
F_i[1] pure of weight i+1. Structure PROVED (re-derived: twist-and-extend of the
previous sequence, punctual kernel of j_! -> j_{!*} with stalks the new inertia
invariants, compact-to-middle sequence with constant kernel, weight bookkeeping).
Inputs: CITED Katz 2.9.4(3) verbatim: "The middle and ! additive convolutions of
K with j_*L_chi[1] sit in an exact sequence of perverse sheaves on A^1: 0 ->
(the constant sheaf (H^{-1}(K) tensor L_chi)^{I(inf)})[1] -> K *_! L_chi[1] ->
K *_mid L_chi[1] -> 0"; CITED Katz 5.5.5.10 for purity of the middle
convolution: "the image of the 'forget supports' map is precisely the weight
= w+1 quotient of the group H^1_c(U, G tensor L_chi(t-x)), which is a priori of
weight <= w+1" (this is Katz's own use of [De-Weil II, 3.2.3]). UNVERIFIED
CITATION (locally): the weight statements "inertia invariants of a pure lisse
sheaf of weight w have weights <= w" (Weil II 1.8.4) and "j_!, R f_!, external
product preserve weights <= w" (Weil II 3.3.1; BBD 5.1.14), "H^q stalks of a
perverse sheaf of weights <= w have weights <= w+q" (BBD 5.1.8), "j_{!*} of pure
is pure" (BBD 5.3.1-5.3.2). The BBD PDF in sources/ has no text layer
(pdftotext returns nothing), so this pass could not quote it; the note reports
visual inspection. These are textbook statements and I have no doubt about
them, but under the brief's rules they are unverified here. Also unverified:
the passage from Katz's geometric statements (k algebraically closed) to the
arithmetic objects over F_p with their Frobenius actions (base change of the
image of a morphism of perverse sheaves; Frobenius on V_i induced from the
inertia invariants at infinity). Standard, but asserted rather than argued in
the note.

A8. Constituent mass (C): c(Phi_T K) <= (2a+2) c(K). PROVED given A2/A3: for a
simple middle extension of rank r, at most |T| r <= a r punctual constituents
from extension by zero, middle convolution of rank <= (a+1) r by 3.3.7, constant
kernel of rank <= r by 2.9.4(3); punctual constituents become rank-one Kummer
translates; constants convolve to zero (3.3.3 1a) or, after a nonempty twist,
become rank-one middle extensions of mass <= a+2. Additivity of c and exactness
of Phi_T give the bound. Re-derived; the pass-9 verifier checks new + punctual +
infinity <= (2a+2) old on every transition.

A9. Pointwise bound (P): |R_i(x,y)| <= c_i p^{i/2} at all finite x. PROVED given
A7's weight facts: middle extensions have no H^0 stalks and H^{-1} stalks (inertia
invariants) of weight <= i; E_i has H^{-1} stalks of weight <= i-1 and H^0 stalks
of weight <= i; total stalk dimension <= mass by subadditivity.

A10. Closing pairing: for y not in A, sum_{x in U_y} of the trace of F_j tensor
L_{d_{Z_1}} tensor L_chi(x-y) is bounded by a d_j p^{(j+1)/2}. PROVED given A3
(geometric irreducibility, rank d_j >= 2 => H^0_c = H^2_c = 0), tame Euler
characteristic on P^1 minus a+2 points (CITED by the note from Katz GKM; standard
Grothendieck-Ogg-Shafarevich) and Weil II 3.3.1 for the weights of H^1_c
(UNVERIFIED CITATION locally, as in A7). The E_j part costs (c_j - d_j)
p^{(j+1)/2} by the trivial pointwise bound; the <= a-1 anchor points x in A \
Z_1 cost (a-1) c_j p^{j/2} by A9; x = y is killed by chi(0) = 0.

A11. Exceptional y in A \ Z_k and the final constant: |I(y)| <= p^{(j+2)/2} by
operator norms; (2a-1)(2a+2)^{k-2} + (a-1)/sqrt p <= 3a(2a+2)^{k-2}. PROVED
(elementary; re-checked).

A12. The claim "this proves Kunisky's Conjecture 1.14 at every fixed degree":
correct as an implication (Section 1.4, item 3) conditional on (T). The claim
that the argument is uniform in the anchor positions: correct, nothing in A1-A11
depends on the positions beyond distinctness.

Summary for (a): PLAUSIBLE BUT UNREVIEWED GAPS. Gaps: (G1) arithmetic
Frobenius structure of the constant kernels V_i and of the punctual kernels,
with weights <= i-1, asserted via Weil II 1.8.4 without argument; (G2) the BBD /
Weil II weight statements could not be re-verified from a local text in this
pass; (G3) the note never states Kunisky's Theorem 1.18 (a = 1 already proved
with the better constant k) or Theorem B.1, and presents a = 1 as a corollary.
No mathematical error found; every combinatorial and elementary step
re-derived and machine-checked.

## 5. Line-by-line classification of the argument of (b)

B1. Word recovery (1) and non-isomorphism (2). PROVED: MC_chi is an equivalence
(Katz 2.9.7), so F_w determines G = ME(F_prefix tensor L_{d_T}) up to
isomorphism; t_v(F_prefix) >= delta >= 1 by A6 and t_v(F_0) = 1, so twisting
makes t_v(G) < 0 exactly at v in T; rank F_w = 1 iff w is empty (d_w >= |w|+1)
fixes the length. Verified by exhaustive backwards recovery and distinct
signatures for every word to depth 14/7/5/4 (a = 1/2/3/4).

B2. Self-duality (3). PROVED given the CITED facts that duality exchanges *_! and
*_* (Katz 2.6, "D(L*_!K) = D(L)*_*D(K)"), commutes with middle extension (2.3.3
"it commutes with duality"), and D(L_chi[1]) = L_chi[1](1) for quadratic chi:
induction gives D(F_w[1]) = F_w[1](|w|+1), i.e. F_w^vee = F_w(|w|), hence real
traces on U. Numerically: the raw Gram diagonals p^{-1} sum_U q_w^2 equal 1 to
within 3*10^{-4} at p = 100049 for all 343 words of length 3 at a = 3, and to
within 2*10^{-5} at p = 1,000,033 for the 49 words of length 2 at a = 3 (Part E),
exactly the prediction of irreducibility plus (3).

B3. Gram bound (4). PROVED given B1-B2 and the weight input of A10: off-diagonal
H^2_c vanishes (non-isomorphic irreducibles), dim H^1_c = a d_u d_v; diagonal
H^2_c is one-dimensional with eigenvalue p^{|w|+1}; entrywise |G_uv - delta_uv|
<= (a d_u d_v + 1_{u=v})/sqrt p; operator norm by domination by the rank-one
matrix a d d^T/sqrt p plus 1/sqrt p. Numerically (Part E): in the four
non-vacuous cases (a = 1, lengths 1..6, and a = 2, m = 2, at p = 100049 and
1,000,033) the actual ||G - I||_op is 0.0049, 0.0103, 0.0033, 0.0045 against
claimed bounds 0.443, 0.819, 0.140, 0.257; ||G - I||_op sqrt p stays in
[1.5, 4.5] while the claimed constant a D_W + 1 is 140 or 247, so (10) is valid
but loose by a factor ~50 there. Max off-diagonal times sqrt p is between 1.5
and 5.8 for every configuration up to 343 words (a = 3, m = 3) and the
entrywise ratio |G_uv - delta_uv| sqrt p /(a d_u d_v + delta_uv) never exceeds
0.20; distinct words are numerically uncorrelated at the predicted rate.

B4. Cell bound (5). PROVED given B3 (scalar quadratic inertia at fresh points
kills H^2_c on the diagonal; a+s+2 punctures). Not re-tested numerically here;
the pass-10 verifier has exact certificates.

B5. Raw-rank recurrence (6). PROVED here by an independent Euler-characteristic
computation: for a tame perverse K with finite singular set S and generic t,
r(K *_! L_chi[1]) = h^1_c(A^1, K tensor L_chi(t-x)) = -chi_c = sum_{s in S}
c_s(K) (Grothendieck-Ogg-Shafarevich on P^1 minus S u {t, inf}, plus the stalks
at S, minus nothing at t where L_chi(t-x) is extended by zero; h^0_c = h^2_c = 0
generically), and at v in S the same computation with t = v gives c_v(K *_! L)
= r' - sum_{s != v} c_s(K) = c_v(K). Extension by zero at T sets c_v = r for v in
T. Additivity extends this from middle extensions and punctual objects to all
raw objects. The note's phrase "the constant infinity kernel restores the
subtracted infinity rank" is a correct description of the same fact.

B6. Word sums (7) and growth (8). PROVED: the per-word identity r' = r +
sum_{v in T}(r - c_v) summed over nonempty T gives the displayed linear
recurrences; my enumeration reproduces (W,R,Q,Z) at every depth checked (a = 2:
Q = 1, 17, 199, 2001, 18687, 167801, ...), and Lambda_a equals the Perron root of
[[A_2, h/2],[a h, h-1]] (Lambda_2 = 8.5311, Lambda_3 = 33.7797, Lambda_4 =
108.264); Q_{m+1}/Q_m decreases toward Lambda_a (8.98 at m = 5 for a = 2).

B7. (9)-(10). PROVED given A7 (weight gap) and B5 (generic rank of E_w is r_w -
d_w): triangle inequality in L^2(U). Inherits gaps G1-G2.

B8. (11)-(12). PLAUSIBLE, main terms re-derived (expansion of the cell indicator
over subsets B' of B gives 2^{-s} sum_{B'}(a + |B'|) D_W/sqrt p <= (a + s/2)
D_W/sqrt p, plus beta eps_W, plus the removal of the s boundary points where
|h_w| <= d_w, giving s D_W/(2p)); the exact bookkeeping of the beta terms was not
re-checked line by line.

B9. Dimension obstruction (13). PROVED (linear algebra); the pass-10 verifier's
integer null vector at p = 13, A = {0,1}, y = 2, m = 3 is a valid witness. This
is a correct and useful negative statement about the approach: unrestricted
coefficient isometry cannot reach m/log p -> infinity.

Summary for (b): PLAUSIBLE BUT UNREVIEWED GAPS, inheriting exactly G1-G2 of (a)
through (9)-(10); B1, B5, B6, B9 are PROVED here independently; the numerical
Gram data give strong independent evidence for B1-B3.

## 6. Verdicts

(a) parallel9, inequality (T) and "proves Conjecture 1.14 at every fixed degree":
PLAUSIBLE BUT UNREVIEWED GAPS (G1, G2, G3 of Section 4). Not REFUTED: no
integer witness in ~10^5 exact necklace and inner-sum checks, including all
non-vacuous (a,k) reachable below p = 2*10^6. Not PROVED under the brief's
standard: the weight inputs are cited, not re-verified from a local text, and
the arithmetic-versus-geometric bookkeeping is asserted. Conditional
statement I am willing to sign: **CONDITIONAL on (i) Katz RLS 2.9.4(3), 2.9.7,
3.3.3, 3.3.6, 3.3.7 (quoted above), (ii) Deligne Weil II 1.8.4 and 3.3.1, (iii) BBD
5.1.8, 5.1.14, 5.3.1, (iv) the compatibility of middle convolution with base
change from F_p to its algebraic closure, inequality (T) holds for all a >= 1,
k >= 3, p = 1 mod 4 (p > a), and Kunisky's Conjecture 1.14 holds at every degree
a.** Independently PROVED here: the combinatorial core (R), (R'), the
pseudoreflection structure, the mass bound (C), and the elementary closing
arithmetic.

(b) parallel10: PLAUSIBLE BUT UNREVIEWED GAPS (same G1-G2 via (9)-(10)); (1), (2),
(6), (7), (8), (13) PROVED; (3), (4), (5) PROVED conditional on the same inputs.

Neither note claims, and neither implies, anything about the signed spectral
aggregate, the extreme eigenvalues, the clique number beyond Kunisky's Theorem
1.17 consequence, or the two-set conjecture; those remain OPEN.

## 7. Remaining obligations

1. Quote Weil II 1.8.4, 3.3.1 and BBD 5.1.8, 5.1.14, 5.3.1 from a text source
   (or OCR the local scan) and write the two-line arithmetic base-change
   argument for V_i and the punctual kernels (closes G1-G2).
2. State explicitly in (a) that a = 1 is Kunisky's Theorem 1.18 with the better
   constant k, and that the constant two-anchor word is his Theorem B.1 (G3).
3. If the sharper Remark-1.16 form (O(p^{k/2}) for non-degenerate words) is
   wanted, one needs the weight analysis of H^k_c of the full k-variable
   arrangement or cancellation in y of the rank-(2a+2)^{k-2} trace function
   y -> I(y); nothing in (a) or (b) addresses it.
4. For (b): the useful depth m <= (1/2 - eps) log p / log Lambda_a and the
   obstruction (13) together show that any spectral use must exploit the
   specific coefficients of the cyclic trace; this is a research problem, not
   a gap in (b).

## 8. Verification

experiments/sigma_referee_2026_09_05.py (standard library + numpy only; no
import of the pass-9/10 code; independent implementations of Katz's rules, of
exact FFT convolution by 15-bit digit planes with asserted rounding residual,
and of the dense exact traces) writes results/sigma_referee_2026_09_05.json.
Full run (default mode, 308.6 s on this machine): 1,417,148 checks, 0
witnesses. Counts: Part A, 18 primes 101..3389, 11 patterns: 1,188 exact
comparisons with (T) and 1,188 with the trivial bound, 198 k = 1 zeros, 198
k = 2 bounds, 36 float/int64/object cross-checks, 3 literal Definition-1.13
sums; Part D, 77,140 middle-convolution transitions (a = 1..4 to depth
14/7/5/4) each checked for rank growth, formula (R), block-count formulas
(282,353 finite and 77,140 infinity identities), t_v >= delta, pseudoreflection,
involution, distinct signature, word recovery and raw rank, plus 19 levels of
recurrence (7) and 3 Perron identities; Part C, 88 pointwise expansions, 69
clique cases, 30 exact trace expansions, 5 checks of S^2 = pI - J; Part B, 22
FFT-vs-dense identities, 528 generic-y and 270 anchor-y inner-sum bounds at
p = 10009, 100049, 1,000,033; Part E, 12 Gram configurations (4 non-vacuous).
A FAST mode (SIGMA_FAST=1, 55 s, 124,640 checks) reproduces every conclusion
at smaller primes. Part F of the JSON summarizes the previous worker's
explore2/3/4 data (slopes, maximal ratios, absence of violations) read from the
scratchpad when present; the JSON also records the a = 1 maxima
|N|/(3*4^{k-2} p^{(k+1)/2}) = 0.166, 0.033, 0.019, 0.004, 0.001, 0.0002 for
k = 3..8 and |N|/(k p^{(k+1)/2}) <= 0.74 against Kunisky's own constant.
