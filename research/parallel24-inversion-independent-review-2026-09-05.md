# Separate-agent review of the unsigned inversion existence theorem

Verdict: **no mathematical gap found** in the stated finite theorem or
its fixed-r,h asymptotic consequence. The inversion identity works for
both congruence classes of odd primes, with the t=0 term included.
The repeated-root character bound and good-pole fraction are valid.
The sixth-moment boundary coefficient is exactly 20.

This gives a bounded unsigned B_h representative in each inversion
family. It does not supply the uniform moment hypothesis on all B_h
sets, or the transported signed estimate for the original set. The
Paley/prize goal remains unproved. This is separate-agent mathematical
review, not independent human refereeing or formal verification.

## 1. Scope and pinned inputs

Read the complete proof and verifier, the recorded result, and the
relevant prior inversion/completion proof. The archived primary
[McDonald–Sahay–Wyman HTML, Lemma 2.1](https://arxiv.org/html/2210.03789v2)
was inspected directly. Its quadratic-character specialization is
applied to nonconstant monic squarefree polynomials here. No fresh
network retrieval, PDF inspection, or review of the proof of Weil's
theorem is claimed in this lane. The sigma note is a pinned transitive
input; the new note supplies its own complete cardinality argument.

Final SHA-256 values, checked after the parent's verifier refresh:

| Input | SHA-256 |
| --- | --- |
| `research/parallel24-inversion-orbit-upper-2026-09-05.md` | `e011ead2f6aa9c73b39458eefad7cfc839803c5b0e62eb77993a5da9567eb552` |
| `experiments/parallel24_inversion_orbit_upper_2026_09_05.py` | `46a26f0b7f8217338e6d355425153fd619f0387a06d8b1fa215800c7b362e5a4` |
| `results/parallel24_inversion_orbit_upper_2026_09_05.json` | `928152d51aeedae125493ccd032e0304dc847032fcb6e231062f65dd3d07f1b6` |
| `research/parallel20-inversion-moments-2026-09-05.md` | `31ec520dcb38decced35eb111128d65383f4a4f6fb2eca545d8bc5406b55922f` |
| `research/sigma-sumproduct-2026-09-05.md` | `1a2c947e55de80867fb76492b6a4aafe030185a611dc52443e025c8597176a4d` |
| `sources/mcdonald-sahay-wyman-2210.03789v2.html` | `3ae1bc2ca60ef412b271e3a03056df70fe1465fda26a89727e1929cfe9e1f622` |

All recorded input hashes match these current files. The refreshed
result includes 573 successful fixed-width accumulation guards.

## 2. Inversion, signs, and the missing row

For z∉C and x≠z, write t=(x−z)^(−1), d=(c−z)^(−1). Then

    χ(t−d)=χ(c−x)χ(x−z)χ(c−z)
           =χ(x−z)χ(x−c)χ(z−c).

There are two factors χ(−1), so they cancel regardless of whether
p≡1 or 3 mod4. Summing gives the stated unsigned coordinate identity.
The map x↦t is a bijection from F_p\{z} to F_p*.

At the omitted row,

    χ(−d)=χ(−1)χ(c−z)=χ(z−c),

so F_(D_z)(0)=F_C(z), while Q_C(z,z)=k. The latter uses z∉C;
it would not equal k for a pole in C. Thus the exact even moment is

    M_(2r)(D_z)=Σ_x Q_C(x,z)^(2r)−k^(2r)+F_C(z)^(2r).

Since |F_C(z)|≤k, the column upper bound follows. The subsequent
two-row sum is allowed to include z∈C because it only enlarges a
nonnegative sum; it does not apply this boundary identity there.

For the transported signs w_d=χ(c−z)=χ(d), the coordinate identity
instead has a remaining factor χ(−1):

    F_(D_z,w)(t)=χ(−1)χ(x−z)F_C(x),
    F_(D_z,w)(0)=χ(−1)k.

Even powers remove that factor and give exactly equation (8) of the
note. No sign or boundary correction has been lost.

## 3. Complete sums and repeated roots

Expand Q_C(x,z)^(2r) as a sum over ordered tuples from C. The x and z
sums of each product are identical real integers, so their product is
the square of the same complete character sum. This proves the exact
identity (4), with no absolute-value inequality at this step.

When every multiplicity is even, the polynomial is a monic square,
and its character is one off its distinct roots and zero on them.
Its sum has absolute value at most p. Every such tuple is obtained
from at least one of A_r pairings of the positions, followed by k
choices for each pair. Thus A_r k^r is a valid upper count; the
possible multiple representations strengthen the upper bound.

Otherwise the number s of odd-multiplicity roots is even and at least
two, since the total degree is 2r. Reduce odd multiplicities to one
and remove roots of positive even multiplicity. At an odd root both
the original and squarefree characters vanish. They can differ only
at the e removed roots, each by at most one. The monic squarefree
polynomial of degree s is not a square. The primary lemma therefore
gives

    |complete sum|≤(s−1)√p+e
                  ≤(s+e−1)√p≤(2r−1)√p.

There is no requirement p>2r in this reduction. The character has
order two, and the field is odd; monicity also avoids any ambiguity
about a constant multiple of a perfect power. Squaring and adding
the even and non-even tuple bounds gives the displayed U_r. Using
k^(2r) as an upper bound for the non-even tuple count is deliberately
coarse but valid.

## 4. Good poles and the positive-fraction quantifier

Cancel the common entries of two different h-element multisets. Their
reciprocal-sum difference is Σ_c ν_c/(c−z), where all surviving ν_c
are nonzero integers of absolute value at most h and their sum is zero.
The hypothesis p>h keeps each ν_c nonzero in F_p. Clearing the
distinct denominators produces a nonzero numerator: evaluation at
any one of the poles c₀ gives
ν_(c₀)Π_(c≠c₀)(c−c₀)≠0. Its top-degree coefficient vanishes because
Σν_c=0, so its degree is at most 2h−2.

There are binom(N_h(k),2) unordered pairs of h-multisets. The union
bound, together with the k excluded poles in C, proves
|G|≥p−Q_h(k). Distinct pairs can have the same polynomial or share
roots; neither possibility invalidates this upper bound on bad poles.

For p>Q_h(k), let H=p−Q_h(k)>0. The sum of nonnegative moments over
G is at most U_r. The minimum is at most U_r/|G|≤U_r/H.
For any real λ>1, if B good poles exceed λU_r/H, then
B≤H/λ. Therefore at least

    |G|−H/λ≥(1−1/λ)H

good poles satisfy the claimed bound. This argument does not assume
independence of the additive-relation and character-moment properties.
It counts poles, including different poles that may produce the same
inverse set. The statement does not require distinct representatives.

## 5. The fixed-order slice and coefficient 20

Keep r≥3 and 2≤h≤floor((r+1)/2) fixed before taking the prime limit.
For k=floor(p^(1/(r+1))), one has k→∞ and k^(r+1)/p→1. Since
N_h(k)∼k^h/h!, the leading term of Q_h(k) is

    [(h−1)/(h!)²] k^(2h).

It follows that Q_h(k)/p tends to zero in the strict interior and
to b_h=(h−1)/(h!)² at 2h=r+1. To check the bound for all h≥2,
b₂=1/4 and

    b_(h+1)/b_h=h/[(h−1)(h+1)²]<1.

Thus b_h≤1/4, and p>Q_h(k) eventually. Moreover k^r/p→0. The
normalized finite bound is exactly

    U_r/[(p−Q_h(k))p k^r]
      =[A_r+(2r−1)²k^r/p]/[1−Q_h(k)/p].

This proves the stated asymptotic and makes its uniformity in C clear:
the displayed expression depends only on p,k,r,h. At r=3,h=2,
A₃=15 and b₂=1/4, giving 15/(3/4)=20. Nothing here permits taking
r or h to grow with p without further error analysis.

## 6. Independent small computations and implementation coverage

Separate plain-Python code, without importing the lane verifier or
using NumPy, checked C={0,2,7} in F_41 and F_43. These cover both
values of χ(−1). For both fields all permitted poles and all inversion
coordinates were checked, including the unsigned and signed boundaries.
All 3^6 ordered tuples were independently expanded. The checks also
tested the squarefree correction at every repeated-root tuple, B₂
status, the finite existence bound, and the fraction for λ=3/2.

| p | χ(−1) | two-row sixth sum, equal to tuple expansion | good B₂ poles | smallest unsigned M₆ on those poles |
| ---: | ---: | ---: | ---: | ---: |
| 41 | 1 | 272346 | 35 | 6054 |
| 43 | −1 | 301878 | 37 | 6720 |

In each field Q₂(3)=33, there are 183 even-multiplicity tuples, and
450 nonsquare tuples have at least one removed even-multiplicity root.
All good poles meet the tested λ=3/2 threshold. At z=1 the unsigned
and transported signed moments differ: (6654,7382) for p=41 and
(6720,7448) for p=43. This is a finite illustration of the distinction,
not an asymptotic obstruction to sign removal.

The parent run records 573 sets, 3,570,540 coordinates for each of
the unsigned and signed transformations, 27 tuple expansions, and
1,149 good-representative bounds. Those counts are not additional
independent runs by this reviewer. The number 573 agrees with the
specified exhaustive small-set slices plus ten larger cases. The
script was read and parsed; its full run was not repeated here.

The array arithmetic uses signed NumPy int64, while scalar bounds use
Python integers. For the recorded sizes its largest conservative
accumulation bound is p²k^6≤1297²·6^6=78485143104, below 2^63−1.
Thus no fixed-width overflow occurs in the recorded checks. During
this review the parent added an explicit guard for this bound and
reran the verifier; this is an implementation-scope clarification,
not a correction to a mathematical formula.

## 7. What remains missing

The unsigned representative estimate controls M_(2r)(D_z,1). The
original moment is transported to M_(2r)(D_z,w), where w is prescribed
by the pole and can have mixed signs. Decomposing its negative part
introduces the moment of a specific subset, which the representative
existence theorem does not bound. Likewise, the earlier completion
argument requires a bound on every B_h completion set of the stated
size, not one selected inverse from each family.

The note correctly identifies this unresolved step and does not claim
that sign removal is impossible. No original-set uniform bound or
Paley implication has been supplied by this existence theorem alone.
Only this new review file was written by this reviewing agent.
