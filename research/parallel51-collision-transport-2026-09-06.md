# Cross-level collision transport and a quartic obstruction

The Paley conjecture and Proximity Prize remain unproved. This note
tracks the actual primitive fibers between subgroup sizes. It gives
an exact compatibility equation and a quartic example where primitive
triple mass disappears at the next level. Neither supplies the remaining
uniform mass bound from [pass 50](parallel50-tower-mass-2026-09-06.md).

## 1. A sharper lower bound retaining the quadratic baseline

Fix dyadic N>=4 and an odd prime p=1 mod N. At each dyadic s>=4 let
c_s(lambda) be the root multiplicities of R_s modulo p, and write

    A_s = sum_lambda c_s(lambda)(c_s(lambda)-1),
    Y_s = sum_(c_s(lambda)>=3) c_s(lambda)(c_s(lambda)-2).

Here A_s is the ordered collision mass of primitive labels; it is
unrelated to the integer A_n in the earlier discriminant note. Put
A_tot=sum_s A_s and Y_tot=sum_s Y_s. In particular A_s>=Y_s.
Let E_s, B_s and T_s have their existing energy meanings, and define
the normalized excess over opposite-pair energy by

    D_s = [E_s-(3s^2-3s)]/s.

The existing primitive identity gives B_s=s^2/4+s A_s. Substituting
this in E_s=2E_(s/2)+6B_s+8T_s proves the exact recursion

    D_s = D_(s/2) + 6A_s + 8T_s/s.

Since E_2=6 and D_2=0, we obtain

    D_N = 6A_tot + 8 sum_s T_s/s,
    E_N >= 3N^2-3N + 6N Y_tot.                         (1)

This strengthens the baseline in pass 50's lower bound. It does not
improve a uniform energy exponent. It also shows exactly why old
energy excess cannot disappear when a new primitive level is simple.
The D_s are nondecreasing, while the Y_s need not be.

## 2. Transport every primitive label to the endpoint

The existing polynomial factorization iterates to

    P_N(Y) = product_(s=4,8,...,N)
               product_(R_s(lambda)=0) [Y-lambda^(N/s)]. (2)

All products retain labeled roots and therefore multiplicity. To check
(2) directly, partition the inverse-pair representatives of nontrivial
N-th roots other than -1 by their exact order s. For such a root h,
its old label is (1+h)^s, and its endpoint label is its (N/s)-th power.
Thus this is a reuse of the earlier primitive/kernel factorization,
not a new integer-polynomial construction.

Define the transported multiplicities

    b_s(beta) = sum_(lambda^(N/s)=beta) c_s(lambda),
    e(beta) = sum_s b_s(beta),  a = 2^N modulo p.

The [kernel identity](kernel-discriminant.md) is

    kappa(beta) = 2e(beta) + 1_(beta=a),
    D_N/4 = sum_beta e(beta)(e(beta)-1) + e(a).          (3)

To separate what can change across levels, put

    U_s = sum_beta b_s(beta)(b_s(beta)-1) - A_s >= 0,
    U = sum_s U_s,
    C = sum_(s<t) sum_beta b_s(beta)b_t(beta),
    H = e(a).

Here U counts additional ordered pairs merged by the power maps;
C counts pairs of roots from different primitive levels sharing an
endpoint label; and H counts roots meeting the distinguished label.
Expanding (3) gives

    D_N/4 = A_tot + U + 2C + H.                        (4)

Comparing (1) and (4) proves the exact compatibility condition

    U + 2C + H = A_tot/2 + 2 sum_s T_s/s.               (5)

In particular U+2C+H>=A_tot/2. A nonzero primitive collision mass
requires additional label merging, coincidences across levels, or a
distinguished-label collision. The locations and sizes of primitive
fibers cannot be assigned independently of their endpoint images.

This is a necessary condition for actual field towers. It does not
exclude every abstract mass profile in pass 50: those profiles did not
specify the transported labels, U, C, H or the T_s. No realization or
nonrealization theorem for them follows from (5) alone. The desired
uniform estimate would require further quantitative control.

## 3. A triple can disappear at the next primitive level in the target range

Take

    p = 2144280833,  N = 256,  g = 231737012.

Trial division through floor(sqrt(p)) certifies p prime; g has exact
order 256, checked by g^256=1 and g^128=-1 modulo p. Also

    1073741824 = N^4/4 < p < N^4 = 4294967296.

The actual primitive multiplicities are:

| Subgroup order s | Primitive multiplicities | Y_s | A_s | E_s | D_s |
|---:|:---|---:|---:|---:|---:|
| 64 | 16 singletons | 0 | 0 | 12096 | 0 |
| 128 | 27 singletons, one doubleton, one tripleton | 3 | 8 | 54912 | 48 |
| 256 | 64 singletons | 0 | 0 | 208128 | 48 |

All smaller levels have zero Y_s, A_s and T_s, and T_128=T_256=0.
Consequently Y_tot=3 while Y_256=0. This refutes the pointwise assertion
Y_tot<=C Y_N for every eligible N,p and any finite constant C, even in
the quartic range. It also refutes pointwise monotonicity of primitive
triple mass there. This finite example does not exclude an eventual
large-N assertion or a bound with an additive error. It does not
refute an energy bound with a quadratic baseline, nor any intermediate
energy, triangle, Paley or prize target.

For a direct triple witness, use h=g^2, K=<h^2> of order 64 and
L=hK. The number x=5438615 has exactly these three K x L representations:

    (1994579121,155140327),
    (1765701062,384018386),
    (559077563,1590641885).

All sums equal x modulo p, and direct membership and pair enumeration
verify completeness. At the next level the primitive R_256 is squarefree.
The old collision nevertheless survives in P_256 through the factor
coming from squared roots of R_128; it is not lost from the full energy.

For this endpoint the exact transport totals are

    A_tot=8,  U=2,  C=1,  H=0,  sum_s T_s/s=0.

Thus both sides of (5) equal 4, and (4) gives D_256=4(8+2+2)=48.
The quadratic baseline is 195840; the actual excess is 12288=256*48.
This explains the apparent disappearance: only the new primitive
polynomial is simple, while the endpoint kernel retains the old roots.

## 4. Verification and next task

The [standard-library checker](../experiments/parallel51_collision_transport.py)
checks the transported labels directly from a common endpoint generator,
the complete kernel identity, every energy recursion, and (1), (4)-(5).
It covers the previous 70 certified splitting cases and the new
order-256 endpoint, for 71 towers, 397 levels and 2,731,208 literal pair
enumerations. All checks passed. The witness prime receives its own
trial-division check. Full witness roots, transported multiplicities and
representations are retained in the
[results](../results/parallel51_collision_transport_2026_09_06.json).

The prime was already a certified triple-fiber prime at order 128 in
pass 48. What is checked here is its extension to an order-256 endpoint
in that endpoint's quartic interval; no new complete prime classification
is claimed. All quantities are exact integers or rational numbers.

The next task remains a uniform upper-tower mass estimate. One cannot
replace that interval by its largest primitive level using monotonicity
or constant-factor domination. The full transported kernel and its
compatibility condition retain the information that such a shortcut
would discard. Equation (5) itself gives no uniform upper bound.
Uniform energy49/20, triangle86/15, period71/72 and absolute exception5/3
remain unchanged. Root completed ordinary proofs and exact checks;
no separate-agent, Lean, external review or novelty claim is made.
No Lean process or Prove2Me submission was changed.
