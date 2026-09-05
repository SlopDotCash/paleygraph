# Prescribed inverse signs: joint moments and the Möbius cocycle

**Status: no uniform sixth-moment bound or useful contraction was
proved.** This pass obtains exact identities for the two actual sign
halves, a quantitative joint average, and a lower-order bound for the
asymmetry of their averages. Their common sixth-moment contribution
retains the original unknown. Repeated Möbius transformations preserve
that contribution through an exact weight cocycle. These are statements
about the actual quadratic character, with no synthetic sign model.

The [pass-24 inverse theorem](parallel24-inversion-orbit-upper-2026-09-05.md)
and [review](parallel24-inversion-independent-review-2026-09-05.md)
give a good unsigned B_h representative of every input. The gap remains
the prescribed weights χ(d), not the unsigned moment of that selected
representative. This note does not repeat the balanced-split estimate
from the [pass-24 local investigation](parallel24-classical-exceptions-2026-09-05.md).

Let p be an odd prime, χ(0)=0, ε=χ(−1), C⊂F_p, n=|C|,
F(x)=Σ_(c∈C)χ(x−c), M_j=Σ_xF(x)^j for even j, and

    Q(x,z)=Σ_(c∈C)χ(x−c)χ(z−c).

The quantitative critical-scale bounds below assume n≥2 and p≥n⁴.
The exact identities themselves do not require that restriction.

## 1. The missing rows cancel separately for the two sign halves

For z∉C, put

    D_z^±={1/(c−z):c∈C, χ(c−z)=±1}.

These are the actual positive and negative coefficient supports after
inversion; since χ(c−z)=χ(1/(c−z)), they are exactly the square and
nonsquare parts of D_z. For x≠z and t=1/(x−z), multiplicativity gives

    F_(D_z^±)(t)=χ(x−z)[Q(x,z)±εF(x)]/2.             (1)

Write m_±=|D_z^±|. At x=z, the expression in brackets divided by two
equals (n±εF(z))/2=m_±. The actual remaining row t=0 has value
±εm_±, which has the same even sixth power. Therefore there is no
leftover boundary error in either identity:

    M_6(D_z^±)=1/64 Σ_x[Q(x,z)±εF(x)]^6.             (2)

In particular, at every permitted pole,

    32[M_6(D_z^+)+M_6(D_z^−)]
      =M_6+15Σ_xF(x)^4Q(x,z)^2
           +15Σ_xF(x)^2Q(x,z)^4+Σ_xQ(x,z)^6
      ≥M_6.                                         (3)

The common term in (3) is exact and has positive sign. Estimating the
two sign halves jointly cannot discard it.

## 2. Complete averages and arithmetic cancellation of their difference

Define formal completed-pole quantities for every z∈F_p by

    J_±(z)=1/64 Σ_x[Q(x,z)±εF(x)]^6,
    T_j=Σ_(x,z) F(x)^(6−j) Q(x,z)^j,  0≤j≤6,
    A_4=Σ_(c∈C)F(c)^4.

Only when z∉C is J_± asserted to equal an inverse-subset moment.
Direct expansion and the complete quadratic-character correlation give

    T_0=pM_6,       T_1=0,
    T_2=p(nM_4−A_4)−M_6,                             (4)
    Σ_z[J_+(z)+J_−(z)]
       =[pM_6+15T_2+15T_4+T_6]/32,                  (5)
    Σ_z[J_+(z)−J_−(z)]
       =ε[20T_3+6T_5]/32.                           (6)

For (4), fix x and use

    Σ_zQ(x,z)^2=p[n−1_(x∈C)]−F(x)^2.

The odd terms in (6) admit a bound independent of M_6. For a tuple
u=(c_1,…,c_6), write K(u)=Σ_xχ(Π_i(x−c_i)), and use the same
notation for shorter tuples. Expansion before either field sum gives

    T_j=Σ_(u∈C^6) K(u)K(c_(7−j),…,c_6),  j≥1.       (7)

There are at most 15n³ full tuples of even multiplicities. Their full
sum has absolute value at most p; all other full tuples have absolute
value at most 5√p. Odd-length three- and five-tuples are nonsquares,
with absolute sums at most 2√p and 4√p respectively. These are the
ordinary individual Weil bounds, including the repeated-root correction
already stated in the pass-24 proof. Thus

    |T_3|≤30p√p n³+10pn⁶,
    |T_5|≤60p√p n³+20pn⁶.                            (8)

Equations (6) and (8) prove the completed-pole estimate

    |Σ_z[J_+(z)−J_−(z)]|≤30p√p n³+10pn⁶.            (9)

To return to actual poles, use

    J_+(z)+J_−(z)≤[M_6+Σ_xQ(x,z)^6]/2.

This is the scalar sixth-power inequality
(a+b)^6+(a−b)^6≤32(a^6+b^6), summed and divided by 64.
Since |Q|≤n and Σ_xQ(x,z)^2≤pn, deleting the n forbidden poles
changes the difference in (9) by at most (nM_6+pn⁶)/2. Consequently

    |Σ_(z∉C)[M_6(D_z^+)−M_6(D_z^−)]|
       ≤30p√p n³+(21/2)pn⁶+(n/2)M_6.                (10)

On n=floor(p^(1/4)), the ordinary bound
M_6≤15pn³+5√p n⁶ makes (10) equal to o(p²n³), uniformly in C.
Thus the **averages** of the two actual halves differ by o(pn³).
This is arithmetic cancellation in their difference. It gives no
estimate for their common contribution in (3) or (5), and does not say
that the halves have comparable moments at every individual pole.

## 3. A joint upper estimate, with its unresolved term displayed

Let

    S=Σ_(z∉C)[M_6(D_z^+)+M_6(D_z^−)].

For p≥n⁴ the following finite bounds hold:

    (p−n)M_6 ≤32S
      ≤(p−n−15)M_6+195p²n³+25pn⁶.                  (11)

The lower bound is (3). For the upper bound, sum (3) over permitted
poles and complete only the three nonnegative mixed terms. Formula
(4) accounts for the −15M_6 term. The standard fourth-moment bound
and its version with the fixed row weights χ(x−c) give

    M_4≤3pn²+3√p n⁴≤6pn²,
    T_4≤(3pn²+3√p n⁴)M_2≤6p²n³.

The squared-complete-sum estimate from pass 24 gives

    T_6≤15p²n³+25pn⁶.

Substitution, dropping the nonpositive −15pA_4 correction, proves
(11). These constants are convenient upper bounds, not claimed sharp.

This also preserves a positive fraction of good poles simultaneously.
Let H=p−Q_h(n)>0, with Q_h and the B_h pole bound as in pass 24,
and define

    U=15p²n³+25pn⁶,
    B=(p−n−15)M_6+195p²n³+25pn⁶.

For any α,β>0 with α+β<1, at least (1−α−β)H poles are B_h and
satisfy both

    M_6(D_z)≤U/(αH),
    M_6(D_z^+)+M_6(D_z^−)≤B/(32βH).                 (12)

This follows by two Markov bounds over the good poles; it does not
assume independence. For example α=β=1/3 leaves at least H/3 poles.
If near balance |F(z)|≤τn is also required, discard at most
(p−n)/(τ²n) further poles, by the exact second moment M_2=pn−n².
At h=2, n=floor(p^(1/4)), and τ=n^(−1/4), this extra loss is
o(p), whereas H∼3p/4; moreover m_±=n/2+O(n^(3/4)).

The unknown M_6 remains explicitly in B. Combining the most favorable
unrestricted-pole average in (11) with M_6≤32S/(p−n) produces only

    M_6≤[1−15/(p−n)]M_6
              +[195p²n³+25pn⁶]/(p−n).               (13)

Although its leading coefficient is strictly below one, absorbing that
term divides by 15/(p−n), yielding merely

    M_6≤13p²n³+(5/3)pn⁶.                            (14)

This is weaker than the existing individual Weil estimate. The
contraction is of order 1/p, whereas its additive error is of order
pn³. Requiring the good B_h poles and using (12) also incurs the
factor (p−n−15)/(βH); at the r=3,h=2 boundary its limiting value is
1/[β(3/4)]>1. Thus the proved joint selection does not close the
prescribed-sign gap. Neither the symmetry of the averaged halves nor
near balance removes their common M_6 term.

## 4. Repeated Möbius transformations retain a single prescribed cocycle

Take a matrix g=(a b; c d) of nonzero determinant Δ, and suppose
cu+d≠0 for each u in the support C. Define

    g(u)=(au+b)/(cu+d),
    w'_(g(u))=w_uχ(cu+d).

For a row x where cx+d≠0, the exact factorization gives

    F_(g(C),w')(g(x))=χ(Δ)χ(cx+d)F_(C,w)(x).         (15)

When c≠0, let z=−d/c. The new missing row a/c has value
χ(Δ)χ(c)Σ_u w_u. Also |Σw'|=|F_(C,w)(z)|. Including these rows
proves the augmented-moment invariance

    M_6(g(C),w')+|Σw'|^6=M_6(C,w)+|Σw|^6.           (16)

For c=0 this is ordinary affine invariance. Under composition the
denominator factors multiply:

    [c_2g_1(u)+d_2][c_1u+d_1]
       =c_(g_2g_1)u+d_(g_2g_1).

Their quadratic characters therefore form exactly the same cocycle;
they are not independent new signs.

In particular, start from all-positive weights and invert at z∉C,
then invert the result at u∉D_z. If u=0, the new set is C−z and
all weights return to +1. Its unsigned moment is the original M_6.
If u≠0, put v=z+1/u. The exact new coordinates and weights are

    1/[1/(c−z)−u]=−1/u−1/[u²(c−v)],
    χ(c−z)χ[1/(c−z)−u]=χ(−u)χ(c−v).                 (17)

Thus two inversions give an affine image of a single inverse at the
new pole v, with precisely its transported signs up to a global sign.
The same conclusion holds for further compositions by (15)–(16).
Choosing additional transforms gives no new independent sign variable
outside this Möbius family. This does not exclude a successful joint
arithmetic selection within that family; such a selection is unproved.

## 5. Verification and scope

The [new exact verifier](../experiments/parallel25_signed_inversion_2026_09_05.py)
and [results](../results/parallel25_signed_inversion_2026_09_05.json)
record 666 actual field/set cases. Coverage comprises every two- and
three-element set in F_5,F_7,F_11,F_13, and six explicit Sidon cases
at (p,n)=(17,2),(19,2),(83,3),(89,3),(257,4),(1297,6).
The last six satisfy the critical-scale inequality p≥n⁴.

- 15,140 separate actual sign-half moment identities;
- 7,570 joint lower bounds and exact sign-half size identities;
- 666 complete even/odd average identities and 3,330 mixed-sum bounds;
- 1,332 completed/actual average-asymmetry bounds;
- 5,964 two-inversion coordinate-and-weight cocycles;
- 28,914 weighted Möbius coordinates and 2,166 augmented invariants;
- two independent ordered-six-tuple expansions for all seven T_j;
- six critical-scale joint upper bounds and simultaneous good-pole checks.

All calculations use Python integers. Radical inequalities are tested
by retaining the sign of their rational part and then squaring, with
no floating acceptance threshold. The finite tests verify the formulas
and implementation; the preceding derivations supply their general
proofs. Input hashes are recorded in the result file.

The strongest new quantitative conclusion here is the lower-order
asymmetry of the two averaged sign-half moments, together with their
exact joint formula. Their unknown common component is not bounded at
the required scale. No uniform classical estimate, improved asymptotic
cancellation exponent, or full proof follows from this pass.
