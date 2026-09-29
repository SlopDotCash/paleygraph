# A smaller exceptional set for the weighted-triangle estimate

For each fixed dyadic n>=4, at most O_c(n^(5/3)/log n) eligible primes
in a quartic interval [c n^4,C n^4] can fail the intermediate estimate
W<<n^(17/3). This improves the pass42 absolute exception bound with
exponent 63/31. It does not bound every prime, establish a proportion of
eligible primes, or prove the full Paley conjecture or Proximity Prize.
No novelty claim is made: the argument combines a classical single-coset
incidence bound with the existing cyclotomic fourth-energy prime average.

## 1. Use additive energy as the sufficient input

For H=-H of order n with n^2<p, let

    f(x)=r_(H-H)(x), E= sum_x f(x)^2, E_*=E-n^2,
    W=sum_(x,y!=0,x!=y) f(x)^2 f(y)^2 f(x-y)^2.

On G=F_p^*/H write a(C)=f(C), b=a^2 and

    A2=sum_C a(C)^2=E_*/n,
    T(u,v)=sum_C b(C)b(uC)b(vC),
    rho(u,v)=#{x in uH:x+1 in vH}, R_max=max_(u,v) rho(u,v).

The exact identities already checked in pass41 give

    W=n sum_(u,v) rho(u,v)T(u,v), sum_(u,v) T(u,v)=A2^3.

Every term is nonnegative. Consequently

    W <= n R_max A2^3 = R_max E_*^3/n^2.                    (1)

The classical single-coset intersection estimate is

    R_max << n^(2/3),                                      (2)

so

    W << E_*^3/n^(4/3).                                    (3)

In particular E<<n^(7/3) implies W<<n^(17/3), with no logarithmic
loss. The imported uniform energy bound E<<n^(49/20) up to logarithms
does not reach this condition: 49/20>7/3. Inserting that larger energy
in (3) gives power 361/60, which is worse than pass41's uniform triangle
power 86/15. Thus this argument improves the exception count, not the
estimate for every prime.

### Exact hypothesis check for (2)

Root visually inspected the published
[Shkredov, Some new inequalities in additive combinatorics, Lemma 2,
equation (18), printed page 197](http://mjcnt.phystech.edu/en/download.php?id=66)
in the pinned PDF used by pass38. Set Gamma=H, Q=H, Q1=uH, Q2=vH.
The cardinality hypothesis is valid here:

    (Q1 x Q2) Delta^(-1)(Q) = Q1 x Q2,
    |(Q1 x Q2) Delta^(-1)(Q)|=n^2=|Q1||Q2||Q|/|H|.

The two size conditions reduce to n^3<<n^5 and n^4<<p^3 and hold
in the stated range, with finitely many bounded n absorbed in the
constant. Its incidence sum is n*rho(u,v), since multiplication by
each element of Q preserves the two cosets. The right side is n^(5/3),
which gives (2) after division by n. This application uses three single
cosets. It does not reuse the invalid union-of-cosets application audited
in pass38. The checker independently verifies this cardinality and
normalization on an actual rich cell of the new quartic witness.

## 2. Apply the fourth-energy prime average

For n dyadic, the complex intrinsic zero-quadruple count is

    T4=3n^2-3n.

The existing [cyclotomic prime-average proof](cyclotomic-prime-average.md),
including its AGM refinement, gives

    sum_(p prime,p=1 mod n) (E(H_p)-T4) log p <= C4(n),
    C4(n)=[n^4-T4]/2 * log[4n^4/(n^4-T4)].                   (4)

Every summand is nonnegative; fixed n has finite support. This is the
ordinary additive energy E, with the zero difference retained. It is not
the third difference energy or the shifted multiplicative excess X.

For any L>1 and threshold U>0, (4) implies

    #{p>=L:p=1 mod n,E(H_p)>T4+U} <= C4(n)/(U log L).        (5)

Set U=n^(7/3) and L=c n^4. For n sufficiently large depending on c,
the interval lies in n^2<p and the right side of (5) is

    O_c(n^(5/3)/log n).                                    (6)

Every remaining eligible prime has E<=T4+n^(7/3)<=4n^(7/3), so (3)
gives W<<n^(17/3). Thus the set that can fail an appropriate fixed
constant in the W bound is contained in the exceptional set of (5).
The two exception definitions in passes42 and43 need not be nested;
(6) improves the bound on their counts for the same triangle target.
The exponent improvement is 63/31-5/3=34/93.

More generally, a threshold U=n^(2+theta), theta>=0, gives
E=O(n^(2+theta)) and W=O(n^(14/3+3theta)) outside an absolute set
of O_c(n^(2-theta)/log n) primes. This is a family of consequences of
the same budget, not a theorem about density in an arithmetic progression.

## 3. What remains pointwise

A particular quartic prime can have E of size n^(49/20) without
violating (4). Its single budget contribution is far smaller than O(n^4).
The prime average cannot exclude that event. The full arbitrary-set Paley
bound and the official prize also require further estimates and bridges;
the weighted-triangle target alone is not their completion criterion.

The [exact checks](../results/parallel43_pointwise_inputs_2026_09_06.json)
verify (1) in four fields, the single-coset incidence normalization, and
the rational exponent arithmetic. They do not prove (2) by enumeration;
that is the explicitly scoped published dependency. Root reviewed the
ordinary derivation and the already existing proof of (4). No separate
agent, Lean, or external peer review was completed for this pass.
