# The ambient elliptic kernel has no fixed norm saving

The Paley conjecture and Proximity Prize remain unproved. Recent tower
passes did not improve the uniform energy exponent. This pass returns
to the direct two-anchor spectral route and rules out an ambient-norm
shortcut there, including on the inversion-odd subspace. The actual
compressed, symmetry-averaged operator remains to be bounded.

## 1. The operator that must be kept

For prime p=1 mod4, let chi be the quadratic character, Q the nonzero
quadratic residues, and

    C={x:chi(x)=chi(x-1)=1},
    L(t)=sum_y chi(y(y-1)(y-t)).

The [existing full-sector reduction](parallel25-spectral-full-operator-2026-09-05.md)
uses K_0=(L(y/x))_(x,y in C), its S_3 conjugacy average A, and

    S_C^2=pI/4+3A/4-3J/2.                              (1)

On the nontrivial S_3 sectors, the desired two-anchor edge is equivalent
to lambda_max(A)<= (2/3+o(1))p. On the trivial sector the recorded
rank-two border must also be retained. This is an intermediate spectral
target; no reduction from proving only this edge to the full two-set
Paley conjecture or the exact prize statements is asserted.

Enlarge K_0 to the genuine ambient operator

    K_Q=(L(y/x))_(x,y in Q).

Then K_0 is its coordinate compression. Bounding K_Q by 2p/3 would
be a convenient sufficient step before averaging. It is false, even
asymptotically. This does not make the smaller averaged bound false.

## 2. Exact Jacobi diagonalization of the ambient operator

All multiplicative characters below are extended by zero at zero,
including the trivial character. For a character psi of F_p^*, write

    J(psi,chi)=sum_(t!=0) psi(t)chi(1-t).

The existing squared-Jacobi convolution identity can be seen directly.
Put f(t)=chi(1-t) on F_p^*. Its multiplicative convolution satisfies

    (f*f)(t)=sum_(y!=0) chi(1-y)chi(1-t/y)=L(t),

because chi(-1)=1. Thus

    sum_(t!=0) L(t)psi(t)=J(psi,chi)^2.                 (2)

Also f(1/t)=chi(t)f(t), giving

    J(psi chi,chi)=J(conjugate(psi),chi)
                 =conjugate(J(psi,chi)).               (3)

The character alpha=psi|_Q is an eigenvector of K_Q with eigenvalue

    lambda_alpha=sum_(t in Q)L(t)psi(t)
                =[J(psi,chi)^2+J(psi chi,chi)^2]/2
                =Re(J(psi,chi)^2).                     (4)

The two extensions psi and psi chi give the same eigenvalue, as
required. The characters of Q form a complete orthogonal basis, so
(4) diagonalizes the ambient operator, not its C compression.

For completeness, the character magnitudes follow from a single
quadratic correlation. At u!=0,

    sum_(x!=0) f(x)f(ux)=p 1_(u=1)-1-chi(u).

Its Mellin transform gives |J(psi,chi)|^2=p for psi!=1,chi. At the two
exceptional characters J=-1, so the constant Q-mode has eigenvalue 1.
Hence ||K_Q||<=p. These identities are elementary and do not require
an assumed random distribution of phases.

Inversion I_Q:v(x)->v(1/x) commutes with K_Q. If alpha is not
self-inverse, alpha and its conjugate share the real eigenvalue (4),
and Im(alpha) is a nonzero real inversion-odd eigenvector. The only
extensions removed by requiring non-self-inverse alpha, apart from
1 and chi already excluded, are the two quartic characters satisfying
psi^2=chi. Inversion-odd here is an ambient condition; it is not the
full sign representation on C.

## 3. A sourced asymptotic obstruction, with an explicit threshold

The primary source [Lu–Zheng–Zheng, Theorem 1.4, equation (1.7)](https://arxiv.org/html/1305.3405v3)
applies with its m=1, k=1 and singleton A_1={chi}. Its delta term
then vanishes. It bounds by 2p^(-1/4) the discrepancy on the unit circle
of J(psi,chi)/sqrt(p), over all psi!=1,chi. This is a full-character
statement; it is not an estimate on an arbitrary smaller character set.
The theorem and these conventions were checked in both the primary
web text and the previously archived v3 HTML. This estimate is imported,
not proved by the finite checker.

Write J(psi,chi)/sqrt(p)=exp(2pi i t). For p>=2^20, take the arc
|t|<=h, where h=2p^(-1/4), interpreted modulo one. Its length is 2h,
so the discrepancy bound gives at least

    (p-3)[2h-2p^(-1/4)]=2(p-3)p^(-1/4)>2

characters in it. After removing the two quartic characters, at least
one non-self-inverse Q-character remains. Its real odd eigenvector
has, by (4),

    lambda/p=cos(4pi t)
             >=1-8pi^2 h^2=1-32pi^2 p^(-1/2).

Therefore

    lambda_max(K_Q restricted to inversion-odd vectors)
       >= p-32pi^2 sqrt(p),                            (5)
    lambda_max(K_Q restricted to inversion-odd vectors)/p -> 1.

The upper side of the limit follows from ||K_Q||<=p. Furthermore,
32pi^2/1024<1/3, so (5) is strictly greater than 2p/3 for every
prime p=1 mod4 with p>=2^20. This rules out the proposed ambient
2p/3 estimate for all sufficiently large eligible primes, rather than
merely exhibiting a small exception. The argument uses no conclusion
about the overlap of these eigenvectors with the support C.

## 4. An exact small inversion-odd eigenvector

At p=89, g=3 is a primitive root. Order Q as x_r=g^(2r), 0<=r<44,
and put

    v_r=(0,1,0,-1)_(r mod4).

Exact integer multiplication verifies

    K_Q v=73v,  I_Q v=-v,
    ||v||^2=22,  v^T K_Q v=1606,
    73>2*89/3.

This vector has ten nonzero coordinates outside C. It is not an
eigenvector certificate for the desired C-supported averaged operator,
and is not a counterexample to the two-anchor conjecture or Paley.
The large-prime statement (5) follows from the source theorem and
the ordinary argument above, not from extrapolating this example.

## 5. Exact checks and next action

The [standard-library checker](../experiments/parallel53_ambient_elliptic.py)
works with cyclic integer coefficient lists for f. It compares the
Legendre trace sum with f*f and verifies the exact autocorrelation
whose Mellin transform supplies the Jacobi magnitudes. Thus the checks
of (2)-(3) use no rounded complex characters. Across 16 prime fields,
2,048 convolution coefficients and 2,048 autocorrelation coefficients
agree exactly. For 12 smaller fields it additionally checks 3,198
entries of the original restricted identity (1), including its -3J/2
boundary term. The integer eigenvector and its support are checked
directly. The [results](../results/parallel53_ambient_elliptic_2026_09_06.json)
also record the rational threshold comparison using pi<22/7.

The useful saving must involve the joint effect of support on C and
the S_3 conjugacy average, or another argument on the actual restricted
operator. Full-character phase equidistribution supplies the ambient
obstruction; it does not bound character weights selected by an
adversarial C-supported vector. There is no new estimate for that
interaction in this pass. This is the next direct spectral task.

The subgroup energy and triangle routes remain available with their
previous unresolved inputs. None of their uniform exponents changes.
The full two-set conjecture, subgroup square-root target and prize
reduction all remain open here. Root completed ordinary derivations
and exact checks; no independent-agent, Lean, external review or novelty
claim is made. No Lean process or Prove2Me submission was changed.
