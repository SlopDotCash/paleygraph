# Weighted matrix closure and the remaining triangle correlation

The full Paley and Proximity Prize goals remain unproved. The following
identities extend the last pass's matrix-square identity to arbitrary
weights. They retain rank corrections that are lost if the compressed
matrices are treated as commuting. The direct norm estimate derived here
is weaker than the current triangle bound; no uniform exponent improves.

Let H=-H have order n in F_p*, let m=(p-1)/n, and write quotient indices
additively. Throughout the projection argument assume n^2<p. Set

    C_ij=#{x in H_i:x+1 in H_j},
    (S_k C)_ij=C_(i-k,j-k),  L(u)=sum_k u_k S_k C.

All vectors below are real. This L(u) is the weighted shift operator,
not the differently named integer period matrix in the older note.
The [coset-convolution algebra](mixed-periods-and-shifted-energy.md) and
[shifted inner products](parallel44-interval-and-matrix-review-2026-09-06.md)
are the existing inputs. No new imported incidence estimate is used.

## 1. Full multiplication and rank corrections

For F_u=sum_i u_i 1_(H_i), its convolution operator on the invariant
subspace, in the orthonormal basis delta_0, n^(-1/2)1_(H_i), is

    M_u = [ 0          sqrt(n) u^T ]
          [ sqrt(n) u  L(u)        ].

The zero coefficient of F_u*F_v is n<u,v>, and its nonzero coset
coefficients are L(u)v. Commutativity of field convolution therefore gives

    L(u)v=L(v)u,
    M_u M_v=n<u,v>I+M_(L(u)v).

Comparing the lower blocks proves the full identity

    L(u)L(v)=n<u,v>I-n u v^T+L(L(u)v).                  (1)

Consequently

    [L(u),L(v)]=n(v u^T-u v^T),                         (2)
    L(u)^2=n||u||_2^2 I-n u u^T+L(L(u)u).               (3)

Only the augmented matrices M_u commute in general. Taking u=v=e_0
in (3) recovers C^2=nI-ne_0e_0^T+L(a), a=C e_0.
This derivation is a consequence of the already established convolution
algebra, not a claim that the algebra or association-scheme method is new.

## 2. Cubic and quartic trace identities

Write S=sum u, Q=||u||_2^2, c=L(u)u, tau=<u,c>, R=||c||_2^2.
The exact row and trace sums are

    L(u)1=nS1-u,    tr L(u)=(n-1)S,
    sum c=nS^2-Q.                                       (4)

Polarizing the preceding pass's Frobenius identity gives

    <L(u),L(v)>_F=n(n-1)(sum u)(sum v)+(p-2n)<u,v>.       (5)

Multiply (3) by L(u), take the trace, and use (4)-(5). The terms
containing SQ cancel, leaving

    tr L(u)^3=n^2(n-1)S^3+(p-3n)tau.                    (6)

Likewise expand the squared Frobenius norm of (3). In the cross term,
u^T L(c)u=c^T L(u)u=R. Substituting sum c=nS^2-Q gives

    tr L(u)^4=n^3(n-1)S^4+n(p-2n)Q^2+(p-4n)R.           (7)

Both identities hold for signed weights. Their validity does not require
the factors p-3n or p-4n to be positive. Direct field convolution also
gives

    sum_(x,y) F_u(x)F_u(y)F_u(-x-y)=n tau,
    sum_x (F_u*F_u)(x)^2=n^2 Q^2+n R.                    (8)

For u=b=a^2, the first expression in (8) is exactly the existing W.
Thus (6) expresses the same triangle quantity through a cubic trace,
including its principal contribution; it does not bound that trace.

## 3. Exact projection onto the shifted matrices

Put alpha=n(n-1), beta=p-2n>0, D=alpha m+beta=n(p-3)+1.
The shifted-matrix Gram matrix is alpha J+beta I. Moreover

    sum_k S_k C=nJ-I.                                   (9)

For i!=j, (9) counts ratios r in the specified nonidentity H-coset,
each determining x uniquely from (x+1)/x=r, giving n. On the diagonal,
omit r=1 and obtain n-1.

For the rank-one matrix U=u u^T, define

    q_k=<U,S_k C>_F=u^T(S_k C)u,  Z=sum q=nS^2-Q.

Its orthogonal projection P(U) onto the real span of {S_k C} is L(w),
where inversion of alpha J+beta I gives

    w_k=q_k/beta-alpha Z/(beta D).

Pythagoras, with q_bar=Z/m, yields the exact identity

    Q^2=Z^2/(mD)+||q-q_bar*1||_2^2/beta
                  +||U-P(U)||_F^2.                     (10)

Since tau=<u,q>, Cauchy-Schwarz now gives

    |tau-SZ/m|^2
      <=(Q-S^2/m) beta [Q^2-Z^2/(mD)].                  (11)

This accounts for the mean and the unused orthogonal component rather
than applying a full Frobenius norm immediately. It is still insufficient
at the current worst-case moment scales. For u=a^2, Q=A4<<n^(8/3) by
the existing weak cubic tail. In the quartic range, the square-root error
in n times (11) is bounded by O(n sqrt(p) A4^(3/2))=O(n^7).
The current uniform W bound has power86/15, and the intermediate target
has power17/3. Therefore this particular insertion of scalar upper bounds
is weaker than the current result. No impossibility claim about using
(1)-(10) more effectively is implied.

The remaining object here is the actual correlation of the centered
weight vector with q, which includes field incidences. Estimating it
requires more than applying (11) with those scalar moments. A generic
function on the quotient, such as the earlier interval example, does not
come equipped with matrices satisfying this full multiplication law.

## 4. Exact checks and scope

The [checker](../experiments/parallel45_weighted_matrix_algebra.py)
constructs actual C by enumerating field elements for (p,n)=(13,2),
(97,8),(257,16),(353,16). It uses e_0,e_1,a,a^2,e_0-e_1 as weights.
All25 ordered weight pairs per field satisfy (1)-(2), covering23,000
matrix entries. Independent matrix products verify40 cubic/quartic
trace identities. Twenty exact rational projections verify orthogonality,
Pythagoras, and (11). Twenty direct field convolutions independently check
the zero coefficient, all nonzero coefficients, and both identities (8).

The [certificate](../results/parallel45_weighted_matrix_algebra_2026_09_06.json)
records these finite checks. The general statements follow from the
ordinary derivations above. Root completed the derivations and checks;
no separate-agent, Lean, external peer review, or full proof is claimed.
