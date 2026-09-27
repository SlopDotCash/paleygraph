# Conditioning and coupling the scalar quotient

All claims concern the stated finite cyclotomic digit model. They are independent of the unresolved Paley norm estimates and proximity prize claims.

## Exact model and inherited erasure lattice

Let N>=4 be a power of two, p an odd prime, g of order2N modulo p, u=g^-1, and f a nonzero integral relation with f(u)=0 modulo p. In R=Z[X]/(X^N+1), let M be multiplication by f, m=(p-1)/2, and F_a[j]=center(a*g^j modulo p). The digit word is D=E(a)=M*F_a/p. Let f*A_f=k*p. The coefficient-zero inverse row w gives w.D=k*F_a[0], hence a=(w.D/k) modulo p.

For erased coordinates U and known coordinates K, write x=D_U and y=D_K. The congruence w_U.x=-w_K.y modulo k yields x=x0+z*L, using the row convention. The integer kernel basis L, its rational inverse T, its index and a Bezout lift are inherited from round21 and rechecked. The scalar step of row L_j is s_j=(w_U.L_j/k) modulo p. Only coordinates with s_j!=0 affect the recovered scalar.

Define q_j=M_U^T*T[:,j]. The old bounds on z_j have center -(x0*T)_j and radii B*||T[:,j]||1 and (m/p)*||q_j||1. Every candidate scalar is eventually re-encoded and compared with all known digits.

## Known digits yield shifted intervals

For every rational lambda on K,

    r_j = q_j - M_K^T*lambda,
    z_j = lambda.y - (x0*T)_j + r_j.F_a/p.

Thus z_j lies in the interval centered at lambda.y-(x0*T)_j with radius (m/p)*||r_j||1. Intersect this interval with the old interval and round the endpoints inward to integers. The centers are generally different. Keeping the old center while shrinking only its radius is unsound; the saved query regression exhibits an actual excluded word under that mistake.

Any rational lambda is a valid proposal. To certify the best possible radius within this formulation, additionally supply eta with

    M_K*eta=0,  ||eta||infinity<=1,
    q_j.eta=||q_j-M_K^T*lambda||1.

For any alternative lambda', q_j.eta=(q_j-M_K^T*lambda').eta<=||q_j-M_K^T*lambda'||1. Equality therefore proves optimality. The floating LP proposes a zero-residual subsystem; exact rational reconstruction supplies lambda. Only explicitly verified duals are called optimal.

These optimality claims concern a fixed functional and the centered cube with the known equations eliminated. They do not prove that a whole candidate box is optimal, that another lattice basis cannot improve it, or that query-specific endpoint optimization cannot improve an interval.

## A changed metric

For each erased unit vector e_i, form its adjugate inverse image A_f*e_i. Apply LLL to the rows consisting of the inverse images of the old lattice basis. Its proposed integer transform S is accepted only when det(S)=+1 or-1 and L_new=S*L_old holds exactly. Recompute the inverse, scalar steps, pullback functionals and known-digit cuts for L_new.

This preserves the entire integer kernel and changes its coordinate geometry. The large consecutive32 case retains13 scalar-visible directions. Its independent-coordinate uniform cap becomes512, but nine visible difference coordinates can still be-1,0,1, so this alone does not prove uniqueness.

## The coupled difference body

Consider two actual words with the same known digits. Put h=(F_a-F_b)/p and let t denote their visible lattice-coordinate difference. Then

    ||h||infinity <= b = (p-1)/p,
    M_K*h=0,
    q_j.h=t_j  for every visible j,
    t is integral.

For any rational direction mu on the visible coordinates and any rational lambda on K,

    r = sum_j mu_j*q_j - M_K^T*lambda

gives the necessary joint inequality

    |mu.t| <= b*||r||1.

If mu is integral, the right side can be rounded down because mu.t is integral. Pair directions e_i+e_j and e_i-e_j already remove many corners of the independently bounded integer box.

At the saved consecutive32 input, the individual difference box contains19683 vectors, including zero. Since the body is centrally symmetric,9841 representatives cover every nonzero vector up to sign. The156 pair cuts leave86 representatives.

## A certificate for a targeted difference

For a nonzero target t, choose a pivot v with t_v!=0. Minimize an L1 residual for

    q0=q_v/t_v

modulo the rows of M_K and the additional functionals

    q_j-(t_j/t_v)*q_v,  j!=v.

If the resulting multiplier on an additional functional is alpha_j, define

    mu_j=-alpha_j for j!=v,
    mu_v=1/t_v+sum_(j!=v) (t_j/t_v)*alpha_j.

Then mu.t=1. The same residual identity gives |mu.t|<=b*||r||1 for any feasible difference. Therefore b*||r||1<1 rules out that target exactly. Optimality of the LP is unnecessary for this exclusion.

If an optimal dual eta is available and nu=||r||1>0 with b*nu>=1, then h=eta/nu is a rational point satisfying M_K*h=0, q_j.h=t_j and ||h||infinity<=b. This certifies feasibility in the continuous relaxation. Such a point need not arise from actual centered scalar orbits; it does not prove codeword ambiguity.

For this input all86 surviving representatives have strict rational separators. Every nonzero visible integer difference is excluded. A zero visible difference gives scalar difference sum_j s_j*t_j=0 modulo p, so the two words are equal. This proves universal unique completion after32 consecutive initial erasures, independently of any query sample.

## Cyclic transport and limitations

The identity X*E(a)=E(a/g) transports an erasure block beginning at any cyclic coordinate to the initial block, with the negacyclic sign at the wrap. Hence the same unique-completion guarantee holds after any32 consecutive cyclic erasures, for every assignment to the remaining32 coordinates. Existence is not promised; empty fibres remain possible. The complete decoder retains the512-combination cap and filters by re-encoding. The joint uniqueness proof does not imply a one-candidate construction algorithm.

No arbitrary32-erasure guarantee is established. Alternating and seeded scattered masks still have enormous certified interval boxes, even after conditioning and metric changes. The joint difference census used here is feasible because the consecutive case has only9841 signed candidates. Scaling joint separation to those other masks requires a different enumeration or covering mechanism. No asymptotic runtime, norm bound, general Paley theorem or ambient proximity decoder follows.
