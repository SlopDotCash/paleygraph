# Root review of the three pass40 lanes

The arguments below were reviewed against the current source notes. The
independent [checker](../experiments/parallel40_independent_verify.py) and
[output](../results/parallel40_independent_verification_2026_09_06.json)
cover the finite arithmetic described here. This is ordinary mathematical
review and exact computation, not Lean certification or external peer review.

## Short shells and reciprocals

The [shell converse](parallel40-shell-converse-2026-09-06.md) is valid.
Rank one modulo p forces divisibility of every (N-1)-minor by p^(N-2),
so T=tau*p/F is integral when Norm(F)=p^(N-1)*tau. Its norm is
p*tau^(N-1), and it vanishes at the unique root where F is nonzero.
The hypothesis V<1 gives 1<=tau<p by the existing norm inequality.
Thus the converse has the correct p-adic valuation, but its cofactor
can be enormous. The scalar reciprocal family must not be confused
with all possible kernel elements.

For p=6700417,n=64, root independently constructed the centered vector
and obtained V=6120237/6700417<1. A 32-by-32 fraction-free Bareiss
determinant, separate from the agent's quadratic norm descent, gives
Norm(F)=p^31*1217. The reciprocal multiplication and its determinant
also check exactly. Since1217 is prime and different from p,
integrality of m*p/F forces1217|m through its norm p*m^32/1217.
The displayed integral choice attains m=1217. This proves minimality
only within that scalar family; a different kernel element of norm641*p
already exists.

The nonprincipality argument is also valid without a class-group
computation. At q=641 the index-ten character expansion bounds every
period by (1+9 sqrt641)/10<47/2. A norm641 evaluation-kernel generator
would, by pass39, force a period >1200/49. The exact gap is97/98.
If the evaluation kernel at p and root2 were principal, its generator
would divide2-X in the integer ring, giving an integral quotient of
norm641 that vanishes at2 modulo641. The generator is invertible
modulo641 because its determinant is p. This contradicts the q bound.
Conjugation transfers the conclusion to the root2^(-1) associated with F.

The optional flat-Gram discussion correctly separates an ideal-class
condition from a norm-of-unit condition. Its converse concerns a flat
integral element with a specified ideal, and does not imply that the
evaluation prime itself is principal. The main finite counterexample
does not need that optional ideal-factorization discussion.

## Actual distinct-edge incidence

The [distinct-edge note](parallel40-distinct-edge-input-2026-09-06.md)
retains the correct normalization n sum m*rho with y-x=1. This is
compatible with the earlier triangle weight after alpha=y,beta=x.
The collision identity counts ordered pairs of distinct triples sharing
a normalization: their unique nonidentity common scale q gives
the product of falling factorials of |D_l intersect qD_l|. Different
levels are disjoint; repeated appearances of one level require those
falling factorials. The correction bound follows from
m-1<=m(m-1)/2, with no injectivity assumption.

Root independently enumerated the actual levels at p=215535361,n=128
from H-H, then directly counted y-x=z through the two smaller level
sets. The result is52736, whereas dropping quotient multiplicity gives
51968. The weighted block is144850944. A separate quotient computation
has image size1438, domain size1440, and collision count4; the level
intersection formula gives the same4. The largest additive count on a
repeated fiber is3, making the correction exactly128*3*4/2=768.
The three displayed normalized solutions were individually checked.
This refutes the actual-level/distinct-edge injectivity candidate in
the quartic window. It does not refute an asymptotic energy inequality.

## Classical structure and the norm convention

The [classical note](parallel40-classical-structure-2026-09-06.md)
uses r_(V-V)(s)<=1 for nonzero s when V is Sidon. Counting mixed
energy and applying Cauchy-Schwarz gives its translation-envelope
bound. The container bounds follow by counting the distinct sums
of V intersect W. The exact character Gram identity then gives the
translation and smoothing inequalities without a probabilistic
independence assumption.

The normalization caveat is essential. Error small relative to ||F_V||2
is impossible for a nonzero translation in the sparse growing regime.
For G_V=F_V/k, however, every shift has squared error <=2p/k, so all
shifts meet the usual larger tolerance epsilon||chi||2 once
k>=2p/[epsilon^2(p-1)]. Root requested this distinction, and the agent
retained the exact ||chi||2^2=p-1 correction. The note therefore does
not claim that standard almost-periodicity fails for this convolution.

Root reproduced all538 Sidon sets of sizes2-4 in F7,F11,F13, including
5916 translation identities and42628 set-pair/smoothing checks. The
source theorem applicability remains scoped to the versioned primary
statements cited in the agent's note. No new uniform classical moment
bound, subgroup period bound, or prize reduction follows.
