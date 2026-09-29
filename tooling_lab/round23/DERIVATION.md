# Exhaustive integer-box covers

This extends the coupled difference model in [round22](../round22/DERIVATION.md). For two actual words agreeing on the known digits, the visible integer difference t lies in a certified symmetric box. A rational direction mu and residual identity give a valid inequality |mu.t|<=r. The new mechanism stores how these inequalities exclude or narrow entire integer boxes.

## Root coverage

Suppose -b_j<=t_j<=b_j. Ignore zero. Central symmetry means it is enough to consider vectors whose first nonzero coordinate is positive. For each j with b_j>0, create the box

    t_i=0 for i<j,
    1<=t_j<=b_j,
    -b_i<=t_i<=b_i for i>j.

These root boxes are disjoint and contain exactly one representative of each nonzero vector up to sign. The checker reconstructs them from the certified individual radii, rather than trusting a supplied list of roots.

## A narrowing step

For a current integer box l_i<=t_i<=u_i and a valid cut |sum a_i*t_i|<=r, the other coordinates have exact contribution range

    L=sum_(i!=j) min(a_i*l_i,a_i*u_i),
    H=sum_(i!=j) max(a_i*l_i,a_i*u_i).

Consequently -r-H<=a_j*t_j<=r-L. Divide by a_j with its sign respected, round the new lower endpoint up and the upper endpoint down, and intersect with the existing integer interval. Each trace entry names the cut and coordinate and stores the resulting interval. A separate implementation reconstructs the implication exactly.

An empty coordinate interval proves the node impossible. A box whose total cut contribution range lies wholly above r or below-r is also impossible. The trace cannot discard an integer point satisfying all the cuts.

## Exact branching and unresolved leaves

When propagation cannot finish a box, split one coordinate at an integer m into [l,m] and [m+1,u]. These children are disjoint and cover every integer point in their parent. The checker derives the two child boxes, verifies both are present, forbids repeated nodes or cycles, and requires every node to be reachable from a root.

The search may propose additional cuts by separating a surviving integer target with round22's rational norm-residual compiler. Every accepted cut is globally valid; even a cut failing to exclude its own target can validly narrow other regions. Discovery budgets limit the work. A budget cutoff produces an explicitly unresolved leaf whose full remaining integer region is retained.

Universal uniqueness is certified only when every root is covered by exclusions and there are no unresolved leaves. It is not inferred from a lack of found codeword pairs. Zero visible difference gives zero scalar difference modulo p, hence equal codewords, as in round22.

## Continuous obstructions and small controls

When the exact norm dual provides a rational continuous point with M_K*h=0, q_j.h=t_j and ||h||infinity<=(p-1)/p, that point witnesses that the retained continuous constraints cannot exclude t. It does not necessarily come from two centered scalar orbits or satisfy the invisible integer-coordinate constraints.

On complete small codebooks, the checker enumerates every actual pair agreeing on the known digits, computes its visible difference, and verifies it remains inside an unresolved leaf. This tests that real ambiguities survive every cut, propagation and partition step. A supplied scalar-p-squared example can be unique as an actual code while the relaxed body still has feasible integer visible directions, exposing the extra information that the cover has omitted.

The covering tree changes the certificate's size, not the mathematical model. Exact branch-and-bound certification has established prior art, recorded in the [audit](prior_art.md). No arbitrary erasure theorem or prize estimate follows automatically from success on consecutive masks.

## Resuming without discarding earlier work

A leaf can remain unresolved when first visited and become contradictory after cuts are learned elsewhere. The refinement prototype reloads its saved input box and narrowing trace, reapplies every currently available cut, and appends any new justified narrowing steps. Earlier excluded or split nodes remain byte-for-byte unchanged. Existing cuts remain an exact prefix of the refined cut list. If more branching is needed, new nodes are appended and the former unresolved leaf becomes their parent.

The refined33-erasure run illustrates this issue sharply: the first run left one singleton unresolved after exhausting its64-cut discovery budget. Later cuts already excluded that singleton, but the first search never revisited it. Resaturation appends the missing implication, producing a complete proof with the same64 cuts and1337 nodes. A separate review checks5620 exact narrowing steps and preservation of all1336 previously finished nodes.

This proves unique completion after33 initial consecutive erasures on the saved p2013265921,N64 input. The scalar transport X*E(a)=E(a/g) extends it to any33 consecutive cyclic erasures. The remaining31 digits may have no completion, but cannot have two. The complete interval decoder retains a uniform12288-combination cap; the uniqueness proof is separate from that candidate-enumeration cost.

The36-erasure cover was also checked completely:64 cuts,1867 nodes and7953 narrowing steps exclude all3690562 nonzero representatives up to sign. Hence any36 consecutive cyclic erasures leave at most one completion from the remaining28 digits on this same finite input. Its complete decoder cap is27648 combinations. At40 erasures, the verified cover retains19 unresolved leaves, and a rational feasible point prevents a uniqueness conclusion from this cover.

## Information carried by the invisible coordinates

Let P_j=(A_f/k)*H_j, where H_j is basis row L_j embedded on the erased digits. These are the inverse images of integer erasure-lattice rows. In all five large tested geometries, exact checks show each P_j is integral and P_j[i]=s_j*g^i modulo p. Consequently invisible rows, with s_j=0, are integer p-multiple vectors. This implication is verified for the saved bases; degenerate scalar-p-squared controls demonstrate that it cannot be assumed for every relation.

For the consecutive33,36,40 bases, the invisible P_j are signed p*e_i vectors on coordinates0..19,0..22,0..26 respectively. At40 erasures write J={0,...,26} and I={27,...,63}. For a fixed visible integer target t, put H(t)=sum_(j visible) t_j*P_j. Every lifted difference with that visible target agrees with H(t) on I; on J it can vary by arbitrary p-multiples. Its scalar difference is delta=H(t)[0] modulo p.

An actual pair with this target exists precisely when some a satisfies

    F_a[i]-F_(a-delta)[i]=H(t)[i] for every i in I.

Necessity follows from the fixed inverse-image coordinates. For sufficiency, actual differences have the same residues delta*g^i modulo p. Their difference from H(t) is therefore supported on J and is an integer combination of the invisible p*e_i rows. Multiplying by f/p yields an integer erasure-lattice difference with the same visible t and zero known digits. This proves that the scalar count below is exactly the number of actual ordered pairs with this target, not just another necessary test.

## Exact orbit intersection

Let m=(p-1)/2. The fixed-coordinate equation above is equivalent to

    F_a[i] in [-m,m] intersect [H(t)[i]-m,H(t)[i]+m].

The allowed interval has p-|H(t)[i]| integer points. Choose its tightest coordinate i0. Since multiplication by g^i0 is a bijection modulo p, enumerating those centered values enumerates every possible a exactly once. Filter the surviving a through the remaining intervals. The producer uses this interval formulation; the independent reviewer enumerates centered F_b[i0], forms a=b+delta, and checks the actual centered differences directly. All large modular products are checked to fit int64; certificate arithmetic otherwise remains exact.

The saved40-erasure continuous target is t=(1,0,...,0), with delta=345549834. Its tightest coordinate is62, where F_a[62] must lie in[978728438,1006632960]. Both implementations exhaust all27904523 candidates and return zero. There are37 fixed coordinates and27 free integer carry coordinates. The small oracle suite checks every nonzero delta and every two-lift assignment on masks{0},{0,2},{0,1,2,3} at p17,41,97:3344 intersections, including empty and nonempty cases, against246928 direct scalar-pair checks.

## An integral gap, beyond a fractional relaxation

Set H[i]=H(t)[i] for i in I and H[i]=center(delta*g^i modulo p) for i in J. The [gap certificate](gap_certificate.json) verifies that H is integral, |H[i]|<=p-1, and H[i]=delta*g^i modulo p. The digit difference D=f*H/p is integral, vanishes on every known digit, satisfies |D[i]|<=2B, and has integer coordinates in the full erasure lattice, with visible target t. Yet the exact orbit count proves no actual pair has that target. Thus even the full integral bounded difference lattice admits a false candidate here; merely restoring invisible integrality does not remove it.

This does not prove uniqueness at40 erasures. It excludes one particular nonzero visible target up to sign and identifies the missing joint scalar-orbit condition. Other unresolved regions remain explicit. A general exact orbit-box solver or certificates grouping many such targets are the next tooling questions.
