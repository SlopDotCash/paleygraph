# A quartic counterexample and exact collision valuations

An actual subgroup inside the quartic window has no diagonal rich mass
but positive off-diagonal collision excess. This strengthens the pass42
counterexample, which was outside that window. It rules out the proposed
constant-factor diagonal domination in the target range, but is not a
Paley counterexample or evidence that the weighted-triangle bound fails.

## 1. The actual subgroup

Take

    p=278177, n=32, H=<160164>.

Trial division certifies p prime, and the generator has exact order 32.
The required window is 262144<=p<=1048576. With a(C)=|(H-1) intersect C|
over nonzero cosets, its positive multiplicities are fifteen values 2
and one value 1. Thus

    E(H)=2976=3n^2-3n,
    sum_C(a(C))_3=0,
    B=E^times((H-1) minus {0})=1963,
    X=B-[2(n-1)^2-(n-1)]=72=X_dist.

The twelve cells with rho>=3 each have rho=3, three pairwise distinct
edge cosets, and T_b=144. Exact triangle enumeration gives

    W=1849344, W_low=1683456, W_high=165888,
    K_quot=8097.

Consequently X_dist<=C*(X-X_dist) fails for every finite C even when
n^4/4<=p<=n^4. Rich off-diagonal cells can also survive minimum additive
energy. This separates the arithmetic parameters; it does not contradict
the [energy-based sufficient triangle condition](parallel43-energy-prime-exceptions-2026-09-06.md).
Here all nonzero difference multiplicities are at most 2, so the simple
bound W<=4E_*^2 already gives W=O(n^4).

The [checker](../experiments/parallel43_pointwise_inputs.py) reproduces the
full rich-cell list by the exact shifted-product-to-triple bijection,
then independently checks the triangle identities using a full-field
rho count and direct differences. The [certificate](../results/parallel43_pointwise_inputs_2026_09_06.json)
contains all cell members. The exploratory search encountered this example
after 82 eligible primes; no all-prime classification or minimality claim
is based on that prefix.

## 2. Distinct entries in the product equation

List the m=n-1 nonzero shifted roots r_a=h_a-1. For every product value z,
let d_z count square pairs r_a^2=z, and let o_z count unordered pairs
of distinct roots with r_a*r_b=z. The ordered product multiplicity is
d_z+2o_z. Subtracting the two trivial pair matchings gives the exact identity

    X = sum_z d_z(d_z-1) +4 sum_z d_z o_z
        +4 sum_z o_z(o_z-1).                               (7)

These are the contributions from two, three, and four distinct shifted
entries. Nontrivial equal products cannot share an entry because every
shifted root is nonzero. At the displayed quartic witness the three
contributions are (0,0,72): all excess uses four distinct entries.
The word square in d_z concerns product pairs. It is different from the
diagonal incidence cells u=1, v=1, or u=v.

## 3. The valuation loss is exactly higher-precision collision mass

For fixed dyadic n, let the nonzero integer mathcal_P_n be the product
over all M nontrivial ordered complex shifted-product differences used in
[pass42](parallel42-excess-equations-2026-09-06.md). At a split prime p,
lift one primitive n-th root compatibly modulo p^r and let H_r be its
order-n subgroup in the units of Z/p^rZ. Define X_r by the same product
count, with the same two trivial pair matchings removed. These rings are
not finite extension fields.

The lifted nonidentity shifted roots are units and remain distinct
modulo p. Products may collide, with exactly the same weights as in (7).
The partition into equal product values refines with precision, so X_r
is nonnegative and nonincreasing. For each nontrivial tuple, its p-adic
valuation equals the number of precisions at which that tuple vanishes.
Summing over the finite tuple family proves

    v_p(mathcal_P_n)=sum_(r>=1) X_r.                         (8)

All factors are nonzero over the complex cyclotomic field, hence also
under its p-adic embedding. Their valuations are finite. Once X_r=0,
all later terms vanish. Formula (8) is stronger bookkeeping than
X_1<=v_p(mathcal_P_n); it supplies no uniform upper bound for X_1.
This is the shifted-product counterpart of the existing additive
[kernel-discriminant lifting identity](kernel-discriminant.md).

Per precision, hashing the n(n-1)/2 unordered products and their weights
computes X_r with O(n^2) modular arithmetic and storage; enumerating all
O(n^4) pairs of products is unnecessary. The checker compares the summed
lifted excess with every independent split-prime norm valuation from
the complete n=4,8,16 factorizations in pass42. It also records lifts at
the new n=32 witness and the known n=64 and n=128 quartic examples.
These are exact local certificates, not conductor-uniform arithmetic
estimates. No new Lean run or Prove2Me submission was needed.
