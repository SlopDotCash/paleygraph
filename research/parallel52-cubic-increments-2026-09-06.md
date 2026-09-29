# Cubic incidence increments along the subgroup tower

The full Paley conjecture and Proximity Prize remain unproved. This
pass gives a local cost for primitive collisions in the increase of
shifted multiplicative energy. The available source bound on that
increase is still too weak. No uniform exponent improves.

## 1. Four fine cells merge into one coarse cell

Fix an odd prime p, dyadic s=2k>=4 dividing p-1, and K=H_k subset H_s.
For a subgroup H_t, partition F_p minus {0,-1} by the power-label pairs

    (z^t,(1+z)^t).

Write C^(t) for their cell sizes and

    X_t = sum_cells (C^(t))_3,
    (c)_3 = c(c-1)(c-2).

By the [existing collision bijection](mixed-periods-and-shifted-energy.md),
X_t is the nontrivial multiplicative energy of (H_t-1) minus {0}, with
the two trivial pair matchings removed. Squaring both labels maps the
k-cells to s-cells. Every parent cell is the disjoint union of at most
four children. If their sizes are c_1,...,c_4, including zeros, and
C=sum_i c_i, then

    (C)_3-sum_i(c_i)_3
      =3 sum_i (c_i)_2(C-c_i)
         +6 sum_(i<j<l)c_i c_j c_l >=0.                (1)

This counts ordered triples using either two or three distinct children.
Summing over all parents proves X_s>=X_k. Put Delta_s=X_s-X_k.
The increases telescope exactly: sum_(M<s<=N)Delta_s=X_N-X_M.

There is also an equality criterion. Delta_s=0 if and only if every
parent cell with at least three points consists of a single nonempty
child cell. A mixed parent with at least three points always has an
ordered triple not contained in one child. Thus unchanged cubic excess
can represent inherited rich cells, not their absence.

## 2. The primitive fibers occupy two equal children

For a parent cell in row zero, z belongs to H_s, so its first fine label
is 1 or -1. For a choice of square root u of its second parent label,
write its four children as

                      second fine label
                        u       -u
    first label  1      a        b
                -1      c        c

The equality in the last row follows from z -> 1/z: for z in H_s\K,
z^k=-1 and ((1+z)/z)^k=-(1+z)^k. Each primitive inverse pair is split
between these two children. Consequently this c is precisely the
primitive root multiplicity c_s(beta) in the preceding notes.

Let S=a+b. The contribution of this block to Delta_s is exactly

    (S+2c)_3-(a)_3-(b)_3-2(c)_3
      =3ab(S-2)+6c^2(c-1)+6cS(S+2c-2).                 (2)

All three terms are nonnegative for nonnegative integer a,b,c. The
last term is at least 6cS: when cS>0, S+2c-2>=1. Summing (2) over
row-zero parents, which form a subset of all parent cells, gives

    Delta_s >= 6Q_s + 6L_s,
    Q_s = sum_beta c_s(beta)^2(c_s(beta)-1),
    L_s = sum_blocks c(a+b).                            (3)

Retain A_s=sum c_s(c_s-1) and the triple mass Y_s=sum_(c_s>=3)c_s(c_s-2).
Since c>=2 whenever c(c-1)>0, Q_s>=2A_s. Counting nonzero sums by their
K-cosets gives two useful identities:

    sum_blocks ab=A_s,       T_s=k L_s.

For the first, the mixed energy is both k^2+s A_s and
k^2+k sum_j a_j a_(j+sigma)=k^2+s sum_blocks ab, where sigma is the
order-two quotient class H_s\K. For the second, each of the two fine
cosets contributes k times its pure representation count times c.

The exact energy recurrence from [pass 51](parallel51-collision-transport-2026-09-06.md)
therefore yields

    D_s-D_k=6A_s+4L_s,
    Delta_s >= (3/2)(D_s-D_k)+3A_s.                     (4)

In particular Delta_s=0 implies A_s=Y_s=T_s=0 and D_s=D_k. The
previous primitive collisions can remain in D_k, as the example below
illustrates. Formula (4) is an increment statement; it is not a new
best global comparison between D_N and X_N.

## 3. Quantitative mass budgets and their current limitation

For c>=3 one has c^2(c-1)>=6c(c-2): after dividing by c the difference
is (c-3)(c-4)>=0. Thus (3) implies

    Delta_s >=6Q_s>=36Y_s.                             (5)

A second useful inequality comes from weighted Cauchy and the degree
budget sum_beta c_s(beta)=s/4:

    Y_s^2 <= (s/4) sum_(c_s>=3)c_s(c_s-2)^2
            <= (s/4)Q_s <= s Delta_s/24.               (6)

Here (c-2)^2<=c(c-1). Consequently, for any dyadic cutoff 2<=M<=N,

    sum_(M<s<=N) sqrt(Y_s)
      <=24^(-1/4) [sum_(M<s<=N) s^(1/3)]^(3/4)
                       (X_N-X_M)^(1/4)
      << N^(1/4)(X_N-X_M)^(1/4).                       (7)

The first line is Holder with exponents 4/3 and 4. The geometric sum
in the bracket is at most N^(1/3)/(1-2^(-1/3)), so no number-of-levels
factor is hidden. Likewise (4)-(5) give the exact summed budgets

    X_N-X_M >=36 sum_(M<s<=N)Y_s,
    X_N-X_M >=(3/2)(D_N-D_M)+3 sum_(M<s<=N)A_s.         (8)

Using (7) in the radical recovery from [pass 50](parallel50-tower-mass-2026-09-06.md)
gives

    E_N << (N/M)E_M+N^2+N^(3/2)(X_N-X_M)^(1/2).        (9)

At the existing cutoff M near N^(80/87)/(1+log N)^(4/29), the base
term is already O(N^(7/3)). Hence an upper-interval increment bound

    X_N-X_M << N^(5/3)                                 (10)

would suffice for the intermediate energy7/3 and triangle17/3 targets.
It is not proved. It is only a sufficient substitute for the direct
square-root mass criterion, not a necessary condition for that criterion
or for the full conjecture.

The already audited [Shkredov theorem, Theorem 6](https://arxiv.org/abs/1504.04522)
gives X_N<<N^2(1+log N) for N^2<p. Its complete statement was reread
in the archived primary text in this pass. Nonnegativity gives only
X_N-X_M<=X_N. Substituting this into (9) returns the known energy
power5/2, with a square-root logarithm, which is weaker than the
already available49/20. The source itself also records that5/2
consequence. Subtracting two coarse upper bounds gives no bound of
the strength (10). Thus this charging argument does not close the gap.

## 4. A positive cubic excess can remain exactly unchanged

For the previously certified quartic endpoint

    p=2144280833,  N=256,  g=231737012,

direct shifted-product counts now give

| s | X_s | Delta_s | Y_s | D_s |
|---:|---:|---:|---:|---:|
| 32 | 0 | 0 | 0 | 0 |
| 64 | 36 | 36 | 0 | 0 |
| 128 | 1260 | 1224 | 3 | 48 |
| 256 | 1260 | 0 | 0 | 48 |

At the last step every rich parent cell is inherited from one child.
This disproves a pointwise claim X_(2k)>=c X_k with any fixed c>1
for every eligible step, including at a quartic endpoint. It does not
exclude an eventual large-size or additive-error version. The example
shows why positive old excess cannot automatically be charged again
at the next level. It does not disprove any target energy estimate.

## 5. Exact verification

The [checker](../experiments/parallel52_cubic_increments.py) reconstructs
all parent cells with at least two points by the exact inverse map

    z=(b-1)/(a-b),  a,b in H_s\{1}, a!=b.

Every point of every repeated cell occurs this way. Cells of size zero
or one contribute no cubic count, nor can they contain a repeated
child. Thus partitioning these recovered points by their fine power
labels suffices to compute both X_s and X_k exactly. The checker
independently counts products in the two shifted subgroups to verify
both totals. Row-zero singleton cells are separately included when
checking (2)-(4).

All 368 distinct step cases passed, covering 71 existing towers and
397 cutoff checks. Reconstruction uses 1,131,420 ordered pair parameters
and produces 427,090 parent cells; 1,424,624 shifted-product pairs are
counted independently. The [results](../results/parallel52_cubic_increments_2026_09_06.json)
retain each step, its positive row blocks, and the unchanged-excess
witness. The fields and generators come from the preceding certified
data. No new prime or complete prime classification is claimed.

Root completed ordinary derivations and exact checks. No independent
agent, Lean, external review or novelty claim is made. No Lean process
or Prove2Me submission was changed. The next task is an actual saving
for the upper-interval collision mass or its increment budget; the
nonnegative accounting identities alone do not supply it. Uniform
exponents and the full goals remain unchanged.
