# What the paired-row energy calculation cannot distinguish

This is an elementary correctness statement for the diagnostic, derived from its coefficient formula. Its historical novelty is unestablished. It concerns an explicitly larger class than quadratic-character matrices, so it does not transfer an example in that class into a counterexample to a Paley conjecture.

## Statement

Let S be a symmetric q-by-q matrix with zero diagonal, off-diagonal entries in `{−1,1}`, zero row sums, and

```
S S^T = q I - J.
```

For a d-subset A of columns, define

```
E_S(A) = sum_x product_(a in A) S_(x,a).
T_(S,d)(C) = sum_(A subset C, |A|=d) E_S(A),  |C|=n.
```

For `d<=min(n,q-n)`, the complete Johnson harmonic L2 energy spectrum of `T_(S,d)` on the uniform n-subset slice depends only on q,n,d. In particular, it is the same for every S satisfying these hypotheses. This statement does not assert equality of pointwise values, value distributions, maxima, or conditional distributions on an arithmetic family.

## Proof from the exact diagnostic

Fix two distinct rows x,y and put `s=S_(x,y)=S_(y,x)`. At columns x,y their pairs are `(0,s)` and `(s,0)`. Off those columns let `N_(a,b)` count pair types with a,b in `{−1,1}`. Zero row sums and the off-diagonal Gram entry give

```
sum N_(a,b) = q-2,
sum a N_(a,b) = -s,
sum b N_(a,b) = -s,
sum ab N_(a,b) = -1.
```

The four counts are therefore

```
N_(1,1)   = (q-3-2s)/4,
N_(-1,-1) = (q-3+2s)/4,
N_(1,-1) = N_(-1,1) = (q-1)/4.
```

Multiplying both entries of every pair by s makes the two zero-site pairs `(0,1),(1,0)` and gives the same normalized type table for every distinct row pair. On equal rows, the types are `(0,0)` once and `(1,1),(-1,-1)` each `(q-1)/2` times.

For `0<=r<=d`, form

```
K_r = sum_(|A|=|B|=d, |A intersection B|=r) E_S(A) E_S(B).
```

After expanding the two row sums, each row pair with entries a_z,b_z contributes

```
[X^r Y^(d-r) Z^(d-r)] product_z (1+a_z*b_z*X+a_z*Y+b_z*Z).
```

Simultaneously negating a_z and b_z leaves this coefficient unchanged, because the total power of Y and Z is `2(d-r)`. Hence every distinct row pair contributes the same coefficient, and every equal row pair contributes the same diagonal coefficient. Thus every K_r depends only on q and d. No field representation or arithmetic labels entered this calculation.

Now let `(C,D)` be a uniformly chosen ordered pair of n-sets at Johnson distance j. For fixed marked d-sets A,B, the probability `A subset C, B subset D` depends only on q,n,j and `r=|A intersection B|`. It follows that

```
E[T_(S,d)(C) T_(S,d)(D)] = sum_r K_r * probability(q,n,j,r,d)
```

is independent of S for every j. The Johnson distance operators have the known Eberlein eigenvalues on harmonic levels0 throughd. Their distance-correlation vector determines the squared norm on each level by the invertible change of basis used in the diagnostic. Therefore those energies are independent of S. This proves the statement.

## Consequence for tool design

The exact L2 compiler can locate the degree at which variation occurs, but it cannot identify additional arithmetic structure absent from these row-pair types. Its speed comes partly from forgetting that structure. A stronger discriminator must retain additional information, such as higher row configurations or a specified arithmetic conditioning rule.

A pair of matrices satisfying the hypotheses but having different distributions of E_S on d-subsets would make this limitation explicit. Different distributions also rule out a simultaneous vertex relabelling, since such a relabelling preserves the multiset of these values. The accompanying prototype is testing that finite witness criterion; this note does not assume that its search succeeds.

The formulas use the established Johnson/Hoeffding decomposition and classical conference-matrix intersection counts. See the [round2 source audit](../../round2/novelty/exchange_novelty.md). This is a derived contract and limitation for the implemented tool, not a claim of a new association-scheme theorem or an answer to either prize.
