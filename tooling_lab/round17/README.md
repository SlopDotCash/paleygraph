# Round17: exact carries, compression and digit realizability

The Paley research record separates Euclidean shell radius from cyclotomic norm defect. [Pass41](../../research/parallel41-shell-arithmetic-recovery-2026-09-06.md) already proves that coefficient centering can increase the norm defect and supplies a real-subfield recovery map. This round builds tools that retain the arithmetic lost at centering, then compresses the norm calculation into a small-digit representation. It does not try to prove the conjecture.

## Centering carries under a unit transformation

Let n=2N be a power of two, p a prime with n dividing p-1, g of order n, and R=Z[X]/(X^N+1). Write F_a for the polynomial whose j-th coefficient is the centered residue of a*g^j. Its field norm has the form p^(N-1)*tau(a). The scalar tau is the norm defect; the squared radius is the coefficient energy divided by p^2.

The element epsilon=1+X+X^-1 is a unit of norm one. Multiplication by epsilon preserves the field norm before centering. It maps F_a modulo p to F_(ma), where m=1+g+g^-1. The [transport tool](norm_carry.py) retains the exact identity

```
epsilon*F_a = F_(ma) + p*C_a,
C_a has coefficients in {-1,0,1}.
```

It also checks the two-step carry identity C_(two steps)=epsilon*C_a+C_(ma), and the exact energy correction caused by the carry. [Independent review](carry_review.json) uses polynomial multiplication/remainders and integer resultants, separately from the producer's coefficient convolution and quadratic norm descent.

The small case exhausts all256 nonzero scalars at p257,n16. Four additional walks contain64,64,64 and32 steps. Across all480 transitions, the independent review computes489 resultants, including the five unit norms. Only the first case is a complete census.

The resulting exact counterexamples identify the operation where the two objectives separate. At p6700417,n64,a837595, a unit step followed by centering uses only **two nonzero carries**. Coefficient energy falls from95318544243171 to78979832004517, while the norm defect rises from1153 to92686337. A bounded order128 walk also contains five such radius-decrease/norm-increase steps at p2013265921, which lies between n^4 and8*n^4. These are finite arithmetic witnesses, not an asymptotic obstruction to every norm-based method.

## Short-relation digit compression

Choose a nonzero integral polynomial f with f(g^-1)=0 modulo p and Norm(f)=k*p. The new [compression tool](norm_compression.py) calculates

```
D = f*F_a/p,
Norm(D) = k*tau(a),
|D_j| <= floor((p-1)*sum|f_j|/(2p)).
```

The division is coefficientwise integral because F_a modulo p is supported only at the split root g^-1. This is a small-coefficient encoding of the same norm defect: compute Norm(D), then divide by the retained cofactor k. Multiplication by a nonzero f is injective in the cyclotomic field, so the representation loses no information when p and f are retained. It is not a new norm identity; the contribution is the executable, checked representation and its measured arithmetic sizes.

At g=2 the specified relation is f=2-X^-1, with k641 at p6700417. The corresponding D is ternary. For the other inputs, a bounded exact lattice reduction in dimension8 or16 supplies f; every relation is checked afterward, and no shortest-relation or minimum-cofactor claim is made.

| p | n | Relation coefficient absolute sum | Digit magnitude bound | Largest original coefficient magnitude bits | Largest digit magnitude bits | Largest original norm bits | Largest digit norm bits |
|---:|---:|---:|---:|---:|---:|---:|---:|
|257|16|3|1|8|1|61|5|
|65537|32|3|1|16|1|257|17|
|6700417|64|3|1|22|1|739|45|
|2013265921|128|13|6|30|3|2055|191|

The two p6700417 walks have the same size bounds. All484 vector records, including repeat scalars across separate walks, pass the [separate compression review](compression_review.json). Bit counts describe magnitudes and exclude signs, fixed relation metadata and the cofactor; they are not total file sizes or measured runtime speedups. At the largest prime the relation cofactor is9985208709332560769028097. Small digits therefore do not imply a small norm defect or a small kernel cofactor.

## Exact membership exposes what the digit relaxation discards

The [realizability checker](realizability.py) constructs an integral adjugate by quadratic norm descent and recovers F=adj(f)*D/k. It checks integral coefficients, their centered range and the common scalar congruences. It handles p dividing k as well: a valid control with f=p has integral centered inverse but fails the scalar congruences. The zero scalar is explicit.

The [separate review](realizability_review.json) uses rational polynomial Euclidean inversion and integer multiplication matrices. It verifies all484 positive records and450 probes, including modifications that remain real cosets. It also exhausts all6561 ternary words at p257,n16, checks every digit norm by a separate resultant, and compares accepted words with direct enumeration of all257 scalars. Ten corrupted certificates or outputs are rejected.

| Small-cube condition | Words admitted, including zero |
|---|---:|
| Digit height and norm divisibility by k=1 |6561|
| Norm belongs to the exact actual norm set {0,1,17} |945|
| Both norm and digit energy belong to an actual observed pair |945|
| Exact inverse realizability |257|

The middle rows have an oracle advantage: they know the complete actual norm or norm-energy set. They still admit688 impossible words. Adding digit energy to this oracle summary eliminates none of them in this example. Every impossible word in the small census fails centering; integrality and scalar congruences alone accept all6561.

For example, D=(-1,-1,-1,0,0,0,0,0) is realized at a207, while D'=(-1,0,0,-1,0,0,-1,0) is not. Both have norm1 and digit energy3. The inverse of D' has coefficients-224 and-191, outside the centered range[-128,128], despite satisfying all scalar congruences. This is a finite separation of proposed summaries, not a proof that every richer norm-based tool fails.

The exact geometric interface is also clear: the image lattice J=(f/p)I has index k in R, and the actual digit words are J intersected with a fundamental parallelotope for fR. This retains both lattice congruences and the section inequalities; [the derivation](DERIVATION.md) does not assume p is coprime to k. The next opportunity is to compile these conditions into a smaller recognizable digit language for sparse relations, instead of storing a dense inverse.

The [source audit](prior_art.md) identifies existing ideal-lattice, weighted-norm and toral-dynamics machinery, and the project's earlier exact sparse compression example. Historical originality is unestablished. This round supplies verified finite tools and witnesses; it supplies neither a uniform shell estimate nor a prize proof.

```sh
/opt/miniconda3/bin/python3 tooling_lab/round17/norm_carry.py
/opt/miniconda3/bin/python3 tooling_lab/round17/carry_review.py
/opt/miniconda3/bin/python3 tooling_lab/round17/norm_compression.py
/opt/miniconda3/bin/python3 tooling_lab/round17/compression_review.py
/opt/miniconda3/bin/python3 tooling_lab/round17/realizability.py
/opt/miniconda3/bin/python3 tooling_lab/round17/realizability_review.py
/opt/miniconda3/bin/python3 tooling_lab/round17/checkpoint.py
```
