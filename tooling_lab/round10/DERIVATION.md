# Decomposition and information audit

Use the round9 conference contract S=Sᵀ, S1=0, S²=qI−J, with signs ±1 off the zero diagonal. Fix C of size n and delete a,c. Let A=C\{a,c}, m=q−n, N=m(m−1), D≤6. Write T_j(B)=Σx e_j(S[x,B]), with negative elementary degrees equal to zero. The numerator of the conditional mean over ordered outside insertion pairs is

```
M_ac = Σx [N eD(A)−2(m−1)R_C(x)e(D−1)(A)
            +(R_C(x)²−(m−1)−1_C(x))e(D−2)(A)].
```

The first [prototype](deletion_pair_means.py) evaluates this for all pairs using three Boolean weights and zero-row corrections. [pair_means_backend.cpp](pair_means_backend.cpp) streams the same algebra over the field. Its total cost is O(q(nD+n²)) with O(q+n²) storage, including the character table. There is no character convolution and no insertion-pair enumeration.

## Independent inward-target identity

Set r=m−D, v=n−D+1,

```
H=(r+1)(r+2)
L=−(2r+5)(n−D)−r−3
J=v(v+1).
```

Then

```
M_ac = H T_D(A)
       −2(r+1)[T_D(C\{a})+T_D(C\{c})] + 2T_D(C)
       +L T_(D−2)(A)
       +2v[T_(D−2)(C\{a})+T_(D−2)(C\{c})]
       +J T_(D−4)(A) + B_ac,
```

where the entire boundary term uses only the induced matrix S[C,C]:

```
B_ac = 2r Σx∈A e(D−2)(S[x,A])
       −2 Σx∈{a,c} e(D−2)(S[x,A])
       −2 Σx∈A (Sxa+Sxc)e(D−3)(S[x,A])
       −2v Σx∈A e(D−4)(S[x,A]).
```

To derive it, put R_A=Σy∈A Sxy and z_A=1_A(x). The elementary recurrence

```
R_A e_j(A)=(j+1)e_(j+1)(A)+(n−2−z_A−j+1)e_(j−1)(A)
```

comes from whether multiplication repeats a selected sign. Use it twice for R_A²e_(D−2), expand R_C=R_A+Sxa+Sxc, and use Sxa²=1−1_{a}(x), Sxc²=1−1_{c}(x). Finally replace (Sxa+Sxc)e_(j−1)(A) by the sum of the two singleton-deletion row targets minus2e_j(A), and SxaSxc e_(j−2)(A) by the full target minus those terms. Collecting gives the displayed identity. The low-degree conventions agree with the separately checked degrees0 through3.

The independent [inward-target backend](inward_targets.cpp) directly accumulates degree D,D−2,D−4 targets for C, all singleton deletions and all pair deletions using counts of positive and negative entries. It does not use the forward weights or weighted Gram algorithm. [decomposition.py](decomposition.py) combines those totals with the induced boundary, giving a complete second arithmetic computation of every large-field pair mean.

## Removing vertex-additive effects

For an edge array F_ac, let R_a=Σc≠a F_ac and E=Σa<c F_ac. Its orthogonal projection away from all arrays u_a+u_c is

```
(PF)_ac = F_ac −(R_a+R_c)/(n−2)+2E/[(n−1)(n−2)].
```

Every projected row sums to zero. This is an ordinary least-squares projection, not a new invariant construction. It eliminates the constant and singleton-deletion terms above:

```
P M = H P[T_D(A)] + L P[T_(D−2)(A)]
      + J P[T_(D−4)(A)] + P B.
```

At D=6, T2(A)=−C(n−2,2) by S²=qI−J, so its projection vanishes. For n≤7, A has fewer than6 points, so T6(A)=0 as well. Thus **every projected mean in that size range is explained by quartic incidence and the induced boundary**. At q17,n7 specifically, H=30,L=−20,J=6,N=90. This applies to the saved round8 twins; their projected separation supplies no new character sum beyond that information.

At n≥8 the high-degree term survives. It is target-related information: summing T6(C without a,c) over all pairs gives C(n−6,2)T6(C). Therefore it must not be treated as an independent explanation of T6(C). At larger tested inputs the normalized high-degree term accounts for the projected mean to within roughly10⁻⁶ in relative Euclidean norm. This is an actual-input observation, not a uniform estimate or a new closure theorem.

## Integer range

Hypothetical Boolean divisions can increase elementary coefficients when the hypothetical signs are not present in a row. A safe bound after two divisions is

```
|e_j after division| ≤ Σt=0..j (t+1)C(64,j−t).
```

At j6 this is92,306,647. At q≤10,000,000,n≤64,D≤6, the conservative row-Phi bound in [controls_verification.json](controls_verification.json) is9,230,676,263,919,979,653,248. The signed128 accumulator bound is below1.48×10³⁰, safely below2¹²⁷. Direct inward targets fit signed64 bits because qC(64,6)<2⁶³. Fractions and final projections use Python integers. These bounds cover accepted inputs, not just the sampled primes.
