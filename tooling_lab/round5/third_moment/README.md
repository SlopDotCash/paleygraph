# Exact third moments beyond six-column sets

This compiler evaluates the third moment of a row-character statistic over **all** n-column subsets from a compressed inventory of row triples. It handles n greater than the statistic's degree, retains repeated-row terms, and uses exact integers and rational probabilities. Its output agrees with independent exhaustive checks, including all 450,978,066 eight-column subsets in each of the two order49 examples.

The million-prime, n31 case uses35 coefficient evaluations, independently of the1,750 input inventory bins or the number of column subsets. The subsequent [symmetry adapter](../trace_invariants/) reduces the same calculation to19 evaluations. These are distributional measurements, not maximum or prize bounds.

## Definition and input contract

Let S be a symmetric conference sign matrix with zero diagonal, `S 1=0` and `S²=q I−J`. For a uniformly chosen n-element column set C, define

```
T_d(C) = sum_x e_d((S_xy)_(y in C)).
```

The public function `global_third_moment(record,n,degree=6)` returns `E[T_d(C)^3]` as an exact rational. Its normalized interface supports d=0,2,4,6. It requires independently justified row-transitivity and sign-complement normalization: each normalized parameter has ordered-row weight q(q−1). This holds for the prime Paley inputs and the separately checked Paley49/Peisert49 inputs used here. Conference identities alone do not guarantee that normalization.

An input record supplies q and the joint counts of `(e01,e02,e12,tau)`, with e01=1 and `tau=sum_y S_0y S_1y S_ty`. Nonnegative type multiplicities, class totals, several trace identities and optional metadata are checked. Those checks are necessary diagnostics, not certificates that an arbitrary supplied inventory is realized by a graph.

## Union coefficients retain the information lost by the n=d shortcut

For an ordered row triple, write its column sign pattern as `(a_y,b_y,c_y)` and put

```
H_y = (1+a_y t1)(1+b_y t2)(1+c_y t3) − 1.
F(z,t1,t2,t3) = product_y (1+z H_y).
A_k = [z^k t1^d t2^d t3^d] F.
```

The z exponent counts distinct columns in the union of the three selected d-sets. Consequently the expected product of the three row statistics is

```
sum_(k=0)^min(n,3d) A_k (n)_k/(q)_k.
```

There is no additional binomial normalization. At n=d only one union size survives, allowing the earlier product-of-signs shortcut. For larger n the different union sizes and the row-boundary zeros must remain separate.

The C++ backend computes the elementary polynomials in the H_y by the exact Newton recurrence. Its power sums are assembled from sign-pattern multiplicities. The small universal coefficient

```
sum_(h=0)^u (−1)^(u−h) binom(u,h) binom(h,i) binom(h,j) binom(h,l)
```

counts three subsets of a u-element set, of sizes i,j,l, whose union is the entire set. This gives `[t1^i t2^j t3^l] H_y^u` after multiplying by sign powers. The backend checks the exact divisions in Newton's identity. All large arithmetic uses Boost `cpp_int`; the declared population limit is10^12 and degree limit6.

Three special histograms handle all-equal rows and the two mutual signs for exactly-two-equal rows. At even degree the latter contributions coincide. The final sum is

```
q * all_equal + 3 q(q−1) * twice_equal
    + q(q−1) * sum_normalized_parameters distinct_row_product.
```

There is no extra factor6 in the distinct-row term.

## Why finitely many synthetic evaluations suffice

At fixed q and mutual signs, all eight nonzero triple-type multiplicities are affine in tau. In the formal logarithm of F, a tau-dependent term must have positive odd degree in each of t1,t2,t3. At most d such factors can reach target degree(d,d,d). Thus **each A_k, and its inclusion-weighted sum, is a polynomial of degree at most d in tau**.

The adapter interpolates from d+1 admissible integral, nonnegative type histograms and checks another node when available. Synthetic inventories need not come from actual graphs. It never uses negative multiplicities. Small q may have too few admissible nodes; those classes use the actual finite inventory instead.

For d6 there are at most four classes times eight evaluations, plus three repeated-row evaluations:35. Twenty-eight distinct-class evaluations determine the polynomial; the four additional nodes are implementation checks. Exact interpolation weights are summed with the actual parameter multiplicities. The [later folding tool](../trace_invariants/) identifies two classes under row permutation and global sign change.

## Evidence and costs

[results.json](results.json) is the canonical final-source experiment record. The two `scale_p*_n*.json` files are retained first-run records from before the final input-metadata checks; they are not final-source validation records. Timings below come from the canonical run and include concurrent machine load; they are not complexity claims.

| q | n | Coefficient evaluations | Compiler seconds | Maximum intermediate bits | Standardized third moment |
|---:|---:|---:|---:|---:|---:|
|101|8|35|1.01|47|−0.0227248|
|1297|8|35|0.32|77|−0.000355452|
|65537|16|35|3.46|146|0.104722|
|1000033|31|35|3.65|181|0.0536039|

The standardized value uses the exact mean and second moment from the frozen round4 compiler; its square is saved exactly, with the sign retained separately. Character-convolution inventory construction is measured separately in [trace_backend](../trace_backend/README.md). Once the inventory is available, this compiler does not construct the q by q matrix.

At q101 the direct inventory requires only20 evaluations and is faster than interpolation. At q1297 it requires66 versus35, and both methods agree exactly. A constant evaluation budget is useful at scale; it is not automatically the fastest choice on small inputs.

The [independent review](../review/README.md) includes156 arbitrary literal-column backend cases with681 union coefficients,80 complete-oracle moment comparisons,227,388 ordered-row normalization checks, small-field fallback checks and rejection of unsupported odd-degree normalization. Complete order49 censuses cover n6,7,8.

Those actual twins also give an information-limit check. Peisert49 has the larger third moment at both n7 and n8, while its maximum is smaller at n7 (61 versus77) and larger at n8 (132 versus116). Third-moment ordering does not order even these two finite extrema.

## Reproduction and originality boundary

```
/opt/miniconda3/bin/python3 tooling_lab/round5/third_moment/run_experiments.py
/opt/miniconda3/bin/python3 tooling_lab/round5/third_moment/third_moment.py \
  tooling_lab/round5/trace_backend/inventory_p1000033.json --n 31
```

The build helper uses `clang++ -O3 -std=c++17 -I/opt/homebrew/include` and the installed Boost headers. Sources, binary and input hashes are recorded in the outputs; independent review records bind their own implementations.

Newton identities, Krawtchouk coefficients, finite-population inclusion probabilities, interpolation and cubic character sums are established tools. The [prior-art audit](../review/README.md) and [preflight ledger](../preflight/README.md) document close literature and earlier local work. The demonstrated result is a specific exact compiler, a proved finite input budget, independently checked scale results and a concrete limitation. Historical uniqueness has not been established.
