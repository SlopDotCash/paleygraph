# An exact obstruction to partition-based root-count certificates

All 16 saved hard pencils are now ruled out for the **entire disjoint-partition cap-certificate family**, independently of the polynomial-discovery heuristic. The census finds the exact maximum joint agreement M, derives the smallest possible root-count cap, and exhibits an actual two-track partition attaining it. Every one of those optimal caps is at least the required agreement threshold, so the strict certificate inequality cannot hold.

This is an obstruction to a specific certificate family. It does not rule out overlapping-support arguments, other decoding methods, the prize statement, or improved structural inequalities. Concavity, interpolation and polynomial root counting are established foundations; this tool packages exact maximum-agreement evidence and a family-level failure diagnosis. No historical-first claim is made.

## Formula and proof

For n≥1, t=k−1, and an integer maximum block size1≤M≤n, the exact minimum over all integer partitions into positive parts r_j≤M is

\[
\min\sum_j\min(r_j,t)
=\left\lfloor\frac nM\right\rfloor\min(M,t)+\min(n\bmod M,t).
\]

In the relevant regime M≥t this is `floor(n/M)*(k−1)+min(n mod M,k−1)`. If M<t, every part contributes its whole size and the minimum is n. The implementation handles both cases and k1.

To prove the formula, sort the parts and consider two nonfull parts a≥b>0. Transfer one unit from b to a. The incremental cost of increasing a is no greater than the saving from decreasing b, because the discrete increments of `min(r,t)` are nonincreasing. Repeat until a reaches M or b becomes zero, and discard zero parts. Repeating leaves only full M-blocks and at most one remainder block, without increasing cost. That partition attains the displayed value.

Now suppose every degree<k polynomial matching a scalar word—or every degree<k polynomial pair jointly matching a pencil—has agreement support size≤M. Every region of any valid piece/track partition then has size≤M. The formula is therefore a lower bound for **every** such partition's root-count cap. A lower bound≥s proves that none can satisfy the required strict cap<s.

An observed candidate support of size M gives a lower bound on maximum agreement, not the upper bound needed here. `cap_barrier` rejects a candidate-only evidence label. The saved obstruction records are backed by complete maximum-agreement enumerations, or by an explicit polynomial-partition upper-bound proof.

## Exact maximum agreement without scanning field scalars

Every scalar word and every vector of words interpolates uniquely on each k-coordinate base. Hence the maximum joint agreement is at least k. Conversely, any polynomial vector with at least k matching coordinates contains one of those bases and must equal its interpolant. Enumerating all k-bases therefore finds the exact maximum for scalar words and for joint polynomial pairs.

`base_census.cpp` generates maximal support masks by barycentric interpolation and independently replays **every mask** using Newton divided differences. The two methods share only basic finite-field arithmetic and the domain. Generation and replay both enumerate bases in lexicographic order, and the driver checks the total against the exact binomial coefficient. This is exhaustive finite computation, not a candidate search or a sample-based upper bound.

The binary ledger contains, for every base and every input sample, three complete support masks: u0, u1, and their joint intersection. It therefore records all joint polynomial supports of size≥k, with repeated supports permitted. Each exact maximum also exports a concrete polynomial or polynomial-pair witness and its maximal support. No scalar enumeration or full q^k polynomial enumeration is used for the hard instances.

For these hard instances n=2k. If a maximum support has size M≥k, its complement has at most k points and can be interpolated. The maximum-support polynomial pair plus the complement interpolant forms a valid partition of sizes M and n−M. This attains the abstract minimum above. Thus the reported caps are the **actual optimal caps over polynomial-track partitions**, not merely a relaxation that might be unattainable. Scalar-word analogues are exported as well.

## Measured hard-case results

| Configuration | n,k,s | Exact joint maxima by sample | Optimal joint partition caps | Outcome |
|---|---|---|---|---|
| F41 linear-middle |20,10,14|11,11,12,11|18,18,17,18|All four blocked|
| F17 small characteristic |16,8,11|9,9,9,10|14,14,14,13|All four blocked|
| F1009 medium characteristic |16,8,11|8,8,8,8|14,14,14,14|All four blocked|
| F2147483713 large characteristic |16,8,11|8,8,8,8|14,14,14,14|All four blocked|

The complete hard-case computation covers 893,464 sample/base combinations: four times C(20,10)=184756 and twelve times C(16,8)=12870. It produces 48 exact maxima and 48 attaining partitions across u0,u1,joint modes. The support ledgers occupy about 10.7MB, including a tiny C++ control ledger.

The component and joint questions are distinct. For example, F17 sample 0 has maximum u0 agreement 13, permitting and attaining cap 10<s11 for a scalar-word certificate, while its joint maximum 9 forces cap 14≥s11. The tool does not infer a joint certificate from successful decoding of an individual component.

The partition formula matches 22140 dynamic-programming checks. The generic-field maximum oracle matches 793 complete polynomial-pair enumerations over F3, F5 and GF9, and the C++ engine separately matches eight tiny prime-field codebook cases. A corrupted support-mask ledger is rejected by the Newton replay. Exact counts and source bindings are in `results.json` and `verification.json`.

## Rate-quarter boundary

At n=4k, let s=2k+delta with 0≤delta≤k−3. The exact necessary maximum agreement for escaping this abstract obstruction is

\[
M\ge 2k-\left\lfloor\frac{\delta+1}{2}\right\rfloor.
\]

For M≤n/3, the minimum cap is at least 3(k−1)≥s. In the remaining relevant range the packed partition has two full blocks and a remainder, so the strict condition becomes `4k−2M<delta+2`, giving the formula. This is necessary for this certificate family; an actual partition need not exist merely because M is large enough.

| n | k | s | Necessary M for this cap family |
|---:|---:|---:|---:|
|1024|256|512|512|
|1024|256|513|511|
|1024|256|566|485|
|1073741824|268435456|592794966|508908885|

The last row is a parameter-level necessary condition at the saved P1 profile, not a maximum-agreement computation on an official input. It gives no official event certificate.

The earlier n1024 examples used k64, rate 1/16. Their exact joint maxima are512 and 342, with caps126 and 189. Those results cannot simply be relabeled rate-quarter results because increasing k enlarges the polynomial class. The tool explicitly rechecks both inputs at k256 using their known polynomial-track partitions. The two-track input still has exact maximum 512: an outside pair has at most510 joint matches. For the three-track input the same reasoning gives only the upper bound765; its exact maximum in the larger degree class remains unresolved. A conditional calculation with M342 at k256 has cap 765≥s566, but M342 is not proved for that larger class.

## Files and reproduction

```sh
clang++ -O3 -std=c++17 tooling_lab/round6/cover_barriers/base_census.cpp -o tooling_lab/round6/cover_barriers/base_census
/opt/miniconda3/bin/python tooling_lab/round6/cover_barriers/run_experiments.py
/opt/miniconda3/bin/python tooling_lab/round6/cover_barriers/verify_results.py
```

To reuse the unchanged support bytes while rerunning the full Newton replay and arithmetic checks, pass `--reuse-ledgers` to `run_experiments.py`. The final saved run used that mode after adding standalone prime validation and witness/provenance checks; it does not claim the ledger was freshly generated in that run.

- `cover_barriers.py`: exact partition minimum, barrier calculation, generic-field k-base oracle, and partition-derived maximum bounds.
- `base_census.cpp`, `base_census`: prime-field generator and separate Newton replay implementation.
- `*.input.txt`: exact domains and saved input words used by the executable.
- `*.supports.bin`: complete masks; little-endian uint32; lexicographic k-bases, then sample order, then u0/u1/joint.
- `*.enumerated.json`, `*.verified.json`, `*.readback.json`: exact maximum summaries from generation and independent-algorithm replay.
- `results.json`: all maximum witnesses, actual cap-attaining partitions, family obstructions and rate comparisons.
- `verify_results.py`, `verification.json`: readback binding original input files, frozen dependencies, executable/source bytes and masks, followed by full Newton replay and witness checks.

The C++ program validates primality by trial division, checks distinct canonical domain elements, and verifies all required inverses. Its supported range is n≤32, p≤2147483713 and at most one million bases. Products fit uint64 because both factors are reduced modulo p; sums are reduced after each addition. The generic Python path supports the inherited GF9 field representation.

The oracle remains exponential in n through C(n,k). The prime-field implementation makes the complete finite audit affordable; it is not a scalable official-size maximum-likelihood decoder. Large examples use a different, explicit root-count upper-bound certificate from their known polynomial partitions.

## Prior-art boundary and the next requirement

The packed-partition proof is an elementary discrete-concavity/majorization argument; [Marshall, Olkin and Arnold, Inequalities: Theory of Majorization and Its Applications](https://link.springer.com/book/10.1007/978-0-387-68276-1) is a standard reference for that framework. The short proof above is supplied directly rather than treating the formula as a fitted numerical pattern.

Maximum-agreement search is the nearest-codeword problem. [Guruswami and Vardy's maximum-likelihood decoding paper](https://www.cs.cmu.edu/~venkatg/pubs/papers/mldrs.pdf) establishes hardness for general Reed–Solomon code families; it is not being applied as a hardness result for the particular subgroup/prize parameter family. The present tool uses elementary exhaustive interpolation, a familiar decoder baseline.

Locally, round2 already enumerates interpolation bases and round4 verifies disjoint polynomial-piece caps. The existing [pencil-count charge probe](/Users/shawwalters/proximityprize/scripts/probes/probe_rate_quarter_p1_pencil_count_charge.py) studies aligned regions, interpolation pencils and per-coordinate scalar votes. The addition here is an exact upper-bound certificate for maximum joint support, followed by optimization over **all possible block sizes**, allowing arbitrarily many pieces.

For the hard inputs, better factor recovery cannot repair the partition cap: even the optimal realizable partition fails. A next approach must change the completeness argument—for example by exploiting overlapping support families—or use information beyond a disjoint partition with a k−1 root charge per block. This conclusion identifies a tooling requirement; it does not prove that any proposed replacement succeeds.
