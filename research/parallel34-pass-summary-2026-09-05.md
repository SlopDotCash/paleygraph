# Pass 34: the centered opposite-free moment target

**The Paley conjecture and the Proximity Prize remain unproved.**
The previous mathematical pass removed repeated coordinates from the
centered moment target. The intervening performance turn was a verified
wait: it inspected a live Lean compiler and found existing build errors,
but did not prove a new mathematical statement. This pass completes a
second reduction, separating opposite pairs without discarding their
Gaussian-scale contribution.

For A=-A subset F_p^*, n=2N, let I_r count ordered distinct zero-sum
words and O_r count those that are also opposite-free. Their correct
centerings are Q_s=I_(2s)-(n)_(2s)/p and
B_s=O_(2s)-2^(2s)(N)_(2s)/p. At 2s<=N, the ordinary proof gives

    Q_s = sum_(t=0)^s [(2s)!/(2t)!]
                       binomial(N-2t,s-t) B_t.

The inverse has alternating signs and cycle-independent-set
coefficients. Both transforms have absolute coefficients bounded by
binomial(s,t) G_s/G_t, where G_s=(2s-1)!! n^s. Thus one-sided bounds
B_t<=K^t G_t imply Q_s<=(1+K)^s G_s. Together with pass33, at
logarithmic depth in the quartic window this gives an equivalent
Gaussian hierarchy up to absolute constants. The required upper
bound on B_t is still missing. No period exponent improves.

A second ordinary argument counts the complement of repeated and
opposite positional pairs. It proves, for every symmetric nonzero set,

    O_(2s) >= [1-binomial(2s,2)(n^(-1/2)+p^(1/s)/n)] E_s.

In particular O_10>=n^6/2 whenever n>=2^35 and p<=n^4. This is a
conditional statement for every eligible field and set, including
subgroups, without asserting existence of quartic subgroups at those
orders. It shows why a raw O(n^5) target is insufficient in large
eligible cases: the principal term must be subtracted. It is not a
counterexample to Paley, the prize, or a centered moment conjecture.

The [proof](parallel34-opposite-pair-transform-2026-09-05.md) includes
both transforms, the hierarchy implications, and the degree-ten bound.
The [verifier](../experiments/parallel34_verify_2026_09_05.py) passes
9,852 exact checks across nine sets and 52 even-order cases, plus cycle
enumeration and coefficient matrices. Separate dynamic programs count
individual-element selections and opposite-class selections. The
recorded run takes about 0.23 seconds and launches no Lean compiler.
These checks validate the finite implementations; separate-author
review and formal verification of the uniform arguments remain open.

The three existing agents still report terminal usage-limit errors.
No replacements, new Lean builds, or hosted proof submissions were
started. The root agent can continue, so the goal remains active.

The next mathematical input must bound the centered opposite-free
aggregate at growing degree. The transform and the lower bound rule
out treating opposite-pair deletion as a negligible error or trying to
replace the centered high-degree target with a small raw count.
