Let $S$ and $P$ be finite sets with a relation $R$, and let integers $q\ge2$, $t\ge1$, and $0\le K<q$ satisfy $|P|+1=qt$. Assume that each element of $S$ is related to at least $(q-1)t$ elements of $P$, while each element of $P$ is related to at most $K$ elements of $S$. Then

$$|S|\le K.$$

Count the related pairs in the two orders to obtain

$$|S|(q-1)t\le |P|K=(qt-1)K.$$

If $|S|\ge K+1$, the proposed lower bound exceeds the upper bound by

$$(K+1)(q-1)t-(qt-1)K=(q-1)t-K(t-1)\ge q-1>0.$$

This is a contradiction. The strict gap includes the endpoint $K=q-1$ and the case $t=1$.

This is the finite counting component of a nonzero-row-projection argument for interleaved linear codes. The proof treats the relation and its degree hypotheses abstractly. It does not formalize finite-field hyperplane counts, quotient witnesses, code maxima, the full interleaving theorem, or either open conjecture.
