# From an exact section to a digit language

Fix N>=2 and Q=2^N+1. Let k be a positive divisor of Q and p=Q/k. The digit-section argument works for these odd moduli without assuming primality. For the field-norm applications we additionally require N to be a power of two and p prime. In that case 2 has exact order 2N modulo p.

## Centered doubling and the endpoint

For a scalar a modulo p, let F_j be the centered residue of a*2^j, j=0,...,N-1, and define F_N=-F_0. Then

```
D_j = (2*F_j-F_(j+1))/p is in {-1,0,1},
(2-X^-1)F = pD in Z[X]/(X^N+1).
```

When a is nonzero, no F_j is zero, because 2 is invertible. A digit+1 occurs exactly when the current centered value is positive and the next is negative; a digit-1 is the opposite sign change. A zero digit preserves the sign. Thus consecutive nonzero digits alternate. The endpoint F_N=-F_0 forces an odd number of sign changes. These conditions are necessary for every nonzero actual digit word.

For Q itself, the map from its Q scalars to D is injective: the recurrence and the endpoint determine F_0 uniquely. A nonempty support of odd cardinality determines precisely two alternating digit words, one for either starting sign. There are 2*sum_(w odd) binom(N,w)=2^N such words. Including zero gives exactly Q words. Necessity plus injection and equal cardinality proves sufficiency over Q.

## The cofactor selects the actual prime section

Telescoping the recurrence gives

```
S(D) = sum_(j=0)^(N-1) 2^(N-1-j)*D_j,
(2^N+1)*F_0 = p*S(D),
F_0 = S(D)/k.
```

An alternating odd-support word is realized over Q, with first coefficient S(D). If k divides S(D), the recurrence with modulus Q shows that every coefficient of its Q-vector is divisible by k. Dividing by k gives the p-vector, with coefficients of magnitude strictly below p/2. Conversely, multiplying any centered p-vector by k gives its centered Q-vector with the same digits. Therefore the exact membership rule is:

1. The zero word represents the zero scalar.
2. A nonzero word has alternating nonzero signs and odd support.
3. Its weighted sum S(D) is divisible by k.

The decoded centered scalar is S(D)/k. No inverse matrix and no resultant is required for this special relation. Each condition is separately necessary; the review saves words satisfying the other two but failing that condition.

## A binary sign representation for independent counting

Let b_j indicate that F_j is negative. Then D_j=b_(j+1)-b_j, with b_N=1-b_0. Conversely every N-bit word b_0,...,b_(N-1) gives one nonzero alternating odd-support digit word. Its scalar numerator is

```
S = 1-(2^(N-1)+1)*b_0 + sum_(j=1)^(N-1) 2^(N-1-j)*b_j.
```

The choices b_0=0 and1 give the consecutive integer intervals[1,2^(N-1)] and[-2^(N-1),-1], respectively. With zero, S ranges once over all integers in[-2^(N-1),2^(N-1)]. As k divides the interval's cardinality Q, every residue modulo k has exactly p words. This explains the total-count check for every residue of the transfer table.

The producer counts alternating signed prefixes, updating the residue by r->2r+d, and records support size. The separate review counts binary sign words using fixed positional weights. Both use O(k*N^2) integer additions over O(k*N) count states; integer bit lengths and memory for saved layers/tables remain additional costs. They compute a distribution, not every digit word. At k641,N32 the producer executes862064 transitions and reaches36025 states at its largest layer. These are instrumentation counts, not a controlled runtime comparison.

The digit sum is also constrained: sum_j D_j=1-2b_0 is+1 or-1. In characteristic2, X^N+1=(X+1)^N when N is a power of two; D(1) is odd, so its field norm is odd. This elementary necessary condition alone does not characterize the section.

## Norms and changes of generator

For N a power of two, Norm(2-X^-1)=2^N+1=kp. Since Norm(F)=p^(N-1)*tau, the exact norm relation remains Norm(D)=k*tau. The transfer table measures support, which equals the coefficient energy of ternary D. It is not the coefficient energy of F or a norm estimate.

If the supplied generator is g=2^e with e odd, the cyclotomic automorphism X->X^e sends F constructed with g to the canonical F constructed with generator2. Suppose it also sends the supplied relation f to s*X^h*(2-X^-1), with s in{+1,-1}. It then sends the supplied D to s*X^h times the canonical digit word. Inverting this signed monomial gives a coefficient permutation with signs. The adapter checks these exact identities; it does not assume that an arbitrary short f has this form.

For the saved N64 relation at p2013265921, generator2 does not have the required order. This adapter explicitly declines that input. The general inverse checker in round17 still works there. A compact realizability interface for that general short relation remains open tooling work.
