# Prior-art and novelty audit, September6,2026

The request asks for methods that have never been invented or tried. A literature search cannot prove that historical negative. This round claims a new implementation and new checked project results, with global novelty unestablished. Renaming established machinery would not meet the request.

| Source inspected | Relevant established ingredient | Boundary for this round |
|---|---|---|
|[Bryant, Symbolic Boolean Manipulation with Ordered Binary-Decision Diagrams, 1992](https://www.cs.cmu.edu/~bryant/pubdir/acmcs92.pdf), sections1.2–1.4|Canonical reduction for a fixed order; strong dependence of representation size on ordering, including linear/exponential examples.|Our diagrams are layered and multi-valued with no skipped levels. We prove the corresponding finite residual statement directly; ordering sensitivity is not novel.|
|[Angluin, Learning Regular Sets from Queries and Counterexamples, 1987](https://homepages.math.uic.edu/~lreyzin/papers/angluin87.pdf), sections1.1–1.2 and2.1|Membership and equivalence are different oracle obligations; observation-table rows use distinguishing suffix experiments.|Our refinement is a fixed-length finite observation-table application. We do not implement or claim the full L* learning guarantee, and lack a general efficient equivalence oracle.|
|[Stanford, Guide to Myhill–Nerode](https://web.stanford.edu/class/cs103/guide_to_myhill_nerode), distinguishing-set explanation and examples|A suffix accepted after only one of two prefixes proves different residual languages.|Expository cross-check only. The finite lower-bound argument is given explicitly in the derivation.|
|[Heuberger, Katti, Prodinger and Ruan, 2005](https://www.math.aau.at/heuberger/publications/pdf/alg1.pdf), already read in round18|Signed-digit transformations based on binary digit differences.|The generator-two language is not advertised as an invention. Round20 adds complete minimized diagrams, arbitrary-coordinate completion queries, and exact certificates for the project's orientations.|

Searches included “Angluin1987 learning regular sets queries counterexamples observation table pdf” and “algebraic number field lattice digit language ordered decision diagram cyclotomic.” The latter also returned [Construction A over number fields](https://arxiv.org/abs/1404.2904) and [rational digit systems over finite fields](https://arxiv.org/abs/1512.07824); those search abstracts do not establish identity with or novelty of this interface. No theorem from those abstracts is imported here.

The Cornell Spring2026 lexer page opened as a short JavaScript shell, so its search excerpt is not used as primary evidence. The date reported by a crawler for the Angluin PDF is not the paper's date; the paper itself is from1987.

The local inputs remain the frozen norm-compression records, round18 signed-digit language/orientation work, and round19 exact codec and obstruction results. The proof research was previously read through pass43 and the pass44 structured-level note. No proof files were rebuilt or altered during this round, and a local negative search is not a proof of historical novelty.
