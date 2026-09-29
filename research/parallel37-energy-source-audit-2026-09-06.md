# Pass 37: unresolved energy-source discrepancy

Status: source verification in progress; no new accepted energy bound or Paley proof. The performance question was handled in [a separate diagnosis](lean-performance-2026-09-06.md).

The current arXiv version of Shkredov's *Some new inequalities in additive combinatorics* is v3 (November 6, 2012). Its PDF page 21 was rendered and visually inspected: Theorem 34 really advertises E(H) ≪ |H|^(22/9) log^(2/3)|H| for |H| ≪ p^(3/5) log^(-6/5)p. This is not an HTML extraction error. The proof uses the three-invariant-set incidence estimate in Lemma 7 on page 6. Its small-subgroup branch on pages 21–23 has size restrictions satisfied in the quartic regime.

This remains unreconciled with the later MRSS paper's 49/20 bound and historical discussion. The existence of the older statement does not by itself justify replacing the project's imported bound. No erratum was located in this pass, and no claim that the older result is false is made. The journal-version links reached through MSP's author index returned HTTP 502; the published text was not obtained.

Next: compare the published version and independently verify Lemma 7, including multiplicities when three unions of multiplicative cosets are converted into pairs of coset ratios. A straightforward reduction to a two-coset incidence estimate can repeat those pairs; that reduction needs justification. This is a concern in our attempted reconstruction, not a demonstrated error in the source.

Conditional arithmetic only: replacing 49/20 by 22/9 in the repeated-six estimate sqrt(E₂E₃), with E₃ ≪ n⁴ log n, would replace n^(129/40) log^(3/5)n by n^(29/9) log^(5/6)n. The power improvement is 1/360. This does not establish the missing general bound. Source T₃ is the project's six-variable energy; source E₃ is the cubic moment of difference multiplicities. They must not be interchanged.

Sources: [arXiv v3 PDF](https://arxiv.org/pdf/1208.2344v3), [MRSS](https://arxiv.org/html/1712.00410), [MSP author index and journal links](https://msp.org/index/ail.php?jpath=cnt&l=S). The two raw Shkredov downloads are retained in `sources/parallel37-energy-source/`; their hashes are in the companion source-state record. No central frontier bound was changed.
