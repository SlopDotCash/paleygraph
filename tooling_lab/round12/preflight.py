#!/usr/bin/env python3
"""Compare the corner optimizer to every literal subset on small real columns."""
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb
from pathlib import Path
import random
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'round11'))
from rank_three_verifier import rank
from capacitated_pinning import optimize

HERE=Path(__file__).resolve().parent


def main():
    rng=random.Random(120917);patterns=[]
    # Every five-column projective multiset over the F2 plane, including zero;
    # sorting is allowed because the target is coordinate-permutation invariant.
    symbols=list(product(range(2),repeat=3))
    from itertools import combinations_with_replacement
    for indices in combinations_with_replacement(range(8),5):
        columns=[symbols[i] for i in indices]
        if rank(columns,2)==3:patterns.append((2,columns))
    # Additional characteristic and multiplicity cases, including large zeros.
    for p in (3,5,7,17):
        for _ in range(30):
            columns=[tuple(rng.randrange(p) for _ in range(3)) for _ in range(8)]
            columns[0]=(0,0,0);columns[1]=columns[2]
            if rank(columns,p)==3:patterns.append((p,columns))
    checks=subsets=0;partial_witnesses=0
    for p,columns in patterns:
        n=len(columns);triples=[t for t in combinations(range(n),3) if rank([columns[i] for i in t],p)==3]
        minima=[comb(n,3)+1]*(n+1)
        for A in range(1<<n):
            value=sum(all(A>>i&1 for i in t) for t in triples);s=A.bit_count();minima[s]=min(minima[s],value);subsets+=1
        for s in range(n+1):
            out=optimize(p,columns,s);w=set(out['worst_agreement_set'])
            assert out['minimum_injective_triples']==minima[s]==sum(set(t)<=w for t in triples)
            assert sum(0<a<len(g['coordinates']) for a,g in zip(out['class_occupancies'],out['projective_classes']))<=1
            checks+=1;partial_witnesses+=out['partial_class'] is not None
    out={'status':'toy_passed','scope':'Complete projective-multiset controls plus seeded finite-field multiplicity cases. Actual polynomial-cluster replay and independent corner optimality verifier pending.',
         'patterns':len(patterns),'threshold_checks':checks,'literal_subsets':subsets,'partial_class_witnesses':partial_witnesses,
         'source_sha256':{'preflight.py':sha256(Path(__file__).read_bytes()).hexdigest(),'capacitated_pinning.py':sha256((HERE/'capacitated_pinning.py').read_bytes()).hexdigest(),
                          '../round11/rank_three_verifier.py':sha256((HERE.parent/'round11/rank_three_verifier.py').read_bytes()).hexdigest()}}
    (HERE/'preflight_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('sha256')}),flush=True)


if __name__=='__main__':main()
