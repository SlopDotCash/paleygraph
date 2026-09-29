#!/usr/bin/env python3
"""Independent Gaussian-elimination and full subset/word review of the collision."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import random

HERE=Path(__file__).resolve().parent


def rank(vectors,p):
    pivots={}
    for vector in vectors:
        v=list(vector)
        for j,b in sorted(pivots.items()):
            a=v[j];v=[(x-a*y)%p for x,y in zip(v,b)]
        pivot=next((j for j,x in enumerate(v) if x),None)
        if pivot is not None:
            inv=pow(v[pivot],-1,p);pivots[pivot]=[x*inv%p for x in v]
            if len(pivots)==3:return 3
    return len(pivots)


def horner(poly,x,p):
    v=0
    for a in poly[::-1]:v=(v*x+a)%p
    return v


def all_subset_profile(columns,p,s):
    n=len(columns);rank_counts=Counter();minima={j:comb(j,3) for j in range(n+1)}
    triples=[t for t in combinations(range(n),3) if rank([columns[i] for i in t],p)<3]
    masks=[sum(1<<i for i in t) for t in triples]
    hist=Counter();witnesses={}
    for mask in range(1<<n):
        indices=[i for i in range(n) if mask>>i&1];j=len(indices)
        rank_counts[j,rank([columns[i] for i in indices],p)]+=1
        informative=comb(j,3)-sum(mask&t==t for t in masks)
        if informative<=minima[j]:minima[j]=informative;witnesses[j]=indices
        if j==s:hist[informative]+=1
    return {'subset_rank_distribution':[[j,r,c] for (j,r),c in sorted(rank_counts.items())],
            'minimum_by_agreement_size':[minima[j] for j in range(n+1)],
            'agreement_histogram':[[a,b] for a,b in sorted(hist.items())],
            'minimum':minima[s],'worst_agreement_indices':witnesses[s],
            'dependent_triples':triples}


def main():
    source=json.loads((HERE/'rank_three_results.json').read_text())
    witnesses=next(a['witness'] for a in source['ablations'] if a['features']=='line_sizes')
    p=source['p'];n=source['n'];s=source['s'];domain=source['domain'];rows=[];controls=[]
    rng=random.Random(11317)
    for w in witnesses:
        f=w['third_polynomial_coefficients'];columns=[(1,x,horner(f,x,p)) for x in domain]
        assert all(rank([columns[i],columns[j]],p)==2 for i,j in combinations(range(n),2))
        profile=all_subset_profile(columns,p,s);assert profile['minimum']==w['minimum_injective_triples']
        assert profile['minimum']==comb(s,3)-sum(set(t)<=set(w['worst_agreement_set']) for t in profile['dependent_triples'])
        weights=Counter(sum(sum(a*v for a,v in zip(coefficients,col))%p!=0 for col in columns) for coefficients in product(range(p),repeat=3))
        # Independently form each maximal line by span tests of its generating pair.
        lines=sorted(set(tuple(t for t in range(n) if rank([columns[i],columns[j],columns[t]],p)==2) for i,j in combinations(range(n),2)))
        assert sorted(map(len,lines))==w['features']['line_size_deck']
        assert sorted(sum(i in t for t in profile['dependent_triples']) for i in range(n))==w['features']['triple_degree_deck']
        # A simple rank-three matroid's rank-two j-subsets lie on unique lines.
        derived=Counter({(0,0):1,(1,1):n})
        for j in range(2,n+1):
            count=sum(comb(len(L),j) for L in lines if len(L)>=j)
            if count:derived[j,2]=count
            if comb(n,j)-count:derived[j,3]=comb(n,j)-count
        assert profile['subset_rank_distribution']==[[j,r,c] for (j,r),c in sorted(derived.items())]
        matrices=[((0,1,0),(0,0,1),(1,0,0)),((1,3,2),(0,1,5),(0,0,1)),((2,0,0),(0,3,0),(0,0,5))]
        for index,M in enumerate(matrices):
            assert rank(M,p)==3
            order=list(range(n));rng.shuffle(order);scales=[rng.randrange(1,p) for _ in range(n)]
            transformed=[tuple(scales[j]*sum(a*v for a,v in zip(row,columns[i]))%p for row in M) for j,i in enumerate(order)]
            before={tuple(sorted(t)) for t in profile['dependent_triples']}
            after={tuple(sorted(order[i] for i in t)) for t in combinations(range(n),3) if rank([transformed[i] for i in t],p)<3}
            assert before==after
            controls.append({'third_polynomial_coefficients':f,'combined_GL3_permutation_scaling_control':index,'all_triples_checked':comb(n,3)})
        rows.append({'third_polynomial_coefficients':f,'maximal_lines':[list(L) for L in lines],
                     'word_weight_distribution':[[a,b] for a,b in sorted(weights.items())],
                     'generalized_hamming_weights':[n-max(map(len,lines)),n-1,n],**profile})
    assert rows[0]['subset_rank_distribution']==rows[1]['subset_rank_distribution']
    assert rows[0]['word_weight_distribution']==rows[1]['word_weight_distribution']
    assert rows[0]['minimum']!=rows[1]['minimum']
    out={'status':'passed','scope':'Two actual F17 polynomial subspaces. Complete subset ranks, pinning spectra and word-weight distributions. No universal novelty claim.',
         'p':p,'n':n,'k':8,'s':s,'domain':domain,'basis_first_two':[[1],[0,1]],'rows':rows,'controls':controls,
         'all_subsets_reviewed':2*(1<<n),'all_codewords_reviewed':2*p**3,
         'same_Tutte_polynomial':True,'same_word_weight_distribution':True,'same_generalized_hamming_weights':True,
         'source_sha256':{'rank_three_review.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'rank_three_results.json':sha256((HERE/'rank_three_results.json').read_bytes()).hexdigest()}}
    (HERE/'rank_three_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','all_subsets_reviewed':out['all_subsets_reviewed'],'all_codewords_reviewed':out['all_codewords_reviewed'],
                      'minima':[r['minimum'] for r in rows],'weight_distribution':rows[0]['word_weight_distribution'],'GL3_controls':len(controls)}),flush=True)


if __name__=='__main__':main()
