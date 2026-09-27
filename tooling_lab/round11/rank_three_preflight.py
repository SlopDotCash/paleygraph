#!/usr/bin/env python3
"""Actual degree<8 polynomial planes with identical pair-rank data."""
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent


def det(a,b,c,p):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p


def cross(a,b,p):return ((a[1]*b[2]-a[2]*b[1])%p,(a[2]*b[0]-a[0]*b[2])%p,(a[0]*b[1]-a[1]*b[0])%p)


def normalize(v,p):
    pivot=next(x for x in v if x);inv=pow(pivot,-1,p)
    return tuple(x*inv%p for x in v)


def build(p,domain,coefficients):
    columns=[(1,x,sum(c*pow(x,j,p) for j,c in enumerate(coefficients))%p) for x in domain]
    n=len(columns);lines=defaultdict(set)
    for i,j in combinations(range(n),2):
        normal=normalize(cross(columns[i],columns[j],p),p);lines[normal].update((i,j))
    blocks=sorted(set(tuple(sorted(v)) for v in lines.values()))
    dependent=[]
    for block in blocks:
        if len(block)>=3:dependent.extend(combinations(block,3))
    assert len(dependent)==len(set(dependent))
    # Independent direct determinant check of the complete triple catalogue.
    direct=[t for t in combinations(range(n),3) if det(*(columns[i] for i in t),p)==0]
    assert sorted(dependent)==direct
    degrees=[sum(i in t for t in direct) for i in range(n)]
    incidence_deck=sorted(tuple(sorted(len(b) for b in blocks if i in b)) for i in range(n))
    return columns,blocks,direct,{'dependent_triples':len(direct),'triple_degree_deck':sorted(degrees),
                                'line_size_deck':sorted(map(len,blocks)),'vertex_line_size_deck':incidence_deck}


def main():
    p=17;domain=list(range(1,p));n=len(domain);s=11
    agreement=[sum(1<<i for i in A) for A in combinations(range(n),s)]
    polynomial_list=[]
    for t in range(2,8):
        f=[0]*(t+1);f[t]=1;polynomial_list.append(f)
        for u in range(2,t):
            for c in range(1,p):
                g=f[:];g[u]=c;polynomial_list.append(g)
    rows=[];fibres={name:defaultdict(list) for name in ('pair_ranks','counts_degrees','line_sizes','vertex_line_sizes')}
    for f in polynomial_list:
        columns,lines,triples,features=build(p,domain,f)
        masks=[sum(1<<i for i in t) for t in triples]
        scores=[sum(A&t==t for t in masks) for A in agreement]
        maximum=max(scores);witness=agreement[scores.index(maximum)]
        row={'third_polynomial_coefficients':f,'features':features,'max_dependent_triples_in_agreement':maximum,
             'minimum_injective_triples':comb(s,3)-maximum,
             'worst_agreement_set':[i for i in range(n) if witness>>i&1]}
        rows.append(row)
        keys={'pair_ranks':(),
              'counts_degrees':(features['dependent_triples'],tuple(features['triple_degree_deck'])),
              'line_sizes':(features['dependent_triples'],tuple(features['triple_degree_deck']),tuple(features['line_size_deck'])),
              'vertex_line_sizes':(features['dependent_triples'],tuple(features['triple_degree_deck']),tuple(features['line_size_deck']),tuple(features['vertex_line_size_deck']))}
        for name,key in keys.items():fibres[name][key].append(row)
    ablations=[]
    for name,groups in fibres.items():
        ambiguous=[v for v in groups.values() if len({r['minimum_injective_triples'] for r in v})>1]
        witness=None
        if ambiguous:
            group=max(ambiguous,key=lambda g:max(r['minimum_injective_triples'] for r in g)-min(r['minimum_injective_triples'] for r in g))
            witness=[min(group,key=lambda r:r['minimum_injective_triples']),max(group,key=lambda r:r['minimum_injective_triples'])]
        ablations.append({'features':name,'fibres':len(groups),'ambiguous_fibres':len(ambiguous),'witness':witness})
    out={'status':'passed','scope':'Complete eleven-set census for a specified246-family of actual polynomial subspaces over F17; not all degree<8 subspaces.',
         'p':p,'n':n,'k':8,'s':s,'domain':domain,'basis_first_two_polynomials':[[1],[0,1]],
         'polynomial_subspaces':len(rows),'agreement_sets_per_subspace':len(agreement),
         'all_pairs_have_rank_two':True,'ablations':ablations,'rows':rows,
         'source_sha256':{'rank_three_preflight.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'rank_three_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','subspaces':len(rows),'agreement_sets_per_subspace':len(agreement),
                      'ablations':[{k:v for k,v in a.items() if k!='witness'} for a in ablations]}),flush=True)


if __name__=='__main__':main()
