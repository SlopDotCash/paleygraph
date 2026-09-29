#!/usr/bin/env python3
"""Independent certificate replay; no constructor or optimizer imports."""
from fractions import Fraction
from itertools import combinations
from math import comb,isqrt


def rank(vectors,p):
    # Row reduction of a 3 by number-of-vectors matrix.
    matrix=[list(row) for row in zip(*vectors)] if vectors else [[],[],[]]
    pivot=0
    for j in range(len(vectors)):
        found=next((i for i in range(pivot,3) if matrix[i][j]%p),None)
        if found is None:continue
        matrix[pivot],matrix[found]=matrix[found],matrix[pivot]
        inv=pow(matrix[pivot][j]%p,-1,p);matrix[pivot]=[v*inv%p for v in matrix[pivot]]
        for i in range(pivot+1,3):
            a=matrix[i][j];matrix[i]=[(x-a*y)%p for x,y in zip(matrix[i],matrix[pivot])]
        pivot+=1
        if pivot==3:break
    return pivot


def validate(c):
    assert c['schema']=='simple_rank_three_pinning_v1'
    data=c['input'];p,k,s=data['p'],data['k'],data['s'];domain=data['domain'];basis=data['basis'];n=len(domain)
    assert type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
    assert type(k) is int and 3<=k<=n<=256 and n<=p and type(s) is int and 0<=s<=n
    assert len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==3 and all(len(b)<=k and all(type(x) is int and 0<=x<p for x in b) for b in basis)
    columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(b))%p for b in basis) for x in domain]
    assert all(rank([columns[i],columns[j]],p)==2 for i,j in combinations(range(n),2))
    assert rank(columns,p)==3
    lines=c['nontrivial_maximal_lines'];assert lines==sorted(lines)
    assert all(L==sorted(set(L)) and len(L)>=3 and all(type(i) is int and 0<=i<n for i in L) for L in lines)
    triple_set=set();pair_lines={}
    for L in lines:
        assert rank([columns[i] for i in L],p)==2
        assert all(rank([columns[L[0]],columns[L[1]],columns[j]],p)==3 for j in range(n) if j not in L)
        for pair in combinations(L,2):assert pair not in pair_lines;pair_lines[pair]=L
        for triple in combinations(L,3):assert triple not in triple_set;triple_set.add(triple)
    direct={t for t in combinations(range(n),3) if rank([columns[i] for i in t],p)==2}
    assert direct==triple_set
    r=c['pinning'];witness=r['worst_known_agreement_set'];assert witness==sorted(set(witness)) and len(witness)==s and all(0<=i<n for i in witness)
    incumbent=sum(set(t)<=set(witness) for t in direct)
    assert r['max_dependent_triples_interval'][0]==incumbent
    tree=r['proof_tree'];position=0;stats={'branch_nodes':0,'bound_leaves':0,'exact_leaves':0,'unresolved_leaves':0}
    frontier=incumbent;line_sets=[set(L) for L in lines];all_points=set(range(n))
    def upper_bound(chosen,available,needed):
        already=sum(len(chosen.intersection(t))==3 for t in direct)
        separate=0;charge={v:0 for v in available}
        for L in line_sets:
            fixed=len(chosen&L);free=L&available;maximum=fixed+min(len(free),needed)
            separate+=comb(maximum,3) if maximum>=3 else 0
            for v in free:
                others=min(len(free-{v}),needed-1)
                charge[v]+=3*fixed*(fixed-1)+3*fixed*others+others*(others-1)
        charged=already+sum(sorted(charge.values(),reverse=True)[:needed])//6
        return min(separate,charged,comb(s,3))
    def visit(chosen,available):
        nonlocal position,frontier
        assert position<len(tree);entry=tree[position];position+=1
        assert type(entry) is int
        needed=s-len(chosen);assert 0<=needed<=len(available)
        if entry==-2:
            assert needed==0 or needed==len(available)
            full=chosen if needed==0 else chosen|available
            assert sum(set(t)<=full for t in direct)<=incumbent
            stats['exact_leaves']+=1
        elif entry in (-1,-3):
            assert 0<needed<len(available)
            bound=upper_bound(chosen,available,needed)
            if entry==-1:assert bound<=incumbent;stats['bound_leaves']+=1
            else:frontier=max(frontier,bound);stats['unresolved_leaves']+=1
        else:
            assert entry in available and 0<needed<len(available)
            stats['branch_nodes']+=1
            visit(chosen|{entry},available-{entry});visit(chosen,available-{entry})
    visit(set(),all_points);assert position==len(tree)==r['visited_nodes']
    assert all(r[k]==v for k,v in stats.items())
    assert r['max_dependent_triples_interval']==[incumbent,frontier]
    interval=[comb(s,3)-frontier,comb(s,3)-incumbent]
    assert r['minimum_injective_triples_interval']==interval
    assert r['status']==('exact' if frontier==incumbent else 'bounded')
    fractions=[Fraction(x,comb(n,3)) for x in interval]
    assert r['uniform_triple_success_probability_interval']==[[x.numerator,x.denominator] for x in fractions]
    return {'p':p,'n':n,'k':k,'s':s,'status':r['status'],'interval':interval,'all_triples_checked':comb(n,3),
            'nontrivial_lines':len(lines),'tree_nodes_replayed':position,**stats}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    print(json.dumps(validate(json.loads(Path(sys.argv[1]).read_text())),indent=2))
