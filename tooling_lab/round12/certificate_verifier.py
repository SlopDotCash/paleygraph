#!/usr/bin/env python3
"""Constructor-free certificate verifier, using row reduction and occupancy search.

The exact corner certificate is checked by enumerating ALL feasible integer
occupancies, with nonnegative partial-count pruning; no corner theorem, NumPy or
subset transform is imported or assumed by that numerical optimality check.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb,isqrt


class VerificationLimit(RuntimeError):pass


def rref(vectors,p):
    rows=[list(v) for v in vectors];r=0
    for j in range(3):
        found=next((i for i in range(r,len(rows)) if rows[i][j]%p),None)
        if found is None:continue
        rows[r],rows[found]=rows[found],rows[r]
        inv=pow(rows[r][j]%p,-1,p);rows[r]=[x*inv%p for x in rows[r]]
        for i in range(len(rows)):
            if i==r:continue
            a=rows[i][j];rows[i]=[(x-a*y)%p for x,y in zip(rows[i],rows[r])]
        r+=1
        if r==3:break
    return tuple(tuple(row) for row in rows[:r])


def inputs(c):
    assert c['schema']=='polynomial_rank_three_pinning_v1';data=c['input']
    p,n,k,s=(data[x] for x in ('p','n','k','s'));domain=data['domain'];basis=data['basis'];origin=data['origin']
    assert type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
    assert type(n) is int and type(k) is int and 3<=k<=n<=min(p,4096) and type(s) is int and 0<=s<=n
    assert len(domain)==n and len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==3 and all(len(poly)<=k and all(type(a) is int and 0<=a<p for a in poly) for poly in [origin,*basis])
    columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(poly))%p for poly in basis) for x in domain]
    span=[]
    for v in columns:
        reduced=rref(span+[v],p)
        if len(reduced)>len(span):span=list(reduced)
        if len(span)==3:break
    assert len(span)==3
    return data,columns


def verify_capacity(data,columns,out,max_states):
    p,n,s=data['p'],data['n'],data['s'];assert out['schema']=='rank_three_capacitated_pinning_v1'
    assert (out['p'],out['n'],out['s'])==(p,n,s)
    zero=[i for i,v in enumerate(columns) if not any(v)];nonzero=[i for i in range(n) if i not in zero]
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in combinations(nonzero,2):
        if len(rref([columns[i],columns[j]],p))==1:parent[find(i)]=find(j)
    groups=defaultdict(list)
    for i in nonzero:groups[find(i)].append(i)
    exported=out['projective_classes'];assert out['zero_coordinates']==zero
    assert sorted(tuple(g['coordinates']) for g in exported)==sorted(tuple(g) for g in groups.values())
    assert all(g['coordinates']==sorted(set(g['coordinates'])) for g in exported)
    for g in exported:
        direction=g['direction'];assert len(direction)==3 and all(type(x) is int and 0<=x<p for x in direction) and any(direction)
        assert next(x for x in direction if x)==1
        assert len(rref([direction,columns[g['coordinates'][0]]],p))==1
    reps=[columns[g['coordinates'][0]] for g in exported];m=len(reps);assert 3<=m<=22
    triples=[list(t) for t in combinations(range(m),3) if len(rref([reps[i] for i in t],p))==3]
    assert out['independent_class_triples']==triples
    capacities=[len(g['coordinates']) for g in exported];target=out['minimum_injective_triples'];assert type(target) is int and 0<=target<=comb(s,3)
    witness=out['worst_agreement_set'];assert witness==sorted(set(witness)) and len(witness)==s and all(type(i) is int and 0<=i<n for i in witness)
    counts=[sum(i in witness for i in g['coordinates']) for g in exported];assert counts==out['class_occupancies']
    assert sum(counts[i]*counts[j]*counts[k] for i,j,k in triples)==target
    # Arbitrary integer occupancies, including multiple partially filled classes.
    # Nonnegative remaining monomials can only increase a partial assignment.
    suffix=[0]*(m+1)
    for i in range(m-1,-1,-1):suffix[i]=suffix[i+1]+capacities[i]
    by_last=[[] for _ in range(m)]
    for i,j,k in triples:by_last[k].append((i,j))
    assigned=[0]*m;visited=0;leaves=0
    def visit(i,needed,value):
        nonlocal visited,leaves
        visited+=1
        if visited>max_states:raise VerificationLimit('occupancy verification budget exhausted; certificate is not verified')
        if value>=target:return
        if i==m:
            if needed==0:raise AssertionError('found an occupancy below the claimed minimum')
            return
        low=max(0,needed-suffix[i+1]);high=min(capacities[i],needed)
        coefficient=sum(assigned[a]*assigned[b] for a,b in by_last[i])
        for a in range(low,high+1):
            assigned[i]=a;visit(i+1,needed-a,value+a*coefficient)
        assigned[i]=0
    for zeros_chosen in range(min(s,len(zero))+1):
        remaining=s-zeros_chosen
        if remaining<=suffix[0]:visit(0,remaining,0)
    prob=Fraction(target,comb(n,3));assert out['uniform_triple_success_probability']==[prob.numerator,prob.denominator]
    partial=[[i,a] for i,a in enumerate(counts) if 0<a<capacities[i]]
    assert len(partial)<=1 and out['partial_class']==(partial[0] if partial else None)
    return {'status':'verified','method':'capacitated_corners','minimum_interval':[target,target],
            'nonzero_classes':m,'all_nonzero_pairs_checked':comb(len(nonzero),2),'all_class_triples_checked':comb(m,3),
            'integer_occupancy_states_checked':visited,'success_probability_lower':[prob.numerator,prob.denominator]}


def verify_pairs(data,columns,out):
    p,n,s=data['p'],data['n'],data['s'];assert n<=1024
    # A separately implemented grouping: store every generating pair, then
    # reconstruct each flat from the union of endpoints of its row-space class.
    pair_groups=defaultdict(list)
    for i in range(n):
        for j in range(i):
            key=rref([columns[i],columns[j]],p);assert len(key)==2
            pair_groups[key].append((j,i))
    lines=[];ordinary=0
    for pairs in pair_groups.values():
        if len(pairs)==1:ordinary+=1;continue
        vertices=sorted({v for pair in pairs for v in pair})
        assert len(pairs)==comb(len(vertices),2)
        lines.append(vertices)
    lines.sort()
    assert lines==out['nontrivial_maximal_lines'] and ordinary==out['ordinary_two_point_lines']
    assert ordinary+sum(comb(len(L),2) for L in lines)==comb(n,2)
    degrees=[0]*n;total=0
    for L in lines:
        total+=comb(len(L),3)
        contribution=(len(L)-1)*(len(L)-2)//2
        for i in L:degrees[i]+=contribution
    b=out['bound'];assert b['total_dependent_triples']==total and b['dependency_vertex_degrees']==degrees
    witness=b['worst_known_agreement_set'];assert witness==sorted(set(witness)) and len(witness)==s and all(type(i) is int and 0<=i<n for i in witness)
    selected=set(witness);attained=0
    for L in lines:
        size=len(selected.intersection(L))
        if size>=3:attained+=size*(size-1)*(size-2)//6
    upper=min(comb(s,3),total,sum(sorted(degrees,reverse=True)[:s])//3)
    assert attained<=upper and b['maximum_dependent_triples_interval']==[attained,upper]
    interval=[comb(s,3)-upper,comb(s,3)-attained]
    assert b['minimum_injective_triples_interval']==interval and b['status']==('exact' if attained==upper else 'bounded')
    probs=[Fraction(v,comb(n,3)) for v in interval]
    assert b['uniform_triple_success_probability_interval']==[[v.numerator,v.denominator] for v in probs]
    return {'status':'verified','method':'simple_pair_degree','minimum_interval':interval,'all_pairs_checked':comb(n,2),
            'all_triple_dependencies_covered':comb(n,3),'literal_triples_enumerated':False,'nontrivial_lines':len(lines),
            'success_probability_lower':[probs[0].numerator,probs[0].denominator]}


def validate(c,max_states=5_000_000):
    data,columns=inputs(c)
    if c['method']=='capacitated_corners':result=verify_capacity(data,columns,c['proof'],max_states)
    else:
        assert c['method']=='simple_pair_degree';result=verify_pairs(data,columns,c['proof'])
    return {'p':data['p'],'n':data['n'],'k':data['k'],'s':data['s'],**result}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    print(json.dumps(validate(json.loads(Path(sys.argv[1]).read_text())),indent=2))
