#!/usr/bin/env python3
"""Exact integer branch certificates for simple rank-three evaluation spaces.

The search is bounded and may return an interval. Its proof tree explicitly marks
unresolved leaves; a time/node limit never turns an incumbent into an optimum.
"""
from collections import defaultdict
from fractions import Fraction
from math import comb,isqrt


def c3(a):return comb(a,3) if a>=3 else 0
def c2(a):return a*(a-1)//2 if a>=2 else 0


def compile_lines(data):
    p,k,s=data['p'],data['k'],data['s'];domain=data['domain'];basis=data['basis'];n=len(domain)
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
            and type(k) is int and 3<=k<=n<=256 and n<=p and type(s) is int and 0<=s<=n):
        raise ValueError('prime field, 3<=k<=n<=256, and valid threshold required')
    if len(set(domain))!=n or any(type(x) is not int or not 0<=x<p for x in domain):raise ValueError('invalid domain')
    if len(basis)!=3 or any(len(b)>k or any(type(x) is not int or not 0<=x<p for x in b) for b in basis):raise ValueError('three degree-bounded polynomials required')
    def evaluate(poly,x):
        a=0
        for v in poly[::-1]:a=(a*x+v)%p
        return a
    columns=[tuple(evaluate(poly,x) for poly in basis) for x in domain]
    groups=defaultdict(set)
    for i in range(n):
        for j in range(i):
            a,b=columns[i],columns[j]
            normal=((a[1]*b[2]-a[2]*b[1])%p,(a[2]*b[0]-a[0]*b[2])%p,(a[0]*b[1]-a[1]*b[0])%p)
            pivot=next((x for x in normal if x),None)
            if pivot is None:raise ValueError('simple evaluation matroid required: all columns nonzero and all pairs independent')
            inv=pow(pivot,-1,p);key=tuple(x*inv%p for x in normal);groups[key].update((i,j))
    if len(groups)==1:raise ValueError('evaluation space must have rank three')
    return columns,sorted(tuple(sorted(L)) for L in groups.values() if len(L)>=3)


def optimize(n,lines,s,node_limit=200_000):
    if type(node_limit) is not int or node_limit<1:raise ValueError('positive integer node limit required')
    masks=[sum(1<<i for i in L) for L in lines];full=(1<<n)-1
    def score(S):return sum(c3((S&L).bit_count()) for L in masks)
    def bounds(S,U,r):
        current=0;line_bound=0;weights=[0]*n
        for L in masks:
            a=(S&L).bit_count();rest=U&L;b=rest.bit_count()
            current+=c3(a);line_bound+=c3(a+min(b,r))
            # Six times a fractional incidence charge: a new triple charges
            # 6/h to each of its h previously undecided members.
            h=min(b-1,r-1)
            charge=6*c2(a)+3*a*max(h,0)+2*c2(max(h,0))
            while rest:
                bit=rest&-rest;weights[bit.bit_length()-1]+=charge;rest-=bit
        eligible=[i for i in range(n) if U>>i&1]
        order=sorted(eligible,key=lambda i:(-weights[i],i))
        charge_bound=current+sum(weights[i] for i in order[:r])//6
        return min(line_bound,charge_bound,c3(s)),order
    # Feasible incumbent by reverse greedy and single-coordinate improvements.
    S=full
    for _ in range(n-s):
        candidates=[i for i in range(n) if S>>i&1]
        i=min(candidates,key=lambda i:(score(S)-score(S^(1<<i)),i));S^=1<<i
    best=score(S);witness=S
    improving=True
    while improving:
        improving=False
        for i in range(n):
            if not S>>i&1:continue
            for j in range(n):
                if S>>j&1:continue
                T=S^(1<<i)^(1<<j);v=score(T)
                if v>best:S=T;witness=T;best=v;improving=True;break
            if improving:break
    tree=[];stats={'branch_nodes':0,'bound_leaves':0,'exact_leaves':0,'unresolved_leaves':0};frontier=0;nodes=0
    def visit(S,U,r):
        nonlocal best,witness,frontier,nodes
        nodes+=1
        if r==0 or r==U.bit_count():
            T=S if r==0 else S|U;value=score(T)
            if value>best:best=value;witness=T
            tree.append(-2);stats['exact_leaves']+=1;return
        upper,order=bounds(S,U,r)
        if upper<=best:tree.append(-1);stats['bound_leaves']+=1;return
        if nodes>=node_limit:
            tree.append(-3);stats['unresolved_leaves']+=1;frontier=max(frontier,upper);return
        i=order[0];bit=1<<i;tree.append(i);stats['branch_nodes']+=1
        visit(S|bit,U^bit,r-1);visit(S,U^bit,r)
    visit(0,full,s)
    upper=max(best,frontier);minimum_lower=c3(s)-upper;minimum_upper=c3(s)-best
    probability=[Fraction(v,comb(n,3)) for v in (minimum_lower,minimum_upper)]
    return {'status':'exact' if upper==best else 'bounded','minimum_injective_triples_interval':[minimum_lower,minimum_upper],
            'uniform_triple_success_probability_interval':[[v.numerator,v.denominator] for v in probability],
            'max_dependent_triples_interval':[best,upper],'worst_known_agreement_set':[i for i in range(n) if witness>>i&1],
            'tree_encoding':'preorder: nonnegative=branch coordinate/include child then exclude child; -1=upper-bound leaf; -2=unique completion; -3=unresolved frontier',
            'proof_tree':tree,'node_limit':node_limit,'visited_nodes':nodes,**stats}


def certify(data,node_limit=200_000):
    columns,lines=compile_lines(data)
    return {'schema':'simple_rank_three_pinning_v1','input':data,'nontrivial_maximal_lines':[list(L) for L in lines],
            'scope':'Declared simple rank-three polynomial evaluation space. No cluster-discovery or coverage claim. Bounds retain unresolved search frontiers.',
            'pinning':optimize(len(columns),lines,data['s'],node_limit)}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv) not in (2,3):raise SystemExit('Usage: rank_three_pinning.py input.json [node_limit]')
    data=json.loads(Path(sys.argv[1]).read_text())
    print(json.dumps(certify(data,int(sys.argv[2]) if len(sys.argv)==3 else 200_000),separators=(',',':')))
