#!/usr/bin/env python3
"""Exact integer-box covers driven by rational coupled cuts."""
from fractions import Fraction as F
from math import ceil,floor
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_erasure import discover,pair


def geometry(c):
    N,f=c['N'],c['relation']; K=c['conditioning']['known_coordinates']; V=c['visible_directions']
    rows=[[F(f[i-a] if i>=a else -f[N+i-a]) for a in range(N)] for i in K]
    Q=[[F(*v) for v in c['centered_coordinate_forms'][j]] for j in V]
    bounds=[floor(2*min(F(*c['coordinate_radii'][j]),F(*c['conditioning']['cuts'][j]['radius']))) for j in V]
    return rows,Q,bounds


def separate(c,t,rows,Q):
    pivot=next(i for i,v in enumerate(t) if v); others=[j for j in range(len(t)) if j!=pivot]
    q0=[v/F(t[pivot]) for v in Q[pivot]]
    extra=[[Q[j][a]-F(t[j],t[pivot])*Q[pivot][a] for a in range(c['N'])] for j in others]
    cut=discover(q0,rows+extra); lam=[F(*v) for v in cut['lambda']]
    mu=[F(0)]*len(t); mu[pivot]=F(1,t[pivot])
    for j,value in zip(others,lam[len(rows):]): mu[j]=-value;mu[pivot]+=F(t[j],t[pivot])*value
    norm=F(*cut['l1']); radius=F(c['p']-1,c['p'])*norm
    assert sum(x*y for x,y in zip(t,mu))==1
    out={'target':t,'pivot':pivot,'certificate':cut,'direction':[pair(v) for v in mu],
         'difference_radius':pair(radius),'separated':radius<1}
    if cut['optimality_certified'] and norm and radius>=1:
        point=[F(*v)/norm for v in cut['dual']]
        assert all(abs(v)<=F(c['p']-1,c['p']) for v in point)
        out['normalized_centered_difference']=[pair(v) for v in point]
    return out


def propagate(box,cuts):
    """Return exact narrowing trace; every trace operation has one cut witness."""
    box=[v[:] for v in box];trace=[]
    changed=True
    while changed:
        changed=False
        for index,(mu,radius) in enumerate(cuts):
            lows=[min(a*lo,a*hi) for a,(lo,hi) in zip(mu,box)]
            highs=[max(a*lo,a*hi) for a,(lo,hi) in zip(mu,box)]
            lower,upper=sum(lows),sum(highs)
            if lower>radius or upper<-radius:
                return box,trace,index
            for j,a in enumerate(mu):
                if not a:continue
                left=(-radius-(upper-highs[j]))/a
                right=(radius-(lower-lows[j]))/a
                lo=max(box[j][0],ceil(min(left,right)));hi=min(box[j][1],floor(max(left,right)))
                if [lo,hi]!=box[j]:
                    trace.append({'cut':index,'coordinate':j,'interval':[lo,hi]});box[j]=[lo,hi]
                    if lo>hi:return box,trace,index
                    changed=True;break
            if changed:break
    return box,trace,None


def cover(c,cut_budget=128,node_budget=20000):
    rows,Q,bounds=geometry(c);cuts=[];encoded=[];targets=set();nodes=[];continuous=[]
    roots=[]
    for j,b in enumerate(bounds):
        if not b:continue
        root=[[0,0] if i<j else [1,b] if i==j else [-v,v] for i,v in enumerate(bounds)]
        roots.append(root)
    stack=[(box,None,None) for box in reversed(roots)];rootids=[]
    # Each stack entry supplies its exact parent partition region.
    while stack:
        initial,parent,side=stack.pop();nodeid=len(nodes)
        node={'input_box':initial};nodes.append(node)
        if parent is None:rootids.append(nodeid)
        else:nodes[parent]['children'][side]=nodeid
        box,trace,contradiction=propagate(initial,cuts);node['trace']=trace
        if contradiction is not None:
            node.update(kind='excluded',witness_cut=contradiction);continue
        if len(nodes)+len(stack)>=node_budget:
            node.update(kind='unresolved',reason='node_budget');continue
        # The closest integer point to zero in the current box is hardest to
        # separate in a centrally symmetric convex body; roots exclude zero.
        target=[lo if lo>0 else hi if hi<0 else 0 for lo,hi in box]
        key=tuple(target)
        if key not in targets and len(encoded)<cut_budget:
            targetcut=separate(c,target,rows,Q);targets.add(key);encoded.append(targetcut)
            if 'normalized_centered_difference' in targetcut:
                continuous.append(len(encoded)-1)
            mu=[F(*v) for v in targetcut['direction']];radius=F(*targetcut['difference_radius'])
            cuts.append((mu,radius))
            narrowed,more,contradiction=propagate(box,cuts);node['trace']+=more;box=narrowed
            if contradiction is not None:
                node.update(kind='excluded',witness_cut=contradiction);continue
        if all(lo==hi for lo,hi in box):
            node.update(kind='unresolved',reason='continuous_or_unseparated_singleton');continue
        if len(nodes)+len(stack)+2>node_budget:
            node.update(kind='unresolved',reason='node_budget');continue
        if len(encoded)>=cut_budget and len(nodes)>=node_budget//2:
            node.update(kind='unresolved',reason='combined_search_budget');continue
        axis=max(range(len(box)),key=lambda j:box[j][1]-box[j][0]);lo,hi=box[axis];middle=(lo+hi)//2
        left=[v[:] for v in box];right=[v[:] for v in box]
        left[axis][1]=middle;right[axis][0]=middle+1
        node.update(kind='split',axis=axis,pivot=middle,children=[None,None])
        stack.append((right,nodeid,1));stack.append((left,nodeid,0))
        if nodeid%100==0:print(f'cover nodes={len(nodes)} cuts={len(encoded)} pending={len(stack)}',flush=True)
    unresolved=[i for i,node in enumerate(nodes) if node['kind']=='unresolved']
    return {'visible_directions':c['visible_directions'],'individual_difference_bounds':bounds,
            'cut_budget':cut_budget,'node_budget':node_budget,'cuts':encoded,'nodes':nodes,'roots':rootids,
            'continuous_witness_cuts':continuous,'unresolved_leaves':unresolved,
            'universal_unique_completion':not unresolved,
            'scope':'Exact sign-symmetric integer-box cover. Unresolved leaves are retained; a continuous point is not an actual word pair.'}
