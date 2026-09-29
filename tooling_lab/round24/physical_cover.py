#!/usr/bin/env python3
"""Exact box propagation on primitive integer physical-coordinate cuts."""
from copy import deepcopy
from fractions import Fraction as F
from math import floor,prod


def propagate(box,cuts):
    box=deepcopy(box);trace=[]
    while True:
        changed=False
        for index,cut in enumerate(cuts):
            a,R=cut['integer_direction'],cut['integer_radius']
            mins=[min(v*lo,v*hi) for v,(lo,hi) in zip(a,box)]
            maxs=[max(v*lo,v*hi) for v,(lo,hi) in zip(a,box)]
            low,high=sum(mins),sum(maxs)
            if low>R or high<-R:return box,trace,index
            for j,v in enumerate(a):
                if not v:continue
                left=-R-(high-maxs[j]);right=R-(low-mins[j])
                if v<0:left,right=right,left
                lo=max(box[j][0],-((-left)//v));hi=min(box[j][1],right//v)
                if [lo,hi]!=box[j]:
                    trace.append({'cut':index,'coordinate':j,'interval':[lo,hi]});box[j]=[lo,hi]
                    if lo>hi:return box,trace,index
                    changed=True;break
            if changed:break
        if not changed:return box,trace,None


def cover(c,constraints,node_budget=4000):
    V=c['visible_directions'];bounds=[floor(2*min(F(*c['coordinate_radii'][j]),F(*c['conditioning']['cuts'][j]['radius']))) for j in V]
    roots=[]
    for pivot,b in enumerate(bounds):
        if b:roots.append([[0,0] if i<pivot else [1,b] if i==pivot else [-x,x] for i,x in enumerate(bounds)])
    if node_budget<len(roots):raise ValueError('budget cannot hold the roots')
    nodes=[];rootids=[];stack=[(box,None,None) for box in reversed(roots)]
    while stack:
        initial,parent,side=stack.pop();nodeid=len(nodes);node={'input_box':initial};nodes.append(node)
        if parent is None:rootids.append(nodeid)
        else:nodes[parent]['children'][side]=nodeid
        box,trace,witness=propagate(initial,constraints['cuts']);node['trace']=trace
        if witness is not None:node.update(kind='excluded',witness_cut=witness);continue
        if all(lo==hi for lo,hi in box):node.update(kind='unresolved',reason='integer_point_in_projected_body');continue
        if len(nodes)+len(stack)+2>node_budget:node.update(kind='unresolved',reason='node_budget');continue
        axis=max(range(len(box)),key=lambda j:box[j][1]-box[j][0]);middle=sum(box[axis])//2
        left=deepcopy(box);right=deepcopy(box);left[axis][1]=middle;right[axis][0]=middle+1
        node.update(kind='split',axis=axis,pivot=middle,children=[None,None])
        stack.extend([(right,nodeid,1),(left,nodeid,0)])
    unresolved=[i for i,n in enumerate(nodes) if n['kind']=='unresolved'];volume=0
    for i in unresolved:
        box=deepcopy(nodes[i]['input_box'])
        for op in nodes[i]['trace']:box[op['coordinate']]=op['interval']
        volume+=prod(hi-lo+1 for lo,hi in box)
    return {'visible_directions':V,'individual_difference_bounds':bounds,'roots':rootids,'nodes':nodes,
            'node_budget':node_budget,'unresolved_leaves':unresolved,'unresolved_integer_points_up_to_sign':volume,
            'universal_unique_completion':not unresolved}
