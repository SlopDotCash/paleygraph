#!/usr/bin/env python3
"""Separate exact matrix, projection, integer-rounding and coverage review."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd,prod
from pathlib import Path
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,rational
from metric_review import check_transform


def ceil_q(v):return -((-v.numerator)//v.denominator)
def floor_q(v):return v.numerator//v.denominator


def verify_constraints(c,constraints):
    prepared=certify(c)
    if 'metric_change' in c:check_transform(c)
    N,p,f,U,V,I=c['N'],c['p'],c['relation'],c['erased'],c['visible_directions'],c['invisible_directions']
    P=[[rational(x) for x in row] for row in constraints['image_rows']]
    assert len(P)==len(U) and all(len(row)==N for row in P)
    M=sp.Matrix([[f[(i-j)%N]*(1 if i>=j else -1) for j in range(N)] for i in range(N)])
    for j,row in enumerate(P):
        H=[0]*N
        for u,x in zip(U,c['basis'][j]):H[u]=p*x
        assert M*sp.Matrix(row)==sp.Matrix(H)
    assert constraints['visible_directions']==V and constraints['invisible_directions']==I
    supports={j:[i for i,x in enumerate(P[j]) if x] for j in I}
    assert constraints['carry_supports']==[{'basis_index':j,'coordinates':supports[j]} for j in I]
    owners={i:[j for j in I if i in supports[j]] for i in range(N)}
    disjoint=all(len(js)<=1 for js in owners.values())
    assert constraints['continuous_projection_complete']==disjoint
    assert constraints['projection_mode']==('disjoint_invisible_supports' if disjoint else 'direct_rows_only')
    expected=[('fixed_coordinate',i) for i in range(N) if not owners[i]]
    if disjoint:expected += [('single_carry_elimination',j,i,h) for j in I for i,h in combinations(supports[j],2)]
    signatures=[];cuts=[]
    for cut in constraints['cuts']:
        if cut['kind']=='fixed_coordinate':
            i=cut['coordinate'];assert type(i)is int and i in range(N) and not owners[i]
            wanted=[P[j][i] for j in V];radius=F(p-1);signatures.append((cut['kind'],i))
        elif cut['kind']=='single_carry_elimination':
            j=cut['invisible_direction'];i,h=cut['coordinates'];assert disjoint and j in I and i<h and i in supports[j] and h in supports[j]
            wanted=[P[v][i]/P[j][i]-P[v][h]/P[j][h] for v in V]
            radius=F(p-1)*(1/abs(P[j][i])+1/abs(P[j][h]));signatures.append((cut['kind'],j,i,h))
        else:raise AssertionError('unknown physical cut type')
        raw=[rational(x) for x in cut['direction']];assert raw==wanted and rational(cut['radius'])==radius
        a,R=cut['integer_direction'],cut['integer_radius'];assert len(a)==len(V) and all(type(v)is int for v in a) and type(R)is int
        if any(raw):
            j=next(j for j,x in enumerate(raw) if x);scale=F(a[j])/raw[j]
            assert scale>0 and all(F(x)==scale*y for x,y in zip(a,raw)) and gcd(*a)==1
            assert R==floor_q(radius*scale)
        else:assert not any(a) and R==floor_q(radius)
        assert R>=0;cuts.append(([F(x) for x in a],F(R)))
    assert signatures==expected
    assert constraints['integer_rounding']=='clear_direction_denominators_then_divide_gcd_and_floor_radius'
    return prepared,cuts,P


def verify(c,constraints,proof):
    prepared,cuts,P=verify_constraints(c,constraints);V=c['visible_directions']
    bounds=[floor_q(2*min(rational(c['coordinate_radii'][j]),rational(c['conditioning']['cuts'][j]['radius']))) for j in V]
    assert proof['visible_directions']==V and proof['individual_difference_bounds']==bounds
    rootboxes=[]
    for pivot,b in enumerate(bounds):
        if b:rootboxes.append([[0,0] if i<pivot else [1,b] if i==pivot else [-x,x] for i,x in enumerate(bounds)])
    ids=proof['roots'];assert len(ids)==len(rootboxes) and len(set(ids))==len(ids)
    nodes=proof['nodes'];stack=list(zip(ids,rootboxes));seen=set();unresolved=[];volume=0;trace_count=0;leaves=[];singletons=[]
    while stack:
        nodeid,expected=stack.pop();assert type(nodeid)is int and nodeid in range(len(nodes)) and nodeid not in seen;seen.add(nodeid)
        node=nodes[nodeid];assert node['input_box']==expected;box=deepcopy(expected)
        for op in node['trace']:
            assert all(lo<=hi for lo,hi in box)
            k,j=op['cut'],op['coordinate'];assert type(k)is int and k in range(len(cuts));a,R=cuts[k]
            assert type(j)is int and j in range(len(V)) and a[j]
            terms=[(a[i]*box[i][0],a[i]*box[i][1]) for i in range(len(V)) if i!=j]
            low=sum(min(v) for v in terms);high=sum(max(v) for v in terms)
            if a[j]>0:L,H=(-R-high)/a[j],(R-low)/a[j]
            else:L,H=(R-low)/a[j],(-R-high)/a[j]
            wanted=[max(box[j][0],ceil_q(L)),min(box[j][1],floor_q(H))]
            assert wanted==op['interval'];box[j]=wanted;trace_count+=1
        kind=node['kind']
        if kind=='excluded':
            k=node['witness_cut'];assert type(k)is int and k in range(len(cuts))
            if all(lo<=hi for lo,hi in box):
                a,R=cuts[k];low=sum(min(v*lo,v*hi) for v,(lo,hi) in zip(a,box));high=sum(max(v*lo,v*hi) for v,(lo,hi) in zip(a,box))
                assert low>R or high<-R
        elif kind=='split':
            assert all(lo<=hi for lo,hi in box);axis,pivot=node['axis'],node['pivot']
            assert type(axis)is int and axis in range(len(V)) and type(pivot)is int and box[axis][0]<=pivot<box[axis][1]
            assert len(node['children'])==2
            left=deepcopy(box);right=deepcopy(box);left[axis][1]=pivot;right[axis][0]=pivot+1
            stack.extend(zip(node['children'],[left,right]))
        elif kind=='unresolved':
            assert all(lo<=hi for lo,hi in box);unresolved.append(nodeid);leaves.append(box)
            volume+=prod(hi-lo+1 for lo,hi in box)
            if all(lo==hi for lo,hi in box):singletons.append([lo for lo,hi in box])
        else:raise AssertionError('unknown covering state')
    assert seen==set(range(len(nodes))) and sorted(unresolved)==proof['unresolved_leaves']
    assert volume==proof['unresolved_integer_points_up_to_sign']
    assert proof['universal_unique_completion']==(not unresolved) and len(nodes)<=proof['node_budget']
    pairs=0;classes=set()
    if c['p']<100:
        words=prepared[2](list(range(c['p'])));K=c['conditioning']['known_coordinates']
        for a in range(c['p']):
            for b in range(a+1,c['p']):
                if any(words[a][j]!=words[b][j] for j in K):continue
                t=[sum(prepared[0][i][j]*(words[b][u]-words[a][u]) for i,u in enumerate(c['erased'])) for j in V]
                assert all(x.denominator==1 for x in t) and any(t)
                if next(x for x in t if x)<0:t=[-x for x in t]
                assert any(all(lo<=x<=hi for x,(lo,hi) in zip(t,leaf)) for leaf in leaves)
                pairs+=1;classes.add(tuple(t))
    return {'inverse_rows_checked':len(P),'cuts_checked':len(cuts),'nodes_checked':len(nodes),'narrowing_steps_checked':trace_count,
            'unresolved_leaves':len(unresolved),'unresolved_integer_points_up_to_sign':volume,'singleton_targets':sorted(singletons),
            'universal_unique_completion':not unresolved,'actual_small_ambiguity_pairs_preserved':pairs,'actual_difference_classes_preserved':len(classes)}


def main():
    inputs={};cases=[];sample=None;gap_checked=False
    for summaryname in ['physical_summary.json','small_summary.json']:
        path=HERE/summaryname;summary=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
        for row in summary['cases']:
            path=HERE/row['artifact'];data=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
            source=HERE/data['source'];inputs[data['source']]=sha256(source.read_bytes()).hexdigest();c=json.loads(source.read_text())['certificate']
            checked=verify(c,data['constraints'],data['cover']);result={'name':data['name'],**checked};cases.append(result)
            print(json.dumps({k:v for k,v in result.items() if k!='singleton_targets'}),flush=True)
            if data['name']=='random32':sample=(c,data)
            if data['name']=='consecutive40':
                path=HERE.parent/'round23'/'gap_certificate.json';gap=json.loads(path.read_text());inputs['../round23/'+path.name]=sha256(path.read_bytes()).hexdigest()
                target=gap['visible_target'];assert len(target)==len(c['visible_directions'])
                assert all(abs(sum(a*t for a,t in zip(cut['integer_direction'],target)))<=cut['integer_radius'] for cut in data['constraints']['cuts'])
                gap_checked=True
    rejected=[]
    for label in ['false_inverse_image','missing_elimination_pair','understated_integer_radius','false_direction']:
        c,data=sample;bad=deepcopy(data['constraints'])
        if label=='false_inverse_image':bad['image_rows'][0][0][0]+=bad['image_rows'][0][0][1]
        if label=='missing_elimination_pair':bad['cuts'].pop()
        if label=='understated_integer_radius':bad['cuts'][0]['integer_radius']-=1
        if label=='false_direction':bad['cuts'][-1]['direction'][0][0]+=bad['cuts'][-1]['direction'][0][1]
        try:verify_constraints(c,bad)
        except AssertionError:rejected.append(label)
        else:raise AssertionError('corrupt physical constraint accepted')
    out={'status':'passed','cases':cases,'integral_gap_survives_physical_projection':gap_checked,
         'corrupt_constraints_rejected':rejected,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['physical_review.py','../round22/conditioned_review.py','../round22/metric_review.py']}}
    (HERE/'physical_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
