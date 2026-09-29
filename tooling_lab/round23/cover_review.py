#!/usr/bin/env python3
"""Independent verification of rational cuts and exhaustive box partitions."""
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
from copy import deepcopy
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,rational
from coupled_review import setup,target_certificate
from metric_review import check_transform


def ceil_q(v):return -((-v.numerator)//v.denominator)
def floor_q(v):return v.numerator//v.denominator


def verify(c,proof):
    prepared=certify(c)
    if 'metric_change' in c:check_transform(c)
    _,V,rows,Q=setup(c)
    bounds=[floor_q(2*min(rational(c['coordinate_radii'][j]),rational(c['conditioning']['cuts'][j]['radius']))) for j in V]
    assert proof['visible_directions']==V and proof['individual_difference_bounds']==bounds
    cuts=[];continuous=[]
    for index,record in enumerate(proof['cuts']):
        norm=target_certificate(c,record,rows,Q);r=F(c['p']-1,c['p'])*norm
        assert r==rational(record['difference_radius']) and record['separated']==(r<1)
        mu=[rational(v) for v in record['direction']];cuts.append((mu,r))
        if 'normalized_centered_difference' in record:
            h=[rational(v) for v in record['normalized_centered_difference']]
            assert len(h)==c['N'] and all(abs(v)<=F(c['p']-1,c['p']) for v in h)
            assert all(sum(a*b for a,b in zip(row,h))==0 for row in rows)
            assert all(sum(a*b for a,b in zip(row,h))==t for row,t in zip(Q,record['target']))
            continuous.append(index)
    assert proof['continuous_witness_cuts']==continuous
    roots=[]
    for pivot,b in enumerate(bounds):
        if b:roots.append([[0,0] if i<pivot else [1,b] if i==pivot else [-x,x] for i,x in enumerate(bounds)])
    nodes=proof['nodes'];ids=proof['roots'];assert len(ids)==len(roots) and len(set(ids))==len(ids)
    seen=set();unresolved=[];trace_count=0;volume=0;leaves=[]
    stack=list(zip(ids,roots))
    while stack:
        nodeid,expected=stack.pop();assert type(nodeid)is int and nodeid in range(len(nodes)) and nodeid not in seen
        seen.add(nodeid);node=nodes[nodeid];assert node['input_box']==expected
        box=deepcopy(expected)
        for op in node['trace']:
            assert all(lo<=hi for lo,hi in box)
            index,j=op['cut'],op['coordinate'];assert type(index)is int and index in range(len(cuts))
            mu,r=cuts[index];assert type(j)is int and j in range(len(box)) and mu[j]
            others=[(mu[i]*box[i][0],mu[i]*box[i][1]) for i in range(len(box)) if i!=j]
            low=sum(min(pair) for pair in others);high=sum(max(pair) for pair in others)
            if mu[j]>0:lower,upper=(-r-high)/mu[j],(r-low)/mu[j]
            else:lower,upper=(r-low)/mu[j],(-r-high)/mu[j]
            interval=[max(box[j][0],ceil_q(lower)),min(box[j][1],floor_q(upper))]
            assert interval==op['interval'];box[j]=interval;trace_count+=1
        kind=node['kind']
        if kind=='excluded':
            index=node['witness_cut'];assert type(index)is int and index in range(len(cuts))
            if all(lo<=hi for lo,hi in box):
                mu,r=cuts[index];lo=sum(min(a*x,a*y) for a,(x,y) in zip(mu,box));hi=sum(max(a*x,a*y) for a,(x,y) in zip(mu,box))
                assert lo>r or hi<-r
        elif kind=='split':
            assert all(lo<=hi for lo,hi in box)
            axis,pivot=node['axis'],node['pivot'];assert type(axis)is int and axis in range(len(box))
            assert type(pivot)is int and box[axis][0]<=pivot<box[axis][1]
            assert len(node['children'])==2
            left=deepcopy(box);right=deepcopy(box);left[axis][1]=pivot;right[axis][0]=pivot+1
            stack.extend(zip(node['children'],[left,right]))
        elif kind=='unresolved':
            assert all(lo<=hi for lo,hi in box)
            unresolved.append(nodeid);count=prod(hi-lo+1 for lo,hi in box);volume+=count;leaves.append(box)
        else:raise AssertionError('unknown node kind')
    assert seen==set(range(len(nodes))) and sorted(unresolved)==proof['unresolved_leaves']
    assert proof['universal_unique_completion']==(not unresolved)
    assert len(nodes)<=proof['node_budget'] and len(cuts)<=proof['cut_budget']
    # Complete codebook check: every actual nonzero visible difference survives
    # in an unresolved leaf. This also rules out false uniqueness on toys.
    pairs=0;classes=set()
    if c['p']<100:
        words=prepared[2](list(range(c['p'])));K=c['conditioning']['known_coordinates']
        for a in range(c['p']):
            for b in range(a+1,c['p']):
                if any(words[a][j]!=words[b][j] for j in K):continue
                t=[sum(prepared[0][i][j]*(words[b][u]-words[a][u]) for i,u in enumerate(c['erased'])) for j in V]
                assert all(v.denominator==1 for v in t) and any(t);pairs+=1
                if next(v for v in t if v)<0:t=[-v for v in t]
                assert any(all(lo<=v<=hi for v,(lo,hi) in zip(t,leaf)) for leaf in leaves)
                classes.add(tuple(t))
    return {'nodes_verified':len(nodes),'cuts_verified':len(cuts),'exact_narrowing_steps':trace_count,
            'unresolved_leaves':len(unresolved),'unresolved_integer_points_up_to_sign':volume,
            'continuous_witnesses_checked':len(continuous),'universal_unique_completion':not unresolved,
            'actual_small_ambiguity_pairs_preserved':pairs,'actual_difference_classes_preserved':len(classes)}


def main():
    inputs={};cases=[];sample=None
    for summaryname in ['cover_summary.json','small_cover_summary.json']:
        path=HERE/summaryname;summary=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
        for row in summary['cases']:
            path=HERE/row['artifact'];data=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
            source=HERE/data['source'];inputs[data['source']]=sha256(source.read_bytes()).hexdigest();c=json.loads(source.read_text())['certificate']
            checked=verify(c,data['cover']);result={'name':data['name'],**checked};cases.append(result);print(json.dumps(result),flush=True)
            if data['name']=='cover32':sample=(c,data['cover'])
    rejected=[]
    for label in ['missing_root','missing_child','false_narrowing','false_uniqueness']:
        c,proof=sample;changed=deepcopy(proof)
        if label=='missing_root':changed['roots'].pop()
        if label=='missing_child':
            node=next(n for n in changed['nodes'] if n['kind']=='split');node['children'].pop()
        if label=='false_narrowing':
            op=next(n['trace'][0] for n in changed['nodes'] if n['trace']);op['interval'][0]+=1
        if label=='false_uniqueness':changed['universal_unique_completion']=not changed['universal_unique_completion']
        try:verify(c,changed)
        except AssertionError:rejected.append(label)
        else:raise AssertionError('corrupt covering accepted: '+label)
    out={'status':'passed','cases':cases,'corrupt_coverings_rejected':rejected,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['cover_review.py','../round22/coupled_review.py','../round22/conditioned_review.py','../round22/metric_review.py']}}
    (HERE/'cover_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
