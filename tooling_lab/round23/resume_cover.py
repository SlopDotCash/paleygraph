#!/usr/bin/env python3
"""Refine saved unresolved boxes while preserving every earlier derivation."""
from copy import deepcopy
from fractions import Fraction as F
from collections import Counter
from hashlib import sha256
from pathlib import Path
import argparse
import json
from box_cover import geometry,separate,propagate
from cover_review import verify

HERE=Path(__file__).resolve().parent


def refine(c,original,extra_cuts=32,extra_nodes=2000):
    verify(c,original)
    p=deepcopy(original);rows,Q,bounds=geometry(c)
    encoded=p['cuts'];cuts=[([F(*v) for v in cut['direction']],F(*cut['difference_radius'])) for cut in encoded]
    targets={tuple(cut['target']) for cut in encoded};nodes=p['nodes'];oldnodes=len(nodes);oldcuts=len(cuts)
    p['cut_budget']=len(encoded)+extra_cuts;p['node_budget']=len(nodes)+extra_nodes
    stack=list(reversed(p['unresolved_leaves']))
    while stack:
        nodeid=stack.pop();node=nodes[nodeid];box=deepcopy(node['input_box'])
        for op in node['trace']:box[op['coordinate']]=op['interval'][:]
        box,trace,contradiction=propagate(box,cuts);node['trace']+=trace
        node.pop('reason',None)
        if contradiction is not None:
            node.update(kind='excluded',witness_cut=contradiction);continue
        target=[lo if lo>0 else hi if hi<0 else 0 for lo,hi in box]
        if tuple(target) not in targets and len(encoded)<p['cut_budget']:
            proposal=separate(c,target,rows,Q);targets.add(tuple(target));encoded.append(proposal)
            if 'normalized_centered_difference' in proposal:p['continuous_witness_cuts'].append(len(encoded)-1)
            cuts.append(([F(*v) for v in proposal['direction']],F(*proposal['difference_radius'])))
            box,trace,contradiction=propagate(box,cuts);node['trace']+=trace
            if contradiction is not None:
                node.update(kind='excluded',witness_cut=contradiction);continue
        if all(lo==hi for lo,hi in box):
            node.update(kind='unresolved',reason='continuous_or_unseparated_singleton');continue
        if len(nodes)+2>p['node_budget']:
            node.update(kind='unresolved',reason='resumed_node_budget');continue
        axis=max(range(len(box)),key=lambda j:box[j][1]-box[j][0]);middle=sum(box[axis])//2
        left=deepcopy(box);right=deepcopy(box);left[axis][1]=middle;right[axis][0]=middle+1
        a=len(nodes);b=a+1
        nodes.extend([{'input_box':left,'trace':[],'kind':'unresolved'}, {'input_box':right,'trace':[],'kind':'unresolved'}])
        node.update(kind='split',axis=axis,pivot=middle,children=[a,b]);stack.extend([b,a])
        if len(nodes)%100<2:print(f'resumed nodes={len(nodes)} cuts={len(cuts)} pending={len(stack)}',flush=True)
    p['unresolved_leaves']=[i for i,n in enumerate(nodes) if n['kind']=='unresolved']
    p['universal_unique_completion']=not p['unresolved_leaves']
    p['refinement']={'old_nodes':oldnodes,'old_cuts':oldcuts,'old_unresolved_leaves':original['unresolved_leaves'],
                     'added_cuts':len(cuts)-oldcuts,'added_nodes':len(nodes)-oldnodes}
    return p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--extra-cuts',type=int,default=32)
    parser.add_argument('--extra-nodes',type=int,default=2000)
    args=parser.parse_args();path=args.source.resolve();old=json.loads(path.read_text())
    source=path.parent/old['source'];c=json.loads(source.read_text())['certificate']
    proof=refine(c,old['cover'],args.extra_cuts,args.extra_nodes)
    name=old['name']+'_refined';output=HERE/(name+'.json')
    out={'name':name,'source':old['source'],'cover':proof,
         'input_sha256':{path.name:sha256(path.read_bytes()).hexdigest(),old['source']:sha256(source.read_bytes()).hexdigest()},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['resume_cover.py','box_cover.py','cover_review.py']}}
    output.write_text(json.dumps(out,separators=(',', ':'))+'\n')
    print(json.dumps({'name':name,**proof['refinement'],'node_kinds':dict(Counter(n['kind'] for n in proof['nodes'])),
                      'unresolved_leaves':len(proof['unresolved_leaves']),'unique_completion':proof['universal_unique_completion']}),flush=True)


if __name__=='__main__':main()
