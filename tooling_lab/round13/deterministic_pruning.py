#!/usr/bin/env python3
"""Execute every query in a checked deterministic coordinate-cover certificate."""
from cover_verifier import validate
from query_cover import evaluate,determinant


def solve(rows,values,p):
    d=determinant(*rows,p)
    if not d:raise ValueError('certified query lost rank')
    inv=pow(d,-1,p);result=[]
    for j in range(3):
        changed=[list(row) for row in rows]
        for i in range(3):changed[i][j]=values[i]
        result.append(determinant(*changed,p)*inv%p)
    return tuple(result)


def prepare(c):
    review=validate(c);data=c['input'];p=data['p']
    return {'data':data,'review':review,
            'columns':[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']],
            'origin_word':[evaluate(data['origin'],x,p) for x in data['domain']],
            'queries':[t for block in c['blocks'] for t in block['queries']]}


def run(prepared,received):
    data=prepared['data'];p,n,k,s=(data[x] for x in ('p','n','k','s'));origin=prepared['origin_word'];columns=prepared['columns']
    if len(received)!=n or any(type(v) is not int or not 0<=v<p for v in received):raise ValueError('canonical received word of certified length required')
    seen=set();outputs={};duplicates=0
    for j,B in enumerate(prepared['queries']):
        params=solve([columns[i] for i in B],[(received[i]-origin[i])%p for i in B],p)
        if params in seen:duplicates+=1;continue
        seen.add(params)
        word=[(a+sum(c*v for c,v in zip(params,column)))%p for a,column in zip(origin,columns)]
        support=[i for i,(a,b) in enumerate(zip(word,received)) if a==b]
        if len(support)<s:continue
        coeffs=[((data['origin'][i] if i<len(data['origin']) else 0)+sum(params[j]*(data['basis'][j][i] if i<len(data['basis'][j]) else 0) for j in range(3)))%p for i in range(k)]
        outputs[params]={'parameters':list(params),'coefficients':coeffs,'codeword':word,
                         'agreement_support':support,'first_query_index':j}
    return {'schema':'deterministic_declared_space_pruning_v1','received':list(received),
            'queries_executed':len(prepared['queries']),'distinct_candidates_evaluated':len(seen),
            'duplicate_candidates':duplicates,'outputs':[outputs[v] for v in sorted(outputs)],
            'scope':'Every nearby member of the certified affine space is returned for this received word. No probability assumption or coverage outside the space.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv)!=3:raise SystemExit('Usage: deterministic_pruning.py cover.certificate.json received.json')
    print(json.dumps(run(prepare(json.loads(Path(sys.argv[1]).read_text())),json.loads(Path(sys.argv[2]).read_text())),separators=(',',':')))
