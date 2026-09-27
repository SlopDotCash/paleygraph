#!/usr/bin/env python3
"""Coordinate pruning inside a supplied affine polynomial space.

The sampling mechanism and union bound are standard. The new lab interface uses
the checked instance certificate to choose an exact finite sampling budget.
"""
from fractions import Fraction
from math import comb
import random
from certificate_verifier import validate


def trials_for(p,probability,failure_bits,max_trials=100_000):
    delta=Fraction(*probability)
    if not (0<delta<=1 and type(failure_bits) is int and 1<=failure_bits<=128):raise ValueError('positive certified probability and 1..128 failure bits required')
    a,b=delta.numerator,delta.denominator
    def enough(t):return p**3*pow(b-a,t)*(1<<failure_bits)<=pow(b,t)
    high=1
    while not enough(high):
        high*=2
        if high>max_trials:
            high=max_trials
            if not enough(high):raise ValueError('certified probability needs more than the allowed trial budget')
            break
    low=0
    while high-low>1:
        mid=(low+high)//2
        if enough(mid):high=mid
        else:low=mid
    assert enough(high) and (high==1 or not enough(high-1))
    return {'trials':high,'failure_probability_at_most':[1,1<<failure_bits],
            'per_candidate_success_lower':[a,b],'candidate_count_upper':p**3,
            'inequality':'p^3 * (denominator-numerator)^trials * 2^failure_bits <= denominator^trials',
            'failure_bits':failure_bits,'model':'Independent uniformly distributed unordered coordinate triples; the bound covers all nearby members of the supplied affine space.'}


def determinant(a,b,c,p):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p


def solve(rows,values,p):
    d=determinant(*rows,p)
    if d==0:return None
    inv=pow(d,-1,p);answer=[]
    for j in range(3):
        changed=[list(row) for row in rows]
        for i in range(3):changed[i][j]=values[i]
        answer.append(determinant(*changed,p)*inv%p)
    return tuple(answer)


def evaluate(poly,x,p):
    v=0
    for a in poly[::-1]:v=(v*x+a)%p
    return v


def prepare(c,failure_bits=40):
    verified=validate(c);data=c['input'];p=data['p']
    columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']]
    origin_word=[evaluate(data['origin'],x,p) for x in data['domain']]
    budget=trials_for(p,verified['success_probability_lower'],failure_bits)
    return {'data':data,'columns':columns,'origin_word':origin_word,'verified':verified,'budget':budget}


def run(prepared,received,queries=None):
    data=prepared['data'];p,n,k,s=(data[x] for x in ('p','n','k','s'));columns=prepared['columns'];origin=prepared['origin_word'];budget=prepared['budget']
    if len(received)!=n or any(type(v) is not int or not 0<=v<p for v in received):raise ValueError('canonical received word of the certified length required')
    generated=queries is None
    if generated:
        rng=random.SystemRandom();queries=[sorted(rng.sample(range(n),3)) for _ in range(budget['trials'])]
    if len(queries)!=budget['trials'] or any(len(B)!=3 or B!=sorted(set(B)) or any(type(i) is not int or not 0<=i<n for i in B) for B in queries):raise ValueError('one canonical triple per certified trial required')
    seen=set();outputs={};singular=0;duplicate=0
    for j,B in enumerate(queries):
        parameters=solve([columns[i] for i in B],[(received[i]-origin[i])%p for i in B],p)
        if parameters is None:singular+=1;continue
        if parameters in seen:duplicate+=1;continue
        seen.add(parameters)
        word=[(a+sum(c*v for c,v in zip(parameters,column)))%p for a,column in zip(origin,columns)]
        support=[i for i,(a,b) in enumerate(zip(word,received)) if a==b]
        if len(support)<s:continue
        coefficients=[((data['origin'][i] if i<len(data['origin']) else 0)+sum(parameters[j]*(data['basis'][j][i] if i<len(data['basis'][j]) else 0) for j in range(3)))%p for i in range(k)]
        outputs[parameters]={'parameters':list(parameters),'coefficients':coefficients,'codeword':word,'agreement_support':support,'first_query_index':j}
    return {'schema':'declared_space_pruning_run_v1','received':list(received),'queries':queries,'sampling_budget':budget,
            'sampling_mode':'system_random' if generated else 'replay_of_provided_queries',
            'outputs':[outputs[x] for x in sorted(outputs)],'singular_queries':singular,'duplicate_candidates':duplicate,'distinct_candidates_evaluated':len(seen),
            'scope':'Output membership and agreement are checked exactly. Completeness probability is conditional on the certified space and independent uniform sampling; no coverage outside that space is asserted.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv) not in (3,4):raise SystemExit('Usage: pruning_sampler.py certificate.json received.json [failure_bits]')
    c=json.loads(Path(sys.argv[1]).read_text());received=json.loads(Path(sys.argv[2]).read_text())
    print(json.dumps(run(prepare(c,int(sys.argv[3]) if len(sys.argv)==4 else 40),received),separators=(',',':')))
