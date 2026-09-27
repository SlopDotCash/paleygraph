#!/usr/bin/env python3
"""Exact joint-cut review and exhaustive coverage of visible differences."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from math import floor
from pathlib import Path
import json
import sympy as sp
from conditioned_review import certify, rational
from metric_review import check_transform

HERE = Path(__file__).resolve().parent


def norm_certificate(q, rows, cut):
    N = len(q); h = len(rows); M = sp.Matrix(h,N,lambda i,j:rows[i][j])
    lam = [rational(v) for v in cut['lambda']]; assert len(lam)==h
    r = [rational(v) for v in cut['residual']]; assert len(r)==N
    assert sp.Matrix(q)-M.T*sp.Matrix(h,1,lam)==sp.Matrix(r)
    value = sum(map(abs,r)); assert value==rational(cut['l1'])
    if cut['optimality_certified']:
        eta=[rational(v) for v in cut['dual']]; assert len(eta)==N and all(abs(v)<=1 for v in eta)
        assert M*sp.Matrix(eta)==sp.zeros(h,1)
        assert sum(x*y for x,y in zip(q,eta))==value
    else: assert cut['dual'] is None
    return value


def setup(c):
    N,f,U=c['N'],c['relation'],c['erased']; K=c['conditioning']['known_coordinates']; V=c['visible_directions']
    M=sp.Matrix([[f[(i-j)%N]*(1 if i>=j else -1) for j in range(N)] for i in range(N)])
    T=sp.Matrix([[rational(v) for v in row] for row in c['inverse']])
    forms=M[U,:].T*T
    return K,V,[[F(v) for v in M[i,:]] for i in K],[[F(v) for v in forms[:,j]] for j in V]


def target_certificate(c, record, rows, Q):
    t=record['target']; pivot=record['pivot']; assert len(t)==len(Q) and t[pivot]!=0
    others=[j for j in range(len(t)) if j!=pivot]; N=c['N']
    q0=[v/F(t[pivot]) for v in Q[pivot]]
    extra=[[Q[j][a]-F(t[j],t[pivot])*Q[pivot][a] for a in range(N)] for j in others]
    norm=norm_certificate(q0,rows+extra,record['certificate'])
    mu=[rational(v) for v in record['direction']]; assert len(mu)==len(t)
    lam=[rational(v) for v in record['certificate']['lambda']]
    assert sum(a*b for a,b in zip(mu,t))==1
    r=[sum(mu[j]*Q[j][a] for j in range(len(t)))-sum(lam[i]*rows[i][a] for i in range(len(rows))) for a in range(N)]
    assert r==[rational(v) for v in record['certificate']['residual']]
    return norm


def main():
    source=HERE/'metric_consecutive32.json'; c=json.loads(source.read_text())['certificate']
    certify(c); check_transform(c); K,V,rows,Q=setup(c)
    path=HERE/'coupled_summary.json'; data=json.loads(path.read_text())
    assert data['visible_directions']==V
    bounds=[floor(2*min(rational(c['coordinate_radii'][j]),rational(c['conditioning']['cuts'][j]['radius']))) for j in V]
    assert bounds==data['individual_difference_bounds']
    tuples=[()]
    for b in bounds: tuples=[t+(x,) for t in tuples for x in range(-b,b+1)]
    nonzero={t for t in tuples if any(t)}
    candidates={t for t in nonzero if next(v for v in t if v)>0}
    assert len(nonzero)==2*len(candidates) and len(candidates)==data['initial_sign_representatives']
    initial=len(candidates); pair_optimal=0
    for item in data['pair_cuts']:
        direction=item['direction']; assert len(direction)==len(V) and all(type(v) is int for v in direction)
        q=[sum(a*Q[j][i] for j,a in enumerate(direction)) for i in range(c['N'])]
        norm=norm_certificate(q,rows,item['certificate']); bound=floor(F(c['p']-1,c['p'])*norm)
        assert bound==item['integer_difference_bound']; before=len(candidates)
        candidates={t for t in candidates if abs(sum(a*b for a,b in zip(t,direction)))<=bound}
        assert before-len(candidates)==item['candidates_removed']
        pair_optimal+=item['certificate']['optimality_certified']
    assert len(candidates)==data['remaining_after_pairs']; after=len(candidates)
    separated=set(); feasible=set(); target_optimal=0
    for r in data['target_separators']:
        t=tuple(r['target']); assert t in candidates and t not in separated
        norm=target_certificate(c,r,rows,Q)
        assert F(c['p']-1,c['p'])*norm<1
        separated.add(t); target_optimal+=r['certificate']['optimality_certified']
    for r in data['continuous_feasible_targets']:
        t=tuple(r['target']); assert t in candidates and t not in separated|feasible
        norm=target_certificate(c,r,rows,Q)
        point=[rational(v) for v in r['normalized_centered_difference']]
        assert len(point)==c['N'] and all(abs(v)<=F(c['p']-1,c['p']) for v in point)
        assert all(sum(a*b for a,b in zip(row,point))==0 for row in rows)
        assert all(sum(a*b for a,b in zip(row,point))==target for row,target in zip(Q,t))
        feasible.add(t); target_optimal+=r['certificate']['optimality_certified']
    unresolved={tuple(t) for t in data['unresolved_targets']}
    assert len(unresolved)==len(data['unresolved_targets']) and not unresolved&(separated|feasible)
    assert separated|feasible|unresolved==candidates
    unique=not feasible and not unresolved
    assert data['universal_unique_completion']==unique
    rejected=[]; sample=data['target_separators'][0]
    for label in ['false_target','false_direction','false_joint_residual']:
        changed=deepcopy(sample)
        if label=='false_target': changed['target'][changed['pivot']]+=1
        if label=='false_direction': changed['direction'][0][0]+=changed['direction'][0][1]
        if label=='false_joint_residual': changed['certificate']['residual'][0][0]+=changed['certificate']['residual'][0][1]
        try: target_certificate(c,changed,rows,Q)
        except AssertionError: rejected.append(label)
        else: raise AssertionError('corrupted joint cut accepted')
    # Deleting one certificate must leave its direction uncovered.
    assert {tuple(r['target']) for r in data['target_separators'][1:]}|feasible|unresolved!=candidates
    rejected.append('missing_coverage_certificate')
    out={'status':'passed','initial_nonzero_differences':len(nonzero),'sign_representatives':initial,
         'pair_cuts_checked':len(data['pair_cuts']),'pair_optimality_certificates':pair_optimal,
         'remaining_after_pair_cuts':after,'target_separators_checked':len(separated),
         'target_optimality_certificates':target_optimal,'continuous_feasible_targets':len(feasible),
         'unresolved_targets':len(unresolved),'universal_unique_completion':unique,
         'corrupt_controls_rejected':rejected,
         'input_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in [source,path]},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['coupled_review.py','conditioned_review.py','metric_review.py']}}
    (HERE/'coupled_review.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out),flush=True)


if __name__=='__main__': main()
