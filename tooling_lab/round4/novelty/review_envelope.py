#!/usr/bin/env python3
"""Independent centered norm bound and rational square-root rounding review."""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import json
import random

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'marked_moments'/'certified_envelope.py'
sys.path.insert(0,str(SOURCE.parent))


def load(path):
    spec=spec_from_file_location(path.stem,path)
    module=module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    candidate=load(SOURCE)
    rng=random.Random(50138)
    values=[Fraction(0),Fraction(1),Fraction(4,9),Fraction(1,10**80),Fraction(10**80),Fraction(9999,10000)]
    values.extend(Fraction(rng.randrange(10**80),rng.randrange(1,10**70)) for _ in range(100))
    rounded=0
    for value in values:
        for digits in (0,1,5,20):
            upper=candidate.rational_sqrt_upper(value,digits)
            assert upper>=0 and upper*upper>=value
            step=Fraction(1,10**digits)
            assert upper==0 or (upper-step)**2<value
            rounded+=1
    ledger_path=HERE.parent/'marked_moments'/'contraction_envelopes.json'
    ledger=json.loads(ledger_path.read_text())
    for filename,digest in ledger['source_sha256'].items():
        assert sha256((SOURCE.parent/filename).read_bytes()).hexdigest()==digest
    cases=[]
    for row in ledger['cases']:
        q=row['q']
        stats=row['vector_statistics']
        # Compute centered norms as fractions first; this rearrangement avoids
        # copying the numerator-product implementation used by the candidate.
        centered2=stats['h2_norm_squared']-Fraction(stats['h2_sum']**2,q)
        centered3=stats['h3_norm_squared']-Fraction(stats['h3_sum']**2,q)
        bound=q*centered2*centered3
        assert bound==Fraction(*row['Q_bound_squared']) and row['Q']**2<=bound
        c=Fraction(*row['Q_coefficient'])
        radius2=c*c*bound
        assert radius2==Fraction(*row['cell_only_envelope_radius_squared'])
        upper=Fraction(*row['envelope_radius_rational_upper'])
        assert upper*upper>=radius2
        base=Fraction(*row['cell_only_term'])
        actual=Fraction(*row['second_moment'])
        assert (actual-base)**2<=radius2
        if base:
            relative=Fraction(*row['relative_radius_rational_upper'])
            assert relative*relative>=radius2/(base*base)
        else:
            assert row['relative_radius_rational_upper'] is None
        if q in (65537,1000033):
            threshold=Fraction(34,10**13) if q==65537 else Fraction(12,10**16)
            assert relative<threshold
        cases.append({'q':q,'n':row['n'],'marks':row['marks'],
                      'relative_radius_rational_upper':row['relative_radius_rational_upper']})
    base_module=load(HERE/'review_conditional_compiler.py')
    direct=0
    for name,S in base_module.matrices().items():
        if len(S)<17:continue
        q=len(S);M=[0,1,2]
        r=candidate.envelope(base_module.load(SOURCE.parent/'conditional_moments.py').matrix_counts(S,M),7)
        h2=[0 if x in M else sum(S[x][M[i]]*S[x][M[j]] for i,j in ((0,1),(0,2),(1,2))) for x in range(q)]
        h3=[0 if x in M else S[x][0]*S[x][1]*S[x][2] for x in range(q)]
        Q=sum(h2[x]*S[x][y]*h3[y] for x in range(q) for y in range(q))
        assert Q==r['Q']
        assert r['vector_statistics']=={'h2_norm_squared':sum(a*a for a in h2),'h3_norm_squared':sum(a*a for a in h3),
                                       'h2_sum':sum(h2),'h3_sum':sum(h3)}
        direct+=1
    out={'status':'passed','date':'2026-09-05','rational_rounding_checks':rounded,
         'saved_envelope_checks':len(cases),'direct_matrix_envelope_checks':direct,'cases':cases,
         'source_sha256':{SOURCE.name:sha256(SOURCE.read_bytes()).hexdigest(),
                          ledger_path.name:sha256(ledger_path.read_bytes()).hexdigest(),
                          Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'scope':'Bound only for the Q-dependent correction to a conditional second moment; not a pointwise target bound.'}
    (HERE/'envelope_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
