#!/usr/bin/env python3
"""Consolidate exact ablation evidence, including direct scalar contractions."""
import datetime,hashlib,json,sys
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.dont_write_bytecode=True;sys.path.insert(0,str(HERE.parent/'marked_moments'))
from contraction_reduction import reduce_from_cells
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def values(u):
    if 0 in u:return 0,0
    a,b,c=u;return a*b+a*c+b*c,a*b*c
def main():
    validation=json.loads((HERE/'direct_validation.json').read_text());rows=[]
    for case in validation['cases']:
        path=HERE/case['witness'];assert sha(path)==case['witness_sha256'];w=json.loads(path.read_text())
        mp=HERE/(case['witness'].replace('witness.json','witness_matrices.json'));assert sha(mp)==case['matrices_sha256'];matrices=json.loads(mp.read_text())['graphs']
        pair=[]
        for item,mat in zip(w['pair'],matrices):
            S=mat['S'];M=item['marks'];q=len(S);h=[values([S[x][m] for m in M]) for x in range(q)]
            Q=sum(h[x][0]*S[x][y]*h[y][1] for x in range(q) for y in range(q))
            reduced={}
            for n in [6,7,8]:
                r=reduce_from_cells(q,M,{tuple(c['pattern']):c['size'] for c in item['record']['cells']},Q,n)
                assert all(r[k]==item['moments'][str(n)][k] for k in ['mean','second_moment','variance']);reduced[n]=r
            pair.append({'graph':item['graph'],'marks':M,'Q_direct_matrix':Q,'direct_Q_pair_terms':q*q,'reduced_moments':reduced})
        assert w['pair'][0]['record']['cells']==w['pair'][1]['record']['cells']
        gaps=[]
        for n in [6,7,8]:
            a,b=(t['reduced_moments'][n] for t in pair)
            assert a['cell_only_term']==b['cell_only_term'] and a['Q_coefficient']==b['Q_coefficient']
            gap=Fraction(*b['second_moment'])-Fraction(*a['second_moment'])
            predicted=Fraction(*a['Q_coefficient'])*(pair[1]['Q_direct_matrix']-pair[0]['Q_direct_matrix'])
            assert gap==predicted
            gaps.append({'n':n,'second_moment_right_minus_left':[gap.numerator,gap.denominator],
                         'Q_coefficient':a['Q_coefficient'],'exact_gap_explained_by_Q':True})
        rows.append({'witness_file':path.name,'q':w['q'],'same_literal_cell_profile':True,'pair':pair,'gaps':gaps,
                    'independently_enumerated_completions':sum(t['completions'] for t in case['all_direct_completions'])})
    out={'status':'cell sizes alone are insufficient; concrete exact ablations and a single-convolution compression',
         'ablation_pairs':rows,'total_direct_completions':sum(r['independently_enumerated_completions'] for r in rows),
         'contraction_validation_file':'contraction_validation.json','contraction_validation_sha256':sha(HERE/'contraction_validation.json'),
         'compiler_source_sha256':sha(HERE.parent/'marked_moments/conditional_moments.py'),
         'reduction_source_sha256':sha(HERE.parent/'marked_moments/contraction_reduction.py'),'source_sha256':sha(__file__)}
    (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    files=[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json']
    (HERE/'manifest.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'total_bytes':sum(f['bytes'] for f in files)},indent=2)+'\n')
    print(json.dumps({'pairs':len(rows),'direct_completions':out['total_direct_completions'],'files':len(files),'bytes':sum(f['bytes'] for f in files)}))
if __name__=='__main__':main()
