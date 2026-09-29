#!/usr/bin/env python3
"""Source binding, exact scalar checks, and preservation of the preceding round."""
from datetime import datetime,timezone
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb,prod,isqrt
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(p,expected):
    assert p.is_file() and digest(p)==expected,str(p)
    BINDINGS.append(str(p.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for x in value:scan(owner,x)
    elif isinstance(value,dict):
        for key,row in value.items():
            if key.endswith('_sha256') and isinstance(row,dict):
                for name,expected in row.items():bind(owner.parent/name,expected)
            else:scan(owner,row)


def main():
    previous=json.loads((HERE.parent/'round8/manifest.json').read_text())
    for name,row in previous['artifacts'].items():bind(HERE.parent/'round8'/name,row['sha256'])
    def read(name):return json.loads((HERE/name).read_text())
    for p in HERE.rglob('*.json'):
        if p.name!='manifest.json':scan(p,json.loads(p.read_text()))
    names=('preflight_results.json','backend_verification.json','interoperability_verification.json',
           'scale_results.json','scale_verification.json','twin_analysis.json','coding_block_results.json')
    for name in names:assert read(name)['status']=='passed',name
    small=read('preflight_results.json');assert len(small['boundary_cases'])==99 and small['literal_distinct_insertion_pair_evaluations']==9618
    backend=read('backend_verification.json');assert len(backend['cases'])==42 and backend['scalar_statistics_compared']==34508
    assert backend['invalid_queries_rejected']==7 and backend['invalid_matrices_rejected']==3
    interop=read('interoperability_verification.json');assert len(interop['cases'])==4 and interop['affine_controls']==21
    for pair in small['pairs']:
        assert not pair['equal_two_swap_histograms'] and pair['equal_two_swap_variances']
    for pair in read('twin_analysis.json')['pairs']:assert not any(pair['equal_features'].values())
    rows=read('scale_results.json')['cases'];checks=read('scale_verification.json')['cases'];assert len(rows)==len(checks)==8
    assert sum(r['full_K_upper_entries_checked'] for r in checks)==314
    assert sum(r['sampled_K_off_diagonal_entries_checked'] for r in checks)==48
    primes=[2013265921,1811939329,469762049]
    for p,root,factors in zip(primes,[31,13,3],[[2,3,5],[2,3],[2,7]]):
        assert all(p%d for d in range(2,isqrt(p)+1))
        assert pow(root,p-1,p)==1 and all(pow(root,(p-1)//f,p)!=1 for f in factors)
    assert prod(primes)>2*10_000_000*(10_000_000-1)*comb(62,5)*comb(62,4)
    for row,check in zip(rows,checks):
        s=row['statistics'];q=s['q'];n=len(s['selected']);m=q-n;N=m*(m-1)
        assert q==check['q'] and n==check['n'] and row['family']==check['family']
        assert s['Q']==check['Q_reconstructed']
        for p,r in zip(check['Q_modular_review']['primes'],check['Q_modular_review']['residues']):assert s['Q']%p==r
        assert check['CRT_modulus']>2*check['Q_bound']
        assert row['ordered_pair_count']==N and row['distinct_final_sets']==N//2
        mean=F(row['target_sum'],N);second=F(row['target_square_sum'],N);var=second-mean**2
        assert mean==F(*row['mean']) and second==F(*row['second_moment']) and var==F(*row['variance'])>=0
        assert var-F(*row['variance_without_Q'])==F(4*s['Q'],N)==F(*row['Q_variance_correction'])
        radius2=F(*row['Q_free_variance_radius_squared'])
        assert F(4*s['Q'],N)**2<=radius2
        if q==6700417:
            assert N//2==22447455617161 and radius2< (F(384,10**13)*var)**2
    coding=read('coding_block_results.json')
    assert coding['complete_monic_polynomials_tested']==5220
    assert sum(n for a,b,n in coding['complete_count_histogram'])==5220
    assert coding['all_basepoints_audit']['intersection_dimension_sum']==2>F(*coding['displayed_lemma_4_1_rhs'])
    assert coding['disjoint_audit']['intersection_dimension_sum']==1
    sources=list(HERE.rglob('*.py'))
    for p in sources:compile(p.read_text(),str(p),'exec')
    for p in HERE.rglob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',p.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('https://','http://','#')):continue
            path=p.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            assert path.exists() or path==HERE/'manifest.json',(p,target)
    assert not list(HERE.glob('*partial.json'))
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
               and p.name not in ('manifest.json','verification.log')}
    out={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Exact completed round with explicitly partial large K readback; no human review, new Lean theorem or historical novelty certification.',
         'previous_round8_artifacts_preserved':len(previous['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round9 verified: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings, all {len(previous['artifacts'])} round8 artifacts preserved.",flush=True)


if __name__=='__main__':main()
