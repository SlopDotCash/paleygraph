#!/usr/bin/env python3
"""Bind final source/input bytes and integrate saved independent round5 checks.

This does not rerun the complete subset censuses or the large convolution.
The separately implemented reviews have already run; their hashes are pinned.
"""
from datetime import datetime,timezone
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
LAB=HERE.parent
BINDINGS=[]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def read(name):return json.loads((HERE/name).read_text())


def bind(path,expected):
    assert path.is_file() and digest(path)==expected,path
    BINDINGS.append(str(path.relative_to(LAB)))


def bind_map(folder,values):
    for name,expected in values.items():bind(HERE/folder/name,expected)


def multiply(a,b):
    out=[Fraction() for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def main():
    sources=list(HERE.rglob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    nested=0
    for folder in ('preflight','trace_backend','trace_invariants'):
        manifest=read(folder+'/manifest.json')
        records=manifest.get('files')
        if records is None:records=[{'path':name,**row} for name,row in manifest['artifacts'].items()]
        for row in records:
            bind(HERE/folder/row['path'],row['sha256']);nested+=1
    previous=json.loads((LAB/'round4/manifest.json').read_text())
    for name,row in previous['artifacts'].items():bind(LAB/'round4'/name,row['sha256'])

    general=read('third_moment/results.json')
    bind_map('third_moment',general['source_sha256'])
    bind(LAB/'round4/marked_moments/conditional_moments.py',general['previous_second_moment_compiler_sha256'])
    for row in general['scale_cases']:
        bind(LAB/row['source_inventory'],row['source_inventory_sha256'])
        bind_map('third_moment',row['source_sha256'])
        mean,second,variance,third=[Fraction(*row[k]) for k in ('mean','second_moment','variance','third_moment')]
        assert variance==second-mean**2
        central=third-3*mean*second+2*mean**3
        assert central==Fraction(*row['central_third_moment'])
        assert central**2/variance**3==Fraction(*row['standardized_third_squared'])
    for row in general['twins']:
        path=HERE/'review'/f'twins_n{row["n"]}_third_oracle.json'
        bind(path,row['oracle_sha256'])
        oracle=next(x for x in json.loads(path.read_text())['records'] if x['graph']==row['graph'])
        assert [row['mean'],row['second_moment'],row['third_moment']]==oracle['moments']
        assert (row['minimum'],row['maximum'],row['subset_count'])==(oracle['minimum'],oracle['maximum'],oracle['subset_count'])
    assert len(general['twins'])==6
    assert all(row['exact_equality'] for row in general['direct_interpolation_comparisons'])

    for name in ('direct_third_oracle.json','twins_n6_third_oracle.json','twins_n7_third_oracle.json',
                 'twins_n8_third_oracle.json','review_union_backend.json','review_global_third.json',
                 'review_piece_discovery.json','review_piece_scale.json','fold_review_derivation.json'):
        bind_map('review',read('review/'+name)['source_sha256'])
    global_review=read('review/review_global_third.json')
    bind_map('review',global_review['oracle_artifact_sha256'])
    assert global_review['compiler_comparisons']==80 and global_review['finite_cases']==40
    assert global_review['ordered_row_triples_normalized_independently']==227388
    assert global_review['odd_degree_cases_rejected']==25
    assert read('review/review_union_backend.json')['coefficient_equalities']==681

    folded=read('trace_invariants/results.json');residue=read('trace_invariants/residue_results.json')
    for record in (folded,residue):bind_map('trace_invariants',record['source_sha256'])
    for name,value in folded['input_sha256'].items():bind(HERE/name,value)
    bind(HERE/'trace_invariants/results.json',residue['baseline_sha256'])
    reference={(r.get('graph',f'Paley{r["q"]}'),r['n']):r for r in general['scale_cases']+general['twins']}
    for record,budget in ((folded,19),(residue,15)):
        assert len(record['comparisons'])==14
        for row in record['comparisons']:
            key=(row.get('graph',f'Paley{row["q"]}'),row['n'])
            assert row['third_moment']==reference[key]['third_moment']
            assert row['coefficient_evaluations']<=budget
            for name,value in row['source_sha256'].items():
                if '/' in name:folder=''
                elif name=='third_moment.py':folder='third_moment'
                else:folder='trace_invariants'
                bind(HERE/folder/name,value)
            assert sum(Fraction(*row[k]) for k in ('all_equal_rows_contribution','exactly_two_equal_rows_contribution','distinct_rows_contribution'))==Fraction(*row['third_moment'])
    universal=[Fraction(1)]
    for _ in range(6):universal=multiply(universal,[0,1,-3,2])
    universal=[c/720 for c in universal]
    assert [Fraction(*c) for c in residue['universal_theta6_union_coefficients']]==universal
    for row in residue['raw_union_high_degree_checks']:
        assert row['theta5_zero'] and row['theta6_equals_universal_polynomial'] and row['seven_statistic_coefficient_rank']==7
    review=read('review/fold_review.json')
    for name,value in review['source_sha256'].items():
        bind(HERE/('review' if name=='fold_review.py' else 'trace_invariants')/name,value)

    piece=read('piece_discovery/results.json');verification=read('piece_discovery/verification.json')
    bind_map('piece_discovery',piece['source_sha256']);bind_map('piece_discovery',verification['source_sha256'])
    bind(HERE/'piece_discovery/results.json',verification['results_sha256'])
    bind(HERE/'piece_discovery/results.json',read('review/review_piece_scale.json')['results_sha256'])
    counts=read('review/review_piece_discovery.json')['counts']
    assert counts['cases']==928 and counts['exact_minimum_cover_sizes_checked']==805
    assert counts['complete_static_list_oracles']==793 and counts['complete_scalar_fiber_oracles']==164
    assert verification['status']=='passed'
    assert len(piece['hard_cosets'])==16
    assert all(not row['summary']['complete_static_list_certified'] for row in piece['hard_cosets'])
    for row in piece['large']:
        result=row['result'];r=result['minimum_polynomial_cover_size_certified']
        assert r in (2,3) and result['complete_static_list_certified'] and result['complete_scalar_fibers_certified']
        assert result['static_certificate']['outside_piece_agreement_cap']==r*(result['k']-1)<result['s']
        assert row['chunk_baseline']['root_count_cap']>=result['s']

    plot=read('plot_metadata.json');bind(HERE/'plot_results.py',plot['source_sha256'])
    for name,value in plot['input_sha256'].items():bind(HERE/name,value)
    broken=[]
    for path in HERE.rglob('*.md'):
        prose=re.sub(r'```.*?```|\\\[.*?\\\]','',path.read_text(),flags=re.S)
        prose=re.sub(r'`[^`\n]*`','',prose)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#','mailto:')):continue
            linked=path.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            if linked==HERE/'manifest.json':continue
            if not linked.exists():broken.append((str(path.relative_to(HERE)),target))
    assert not broken,broken
    artifacts={str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size}
               for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
               and path.name!='verification.log' and path!=HERE/'manifest.json'}
    output={'created_utc':datetime.now(timezone.utc).isoformat(),
            'status':'Final source/input bindings and exact saved independent evidence passed; no prize or historical originality claim.',
            'verification_kind':'Integration audit; independent large censuses and review scripts completed separately.',
            'source_and_input_bindings':len(BINDINGS),'nested_snapshot_files':nested,
            'previous_round4_artifacts_preserved':len(previous['artifacts']),
            'python_sources_compiled':len(sources),'general_compiler_oracle_comparisons':80,
            'refined_compiler_existing_case_comparisons':28,'independent_piece_cases':928,
            'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(output,indent=2)+'\n')
    print(f'Round5 verification passed: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings, 28 refined exact comparisons; round4 preserved.',flush=True)


if __name__=='__main__':main()
