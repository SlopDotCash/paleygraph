#!/usr/bin/env python3
"""Summarize verified exact count artifacts and freeze a local manifest."""
import datetime,hashlib,json,platform,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    validation=json.loads((HERE/'validation.json').read_text());rows=[]
    assert validation['source_sha256']==sha(HERE/'verify_counts.py')
    for v in validation['cases']:
        path=HERE/v['file'];assert v['sha256']==sha(path);d=json.loads(path.read_text())
        rows.append({'file':path.name,'sha256':sha(path),'p':d['p'],'marks':d['marks'],'cells':len(d['cells']),
             'nonzero_relations':len(d['relations']),**d['metadata'],
             'wall_seconds':d['provenance']['wall_seconds'],
             'runner_lifetime_max_child_rss_bytes':d['provenance']['children_peak_rss_bytes'],
             'independently_recomputed_all_relations':v['full_relation_counts_independently_recomputed'],
             'independent_sampled_rows':v['independent_sampled_rows']})
    out={'status':'exact marked-cell input layer implemented and validated; conditional moments computed elsewhere',
         'cases':rows,'total_independent_direct_pairs':sum(v['independent_direct_ordered_pairs'] for v in validation['cases']),
         'total_independent_sampled_row_terms':sum(v['independent_sampled_row_terms'] for v in validation['cases']),
         'affine_controls':validation['affine_covariance_controls'],
         'rejection_guard_tests':len(validation['rejection_guards']),
         'runtime':{'python':sys.executable,'version':platform.python_version(),'machine':platform.machine(),
             'compiler':subprocess.check_output(['/usr/bin/clang++','--version'],text=True).splitlines()[0]},
         'source_sha256':sha(__file__)}
    (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    files=[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json']
    (HERE/'manifest.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'files':files,'total_bytes':sum(p['bytes'] for p in files),'scope':'round4/marked_counts only'},indent=2)+'\n')
    print(json.dumps({'files':len(files),'cases':len(rows),'total_bytes':sum(p['bytes'] for p in files),'source_and_artifact_hashes':'passed'}))
if __name__=='__main__':main()
