#!/usr/bin/env python3
"""Consolidate exact backend results and freeze the local manifest."""
import datetime,hashlib,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    validation=json.loads((HERE/'validation.json').read_text());assert validation['source_sha256']==sha(HERE/'verify_inventory.py');rows=[]
    for row in validation['cases']:
        path=HERE/row['inventory_file'];assert sha(path)==row['inventory_sha256'];data=json.loads(path.read_text());p=data['p']
        assert data['total_power_sums'][:3]==[p-2,2,p*p-2*p-3]
        for cls in data['edge_power_sums']:
            assert cls['powers'][0]==(p-5)//4 if cls['edges']==[1,1,1] else cls['powers'][0]==(p-1)//4
        rows.append({'p':p,'file':path.name,'sha256':sha(path),'joint_bins':len(data['records']),
            'edge_classes':len(data['edge_power_sums']),'n6_third_moment':data['n6_global_third_moment']['third_moment'],
            'n6_third_moment_display':data['n6_global_third_moment']['third_moment_display'],
            'sixth_power_sum_bits':abs(data['total_power_sums'][6]).bit_length(),
            'global_third_moment_numerator_bits':abs(data['n6_global_third_moment']['sum_T6_cubed']).bit_length(),
            'ntt_length':data['metadata']['ntt_length'],'transforms':3,'backend_seconds':data['metadata']['backend_seconds'],
            'wall_seconds':data['provenance']['wall_seconds'],'maximum_child_rss_bytes':data['provenance']['maximum_child_rss_bytes'],
            'full_array_independent_direct_verification':row['full_tau_array_independently_recomputed'],
            'direct_rows':row['direct_tau_rows'],'complete_involution_parameters':row['reflection_and_signed_inversion_parameters_checked']})
    (HERE/'summary.json').write_text(json.dumps({'status':'exact prime-Paley trace adapter and n6 third moment validated',
        'cases':rows,'total_direct_tau_integer_terms':sum(r['direct_integer_terms'] for r in validation['cases']),
        'low_moment_and_edge_class_count_controls':'all passed','compiler':subprocess.check_output(['/usr/bin/clang++','--version'],text=True).splitlines()[0],
        'source_sha256':sha(__file__)},indent=2)+'\n')
    files=[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json']
    (HERE/'manifest.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'total_bytes':sum(f['bytes'] for f in files)},indent=2)+'\n')
    print(json.dumps({'files':len(files),'bytes':sum(f['bytes'] for f in files),'cases':len(rows),'exact_controls':'passed'}))
if __name__=='__main__':main()
