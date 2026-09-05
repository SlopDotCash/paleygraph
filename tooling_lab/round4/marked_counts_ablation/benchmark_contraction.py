#!/usr/bin/env python3
"""Interleave existing full-count and single-contraction backends; exact outputs checked."""
import hashlib,json,statistics,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    rows=[]
    for p in [65537,1000033]:
        timings={'full_relations':[],'single_Q':[]}
        for trial in range(3):
            for name,binary in [('full_relations',HERE.parent/'marked_counts/marked_counts'),('single_Q',HERE/'single_contraction')]:
                start=time.monotonic();out=json.loads(subprocess.check_output([str(binary),str(p),'0','1','2'],text=True));wall=time.monotonic()-start
                if name=='single_Q':
                    reference=json.loads((HERE/f'contraction_p{p}_marks_0_1_2.json').read_text());assert out['Q']==reference['Q'] and out['cells']==reference['cells']
                timings[name].append({'trial':trial,'wall_seconds':wall,'backend_seconds':out['metadata'].get('backend_elapsed_seconds',out['metadata'].get('backend_seconds'))})
        a=statistics.median(t['wall_seconds'] for t in timings['full_relations']);b=statistics.median(t['wall_seconds'] for t in timings['single_Q'])
        rows.append({'p':p,'measurements':timings,'median_full_wall_seconds':a,'median_single_wall_seconds':b,'median_wall_speedup':a/b})
        print(json.dumps(rows[-1]),flush=True)
    (HERE/'contraction_benchmark.json').write_text(json.dumps({'status':'interleaved finite runtime comparison; shared machine load not isolated',
       'cases':rows,'source_sha256':sha(__file__),'full_backend_sha256':sha(HERE.parent/'marked_counts/marked_counts'),
       'single_backend_sha256':sha(HERE/'single_contraction')},indent=2)+'\n')
if __name__=='__main__':main()
