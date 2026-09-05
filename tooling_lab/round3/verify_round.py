#!/usr/bin/env python3
"""Verify independent round3 evidence and bind the current artifact snapshot."""
from datetime import datetime,timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def digest(p):return sha256(p.read_bytes()).hexdigest()


def main():
    for p in HERE.rglob('*.py'):compile(p.read_text(),str(p),'exec')
    # Check the completed phase lane's saved bytes before replaying diagnostics.
    phase=HERE/'phase_transport'
    lane_manifest=json.loads((phase/'manifest.json').read_text())
    for row in lane_manifest['files']:
        assert digest(phase/row['path'])==row['sha256'],row['path']
    env=os.environ.copy();env.update({'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    phase_runtime=json.loads((phase/'results.json').read_text())['runtime']['python']
    commands=[(sys.executable,'pair_type_twins/controls.py'),(sys.executable,'novelty/review_twins.py'),
              (phase_runtime,'phase_transport/verify_phase.py'),(sys.executable,'phase_transport/review_phase_root.py')]
    for runtime,path in commands:
        print('Running '+path,flush=True)
        subprocess.run([runtime,str(HERE/path)],cwd=HERE,env=env,check=True)
    for result,source in [('phase_transport/results.json','phase_transport/ambiguity_lab.py'),
                          ('phase_transport/refinement.json','phase_transport/refine_phase.py'),
                          ('phase_transport/validation.json','phase_transport/verify_phase.py')]:
        assert json.loads((HERE/result).read_text())['source_sha256']==digest(HERE/source)
    review=json.loads((HERE/'novelty/twins_review.json').read_text());assert review['passed']
    assert review['results_sha256']==digest(HERE/'pair_type_twins/results.json')
    for name,h in review['candidate_source_sha256'].items():assert h==digest(HERE/'pair_type_twins'/name)
    broken=[]
    for p in HERE.rglob('*.md'):
        prose=re.sub(r'```.*?```|\\\[.*?\\\]', '', p.read_text(), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('https://','http://','#','mailto:')):continue
            if not (p.parent/target.split('#')[0]).exists():broken.append((str(p.relative_to(HERE)),target))
    assert not broken,broken
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name not in ['verification.log'] and p!=HERE/'manifest.json'}
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),
              'status':'independent finite exact/numerical checks passed; no prize, bound or historical novelty certification',
              'phase_snapshot_files_checked':len(lane_manifest['files']),
              'independent_programs_replayed':[p for _,p in commands],
              'markdown_local_links_passed':True,'artifacts':artifacts,
              'note':'This current snapshot supersedes nested lane snapshot byte listings for regenerated diagnostic files.'}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Round3 verification passed; '+str(len(artifacts))+' artifacts pinned.',flush=True)


if __name__=='__main__':main()
