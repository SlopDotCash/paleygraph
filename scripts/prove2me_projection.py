#!/usr/bin/env python3
"""Idempotent private publishing and verification of the projection count lemma."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import urllib.error
import urllib.parse
import urllib.request
from uuid import uuid4
import prove2me as api

ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / 'results/prove2me_projection_server_2026_09_05.json'
PAYLOAD = ROOT / 'prove2me/projection-survival-problem.json'
PROOF = ROOT / 'prove2me/Sol_paley_projection_survival_count.lean'
EXPLANATION = ROOT / 'prove2me/Sol_paley_projection_survival_count.md'

def save(state):
    state['updated_at_utc'] = datetime.now(timezone.utc).isoformat()
    STATE.write_text(json.dumps(state, indent=2) + '\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['publish', 'verify', 'poll'])
    args = parser.parse_args()
    token, version = api.authenticate()
    payload = json.loads(PAYLOAD.read_text())
    assert payload['private'] is True and payload['env'] == api.ENV
    state = json.loads(STATE.read_text()) if STATE.exists() else {
        'platform_version': version, 'private': True, 'environment': api.ENV,
        'theorem_name': payload['problems'][0]['theorem_name'],
        'proof_sha256': sha256(PROOF.read_bytes()).hexdigest(),
        'payload_sha256': sha256(PAYLOAD.read_bytes()).hexdigest(),
        'scope': 'Finite incidence counting lemma; full code theorem and conjectures remain unproved.',
    }
    assert state['proof_sha256'] == sha256(PROOF.read_bytes()).hexdigest()
    assert state['payload_sha256'] == sha256(PAYLOAD.read_bytes()).hexdigest()
    if args.action == 'publish' and not state.get('job_id') and not state.get('theorem_id'):
        query = urllib.parse.urlencode({'theorem_name': state['theorem_name'], 'env': api.ENV})
        matches = api.request('/theorems?' + query, token)
        if matches.get('total', 0):
            raise RuntimeError('Existing name: inspect its statement before reusing.')
        local = json.loads((ROOT / 'results/prove2me_projection_local_2026_09_05.json').read_text())
        assert local['local_proof_exit_code'] == 0 and not local['sorryAx']
        assert local['proof_sha256'] == state['proof_sha256']
        response = api.request('/submit-problem', token, payload)
        state['publish_response'] = response
        if response.get('jobs'):
            assert len(response['jobs']) == 1 and not response.get('errors')
            state['job_id'] = response['jobs'][0]['job_id']
        save(state)
    if state.get('job_id'):
        job = api.request('/publish-jobs/' + state['job_id'], token)
        state['publish_job'] = job
        if job.get('status') == 'PUBLISHED':
            assert job.get('visibility') == 'private'
            assert job.get('formal_statement') == payload['problems'][0]['formal_statement']
            state['theorem_id'] = job['theorem_id']
        save(state)
    if args.action == 'verify' and not state.get('submission_id'):
        if not state.get('theorem_id'):
            raise RuntimeError('Wait for the private statement to finish compiling.')
        boundary = 'PaleyProof' + uuid4().hex
        parts = []
        for name, value in [('theorem_id', state['theorem_id']), ('proof_type', 'prove'),
                            ('explanation', EXPLANATION.read_text())]:
            parts.append(('--' + boundary + '\r\nContent-Disposition: form-data; name="' + name
                          + '"\r\n\r\n' + value + '\r\n').encode())
        parts.append(('--' + boundary + '\r\nContent-Disposition: form-data; name="file"; '
                      'filename="solution.lean"\r\nContent-Type: text/plain\r\n\r\n').encode()
                     + PROOF.read_bytes() + b'\r\n')
        parts.append(('--' + boundary + '--\r\n').encode())
        req = urllib.request.Request(api.BASE + '/verify', data=b''.join(parts),
            headers={'Authorization': 'Bearer ' + token,
                     'Content-Type': 'multipart/form-data; boundary=' + boundary})
        with api.OPENER.open(req, timeout=25) as response:
            result = json.load(response)
        state['verify_response'] = result
        state['submission_id'] = result['submission_id']
        save(state)
    if state.get('submission_id'):
        state['verdict'] = api.request('/verify?submission_id=' + state['submission_id'], token)
        if state['verdict'].get('status') != 'PENDING':
            theorem = api.request('/theorems/' + state['theorem_id'], token)
            state['theorem_readback'] = {k: theorem.get(k) for k in
                ['theorem_id', 'theorem_name', 'status', 'mathlib_rev', 'formal_statement', 'audits']}
        save(state)
    print(json.dumps({'job_id': state.get('job_id'),
          'publish_status': state.get('publish_job', {}).get('status'),
          'theorem_id': state.get('theorem_id'), 'submission_id': state.get('submission_id'),
          'verdict': state.get('verdict', {}).get('status'), 'private': state['private']}, indent=2))

if __name__ == '__main__':
    try:
        main()
    except urllib.error.HTTPError as e:
        print('Prove2Me request failed: HTTP ' + str(e.code), file=sys.stderr)
        sys.exit(1)
    except (OSError, ValueError, KeyError, RuntimeError, AssertionError):
        print('Stopped: inspect the nonsecret saved state and local proof before retrying.', file=sys.stderr)
        sys.exit(1)
