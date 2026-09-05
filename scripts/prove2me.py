#!/usr/bin/env python3
"""Project connection to the shared Prove2Me workspace; never log secrets."""
import argparse
from datetime import datetime, timezone
import getpass
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = Path.home() / 'prove2me_workspace'
BASE = 'https://prove2.me/api/v1'
ENV = 'c5ea00351c28e24afc9f0f84379aa41082b1188f'
PROPOSAL = 'b5e7121a-eb3c-48f1-a495-29fffd24e71d'

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise RuntimeError('Redirect refused.')

OPENER = urllib.request.build_opener(NoRedirect)

def request(path, token=None, data=None, method=None):
    if not path.startswith('/') or '://' in path:
        raise ValueError('Expected a relative API path.')
    headers = {'Accept': 'application/json'}
    if token:
        headers['Authorization'] = 'Bearer ' + token
    if data is not None:
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(BASE + path, headers=headers,
        data=None if data is None else json.dumps(data).encode(), method=method)
    with OPENER.open(req, timeout=25) as response:
        return json.load(response)

def authenticate(key=None):
    if key is None:
        credentials = WORKSPACE / 'credentials.json'
        if credentials.stat().st_mode & 0o077:
            raise RuntimeError('Credential file permissions must be 600.')
        key = json.loads(credentials.read_text())['api_key']
    if not isinstance(key, str) or not key.startswith('p2m_'):
        raise ValueError('Invalid credential format.')
    session = request('/agent/refresh', data={'api_key': key})
    expected = re.search(r'^  version: "([^"]+)"',
                         (WORKSPACE / 'SKILL.md').read_text(), re.MULTILINE)
    if expected is None or session.get('version') != expected.group(1):
        raise RuntimeError('Platform version changed; refresh the official skill first.')
    return session['access_token'], session['version']

def save_credential(key):
    # The shared, gitignored authentication store is outside this project.
    destination = WORKSPACE / 'credentials.json'
    fd, temporary = tempfile.mkstemp(prefix='.credentials-', dir=WORKSPACE)
    try:
        os.fchmod(fd, 0o600)
        with os.fdopen(fd, 'w') as f:
            json.dump({'api_key': key}, f)
            f.write('\n')
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

def check_connection(token, version):
    profile = request('/me', token)
    environments = request('/environments', token)
    if ENV not in json.dumps(environments):
        raise RuntimeError('Pinned Lean environment unavailable.')
    proposal = request('/mission-proposals/' + PROPOSAL, token)
    milestones = request('/mission-proposals/' + PROPOSAL + '/milestones', token)
    previous = json.loads((ROOT / 'results/sigma_p2m_2026_09_05.json').read_text())
    verdicts = {}
    for name, old in previous['verdicts'].items():
        sid = old['submission_id']
        d = request('/verify?submission_id=' + sid, token)
        verdicts[name] = {'submission_id': sid, 'theorem_id': d.get('theorem_id'),
                          'status': d.get('status')}
    goal_id = '973dac7a-fac5-4e8e-961a-936f59818faa'
    goal = request('/theorems/' + goal_id, token)
    report = {
        'checked_at_utc': datetime.now(timezone.utc).isoformat(),
        'authenticated': True, 'platform_version': version,
        'account_user_id': profile.get('user_id', profile.get('id')),
        'environment': ENV, 'supported_environments': environments,
        'lean_workspace': str(WORKSPACE),
        'proposal': {k: proposal.get(k) for k in
            ('id', 'name', 'status', 'visibility', 'env', 'main_item_id')},
        'proposal_items': proposal.get('items', []),
        'milestones': milestones,
        'verified_existing_submissions': verdicts,
        'goal': {'theorem_id': goal_id, 'theorem_name': goal.get('theorem_name'),
                 'status': goal.get('status')},
        'scope': 'Authenticated readback. No new proof or public release.',
    }
    out = ROOT / 'results/prove2me_connection_2026_09_05.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'authenticated': True, 'platform_version': version,
          'pinned_environment_available': True,
          'proposal_status': proposal.get('status'),
          'proposal_visibility': proposal.get('visibility'),
          'existing_proof_verdicts': {k: v['status'] for k, v in verdicts.items()},
          'full_goal_status': goal.get('status'), 'report': str(out)}, indent=2))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--login', action='store_true', help='Read a replacement API key with echo disabled.')
    args = parser.parse_args()
    if args.login:
        if not sys.stdin.isatty():
            raise RuntimeError('Login requires a terminal with echo disabled.')
        key = getpass.getpass('Prove2Me API key (hidden): ').strip().replace('\\_', '_')
        token, version = authenticate(key)
        request('/me', token)
        save_credential(key)
        del key
    else:
        token, version = authenticate()
    check_connection(token, version)

if __name__ == '__main__':
    try:
        main()
    except urllib.error.HTTPError as e:
        print('Prove2Me request failed: HTTP ' + str(e.code), file=sys.stderr)
        sys.exit(1)
    except (OSError, ValueError, KeyError, RuntimeError):
        print('Connection failed; check authentication, network, supported environment, or API version.', file=sys.stderr)
        sys.exit(1)
