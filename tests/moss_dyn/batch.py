#!/usr/bin/env python3
"""batch.py CASES.json [P=2]: CASES.json = [[name, site, binary, {model_settings}, [extra args]], ...]"""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
H = os.path.dirname(os.path.abspath(__file__))
cases = json.load(open(sys.argv[1])); P = int(sys.argv[2]) if len(sys.argv) > 2 else 2
def run(c):
    name, site, b, ms = c[:4]; extra = c[4] if len(c) > 4 else []
    b = os.path.expanduser(b)
    p = subprocess.run([sys.executable, os.path.join(H, 'run_case.py'), name, site, b, json.dumps(ms)] + extra,
                       capture_output=True, text=True)
    print(name, p.returncode, (p.stdout + p.stderr).strip()[-300:], flush=True)
with ThreadPoolExecutor(P) as ex: list(ex.map(run, cases))
print('BATCH_DONE', flush=True)
