#!/usr/bin/env python
"""Regenerate reference outputs for the RTS regression tests.

  python tests/rts/make_reference.py --baseline-binary /path/to/RTS-branch/dvmdostem \
                                     --rts-binary ./dvmdostem
The baseline reference should be produced by the *unmodified* upstream binary
(RTS branch) with thermokarst off, so the test proves no behavior change.
"""
import argparse, os, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rtslib as r

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference')

def keep(outdir, dest):
    if os.path.exists(dest): shutil.rmtree(dest)
    os.makedirs(dest)
    for f in os.listdir(outdir):
        if f.endswith('.nc') and not f.startswith('restart') and f != 'run_status.nc':
            shutil.copy(os.path.join(outdir, f), dest)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--baseline-binary', default=os.path.join(r.REPO, 'dvmdostem'))
    ap.add_argument('--rts-binary', default=os.path.join(r.REPO, 'dvmdostem'))
    ap.add_argument('--work', default='/tmp/rts_ref_work')
    a = ap.parse_args()
    keep(r.run(a.baseline_binary, a.work + '/baseline', False), os.path.join(REF, 'baseline'))
    keep(r.run(a.rts_binary, a.work + '/rts', True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS}),
         os.path.join(REF, 'rts'))
    print('reference outputs written to', REF)
