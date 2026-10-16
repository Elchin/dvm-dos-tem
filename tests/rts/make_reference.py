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
    ap.add_argument('--legacy-m2d', action='store_true',
                    help='write flag-off (legacy monthly->daily) refs to reference/legacy_m2d')
    ap.add_argument('--only', nargs='*', default=None, help='subset of: baseline rts rts_ice')
    ap.add_argument('--v02', action='store_true',
                    help='write v0.2-default (floor on) rts/rts_ice refs to reference/v02; without it the '
                         'rts/rts_ice refs are the v0.1 ones (written with r.V01_OFF)')
    a = ap.parse_args()
    want = lambda n: a.only is None or n in a.only
    if a.legacy_m2d:
        REF = os.path.join(REF, 'legacy_m2d')
        r.ICE_PARAMS = {**r.LEGACY_ICE_PARAMS, **r.LEGACY_M2D}
        _run = r.run
        r.run = lambda b, w, on, **kw: _run(b, w, on, **{**kw, 'extra': {**(kw.get('extra') or {}), **r.LEGACY_M2D}})
    v01 = {} if a.v02 else r.V01_OFF
    if a.v02:
        REF = os.path.join(REF, 'v02')
    if want('baseline') and not a.v02:
        keep(r.run(a.baseline_binary, a.work + '/baseline', False), os.path.join(REF, 'baseline'))
    if want('rts'):
        keep(r.run(a.rts_binary, a.work + '/rts', True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS}, extra=v01),
             os.path.join(REF, 'rts'))
    if want('rts_ice'):
      keep(r.run(a.rts_binary, a.work + '/rts_ice', True,
               outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS}, extra={**r.ICE_PARAMS, **v01}),
         os.path.join(REF, 'rts_ice'))
    print('reference outputs written to', REF)
