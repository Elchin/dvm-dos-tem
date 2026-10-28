#!/usr/bin/env python3
"""Phase 2 driver: per-cell cohort-response runs for the two-case RTS design.

Per cell (Gydan 10x10, CMT70, NEX-GDDP-CMIP6 MRI-ESM2-0 SSP5-8.5 climate,
Toolik RCP8.5 CO2), from a stage-1 restart-sp.nc written by run_stage1.py:
  stage2a/   historical transient (RTS off), -t HIST, shared by all runs
  base_fNN/  stage 2b baseline (thermokarst on, initiation off): the A=0 column
  yYYYY_fNN/ stage 2b with ONE RTS initiation in year YYYY (Aug, jday 213; no
             repeat slumps), trigger file written by rts_trigger_file.py
             from the cell's projected-explicit-thermokarst-rts.nc
  notrig_fNN/ (--checks) rts mode with an all-zero mask: must equal base_fNN bit-for-bit
fNN = thermokarst_folr x10 (0.8 main, 0.5 sensitivity). The cohort runs are
blended offline (phase 3): case A Jia et al. 2026 (Sci. Data
doi:10.1038/s41597-026-08284-y) cell-mean susceptibility, gt050 sensitivity;
case B Nesterova et al. 2025 counts x 1.44 ha (ARTS) with the SSP5 gt050 ratio
applied to the initiation RATE. N switches (n_nonvasc_extract, n_extract_limiter,
soil_giveup_diag) are left at their defaults (off). Outputs are reduced to what
the blend and closure checks need (no layered LAYERDZ/LWCLAYER/TLAYER).
Resumable: runs with timings.json are skipped. NBAL_DIAG=1 -> nbal.txt per run.
"""
import argparse, json, os, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'tests', 'rts')); sys.path.insert(0, os.path.join(REPO, 'scripts'))
YEARS = tuple(range(2021, 2100, 10))
DROP = ('LAYERDZ', 'LWCLAYER', 'TLAYER')
TRIG = 'projected-explicit-thermokarst-rts.nc'


def outputs():
    import rtslib as r
    o = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS, **r.CH4_OUTPUTS, **r.REGROW_OUTPUTS,
         'AVLN': 'm', 'ORGN': 'm', 'VEGNTOT': 'm', 'ALD': 'y'}
    return {k: v for k, v in o.items() if k not in DROP}


def run_name(kind, year, folr):
    f = 'f%02d' % round(folr * 10)
    return {'base': f'base_{f}', 'notrig': f'notrig_{f}'}.get(kind, f'y{year}_{f}')


def plan(cells, years, folrs, checks):
    runs = []
    for c in cells:
        for f in folrs:
            runs.append((c, 'base', None, f))
            runs += [(c, 'rts', y, f) for y in years]
            if checks: runs.append((c, 'notrig', None, f))
    return runs


def _wrap(w, binary):
    nb = os.path.join(w, 'nbal.txt'); open(nb, 'w').close()
    wr = os.path.join(w, 'wrap.sh')
    open(wr, 'w').write(f'#!/bin/bash\nNBAL_DIAG=1 {binary} "$@" 2>>{nb}\n'); os.chmod(wr, 0o755)
    return wr


def _overrides(data, w, spec, folr=None, trig=None):
    ov = json.load(open(os.path.join(data, 'two_stage_overrides.json')))
    for k, v in ov['IO'].items():
        if isinstance(v, str) and v.endswith('.nc'): ov['IO'][k] = os.path.join(data, v)
    ov['IO']['output_spec_file'] = spec
    ov['IO']['parameter_dir'] = os.path.join(REPO, 'parameters') + '/'
    for k in ('n_nonvasc_extract', 'n_extract_limiter', 'soil_giveup_diag'):
        assert k not in ov.get('model_settings', {}), k   # keep the defaults (off)
    ov['model_settings']['thermokarst_regrow'] = True     # v0.3 default, explicit as in spatial10x10_regrow
    if folr is not None: ov['model_settings']['thermokarst_folr'] = folr
    if trig: ov['IO']['proj_exp_thermokarst_file'] = trig
    p = os.path.join(w, 'overrides.json'); json.dump(ov, open(p, 'w'), indent=1)
    return p


def _call(w, args):
    c0 = os.times(); t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'run_rts_two_stage.py')] + args,
                       capture_output=True, text=True)
    c1 = os.times()
    open(os.path.join(w, 'run.log'), 'w').write(p.stdout + p.stderr)
    if p.returncode: return None
    tp = os.path.join(w, 'timings.json'); t = json.load(open(tp))
    t['wall_total_s'] = time.time() - t0
    t['cpu_children_s'] = (c1.children_user - c0.children_user) + (c1.children_system - c0.children_system)
    json.dump(t, open(tp, 'w'), indent=2)
    return t


def stage2a(args):
    cell, inputs, stage1, out, binary, hist = args
    w = os.path.join(out, cell, 'stage2a'); tp = os.path.join(w, 'timings.json')
    if os.path.exists(tp) and os.path.exists(os.path.join(w, 'stage2a_hist', 'restart-tr.nc')):
        return cell, 'skip', json.load(open(tp))
    os.makedirs(w, exist_ok=True)
    spec = os.path.join(out, cell, 'output_spec.csv')
    import rtslib as r; r.write_outspec(spec, outputs())
    ovp = _overrides(os.path.join(inputs, cell), w, spec)
    t = _call(w, ['--binary', _wrap(w, binary), '--workdir', w, '--config-overrides', ovp,
                  '--stage1-restart', os.path.join(stage1, cell, 'stage1', 'restart-sp.nc'),
                  '--hist', str(hist), '--stop-after', 'stage2a'])
    return cell, ('ok' if t else 'FAIL'), t or {}


def stage2b(args):
    (cell, kind, year, folr), inputs, stage1, out, binary, sc, ch4 = args
    from rts_sus2case import rts_trigger_file as rtf
    name = run_name(kind, year, folr); w = os.path.join(out, cell, name); tp = os.path.join(w, 'timings.json')
    if os.path.exists(tp):
        return cell, name, 'skip', json.load(open(tp))
    os.makedirs(w, exist_ok=True); data = os.path.join(inputs, cell)
    trig = None
    if kind in ('rts', 'notrig'):
        trig = os.path.join(w, TRIG); rtf.write_trigger(os.path.join(data, TRIG), trig, year)
    ovp = _overrides(data, w, os.path.join(out, cell, 'output_spec.csv'), folr, trig)
    a = ['--binary', _wrap(w, binary), '--workdir', w, '--config-overrides', ovp,
         '--stage1-restart', os.path.join(stage1, cell, 'stage1', 'restart-sp.nc'),
         '--stage2a-restart', os.path.join(out, cell, 'stage2a', 'stage2a_hist', 'restart-tr.nc'),
         '--sc', str(sc), '--mode', 'baseline' if kind == 'base' else 'rts']
    if ch4: a.append('--ch4')
    t = _call(w, a)
    if t is not None:
        t.update(kind=kind, year=year, folr=folr)
        json.dump(t, open(tp, 'w'), indent=2)
    return cell, name, ('ok' if t else 'FAIL'), t or {}


def summarize(vals):
    import statistics as st
    v = sorted(x for x in vals if x is not None)
    return dict(n=len(v), min=v[0], median=st.median(v), max=v[-1], sum=sum(v)) if v else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', required=True); ap.add_argument('--stage1', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--binary', required=True)
    ap.add_argument('--nproc', type=int, default=16); ap.add_argument('--cells', default=None)
    ap.add_argument('--years', default=','.join(map(str, YEARS)))
    ap.add_argument('--folr', default='0.8,0.5')
    ap.add_argument('--hist', type=int, default=66, help='1950-2015'); ap.add_argument('--sc', type=int, default=85, help='2016-2100')
    ap.add_argument('--no-ch4', action='store_true', help='drop the RTS CH4 diagnostic (on in spatial10x10_regrow)')
    ap.add_argument('--checks', action='store_true', help='also run notrig_fNN (all-zero mask) per cell/folr')
    a = ap.parse_args()
    inputs, stage1, out, binary = (os.path.abspath(os.path.expanduser(x)) for x in (a.inputs, a.stage1, a.out, a.binary))
    cells = a.cells.split(',') if a.cells else [c['cell'] for c in json.load(open(os.path.join(inputs, 'grid.json')))['cells']]
    years = [int(y) for y in a.years.split(',') if y]; folrs = [float(f) for f in a.folr.split(',')]
    os.makedirs(out, exist_ok=True); log_p = os.path.join(out, 'status.log')
    def log(m):
        line = time.strftime('%Y-%m-%d %H:%M:%S UTC ', time.gmtime()) + m
        open(log_p, 'a').write(line + '\n'); print(line, flush=True)
    runs = plan(cells, years, folrs, a.checks)
    log(f'start: {len(cells)} cells, stage2a x{len(cells)}, stage2b x{len(runs)}, nproc {a.nproc}, '
        f'stage1 {stage1}, binary {binary}, ch4 {not a.no_ch4}')
    t0 = time.time(); T = {'stage2a': {}, 'stage2b': {}}
    with ProcessPoolExecutor(a.nproc) as ex:
        for c, st_, t in ex.map(stage2a, [(c, inputs, stage1, out, binary, a.hist) for c in cells]):
            T['stage2a'][c] = t.get('stage2a'); log(f'{c} stage2a {st_}')
        t1 = time.time()
        ok2a = {c for c in cells if T['stage2a'].get(c) is not None or os.path.exists(os.path.join(out, c, 'stage2a', 'stage2a_hist', 'restart-tr.nc'))}
        todo = [r for r in runs if r[0] in ok2a]
        nfail = 0
        for k, (c, name, st_, t) in enumerate(ex.map(stage2b, [(r, inputs, stage1, out, binary, a.sc, not a.no_ch4) for r in todo], chunksize=4)):
            T['stage2b'][f'{c}/{name}'] = dict(s=t.get('stage2b'), kind=t.get('kind'))
            nfail += st_ == 'FAIL'
            if st_ == 'FAIL' or k % 50 == 0: log(f'{c} {name} {st_} ({k + 1}/{len(todo)})')
    s2b = T['stage2b']
    summ = dict(cells=len(cells), runs=len(runs), failed=nfail, nproc=a.nproc,
                wall_stage2a_s=t1 - t0, wall_stage2b_s=time.time() - t1, wall_total_s=time.time() - t0,
                stage2a_s=summarize(T['stage2a'].values()),
                stage2b_baseline_s=summarize([v['s'] for v in s2b.values() if v['kind'] == 'base']),
                stage2b_rts_s=summarize([v['s'] for v in s2b.values() if v['kind'] == 'rts']),
                stage2b_notrig_s=summarize([v['s'] for v in s2b.values() if v['kind'] == 'notrig']))
    json.dump(dict(summary=summ, **T), open(os.path.join(out, 'timings.json'), 'w'), indent=1)
    log('end: ' + json.dumps({k: v for k, v in summ.items() if not isinstance(v, dict)}))


if __name__ == '__main__':
    main()
