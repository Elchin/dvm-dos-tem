#!/usr/bin/env python3
"""Stage 1 (pre-run/equilibrium/spin-up, RTS off) for every Gydan 10x10 land cell.

Phase 1 of the two-case RTS susceptibility experiment (case A: Jia et al. 2026,
Sci. Data, doi:10.1038/s41597-026-08284-y; case B: Nesterova et al. 2025 West
Siberia RTS inventory). Stage 1 is shared by the baseline and all
single-start-year RTS cohort runs of a cell, so it is run once per cell here
with `scripts/run_rts_two_stage.py --stop-after stage1` (-p 100 -e 1000 -s 200).

Every cell gets its own workdir <out>/<cell>/ with overrides.json (absolute input
paths of inputs/<cell>/two_stage_overrides.json), output_spec.csv, a wrapper that
sets NBAL_DIAG=1 (monthly N balance to nbal.txt, diagnostic only), run.log and the
driver's timings.json (+ CPU seconds of the child processes). All model switches
stay at the code defaults (n_nonvasc_extract, n_extract_limiter,
soil_giveup_diag off). A summary of per-cell times goes to <out>/timings.json.
Resumable: cells whose timings.json and stage1/restart-sp.nc exist are skipped.

usage: run_stage1.py --inputs DIR --out DIR --binary BIN [--nproc N] [--cells c00,c01]
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'tests', 'rts'))


def outputs():
    import rtslib as r
    return {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS, **r.CH4_OUTPUTS,
            **r.REGROW_OUTPUTS, 'AVLN': 'm', 'ORGN': 'm', 'VEGNTOT': 'm', 'ALD': 'y'}


def make_overrides(data, w, repo=REPO):
    """overrides.json for one cell: the cell's two_stage_overrides.json with absolute IO paths."""
    import rtslib as r
    sf = os.path.join(w, 'output_spec.csv'); r.write_outspec(sf, outputs())
    ov = json.load(open(os.path.join(data, 'two_stage_overrides.json')))
    for k, v in ov['IO'].items():
        if isinstance(v, str) and v.endswith('.nc'): ov['IO'][k] = os.path.join(data, v)
    ov['IO']['output_spec_file'] = sf
    ov['IO']['parameter_dir'] = os.path.join(repo, 'parameters') + '/'
    p = os.path.join(w, 'overrides.json'); json.dump(ov, open(p, 'w'), indent=1)
    return p


def done(w):
    return os.path.exists(os.path.join(w, 'timings.json')) and \
        os.path.exists(os.path.join(w, 'stage1', 'restart-sp.nc'))


def run_cell(args):
    cell, inputs, out, binary, pr, eq, sp = args
    w = os.path.join(out, cell)
    if done(w):
        return cell, 'skip', json.load(open(os.path.join(w, 'timings.json')))
    os.makedirs(w, exist_ok=True)
    ovp = make_overrides(os.path.join(inputs, cell), w)
    nb = os.path.join(w, 'nbal.txt'); open(nb, 'w').close()
    wr = os.path.join(w, 'wrap.sh')
    open(wr, 'w').write(f'#!/bin/bash\nNBAL_DIAG=1 {binary} "$@" 2>>{nb}\n'); os.chmod(wr, 0o755)
    c0 = os.times(); t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'run_rts_two_stage.py'),
                        '--binary', wr, '--workdir', w, '--config-overrides', ovp,
                        '--pr', str(pr), '--eq', str(eq), '--sp', str(sp), '--stop-after', 'stage1'],
                       capture_output=True, text=True)
    c1 = os.times()
    open(os.path.join(w, 'run.log'), 'w').write(p.stdout + p.stderr)
    if p.returncode:
        return cell, 'fail rc=%d' % p.returncode, {}
    tp = os.path.join(w, 'timings.json'); t = json.load(open(tp))
    t['wall_total_s'] = time.time() - t0
    t['cpu_children_s'] = (c1.children_user - c0.children_user) + (c1.children_system - c0.children_system)
    t['binary'] = binary
    json.dump(t, open(tp, 'w'), indent=2)
    return cell, 'ok', t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', required=True, help='spatial10x10/inputs (grid.json + cNN/)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--binary', required=True)
    ap.add_argument('--nproc', type=int, default=os.cpu_count())
    ap.add_argument('--cells', default=None, help='comma list (default: all cells in grid.json)')
    ap.add_argument('--pr', type=int, default=100); ap.add_argument('--eq', type=int, default=1000)
    ap.add_argument('--sp', type=int, default=200)
    a = ap.parse_args()
    cells = a.cells.split(',') if a.cells else [c['cell'] for c in json.load(open(os.path.join(a.inputs, 'grid.json')))['cells']]
    os.makedirs(a.out, exist_ok=True)
    st = os.path.join(a.out, 'status.log')
    log = lambda m: (open(st, 'a').write(time.strftime('%Y-%m-%d %H:%M:%S %Z ') + m + '\n'), print(m, flush=True))
    log(f'start: {len(cells)} cells, nproc {a.nproc}, -p {a.pr} -e {a.eq} -s {a.sp}, binary {a.binary}')
    T0 = time.time(); res = {}
    jobs = [(c, a.inputs, a.out, a.binary, a.pr, a.eq, a.sp) for c in cells]
    with ProcessPoolExecutor(a.nproc) as ex:
        for cell, status, t in ex.map(run_cell, jobs):
            res[cell] = dict(status=status, **t); log(f'{cell} {status} stage1={t.get("stage1", float("nan")):.1f}s')
    wall = time.time() - T0
    agg = os.path.join(a.out, 'timings.json')
    old = json.load(open(agg)) if os.path.exists(agg) else {}
    cellsd = {**old.get('cells', {}), **res}
    s1 = sorted(v['stage1'] for v in cellsd.values() if 'stage1' in v)
    json.dump(dict(cells=cellsd, nproc=a.nproc, last_invocation_wall_s=wall,
                   stage1_s=dict(n=len(s1), min=s1[0] if s1 else None, median=s1[len(s1) // 2] if s1 else None,
                                 max=s1[-1] if s1 else None, sum=sum(s1))), open(agg, 'w'), indent=1)
    log(f'end: wall {wall:.0f}s, ok {sum(v["status"] in ("ok", "skip") for v in res.values())}/{len(cells)}')


if __name__ == '__main__':
    main()
