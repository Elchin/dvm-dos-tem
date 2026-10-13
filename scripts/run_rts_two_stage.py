#!/usr/bin/env python3
"""Two-stage dvmdostem workflow for RTS (retrogressive thaw slump) experiments.

Mirrors the restart-based approach used on the original thermokarst (RTS)
branch (spin up once, then `--restart`-style continuation runs), using the
standard `IO.restart_from` mechanism and one output directory per step.

  Stage 1  (original model, RTS/thermokarst disabled)
           pre-run + equilibrium + spin-up:  -p P -e E -s S   (default 100/1000/200)
           -> stage1/restart-sp.nc
  Stage 2a historical transient (RTS disabled), restart from stage 1
           -t H                                 -> stage2a_hist/restart-tr.nc
  Stage 2b scenario stage, restart from 2a, thermokarst module ON:
           * ground-ice (excess-ice) initialization in the first month
           * short spin-up: the scenario years before the RTS initiation year
             (no initiation is forced before then by the projected
             thermokarst file), so the thermal state adjusts to the ice profile
           * RTS initiation in the year flagged in proj_exp_thermokarst_file
           -n N                                 -> stage2b_rts/
  `--mode baseline` runs the same stage 2b with `thermokarst_initiation=false`
  (identical ice initialization / CH4 diagnostics, no slump), for blending.

Only the standard library is used. Example (synthetic one-cell case):
  python3 scripts/run_rts_two_stage.py --workdir /tmp/rts2 \
      --config-overrides testing-data/rts-synthetic/two_stage_overrides.json
"""
import argparse, json, os, re, shutil, subprocess, sys, time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

def load_js(path):
    out = []
    for line in open(path).read().splitlines():
        # strip // comments outside of strings
        res, instr = [], False
        i = 0
        while i < len(line):
            c = line[i]
            if c == '"' and (i == 0 or line[i - 1] != '\\'):
                instr = not instr
            if not instr and line.startswith('//', i):
                break
            res.append(c); i += 1
        out.append(''.join(res))
    return json.loads(re.sub(r',(\s*[}\]])', r'\1', '\n'.join(out)))

def deep_update(d, u):
    for k, v in u.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            deep_update(d[k], v)
        else:
            d[k] = v
    return d

def abspath_io(cfg, base_dir):
    for k, v in cfg['IO'].items():
        if isinstance(v, str) and v and not os.path.isabs(v) and (k.endswith('_file') or k.endswith('_dir')):
            cfg['IO'][k] = os.path.join(base_dir, v) + ('/' if v.endswith('/') else '')

def run_step(binary, cfg, outdir, years, log_level, name):
    os.makedirs(outdir, exist_ok=True)
    cfg['IO']['output_dir'] = outdir.rstrip('/') + '/'
    path = os.path.join(outdir, 'config.js')
    json.dump(cfg, open(path, 'w'), indent=2)
    cmd = [binary, '-f', path, '-l', log_level] + years
    t = time.time()
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    dt = time.time() - t
    print('[%s] %s  (%.1f s)' % (name, ' '.join(cmd[3:]), dt), flush=True)
    if p.returncode != 0:
        sys.stderr.write(p.stdout[-3000:] + p.stderr[-3000:])
        raise SystemExit('step %s failed (rc=%d)' % (name, p.returncode))
    return dt

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default=os.path.join(REPO, 'dvmdostem'))
    ap.add_argument('--config', default=os.path.join(REPO, 'config', 'config.js'), help='base config.js')
    ap.add_argument('--config-overrides', default=None, help='JSON file deep-merged into the base config')
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--pr', type=int, default=100); ap.add_argument('--eq', type=int, default=1000)
    ap.add_argument('--sp', type=int, default=200)
    ap.add_argument('--hist', type=int, default=115, help='historical transient years (stage 2a)')
    ap.add_argument('--sc', type=int, default=30, help='scenario years incl. short spin-up (stage 2b)')
    ap.add_argument('--mode', choices=['rts', 'baseline'], default='rts')
    ap.add_argument('--ch4', action='store_true', help='enable reduced-form CH4 diagnostic (rts-ch4 only)')
    ap.add_argument('--stage1-restart', default=None, help='reuse an existing stage-1 restart-sp.nc')
    ap.add_argument('--log-level', default='fatal')
    a = ap.parse_args()

    cfg0 = load_js(a.config)
    abspath_io(cfg0, REPO)
    if a.config_overrides:
        ov = json.load(open(a.config_overrides))
        base = os.path.dirname(os.path.abspath(a.config_overrides))
        for k, v in ov.get('IO', {}).items():
            if isinstance(v, str) and v and not os.path.isabs(v) and (k.endswith('_file') or k.endswith('_dir')):
                ov['IO'][k] = os.path.join(base, v) + ('/' if v.endswith('/') else '')
        deep_update(cfg0, ov)
    ms = cfg0.setdefault('model_settings', {})
    w = os.path.abspath(a.workdir)
    os.makedirs(w, exist_ok=True)
    times = {}

    # ---- Stage 1: original model behaviour (RTS off)
    if a.stage1_restart:
        sp_restart = os.path.abspath(a.stage1_restart)
    else:
        c = json.loads(json.dumps(cfg0)); c['IO']['restart_from'] = ''
        c['model_settings']['thermokarst'] = False
        times['stage1'] = run_step(a.binary, c, os.path.join(w, 'stage1'),
                                   ['-p', str(a.pr), '-e', str(a.eq), '-s', str(a.sp), '-t', '0', '-n', '0'],
                                   a.log_level, 'stage1 pr/eq/sp')
        sp_restart = os.path.join(w, 'stage1', 'restart-sp.nc')

    # ---- Stage 2a: historical transient, RTS off
    c = json.loads(json.dumps(cfg0)); c['IO']['restart_from'] = sp_restart
    c['model_settings']['thermokarst'] = False
    times['stage2a'] = run_step(a.binary, c, os.path.join(w, 'stage2a_hist'),
                                ['-p', '0', '-e', '0', '-s', '0', '-t', str(a.hist), '-n', '0'],
                                a.log_level, 'stage2a historical')

    # ---- Stage 2b: ice init + short spin-up + RTS (scenario stage)
    c = json.loads(json.dumps(cfg0))
    c['IO']['restart_from'] = os.path.join(w, 'stage2a_hist', 'restart-tr.nc')
    c['model_settings']['thermokarst'] = True
    c['model_settings']['thermokarst_initiation'] = (a.mode == 'rts')
    if a.ch4: c['model_settings']['thermokarst_ch4'] = True
    c['stage_settings']['sc']['dsb'] = True   # thermokarst runs in the disturbance step
    c['IO']['output_nc_sc'] = 1
    times['stage2b'] = run_step(a.binary, c, os.path.join(w, 'stage2b_' + a.mode),
                                ['-p', '0', '-e', '0', '-s', '0', '-t', '0', '-n', str(a.sc)],
                                a.log_level, 'stage2b %s' % a.mode)
    json.dump(times, open(os.path.join(w, 'timings.json'), 'w'), indent=2)
    print('timings (s):', json.dumps({k: round(v, 1) for k, v in times.items()}))

if __name__ == '__main__':
    main()
