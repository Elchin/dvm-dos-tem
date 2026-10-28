#!/usr/bin/env python3
"""Harness test (v5o) for closure.py's N-check logic: the N tolerance is enforced only when the
N1-N3 conservation fixes are all on; otherwise the residual is a warning naming the three bugs.
Synthetic run dirs (meta.json + main/config/config.js); closure() is stubbed, no model run.
usage: python3 tests/regression/test_closure.py   (exit 0 = all tests pass)"""
import json, os, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import closure  # noqa: E402

def mkrun(root, name, ms):
    d = Path(root) / name; (d / "main/config").mkdir(parents=True)
    (d / "meta.json").write_text(json.dumps(dict(case=name, status="OK", main=dict(run_status=[100], fail_logs=[],
        soil_safetynet_events=0, safetynet_csv_rows=0, veg_giveups=0))))
    cfg = {"general": {}, "model_settings": dict(ms)}
    (d / "main/config/config.js").write_text("// comment line\n" + json.dumps(cfg, indent=1))
    return d

class T(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = self.tmp.name
        self._c = closure.closure; self.nres = 5.0     # large N residual (> ntol 1.0)
        closure.closure = lambda rd: dict(n_years=75, W_resid_pctP=0.0, C_resid_max_abs_annual=0.0,
                                          N_resid_mean_gm2yr=self.nres, N_resid_max_abs_annual=self.nres)
        os.environ.pop("REG_N_DEFAULTS", None)
    def tearDown(self): closure.closure = self._c; self.tmp.cleanup(); os.environ.pop("REG_N_DEFAULTS", None)

    def test_defaults_v5o_enforced(self):          # no keys + v5o binary defaults -> enforced -> FAIL
        r = closure.evaluate(mkrun(self.root, "a", {}))
        self.assertEqual(r["N_check"], "enforced"); self.assertFalse(r["closure_pass"])
        self.assertTrue(any("N |mean resid|" in f for f in r["closure_fail"]))
    def test_explicit_on_small_residual_passes(self):
        self.nres = 1e-4
        r = closure.evaluate(mkrun(self.root, "b", dict(n_rootless_fix=2, n_som_move_fix=1, n_nonvasc_extract=1)))
        self.assertEqual(r["N_check"], "enforced"); self.assertTrue(r["closure_pass"]); self.assertEqual(r["closure_warn"], [])
    def test_rootless_1_counts_as_on(self):
        r = closure.evaluate(mkrun(self.root, "c", dict(n_rootless_fix=1)))
        self.assertEqual(r["N_check"], "enforced")
    def test_any_fix_off_is_warning_only(self):
        for k in ("n_rootless_fix", "n_som_move_fix", "n_nonvasc_extract"):
            r = closure.evaluate(mkrun(self.root, "d_" + k, {k: 0}))
            self.assertEqual(r["N_check"], "warn", k); self.assertTrue(r["closure_pass"], k)
            self.assertEqual(len(r["closure_warn"]), 1)
            for bug in ("N1 n_rootless_fix", "N2 n_som_move_fix", "N3 n_nonvasc_extract"): self.assertIn(bug, r["closure_warn"][0])
    def test_legacy_binary_defaults(self):          # no keys with a pre-v5o binary -> legacy -> warning
        r = closure.evaluate(mkrun(self.root, "e", {}), n_defaults="legacy"); self.assertEqual(r["N_check"], "warn")
        os.environ["REG_N_DEFAULTS"] = "legacy"
        r = closure.evaluate(mkrun(self.root, "f", {})); self.assertEqual(r["N_check"], "warn")
        r = closure.evaluate(mkrun(self.root, "g", dict(n_rootless_fix=2, n_som_move_fix=1, n_nonvasc_extract=1)))
        self.assertEqual(r["N_check"], "enforced")   # explicit keys override the binary defaults
    def test_water_carbon_still_fail(self):         # other checks unchanged when N is warn-only
        closure.closure = lambda rd: dict(n_years=75, W_resid_pctP=2.0, C_resid_max_abs_annual=0.0,
                                          N_resid_mean_gm2yr=5.0, N_resid_max_abs_annual=5.0)
        r = closure.evaluate(mkrun(self.root, "h", dict(n_som_move_fix=0)))
        self.assertFalse(r["closure_pass"]); self.assertTrue(any("water" in f for f in r["closure_fail"]))

if __name__ == "__main__":
    unittest.main(verbosity=2)
