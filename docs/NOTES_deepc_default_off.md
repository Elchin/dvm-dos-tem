# deep_carbon binary default flipped to 0 (prod/deepc-default-off)

Date: 2026-11-20. Base: fix/snow-thermal-closure 40279738.
Decision (E. Jafarov, verbatim): "Flip the deep_carbon binary default to 0".

## What changed
- `src/ModelData.cpp`: `model_settings.get("deep_carbon", 1)` changed to `get("deep_carbon", 0)`. The no-config
  `ModelData()` constructor default also changed from 1 to 0. The FATAL message now names 0 as the default.
- `include/ModelData.h`, `src/Cohort.cpp`: comments only.
- `Soil_Bgc::deep_carbon` member default was already 0, so it is unchanged.
- `tests/unit/test_deepc.cpp`: the default test now asserts 0 for these cases: no-arg, key absent, and empty
  `model_settings`. A new test case asserts that `"deep_carbon": 1` in the config still gives 1 and that the
  value reaches Soil_Bgc.

## What did not change
- `config/config.js` (repo template) still sets `"deep_carbon": 1` explicitly. Runs made from the template keep
  deep C on. Its inline comment "default 1" now describes the template value, not the binary default.
- Configs without the key now run the legacy v5q burial path. This covers the v5q-gpp-fleaf2 production configs
  TVC S22, HPC S21H, T125 and S23. That path is bit-identical to v5q-gpp-fleaf2, so the identity check needs no
  explicit pin.
- The CMT06 X3 fl0.03 calibration (`~/cal-cmt06/j4/runs/R0_X3_fl03/ov.json`) sets `"deep_carbon": 1` explicitly,
  so it is not affected.

## Tests
VM conda gxx 13.4, C++20, Fireweed 31051b00 / Burnup 6d56e29b sources in external/ (untracked).
`[doctest] test cases: 477 | 477 passed | 0 failed`. `assertions: 3173006 | 3172998 passed | 8 failed`.
The 8 failures are the known xfail assertions (climate svp 2, richards_water_energy 2, seb 2, stefan_water 2).
