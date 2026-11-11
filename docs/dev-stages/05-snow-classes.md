# Stage 5: snow-class scheme (Sturm et al. 1995 classes as in SiBCASA, Schaefer et al. 2009)

[Back to the stage index](README.md)

## Purpose and physics summary
An optional snowpack structure in place of dvmdostem's single new-snow density / max density: per-class new-snow density from air
temperature and wind (Schaefer et al. 2009 eqs 4-7), Anderson 1976 / CLM4 compaction, a two-layer stratification (eqs 1-3) with a
SWE- and heat-conserving remap, and a fixed depth-hoar conductivity (tundra 0.18, taiga 0.072 W/m/K). Approved by Elchin as a
default-off option with tests. Classes come from a CMT map in the new file `parameters/snow_classes.json` (taiga for forest/bog;
tundra for tundra/fen/bare): TVC (CMT51) is tundra (wind term on), HPC (CMT01) is taiga (no wind term).

## Switches
| key | default | meaning |
|---|---|---|
| `snow_scheme` | **0** (legacy, bit-identical) | 1 = snow classes; anything else FATAL |
| `snow_class_file` | `""` -> `<parameter_dir>/snow_classes.json` | class table + CMT -> class map |
| `snow_class_override` | `""` (use the map) | force one class for all cohorts |
| `snow_class_wind` | -1 (daily forcing wind when the file has it, else 0 m/s with one warning) | >= 0 = constant wind (m/s) for new-snow density |
| `snow_diag` | 0 | 1 writes `output/snow_daily.csv` (diagnostic only) |

Under `snow_scheme = 1` the CMT parameters `snwdennew` and `snwdenmax` are inactive, and snow `maxice` is RHO_CAP*dz (`2ffa807f`).
Tunable class parameters in the JSON: fbotmax, Dmin, Dmax, rho_topref, rho_botref, rho_bulkobs, sigma_bulkobs, k_dh, wind.
Densification constants are hard-coded.

## How to run
```jsonc
"model_settings": { "snow_scheme": 1, "snow_diag": 1 }
```
Do not edit `parameters/snow_classes.json` in place for experiments; point `snow_class_file` at a copy.

## Tests
- `tests/unit/test_snowclass.cpp`, suite `snow_classes` (12 cases, 33,156 assertions): defaults, hand values for eqs 4-7, eqs 1-3 vs
  Table 2, eqs 1-3 with the Fig 4/5 values reproduce the figure curves, compaction, SWE + heat conserving remap, Ground
  stratify/compaction/melt, class-file errors/override, maxice bound, column test. Run: `tests/unit/dvmdostem_unittest -ts=snow_classes`.
- Bit-identity: scheme 0 / key absent is data-identical (8/8) to the `enthalpy/ssf-rdt` runs for both solvers; `snow_diag 1` adds only
  `snow_daily.csv`.

## Key results (snowclass report, 2026-11-10)
- North Slope, CRU monthly (no wind): April bulk density rises from ~209 kg/m3 (pinned at CMT `snwdenmax`) to 240-290 (obs 202-301) at the
  same SWE; DJFM bulk k 0.067 -> 0.089-0.113 W/m/K.
- Nicolsky 2017 gated checks (23 sites, 144 gates): legacy 30 -> 36, new solver 47 -> 71; network MAGT bias new solver +3.6 -> +2.6 C.
- With real daily wind (NASA POWER): Barrow DJF T0.10 bias +14.7 -> +5.1 C, Deadhorse +8.4 -> +3.1 C (legacy).
- TVC/HPC daily forcing: close to neutral; legacy obs_regression 21/21 -> 19/21 (HPC DJF T4cm +3.7 -> +5.4 C with the taiga class),
  new solver 7/21 -> 7/21. Energy closure on new-solver runs <= 4.6e-8 W/m2.
- Bug found and fixed (`2ffa807f`): dense snow-class layers were clipped to `snwdenmax*dz` in `Stefan::cleanSnowSoilLayers`, deleting SWE.

## Known open issues
- The paper's Figs 4/5 "tundra" curves are not drawn with its Table 2 parameters; the model ships Table 2. The Fig 4/5 sensitivity
  (no code change) is better under CRU but mixed to worse under real daily wind.
- k_dh 0.18 for tundra is high vs TVC depth-hoar measurements (0.05 needle-probe, 0.09 indurated; Dutch et al. 2022). k_dh 0.05 gives
  DJFM bulk k 0.080 = TVC measured median. Candidate Stage E tunable (Calibration-expert).
- No blowing-snow redistribution; the CRU North Slope inputs have no wind; no per-gridcell class map yet.
- Any adoption needs Calibration-expert's per-CMT recalibration; the default stays 0.

## Commits in this stage
| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `fda53090` | 2026-11-10 | Adds model_settings.snow_scheme (default 0 = legacy, bit-identical): Sturm et al. 1995 snow classes as in SiBCASA (Schaefer et al. 2009) - new-snow density, compaction, two-layer stratification, fixed depth-hoar k; new parameters/snow_classes.json; snow_diag csv; unit tests test_snowclass.cpp. | 14 files: `src/` (6), `include/` (5), `Makefile`, `parameters/snow_classes.json`, `tests/unit/test_snowclass.cpp` | n/a (snowclass report) |
| 2 | `2ffa807f` | 2026-11-10 | snow_scheme 1 fix: snow maxice = RHO_CAP*dz, because snwdenmax*dz clipped dense snow-class layers in Stefan::cleanSnowSoilLayers and deleted SWE; unit tests for the maxice bound and the paper's Figs 4-5 values. | `src/Ground.cpp`, `tests/unit/test_snowclass.cpp` | n/a (snowclass report) |
