# Development stages on the fork (TEM-dev branches)

This folder documents, stage by stage, every commit that branch **`fix/snow-empty-layers`** adds on top of the production line
`v5q-gpp-fleaf2` (`29f26b5d`). One README per stage; each stage README has a commit table listing every commit of that stage
(short hash, date, plain summary, files touched, design-doc section). Branch code tip documented here: `6c0e7849`
(60 commits since `29f26b5d`).

Rules that apply to every stage: new options are model_settings switches that default to the legacy behaviour and are bit-identical when
off (exceptions: `deep_carbon` defaults to 1, and `thermal_substeps` defaults to 4 under `thermal_solver = 1`, both decided by Elchin);
repo `parameters/` files are never edited in place; test runs use n <= 75 and <= 4 cores on calib-elchin.

## Stages on this branch
| # | stage README | commits | switches added (default) |
|---|---|---|---|
| 1 | [Upstream master sync + deep soil carbon](01-upstream-sync-deep-carbon.md) | 35 (3 ours, 32 imported upstream) | `deep_carbon` (default **1**) |
| 2 | [Energy-closure diagnostic (legacy thermal step)](02-energy-closure-diagnostic.md) | 3 (3 ours) | `energy_diag` (0) |
| 3 | [Enthalpy thermal solver (soil + snow thermal)](03-enthalpy-thermal-solver.md) | 18 (18 ours) | `thermal_solver` (0), `thermal_sfcc` (1), `thermal_substeps` (4 under `thermal_solver` 1), `thermal_face_mode` (0), `thermal_nl_solver` (0), `geothermal_flux_mode` (0), `geothermal_flux` (0.0), `snow_empty_cleanup` (0) |
| 4 | [Richards step-halving fix + water/energy tests](04-richards-soil-hydrology.md) | 4 (4 ours) | `richards_dt_fix` (0) |

Coverage: 60 of 60 commits in `29f26b5d..6c0e7849` appear in exactly one stage README; missing 0, duplicated 0, unknown 0 (checked by script when this file was generated).

## How the stages stack
```
v5q-gpp-fleaf2 @ 29f26b5d  (production; v5 series, see tests/README.md and docs/VERSION_v5q_TVC_HPC_MIP_S22_S21H.md)
 └─ Stage 1  upstream master sync (8d224ec7) + deep C (527a3046, f8f8268f)        branch wetland_deep_C
     └─ Stage 2  energy_diag (517ab675, 1dca209e, a6c718d4)                         branch wetland_deep_C
         └─ Stage 3  enthalpy solver: B -> C -> C2/C2b -> D -> sub-step fix -> default 4 (f8617974)
             └─ Stage 4  Richards tests + richards_dt_fix -> enthalpy/ssf-rdt @ 5c08f890 (Stage E baseline code)
                 └─ Stage 3  melt-out closure fix 67cc8a1e                       (fix/meltout-closure)
                     └─ Stage 3  snow_empty_cleanup d042da83, 6c0e7849   <- this branch
```
Commits are grouped by topic, not strictly by time: Stage 3 and Stage 4 interleave on 2026-11-10 (`f8617974` Stage 3, then
`582cfc3b`...`5c08f890` Stage 4, then the Stage 3 snow fixes on the tip branches).

## Branch map (all on github.com/Elchin/dvm-dos-tem)
| branch | code tip | contains | stage READMEs |
|---|---|---|---|
| `enthalpy/ssf-rdt` | `5c08f890` | Stages 1-4 | (no docs commit on this branch) |
| `fix/meltout-closure` | `67cc8a1e` | Stages 1-4 + melt-out fix | 1, 2, 3, 4 |
| `fix/snow-empty-layers` | `6c0e7849` | Stages 1-4 + melt-out fix + `snow_empty_cleanup` | 1, 2, 3, 4 |
| `feature/snow-classes` | `2ffa807f` | Stages 1-4 + Stage 5 | 1-5 |
| `fix/snow-empty-layers-snw` | `d4f86733` | everything above (most complete) | 1-5 |

Each branch carries only the READMEs and commit rows for commits it actually contains.

## Configurations
- **Production legacy** (T125/S23 calibration valid): `thermal_solver = 0` (default), `richards_dt_fix = 0` (default), `deep_carbon = 0`
  (the binary default is 1), `snow_empty_cleanup = 0`.
- **Stage E baseline (new solver):** `fix/snow-empty-layers-snw` with `thermal_solver = 1`, `thermal_sfcc = 1`, `thermal_substeps` 4, `richards_dt_fix = 1`, `snow_empty_cleanup = 1`, `snow_scheme = 0` (Elchin, 2026-11-10; Calibration-expert switches after their verification checks pass). On `enthalpy/ssf-rdt` alone the baseline is the same without `snow_empty_cleanup`.

## Running the tests
```bash
./build-conda.sh -j8 unittest                 # or: make unittest ... (see tests/README.md)
tests/unit/dvmdostem_unittest -ts=enthalpy_integration   # one suite; -lts lists suites
tests/regression/run_regression.sh -j 4 BIN_A BIN_B CASESET OUTDIR
python3 tests/regression/closure.py RUNDIR --etol-enth 1e-6
```
Suites per stage are listed in each stage README. Harness details: [tests/README.md](../../tests/README.md).

## Related documents
- Enthalpy solver design and measurements: [`tools/enthalpy_bench/enthalpy_solver_design.md`](../../tools/enthalpy_bench/enthalpy_solver_design.md)
  (§15a-§15g). §15h (empty snow-layer cleanup) is in [Stage 3](03-enthalpy-thermal-solver.md).
- Benchmark tools: [`tools/enthalpy_bench/README.md`](../../tools/enthalpy_bench/README.md).
- Production (v5 series) switches and references: [`tests/README.md`](../../tests/README.md), [`docs/daily-mode-REFERENCES.md`](../daily-mode-REFERENCES.md).

## Work that is not on these branches
- Validation suites (local branches on calib-elchin, not pushed): `tests/northslope-valid` (Nicolsky et al. 2017 North Slope checks),
  `tests/unfrozen-valid` (Nicolsky et al. 2007 unfrozen water), `tests/snow-checks` (North Slope snow observations, daily forcing),
  and the Happy Valley / R&O 2000 benchmark on `tests/richards-water-energy`.
- Not started: wetland merge (phase 3), Reburn with the fire switch, RTS; recalibration (Calibration-expert, Stage E).
