# Stage 1: upstream master sync and deep soil carbon (`wetland_deep_C`, phases 1-2)

[Back to the stage index](README.md)

## Purpose
Bring the fork's production line (`v5q-gpp-fleaf2` @ `29f26b5d`) up to date with upstream `master` (`c0fa75cb`), then add H. Genet's
`deep_carbon_v2` burial scheme as a switch. This is the first part of the "wetland + deep C" merge plan; the wetland merge itself
(phase 3) has **not started** and is not on this branch.

## Physics summary (deep C)
With `deep_carbon = 1`, `Soil_Bgc::deepCarbonTransfer()` (from `f8f8268f`):
- fibric -> humic: per-pool fractions (rawc 0.5, soma/sompr/somcr 1.0) of positive increments, spread over the whole humic horizon;
- humic -> mineral: per-pool fractions (0.5/1/1/1) of positive increments into the top 5.2 m;
- `initMslayerCarbon`: mineral SOM initialised to 5 m (legacy 2 m).

Conservation: the fork's `a59ce891` N fix (per-layer N at the donor C:N) is used under `n_som_move_fix = 1`; receiver weights are
normalised by the mineral thickness actually present; no mineral layer means no deep -> mineral transfer; `d2m_transfer_gate` 0/1/2 is
honoured (gate 1: f = eligible thickness / 5.2 m). With `deep_carbon = 0` the v5q burial code path runs textually unchanged.
The `rts-deepc` extra fixes (02d691cb, be82282d, ac022b55, 4b4f07c5, c5f4c45e) were not ported because v5q equivalents exist.

## Switches
| key (`model_settings`) | default | values | meaning |
|---|---|---|---|
| `deep_carbon` | **1** (ModelData and `config/config.js`) | 0 / 1, anything else FATAL | 1 = deep_carbon_v2 burial scheme; 0 = legacy v5q burial, bit-identical |

`deep_carbon` is the one stage switch that is **on by default**. Elchin decided this on 2026-11-09 ("Keep the binary default
`deep_carbon` = 1"). The `Soil_Bgc::deep_carbon` member itself defaults to 0, so hand-built objects in unit tests stay legacy.
To reproduce production (T125/S23) set `"deep_carbon": 0`.

## How to run
```jsonc
"model_settings": { "deep_carbon": 1 }   // 0 = production v5q burial
```

## Tests
- Unit: `tests/unit/test_deepc.cpp`, doctest suite `deepc` (6 cases, 242 assertions): parsing/default/FATAL, legacy untouched, per-pool
  hand check, C/N conservation over the gates, nfeed, `n_som_move_fix`, no-mineral column.
  Run: `make unittest` then `tests/unit/dvmdostem_unittest -ts=deepc`.
- Regression: `tests/regression/run_regression.sh` (see [tests/README.md](../../tests/README.md)), 20-yr spin-up + 75-yr main, 8 TVC/HPC cases.

## Key results (merge plan §8, 2026-11-09, calib-elchin)
- Phase 1 (`8d224ec7`): 26 files, no `src/` changes; v5q base -> phase 1 **PASS, bit-identical 8/8** on calib8 and prod
  (compare.py 3952 / 3976 variables identical, 0 differ).
- Phase 2 with `deep_carbon = 0`: **PASS, bit-identical 8/8** on calib8 and prod.
- Phase 2 with `deep_carbon = 1`, prod: PASS. C closure max 0.0024-0.0059 gC/m2/yr; N mean <= 1e-4, max <= 2e-4 gN/m2/yr;
  water -0.06...-0.53 % P; HPC nfeed safety-net events 53 -> 113 and 49 -> 105.
- Deep C effect (prod, last 10 yr of the 75-yr main run, deep_carbon 1 - 0): HPC nfeed GPP -124/-148 gC/m2/yr (-49/-47 %) and
  NEE -7.5 -> +6.4 and -6.1 -> +10.6 (sink -> source); TVC nfeed SHLWC -18 %, DEEPC +8 %, GPP about -1 %.
- Unit tests after phase 2: 275 cases; the same 6 known xfail assertions as the base.

## Known open issues
- Deep C is a structural soil C change: it **invalidates the T125/S23 calibration and restarts**. Recalibration belongs to
  Calibration-expert's Stage E.
- The burial fractions (0.5/1/1/1) are hard-coded upstream values marked "sensitivity analysis required".
- With N feedback HPC becomes strongly N-limited under deep C (GPP about -48 %).
- The regression starts from legacy restarts (95 yr of deep C dynamics, not an equilibrium); `initMslayerCarbon` at 5 m is not covered
  by a run test.

## Commits in this stage
The stage contains our three commits plus the 32 upstream commits that the two merges bring in (upstream authors;
they are listed so that every commit on the branch is documented).

### Our commits

| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `8d224ec7` | 2026-11-09 | Phase 1 master sync: merges upstream master c0fa75cb into wetland_deep_C (on top of production v5q-gpp-fleaf2). No src/ changes; brings upstream CI, docs and parameters/ updates (CMT13-15 added; CMT65/66/70/71 changed). | merge, vs first parent: 26 files: `parameters/` (8), `.github/workflows/` (6), `pyddt/tests/` (5), `docs_src/sphinx/` (3), `CITATION.cff`, `HOWTO_RELEASE.md`, `pyddt/src/pyddt/calibration/calibration_targets.py`, `scripts/ci/ensure-dev-image.sh` | merge plan §8 (Phase 1) |
| 2 | `527a3046` | 2026-11-09 | Phase 2a: merges upstream deep_carbon_v2 7f607479 for history only; src/Soil_Bgc.cpp kept as ours, so the merge changes no behaviour (empty diff vs first parent). | none (history-only merge, empty diff vs first parent) | merge plan §8 (Phase 2a) |
| 3 | `f8f8268f` | 2026-11-09 | Phase 2b: re-implements the deep_carbon_v2 burial scheme behind model_settings.deep_carbon (default 1; 0 = v5q burial, bit-identical); unit tests test_deepc.cpp. | `config/config.js`, `include/ModelData.h`, `include/Soil_Bgc.h`, `src/Cohort.cpp`, `src/ModelData.cpp`, `src/Soil_Bgc.cpp`, `tests/unit/test_deepc.cpp` | merge plan §8 (Phase 2b) |

### Imported upstream commits (32; not authored on the fork)

| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `b8f142f6` | 2026-01-17 | Upstream (hgenet): 1- Increasing soil column depth for carbon and nitrogen computation 2- Increasing MINEC initial value accordingly 3- Including layer depth in redistributing ORGN from shalow to deep horizons. [history only: 527a3046 kept the fork Soil_Bgc.cpp; the scheme is re-implemented behind deep_carbon in f8f8268f] | `src/Soil_Bgc.cpp` | imported via 527a3046 (deep_carbon_v2) |
| 2 | `a3afec59` | 2026-01-17 | Upstream (hgenet): Change quantity and redistribution of deep to mineral SOC [history only: 527a3046 kept the fork Soil_Bgc.cpp; the scheme is re-implemented behind deep_carbon in f8f8268f] | `src/Soil_Bgc.cpp` | imported via 527a3046 (deep_carbon_v2) |
| 3 | `d2b679a2` | 2026-02-13 | Upstream (hgenet): Changing the pools and burying rates between shallow, deep and mineral horizons. [history only: 527a3046 kept the fork Soil_Bgc.cpp; the scheme is re-implemented behind deep_carbon in f8f8268f] | `src/Soil_Bgc.cpp` | imported via 527a3046 (deep_carbon_v2) |
| 4 | `2c1e2b78` | 2026-02-13 | Upstream (hgenet): minor - indicating which parameters could be included in a sensitivity analysis evaluating the impact of changes in the burial scheme on SOC vertical distribution, RH and ALD [history only: 527a3046 kept the fork Soil_Bgc.cpp; the scheme is re-implemented behind deep_carbon in f8f8268f] | `src/Soil_Bgc.cpp` | imported via 527a3046 (deep_carbon_v2) |
| 5 | `7f607479` | 2026-02-16 | Upstream (hgenet): Renaming variables and parameters [history only: 527a3046 kept the fork Soil_Bgc.cpp; the scheme is re-implemented behind deep_carbon in f8f8268f] | `src/Soil_Bgc.cpp` | imported via 527a3046 (deep_carbon_v2) |
| 6 | `efa2118a` | 2026-10-03 | Upstream (vbriones): Updated calibration targets with all sites found within parameter files. Made minor changes primarily to calparbgc and bgcsoil to update recent model runs conducted 9/2026 | 9 files: `parameters/` (8), `pyddt/src/pyddt/calibration/calibration_targets.py` | imported via 8d224ec7 (upstream master) |
| 7 | `7a4b2173` | 2026-10-06 | Upstream (Tobey Carman): standardize PFT names, CMTs 13, 14, 15 | `parameters/cmt_bgcvegetation.txt`, `parameters/cmt_firepar.txt` | imported via 8d224ec7 (upstream master) |
| 8 | `fabacc73` | 2026-10-06 | Upstream (Tobey Carman): standardize CMT names, comments, 13, 14, 15 | `parameters/cmt_bgcsoil.txt`, `parameters/cmt_calparbgc.txt`, `parameters/cmt_dimground.txt`, `parameters/cmt_firepar.txt` | imported via 8d224ec7 (upstream master) |
| 9 | `fbe5834a` | 2026-10-06 | Upstream (Tobey Carman): standardize names in calibration targets; match to parameters | `pyddt/src/pyddt/calibration/calibration_targets.py` | imported via 8d224ec7 (upstream master) |
| 10 | `54a1b49b` | 2026-10-06 | Upstream (vbriones): Added CMT15 to dimvegetation.txt and added // to space between CMT12 and CMT13 in bgcsoil | `parameters/cmt_bgcsoil.txt`, `parameters/cmt_dimvegetation.txt` | imported via 8d224ec7 (upstream master) |
| 11 | `2f273fa7` | 2026-10-07 | Upstream (tobey): Merge pull request #842 from uaf-arctic-eco-modeling/paramupdates | merge, vs first parent: 9 files: `parameters/` (8), `pyddt/src/pyddt/calibration/calibration_targets.py` | imported via 8d224ec7 (upstream master) |
| 12 | `387b3c62` | 2026-10-08 | Upstream (Tobey Carman): update tests for recent parameter updates: | `pyddt/tests/doctests/doctests_MadsTEMDriver.rst`, `pyddt/tests/doctests/doctests_MadsTEMDriver_alltypestargets.rst`, `pyddt/tests/doctests/doctests_MadsTEMDriver_onlysoiltarget.rst`, `pyddt/tests/doctests/doctests_param_util.md` | imported via 8d224ec7 (upstream master) |
| 13 | `ff67f8a3` | 2026-10-08 | Upstream (Tobey Carman): Fix path, values to testing data for run mask test. | `pyddt/tests/doctests/util.runmask.md` | imported via 8d224ec7 (upstream master) |
| 14 | `2b67be7c` | 2026-10-08 | Upstream (Tobey Carman): Version bump docs, update citation | `CITATION.cff`, `docs_src/sphinx/source/conf.py` | imported via 8d224ec7 (upstream master) |
| 15 | `5c1fbdd5` | 2026-10-09 | Upstream (Tobey Carman): add new github workflow for doctests | `.github/workflows/doctests.yml` | imported via 8d224ec7 (upstream master) |
| 16 | `c84304d0` | 2026-10-09 | Upstream (Tobey Carman): Improve CI efficiency... | `.github/workflows/checkout_build_run.yml`, `.github/workflows/doctests.yml`, `scripts/ci/ensure-dev-image.sh` | imported via 8d224ec7 (upstream master) |
| 17 | `ba7c239f` | 2026-10-09 | Upstream (Tobey Carman): fix bug with GIU, UID in ensure-dev-image.sh | `scripts/ci/ensure-dev-image.sh` | imported via 8d224ec7 (upstream master) |
| 18 | `79b36c7d` | 2026-10-09 | Upstream (Tobey Carman): Update comment on doctest workflow. | `.github/workflows/doctests.yml` | imported via 8d224ec7 (upstream master) |
| 19 | `faba6116` | 2026-10-09 | Upstream (Tobey Carman): remove unused (not comitted yet) test | `.github/workflows/doctests.yml` | imported via 8d224ec7 (upstream master) |
| 20 | `522d5b56` | 2026-10-09 | Upstream (Tobey Carman): remove qcal test which fails on first pass in a clean environment. | `.github/workflows/doctests.yml` | imported via 8d224ec7 (upstream master) |
| 21 | `b85e613b` | 2026-10-12 | Upstream (Tobey Carman): parameterize the ensure dev image script for future support of other images. | `scripts/ci/ensure-dev-image.sh` | imported via 8d224ec7 (upstream master) |
| 22 | `7661c3fa` | 2026-10-12 | Upstream (Tobey Carman): update the testing workflows to use better caching. | `.github/workflows/checkout_build_run.yml`, `.github/workflows/doctests.yml` | imported via 8d224ec7 (upstream master) |
| 23 | `e3a17682` | 2026-10-12 | Upstream (Tobey Carman): re-factor tests/workflow stuff to be more composable. | `.github/workflows/build_images.yml`, `.github/workflows/compile.yml`, `.github/workflows/smoke_test.yml` | imported via 8d224ec7 (upstream master) |
| 24 | `6292e32b` | 2026-10-12 | Upstream (Tobey Carman): fix env variable export in compile step. | `.github/workflows/compile.yml` | imported via 8d224ec7 (upstream master) |
| 25 | `1347aa42` | 2026-10-13 | Upstream (Tobey Carman): update branch triggers for github workflows. | `.github/workflows/build_images.yml`, `.github/workflows/compile.yml`, `.github/workflows/doctests.yml`, `.github/workflows/smoke_test.yml` | imported via 8d224ec7 (upstream master) |
| 26 | `6239d2c6` | 2026-10-13 | Upstream (Tobey Carman): update docs section about branching. | `docs_src/sphinx/source/software_development_info.rst` | imported via 8d224ec7 (upstream master) |
| 27 | `eecf1723` | 2026-10-13 | Upstream (Tobey Carman): update release info: new worflow, update docs, howto file. | `.github/workflows/release.yml`, `HOWTO_RELEASE.md`, `docs_src/sphinx/source/software_development_info.rst` | imported via 8d224ec7 (upstream master) |
| 28 | `9f4a570e` | 2026-10-13 | Upstream (Tobey Carman): update docs: examples, with respect to develop branch. | `docs_src/sphinx/source/examples_and_tutorials/make_contribution.rst` | imported via 8d224ec7 (upstream master) |
| 29 | `d47cb1d8` | 2026-10-13 | Upstream (Tobey Carman): update docs: section on continuous integration | `docs_src/sphinx/source/software_development_info.rst` | imported via 8d224ec7 (upstream master) |
| 30 | `c96558fe` | 2026-10-13 | Upstream (Tobey Carman): add warning about out of date slides in docs, a few minor fixes. | `docs_src/sphinx/source/examples_and_tutorials/make_contribution.rst` | imported via 8d224ec7 (upstream master) |
| 31 | `920d42ef` | 2026-10-13 | Upstream (tobey): Merge pull request #846 from uaf-arctic-eco-modeling/refactor-ci | merge, vs first parent: `.github/workflows/build_images.yml`, `.github/workflows/checkout_build_run.yml`, `.github/workflows/compile.yml`, `.github/workflows/doctests.yml`, `.github/workflows/smoke_test.yml`, `scripts/ci/ensure-dev-image.sh` | imported via 8d224ec7 (upstream master) |
| 32 | `c0fa75cb` | 2026-10-13 | Upstream (tobey): Merge pull request #847 from uaf-arctic-eco-modeling/refactor-ci-2nd-pass | merge, vs first parent: `.github/workflows/build_images.yml`, `.github/workflows/compile.yml`, `.github/workflows/doctests.yml`, `.github/workflows/release.yml`, `.github/workflows/smoke_test.yml`, `HOWTO_RELEASE.md`, `docs_src/sphinx/source/examples_and_tutorials/make_contribution.rst`, `docs_src/sphinx/source/software_development_info.rst` | imported via 8d224ec7 (upstream master) |
