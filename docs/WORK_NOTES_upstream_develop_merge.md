# Upstream `develop` merge into the fork stack (TDDD-7, 2026-11-21)

Brief (E. Jafarov): "Let's do the upstream merge on my fork first". Local only: nothing pushed, no PR, no production change, no GCS writes, no restart reuse.

## Result
- **Branch** `merge/upstream-develop-20261121` (calib-elchin worktree `~/tem-dev/upmerge/wt`), off `fix/restart-state-carry` 66cfec02.
  - `27c723f4` is the merge commit (parents 66cfec02 and upstream develop c7a9fb27). It is a real merge, not a rebase.
  - `13b2f1a9` is a follow-up tooling fix in `setup_working_directory.py`.
  - The notes commit sits on top of these (see `git log`).
- **Identity, all switches off: bit-identical to 66cfec02.** 0 differing variables in 8,440 compared (2,256 files). Covers the harness `prod` set (spin20 + n75) and the full TVC S22 / HPC S21H production mains (spin20 + n150). Determinism B==R holds everywhere.
- **Unit suite:** 482/482 test cases pass. The 8 failing assertions are the same 8 known ones (test_climate:47/52, test_richards_water_energy:198/215, test_seb:177/178, test_stefan_water:142/271).
- **Cost:** 0.681 core-h, under the 1.0 cap (ledger `~/tem-dev/upmerge/LEDGER.txt`).

## 1. Upstream inventory (`66cfec02..c7a9fb27`)
79 commits, 64 of them non-merge (Ruth Rutter 37, Tobey Carman 27). The merge base is c0fa75cb (= upstream master).

| theme | commits (examples) | author | touches model? |
|---|---|---|---|
| pyddt refactor/tooling: util.config, util.test, CMT helpers, outspec helpers, SUCCESS reporting, input checks incl. projected CO2, sim-comparison report moved into pyddt, usage interface | 6567da22, a1a69e78, b4bc8c39, 7b64fc3b, 807ca4ac, bc1715dd, 8bfedde9, 0aebce52, 4ba02357, e36266be, 2fe472d6, 97efe937, 77695b08 … | Rutter | no (Python) |
| pyddt preflight checks, run monitor, agent skill file, region resume, flexible `setup_working_directory` (swd) input matching, config CLI + pytest, param colormaps/tests | 0f4ae3b6, 0fae56b0, 9cea2022, 7eb1eab9, d0de277e, 8e2fc5b0, fc06fc1d, b6ac8173, 4b6f365a, c38ea5f1 … | Carman (some Rutter) | no (Python) |
| **Output names and descriptions** in `config/output_spec.csv` + Runner.cpp keys | ee195125, 985a83e0, **fead82c2** (remove underscores), **888ad6ea** (NRESORB→NRESORPTION), f255c658, 186796e7, 4c003214 (MOSSDEATHC→LFNVC refs), 1f8d82f5, a90821de | Rutter | **output names** |
| Ground layer-array accessor `get_variable_bylayer<T>` (+ geometry, z), `Ground.ipp` | 5ceac73c, 0e800b5f, 265e8042 | Rutter | header only; not called anywhere yet |
| `TEMUtilityFunctions.h` drops `#include <netcdf.h>` | 4837788f | Rutter | build only |
| **Per-cell time-limit fix** (`md.cell_stime = time(0)` in the Runner ctor) | **76e91ed0** | Carman | wall-clock guard only |
| Docker hpcdev container, image version check, docker-compose, build wrapper | 0ef33b47, b8fc5b1d … | Rutter | no |
| git-lfs filter removed from `.gitattributes`; testing-data README | 57abb7ae, 3cccc751, f1c793a2 | Carman | no |
| New demo data `demo-data/crujra-v2.5_cmip6_EML_4x3` (17 nc) | 7e6a9a9e | Carman | new inputs, demo only |

The full list is in VM `~/tem-dev/upmerge/commits.txt`.

**Changes to `parameters/`, physics, or model defaults: none.**
- No commit touches `parameters/`, `config/config.js`, `ModelData`, `Cohort`, `Soil_Bgc`, or any physics source file.
- The only C++ changes are the Runner output keys, the timer reset, a new Ground template, and one `#include`.

## 2. Conflict resolutions (4 files)
1. **`include/Ground.h`** (additive on both sides): kept our `get/set_ch4_bylayer` and upstream's `get_variable_bylayer<T>` declaration.
2. **`src/Runner.cpp`**:
   - Includes: the union (`<ctime>`, `<fstream>`, `<iostream>`).
   - Constructor: kept our dangling-pointer fix, `Cohort(y, x, &this->md)`. Upstream still passes `&mdldata`, the by-value parameter (the ASan stack-use-after-scope bug). Added upstream's `this->md.cell_stime = time(0)` before it.
3. **`config/output_spec.csv`**: row-keyed 3-way merge (script `/tmp/csvmerge.py`, copied in `~/tem-dev/upmerge/`):
   - Kept our row order and all 53 fork-only rows (CH4*, TK*, DS*, BNF*, EVAP_*, ALD_FLAG/FRONT …).
   - Kept our ALD/TRANSPIRATION edits.
   - Took upstream's description text for 96 rows and the units for the 4 Q* rows (mm/time).
   - **Upstream changed no Yearly/Monthly/Daily/PFT/Comp/Layer flag**, so the default output set is unchanged.
   - Applied the 13 renames.
   - Dropped `MOSSDEATHC`. It had no C++ writer on either side.
4. **`pyddt/src/pyddt/util/setup_working_directory.py`**:
   - Took upstream's `pick_input_file()` discovery.
   - Added our CH4 inputs (`ch4_file`, `proj_ch4_file`) through it, with local default names. They are deliberately not in `DEFAULT_INPUT_FILE_NAMES`, because pyddt `input.py` and `preflight` treat that dict as the *required* inputs, and most input sets have no CH4.
   - Added the CH4 keys to the `--copy-inputs` relpath loop.
   - Kept our explicit IO `sort_order`. Upstream's dict-derived list lacks the ch4 and Fireweed keys, so `sort_order.index()` would raise.

The merge also auto-merged `pyddt/.../calibration_targets.py`, `TEMUtilityFunctions.h` and `docker-compose.yml` without conflicts.

### Fork-side additions in the merge
- **Legacy output-name aliases (Runner.cpp, 13 outputs).** Each renamed output first looks up the new name (NRESORPTION, SOC0TO30CM, TSOIL30CM, VWC30CM, …), then the legacy name (NRESORB, SOC0_30cm, TSOIL_30cm, VWC_30cm, …).
  - Why: the output file is created for every spec row in `ModelData::create_netCDF_output_files` but is only filled by the matching `find()` in Runner.
  - Without the alias, the production S22/S21H specs (which request `NRESORB` y/m) would silently get a **fill-only NRESORB file**.
  - With the alias, `NRESORB_monthly_sc.nc` is written and is bit-identical to 66cfec02 (verified in U1/U2).
  - Caveat: a spec that lists both names gets data only under the new name.
- **`13b2f1a9` swd fixes:**
  - (a) Pre-existing fork bug, there since the RTS port 1d7d5fe2: `sort_order` lacked `hist/proj_exp_thermokarst_file`, so swd raised `ValueError` for every config. I added the keys, and unknown keys now sort to the end.
  - (b) Upstream's fuzzy matcher takes the alphabetically first match, so `original-vegetation.nc` beat `vegetation.nc` in the Toolik demo set. It now prefers the exact default filename when present, which restores pre-merge behaviour for standard-named sets.
  - Tested on Toolik 10x10 (plain and `--copy-inputs`) and EML 4x3. pyddt `test_config.py` + `test_util.py`: 17 passed.

Our fixes survive unchanged: no conflict touched them, and identity holds. These are snow thermal closure, Richards fixes A/B, restart-state carry, the deep_carbon binary default 0, moss, run-on, RTS, and the wetland/enthalpy switches.

## 3. Upstream behaviour/default changes and how each was handled
| change | effect | handling |
|---|---|---|
| Output renames (888ad6ea, fead82c2) | The new names are written. Without a fix, old-name specs would get empty files. | **Compat alias**: both names work, so no switch is needed. Our tooling and docs that use old names keep working. Optional later clean-up: migrate the production specs to the new names. |
| Output descriptions/units text (96 rows) | NetCDF `long_name`/units attributes change for files made from the repo template. The production specs carry their own text. | Accepted. Metadata only. |
| Cell time-limit fix (76e91ed0) | Only when `model_settings.cell_timelimit > 0`: the limit is now per cell instead of from program start. Before, later cells in a long regional run could be aborted early. | Accepted as a bug fix, no switch. Output-neutral whenever the limit isn't hit (all our runs). **Decision for Elchin:** keep it as an unconditional fix (recommended), like the restart fix. |
| `MossDeathC` calibration target now maps to output `LFNVC` (4c003214, calibration_targets.py) | Calibration tooling would compare the MossDeathC target to non-vascular litterfall C instead of the (unwritten) MOSSDEATHC. | Accepted as upstream intends. No model effect. **Check if** your calibration uses MossDeathC targets (cal-cmt06 scripts weren't changed here). |
| swd fuzzy input matching (d0de277e, 8e2fc5b0, fc06fc1d) | Can pick a different file in non-standard directories. | Exact default name now preferred (13b2f1a9). Fuzzy matching only as fallback. |
| `.gitattributes` LFS filter removed (57abb7ae) | testing-data is no longer LFS-filtered in new clones. | Accepted (repo hygiene). |
| `TEMUtilityFunctions.h` drops `<netcdf.h>` | Could break files that relied on the transitive include. | Our tree builds cleanly with `-Werror`. |

No upstream physics change was found, so no new default-off switch was needed.

## 4. Build and tests (calib-elchin, conda gxx 13.4, -std=c++20, -j4; VM idle)
- The first build failed because the new worktree had no `external/Fireweed`, `external/Burnup`. These are untracked, git-ignored clones, not submodules. I copied them from `~/tem-dev/restartfix/wt/external` (`diff -rq` identical) and the build then succeeded.
- `dvmdostem` md5 e9fe5eff… (copy `~/tem-dev/upmerge/bins/binM`). Baseline `bins/binA` = 66cfec02 binary md5 972b5705.
- `make unittest`: **482/482 cases**; assertions 3,173,162, of which 8 fail = the same 8 known ones as at 66cfec02.

## 5. Identity check vs 66cfec02 (all switches off; protocol of NOTES_identity_check.md V3/V4)
- Driver: `~/tem-dev/upmerge/runV.sh`. It runs `tests/regression/run_regression.sh -j 4` (≤4 cores) with A=binA, B=binM and R=repeat of B. Outputs go to the new dir `~/tem-dev/upmerge/out/`. Restart0 is read-only from the case sets; nothing is reused or overwritten.

| run | case set | cases | files / vars compared (A vs B) | differing | A==B | B==R | core-h |
|---|---|---|---|---|---|---|---|
| U1 | harness `prod` (`~/tem-dev/regress/cs/prod`), spin20 + main n75 | 8 | 1,088 / 4,136 | **0** | all True | all True | 0.203 |
| U2 | `cs_full` S22 TVC / S21H HPC production copies, spin20 + main n150 (1950–2099) | 8 | 1,168 / 4,304 | **0** | all True | all True | 0.478 |

- U2's harness "FAIL" comes only from its closure gates, as in V4:
  - the production spec does not write the N-closure outputs;
  - HPC_Warming_nonfeed has a −1.16 % P water residual, identical in A.
- **Differences:** none, so there are no magnitudes or source commits to report.
- **Not covered:** MPI build, regional `cell_timelimit` > 0 runs, and a spec using the *new* names (the template has them; that is a string lookup).

Plots:
- `docs/upmerge_identity.png`: variables compared per case and stage, all 0 diff.
- `docs/upmerge_cs_full_series.png`: TVC and HPC Warming_nfeed 1950–2099 ALD, SOC, and the legacy-named NRESORB, A vs B. max|B−A| = 0. NRESORB is filled, not fill-only.

## 6. Notes for the production merge
- This branch is now the most up-to-date candidate for production. Production behaviour is unchanged when the configs pin the legacy switches, as established in §8.9 and §8.9.1 of the merge plan.
- Upstream `RTS` and `cmt_change` are **not** included (14/15 conflicting files, parameter changes). Leave them to the RTS expert.
- Housekeeping (no deletions made): swd test dirs `~/tem-dev/upmerge/swdtest*` and the private venv `~/tem-dev/upmerge/venv` (commentjson + pytest) can be removed later.
