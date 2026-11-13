# Fire C/N balance fixes (classic fire): `fire_n_redep_fix`, `fire_c_bookkeep_fix`

Branch `wetland/reburn-switch`, on top of 5d86bad8 (Reburn merge behind `model_settings.fire`).
Requested by Elchin on 2026-11-12: "fix the fire N re-deposition and C gaps behind default-off switches".
Background: the C/N balance residuals found in the Reburn smoke runs (docs/WORK_NOTES_reburn.md).

## Switches (both default 0 = legacy, bit-identical to 5d86bad8)

| switch (`model_settings`) | default | what 1 does |
|---|---|---|
| `fire_n_redep_fix` | 0 | The fire-emitted N re-deposition really enters the soil. Each month `FirData.fire_a2soi.orgn / 12` (output `BURNAIR2SOILN`, an annual rate = emitted veg+soil N / FRI) is added to **organic N (ORGN) of the top soil layer** (`del_sois.orgn[0] += m_a2soi.orgninput`, nfeed only). It stops when Cohort zeroes `fire_a2soi.orgn` one FRI (`mthsdist >= FRI*12`) after the fire. |
| `fire_c_bookkeep_fix` | 0 | Standing-dead C/N bookkeeping (Vegetation_Bgc + WildFire): (i) the monthly fall `d2wdebrisc/n` is 0 when nothing falls; (ii) at the end of the 9-yr window (`mthsdist >= 108`) the remaining standing dead falls into woody debris (DWD) instead of being deleted; (iii) if `deadc0/deadn0` is unknown (they are not in the restart file) the remainder falls linearly over the rest of the window, `deadc/(108 - mthsdist)` per month; (iv) old standing dead burned by a later fire is added to the fire emissions (`BURNVEG2AIRC/N`, its unretained part used to vanish). The legacy fall rate `deadc0/108` per month is kept. |

Enable (config/config.js or a case config):

```json
"model_settings": { "fire": 1, "fire_n_redep_fix": 1, "fire_c_bookkeep_fix": 1 }
```

Both are integers 0/1; anything else is FATAL at start-up (like the other fork switches). They only do
something when fire is on (`fire = 1`) and a fire has happened (or a restart carries fire state).
To also close the fire-month C gap, use them together with the existing layer switches
`layer_n_remap_fix = 1, layer_bookkeep_fix = 1` (see "Results").

## The bugs (legacy, switches 0)

1. **N re-deposition never applied.** `Soil_Bgc::prepareIntegration` set `orgninput` only
   `if (cd->yrsdist < cd->fri)`. `CohortData::fri` is never assigned anywhere (not initialised either),
   so the N is reported in `BURNAIR2SOILN` but never added: in the classic smoke case about −11.5 g N/m²/yr
   of the N balance. Now `Soil_Bgc::setFireNRedepInput()`; switch 0 keeps the old gate exactly.
2. **Standing-dead C/N created after the fall window.** In `Vegetation_Bgc::prepareIntegration` the fall
   `bd->m_v2soi.d2wdebrisc` was only set inside the window; after it, it kept the last month's value
   (`deadc0/108`), and `Soil_Bgc` added it to DWD every month until the next fire: +8.5 g C/yr in the smoke case.
3. **Standing dead deleted.** At `mthsdist >= 108` the remaining `deadc/deadn` was set to 0 (lost). After a
   restart `deadc0 = 0` (not in the restart file), so nothing fell and the whole pool (72 g C in the smoke case)
   vanished 108 months into the main stage.
4. **Old standing dead burned but not emitted.** `burnVegetation()` burns all old `deadc/deadn` (`comb_deadc`),
   counts it in the retained share, but `fire_v2a.orgc = comb_vegc - reta_vegc` leaves it out of the emissions.
   (Zero in the FRI-15 smoke case, because the standing dead is gone after 9 years; it matters for FRI < 9 or
   re-burns.)

Not changed (found, see "Decisions"): the remaining fire-month C gap of classic fire (3.8 g in the first fire) is
the dynamic soil-layer bookkeeping that `layer_bookkeep_fix` already fixes; dead-root N added to `somcr` (a C pool)
in `burnVegetation()` (≤0.006 g); DWD N decomposed with no destination (`wdebrisn -= rhwdeb*N/C`, Soil_Bgc);
process (Reburn) fire burns soil C with `BURNSOIL2AIRN = 0`.

Files: `include/ModelData.h`, `src/ModelData.cpp` (parse/default/FATAL), `src/Cohort.cpp` (pass-through),
`include/Soil_Bgc.h`, `src/Soil_Bgc.cpp` (`setFireNRedepInput`), `include/Vegetation_Bgc.h`,
`src/Vegetation_Bgc.cpp`, `include/WildFire.h`, `src/WildFire.cpp`, `config/config.js` (keys at 0),
`tests/unit/test_fire_balance.cpp`. No input or `parameters/` change.

## Tests

### Bit identity with both switches 0 (binary of this commit vs b4618e85 = 5d86bad8)

`scripts/identdir.py` on every NetCDF output and restart file (all variables, data bytes; only
`total_runtime` and `config_record.js` ignored), spin-up (20 yr) and main (75 yr):

| set | stage dirs | result |
|---|---|---|
| 16 regression cases (HPC/TVC × warming/no-warming × nfeed/nonfeed, in two settings `leg` and `s4_f1`, 8 + 8), fire off | 32 | identical (65–93 files each, 0 differences) |
| classic explicit fire E1 (fire on) | 2 | identical (78 + 105 files) |
| classic FRI fire F1c (fire on, 4 main-stage fires) | 2 | identical (78 + 105) |
| process fire F1p, fire 0 with a fire request F0g, F1c with monthly DEADC/DWDC | 6 | identical |

### Unit tests (`tests/unit/test_fire_balance.cpp`, suite `fire_balance`, 5 cases, 45 assertions, all pass)

- switch defaults (ModelData(), ModelData(Json), Soil_Bgc, Vegetation_Bgc) and FATAL on 2, −1, 0.5, "1";
- fire-month N balance: with the fix the full `BURNAIR2SOILN/12` lands in top-layer ORGN and the
  soil N balance closes (dN − internal fluxes = input, 1e-12); legacy: not applied; 12 months = annual rate;
- standing-dead C over 168 months after a fire: fix conserves DEADC + fall to DWD exactly, nothing falls after
  the window; legacy keeps adding `d0/108` per month after it (C created);
- after a restart (deadc0 unknown): legacy deletes it at month 108, the fix lets it fall (`d/(108-m)`),
  remainder falls at once past the window;
- DWD side: `Soil_Bgc::deltastate` adds the fall to woody debris.

Full suite: 395/395 test cases pass (`build.sh ~/tem-dev/rb unittest`); the 8 failing assertions are the
known ones inside may-fail cases (climate, Richards, SEB, Stefan), unchanged.

### Balance gaps, TVC smoke cases (CMT51 inputs, 20-yr spin + 75-yr main, FRI 15, July, severity 3)

Annual residuals of the main stage. C: ΔC(VEGC+DEADC+DWDC+SOC) − (NPP − RHSOM − RHDWD − fire C).
N: ΔN(VEGNTOT+ORGN+AVLN+DEADN+DWDN) − (NINPUT − NLOST − fire N + re-deposition due), where the
re-deposition due in a month is the previous month's `BURNAIR2SOILN`/12 (Dsb runs after Bgc). With
`fire_n_redep_fix = 0` the N residual is therefore ≈ −re-deposition (N reported, not applied).

| case (classic fire) | C fire yrs 15/30/45/60 | C max other yrs | C sum 75 yr | N fire yrs | N max all yrs | N sum 75 yr |
|---|---|---|---|---|---|---|
| fixes off (5d86bad8 behaviour) | 3.80 / 6.13 / 2.39 / 1.69 | **72.0** | +23.2 | −1.10 / −0.16 / −1.11 / −1.33 | **15.3** | **−208.7** |
| `fire_n_redep_fix` only | 3.80 / 6.19 / 2.57 / 1.95 | 72.0 | +26.5 | −0.06 / 0.10 / 0.03 / 0.03 | 3.99 | −5.4 |
| `fire_c_bookkeep_fix` only | 3.80 / 1.16 / 0.73 / 0.54 | 0.005 | +6.2 | −1.24 / −0.36 / −1.06 / −1.33 | 11.7 | −209.0 |
| both fixes | 3.80 / 1.22 / 0.86 / 0.69 | **0.005** | +6.6 | −0.16 / −0.09 / −0.03 / −0.02 | **0.36** | **−7.7** |
| both fixes + `layer_n_remap_fix`/`layer_bookkeep_fix` 1 | **0.007 / 0.002 / 0.001 / 0.000** | 0.006 | — | −0.16 / −0.09 / −0.02 / −0.02 | 0.36 | −7.5 |
| layer fixes only (fire fixes 0) | 0.007 / 4.82 / 1.27 / 0.64 | 72.0 | — | −1.10 / −0.17 / −1.10 / −1.31 | 15.3 | −208.2 |

| case (process / Reburn fire) | C fire yrs | N fire yrs | N max | N sum |
|---|---|---|---|---|
| fixes off | 5.18 / 5.13 / 5.04 / 5.04 | −0.05 / −0.32 / −0.43 / −0.25 | 0.69 | −23.0 (= −re-dep. 22.9) |
| both fixes | 5.18 / 5.13 / 5.04 / 5.04 | −6.17 / −0.51 / −9.25 / −1.51 | 9.25 | −17.4 |
| both fixes + layer fixes | −0.001 / 0.009 / 0.000 / 0.009 | −6.17 / −0.41 / −9.01 / −1.27 | 9.01 | −16.9 |
| fire off (F0g) | ≤0.005 | — | 0.0004 | 0.0 |

Reading: with both fire fixes the classic C gap outside the fire month is gone (72 → 0.005 g) and the N gap
drops from −209 to −7.7 g N over 75 yr. The remaining classic fire-month C gap (≈ retained + dead-root C,
3.8 g in the first fire) closes only with the existing layer switches, so it is a soil-layer bookkeeping issue,
not standing dead. The remaining classic N residual (≈ −0.06 to −0.1 g/month in the growing season while DWD
decomposes; −7.5 g over 75 yr) is DWD N lost on decomposition (Soil_Bgc, not fire). For process fire, the N fix
moves the gap into fire years: process fire burns soil C (497 g in 1965) with `BURNSOIL2AIRN = 0`, and the
re-deposited N in the burned top layer disappears with it (ORGN ends the fire month at the same value as with the
fix off). That is the Reburn process soil-burn N accounting, not changed here.

Figures (`docs/figures/fire_balance/`, numbers in `fire_balance_summary.json`):

- `c_balance_gap.png`: annual C gap per year, fixes off vs on (classic and process).
- `n_balance_gap.png`: annual N gap per year, fixes off vs on.
- `pools_on_off.png`: VEGC, DEADC, DWDC, SOC, ORGN, AVLN, classic fire with fixes off vs on, fire off for reference.

![C balance gap](figures/fire_balance/c_balance_gap.png)
![N balance gap](figures/fire_balance/n_balance_gap.png)
![Pools on vs off](figures/fire_balance/pools_on_off.png)

Effect on pools (classic fire, 2024): ORGN +255 g N/m² (2582 → 2838; saw-tooth: +≈190 g N re-deposited over each
FRI, ≈200 g N emitted by the next fire), VEGC +3.8 g C (18.8 → 22.6), DWDC 14.3 → 8.0, SOC +3 g C; DEADC no longer
sits at 72 g for 9 years after the restart (it now falls into DWD).

Runs: `~/tem-dev/wet-work/runs/FB1_{leg,s4_f1}/*` and `~/tem-dev/wet-work/rbsmoke/runs/{E1,F1c,F1p,F0g,F1cdbg}_fb0`,
`F1c_{on,N,C,L,onL}`, `F1p_{on,L,onL}`, `F1cdbg_on`; binary `wet-work/bin/dvmdostem_fb2`; analysis `/tmp/fbana.py`
on the VM (copy in `tests/regression/reburn/analyze_fire_balance.py`). About 0.7 core-h in total.

## Decisions for Elchin

1. **Restart and the re-deposition window** (fixed, see "Fire state in the restart file" below). `mthsdist` (and `deadc0/deadn0`) are not in the restart file;
   `mthsdist` restarts at 0, so with `fire_n_redep_fix = 1` a fire late in the spin-up deposits for up to one more
   full FRI in the main stage (here 5 yr in spin + 15 yr in main of `emitted/FRI`, ≈ 1/3 too much). The balance
   still closes (the input is what `BURNAIR2SOILN` reports), but the "return all emitted N in one FRI" total does
   not. Fix options: write `mthsdist`, `deadc0`, `deadn0` to the restart (restart format change) or use `yrsdist`.
2. **How much N to return.** The TEM design returns 100 % of the fire-emitted N over one FRI (closed-N-cycle
   assumption, comment in WildFire.cpp). With the fix that is +≈13.6 g N/m²/yr in this case and ORGN +255 g by 2024.
   You may want a fraction (a parameter) instead; not added (no parameter changes in this task).
3. **Use with the layer fixes.** Classic fire C closes fully only with `layer_n_remap_fix = 1, layer_bookkeep_fix = 1`
   as well; the same layer switches also close the +5 g/fire C gap of process fire. Consider recommending the four
   together for fire runs.
4. **Not fixed, separate items if you want them:** DWD N decomposed with no destination (≈ −7.5 g N / 75 yr here)
   and dead-root N added to `somcr` instead of `orgn` in `burnVegetation()` (≈0.006 g per fire) — both fixed behind
   switches, see "N leak fixes" below; `BURNVEG2DEADC/N`
   output also counts PFTs with no cover (115 vs 77 g in 1965, output only); process fire soil N not burned/reported
   (`BURNSOIL2AIRN = 0`), which with the N fix shows up as −6 to −9 g N in process-fire years;
   `CohortData::fri` is read uninitialised in the legacy path (kept for bit identity).

## Fire state in the restart file (2026-11-12, follow-up)

Elchin: "Add the fire state to the restart file behind the switch, test it, then push". This fixes the
restart issue from decision 1 above.

**What is written.** When `fire_c_bookkeep_fix = 1` or `fire_n_redep_fix = 1`, the restart file (`restart-*.nc`)
gets three new variables. With both switches 0 nothing new is written or read, so the restart files are unchanged.

| variable | dims | type | content |
|---|---|---|---|
| `fire_mthsdist` | (Y, X) | int, `_FillValue` −9999 | months since the last fire (`CohortData::mthsdist`): the 108-month standing-dead fall window and the FRI×12-month re-deposition window |
| `fire_deadc0` | (Y, X, pft) | double, `_FillValue` −9999 | standing-dead C right after the last fire (sets the fall rate `deadc0/108`) |
| `fire_deadn0` | (Y, X, pft) | double, `_FillValue` −9999 | standing-dead N right after the last fire |

The rest of the fire state was already in the restart: `firea2sorgn` (re-deposition rate), `deadc`/`deadn`,
`yrsdist`. `d2wdebrisc/n` are recomputed every month with `fire_c_bookkeep_fix = 1`. `comb_deadc` lives only inside one
burn, so neither needs saving.
Code: `RestartData::{read,write}_px_fire_vars` (optional variables, added to an existing file the same way as
`d2m_ald_prevyr`; serial builds only), `Cohort::set_fire_state_from_restartdata()` (called last in
`set_state_from_restartdata`) and the writes in `Cohort::set_restartdata_from_state`.

**Old restart files** (no `fire_*` variables) still load. With a switch on, the model prints once
`WARNING: restart file has no fire state (...). Months since fire and the standing-dead fall start as before (legacy
behaviour).` and keeps the previous behaviour (the month counter starts at 0 and the fall uses `deadc/(108-m)`). The warning
showed up in the spin-up of the test case, which starts from the case's old `restart0.nc`.

**Tests**

- Bit identity, both switches 0, this commit vs 19af750f (binary `dvmdostem_fb2`), all NetCDF outputs and restart
  files: 16 regression cases (8 `leg` + 8 `s4_f1`, spin-up + main = 32 stage folders) and the fire cases E1
  (classic explicit), F1c (classic FRI), F1p (process), F0g, F1cdbg (10 stage folders): **42/42 identical,
  0 differences** (restart files have no new variables).
- Restart continuity, both fixes on. Case: TVC explicit classic fire 1965/1980/1995/2010. The main stage was run in
  one piece (75 yr) and again split at 1985, 5 yr after the 1980 fire (inside both the standing-dead window and the
  re-deposition window). Part A covered 1950–1984. Part B covered 1985–2024 from A's `restart-sc.nc`, with the projected
  climate and the explicit-fire file sliced by 35 yr (CO2 follows the calendar year of the climate file). Part B was
  compared with 1985–2024 of the continuous run (`compare_split.py`):

  | run, split vs continuous | DEADC max abs diff | DWDC | DEADN | BURNAIR2SOILN | DEADC in the first year after the split |
  |---|---|---|---|---|---|
  | new (fire state in restart) | **0.00057 g C** (3e-5 rel.) | 0.0018 | 4e-5 | 2e-8 | **0.0** |
  | old, 19af750f, fixes on | 5.78 g C (34 %) | 4.81 | 0.33 | 3.5e-5 | 1.39 |
  | fixes off (legacy) | 13.0 g C (77 %) | 20.4 | 0.67 | 5.5e-4 | 2.83 |

  The remaining differences are not fire-related. They appear the same in all three split runs (for example soil
  temperature up to 0.75 K in the first month, SOC 0.013 g, NPP 0.0013 g, ORGN 0.0019 g), because the dvmdostem
  restart does not restore the full thermal and snow state. That is why the split is not bit-identical: 5 of 101 output variables are identical
  in each split run, and 30–31 of 64–67 restart variables at the end. With the new restart, the fire pools follow
  the continuous run to that same noise level. With the old one they jump: the standing dead restarts a new 9-yr fall.
- Spin-up → main overshoot (classic FRI 15, spin-up fire in July of spin-up year 15, both fixes on), N re-deposited
  from the spin-up fire until the next fire (`redep_overshoot.py`):

  | binary | emitted N | re-deposited N | ratio | months depositing | all re-deposition in the main stage |
  |---|---|---|---|---|---|
  | 19af750f | 170.33 | 221.43 | **1.30** | 234 | 945.7 g N |
  | new | 170.33 | 170.33 | **1.000** | 180 (= FRI × 12) | 734.9 g N |

- Unit tests (suite `fire_balance`, now 8 cases, 154 assertions, all pass), with three new cases: the fire state is
  off by default, `reinitValue()` clears it, and reading a file without the variables leaves them MISSING and adds
  nothing; a round trip of `fire_mthsdist`, `fire_deadc0`, `fire_deadn0` for two pixels plus a never-written pixel
  (fill values); and the standing-dead fall resuming after a restart month for month exactly with the restored state
  (a legacy restart starts a new schedule). Full suite: 398/398 test cases pass (the 8 known may-fail assertions are
  unchanged).

![restart continuity](figures/fire_balance/restart_continuity.png)

`restart_continuity.png`, left: standing-dead C after the 1980 fire, continuous vs restarted in 1985, new vs old
restart. Right: cumulative N re-deposited from the spin-up fire, old (1.30 × emitted) vs new (1.00 × emitted). Numbers
are in `split_summary.json` and `overshoot.json`. Scripts are in `tests/regression/reburn/` (`split_run.py`,
`compare_split.py`, `redep_overshoot.py`, `plot_restart_continuity.py`). Runs (on the VM):
`rbsmoke/runs/E1on_fb{2,3}{,_split}`, `E1_fb3{,_split}`, `F1c_on_fb3`, `runs/FB3_*`. About 0.55 core-h.

Decision for Elchin: with explicit fire the re-deposition rate is `emitted / fri` with `fri` from `fri-fire.nc`. In
this test case `fri = 2000` yr, so only about 1 % of the emitted N comes back in 20 yr. That is unchanged here.

## N leak fixes: `dwd_n_leak_fix`, `deadroot_n_route_fix` (2026-11-12, follow-up; local, not pushed)

Two more `model_settings` switches, both **default 0 = legacy, bit-identical to 06c41688** (integers 0/1 only,
anything else FATAL, parsed in `ModelData`, defaults in both constructors). Pure mass-bookkeeping fixes, no
published source.

| switch | legacy (0) | fix (1) | code |
|---|---|---|---|
| `dwd_n_leak_fix` | woody-debris decomposition removes `RHWDEB × WDEBRISN/WDEBRISC` from `WDEBRISN` with no destination: the N is lost every month DWD decomposes | the same N goes to the organic N (`ORGN`) of the top soil layer (where litter N and the fire re-deposition go; it is then mineralised by the normal SOM path) | `Soil_Bgc::deltastate()` (wood debris block) |
| `deadroot_n_route_fix` | `WildFire::burnVegetation()` adds the N of roots killed (not burned) by fire to `m_sois.somcr`, a **C** pool: the N disappears and the same mass of C appears | the N goes to `m_sois.orgn` by root fraction (the C still goes to `somcr`) | `WildFire::deadRootNPool()` used in `burnVegetation()` |

### Tests

- **Bit identity, new switches 0** (`dvmdostem_fb4` = this commit vs `dvmdostem_fb3` = 06c41688, `scripts/identdir.py`, all
  NetCDF outputs and restart files): **46/46 stage folders identical** — 16 regression cases (`FB4_{leg,s4_f1}` vs
  `FB3_*`, spin + main = 32), the 5 classic fire cases `E1, F1c, F1p, F0g, F1cdbg` (10), and `F1c_on`, `E1on`
  (fire fixes 1, leak fixes 0; 4).
- **Unit tests** (suite `fire_balance`, now 11 cases, 190 assertions, all pass), three new cases: the switches are 0
  by default and 0/1 only; one `Soil_Bgc::deltastate()` month with decomposing DWD, where the legacy soil + DWD N
  balance is short by exactly `RHWDEB × N/C` and closes to round-off with the fix (the N lands in layer-0 ORGN,
  other layers and AVLN unchanged, no DWD → no change); and the dead-root pool (`somcr` legacy, `orgn` fix) with the
  N/C effect of distributing a dead-root N amount by root fraction. Full suite: 401/401 test cases pass (the 8
  known may-fail assertions are unchanged).
- **Balance gaps with the switches on** (TVC smoke cases, 20-yr spin + 75-yr main, annual
  `ΔN − (NINPUT − NLOST − fire N + re-deposition)`, all other fixes as named; `analyze_n_leaks.py`,
  `n_leak_summary.json`):

| case (all with both fire fixes) | N gap sum before | after | max \|annual N gap\| before | after | fire-year N gaps before → after |
|---|---|---|---|---|---|
| classic FRI 15 + layer fixes (`F1c_onL` → `F1c_onLK`) | **−7.665 g N** | **0.000** | 0.420 | 0.0002 | −0.14/−0.09/−0.02/−0.01 → 0.000 |
| explicit classic E1 (`E1on` → `E1onL`) | **−7.798** | **0.000** | 0.420 | 0.0002 | −0.14/−0.09/−0.03/−0.02 → 0.000 |
| process (Reburn) + layer fixes (`F1p_onL` → `F1p_onLK`) | −15.43 | −15.41 | 9.01 | 9.01 | unchanged (−4.7/−0.4/−9.0/−1.3) |

  One switch at a time (classic `F1c_onL`): `dwd_n_leak_fix` alone closes 7.658 of the 7.665 g (−0.007 left, all
  in fire years); `deadroot_n_route_fix` alone closes 0.007 g N (0.0056 g in the first fire) and removes the same
  spurious C from the fire-year C gap (0.0069 → 0.0011 g C in 1966). Effect on pools by 2024 (classic): total N
  +8.2 g N/m² (2791.4 → 2799.5, in ORGN), SOC +0.4 g C, VEGC +0.1 g C, DWDN unchanged within 0.002.

![N leak gaps](figures/fire_balance/n_leak_gaps.png)

`n_leak_gaps.png`: annual N (left) and C (right) gap per year, leak fixes off vs on, classic, process and explicit
fire. Runs (on the VM): `rbsmoke/runs/{F1c_onL,F1c_onLK,F1c_onLD,F1c_onLR,F1p_on,F1p_onL,F1p_onLK,E1on,E1onL,F1c_on}_fb4`,
`runs/FB4_*`; cases `rbsmoke/cs/{F1c_onLK,F1c_onLD,F1c_onLR,F1p_onLK,E1onL}`. About 0.7 core-h (31 runs ≈ 0.5, builds and
unit tests ≈ 0.2).

Decisions for Elchin:

1. Destination of the DWD N: top-layer ORGN (chosen: same place as litter N and the fire re-deposition, robust with
   every `n_rootless_fix` setting). The alternative is top-layer AVLN (treat it as mineralised at once); it would
   make the N available to plants and leaching a little sooner.
2. The process-fire (Reburn) N gap of −15.4 g N over 75 yr, all in fire years, is not one of these two leaks: process
   fire burns soil C with `BURNSOIL2AIRN = 0`, so the N of the burned soil is not booked (item 4 of the first decision
   list). Fixing it is a separate change in the Reburn soil-burn accounting.
3. The E1 C gap in fire years (+3.8 g C in the first fire) remains because `E1on` runs without
   `layer_n_remap_fix`/`layer_bookkeep_fix` (see "Use with the layer fixes" above), not because of these switches.

## Reproducing the tests: scripts and tools in the repo (2026-11-13, follow-up; local, not pushed)

The test scripts and tools are now all in the repo; the case data stays on calib-elchin for now
(`~/tem-dev/srd-work/cs/leg`, the 8 calibrated regression cases, 144 MB; everything else is generated from it).

- `tests/regression/`: `identdir.py` (two output directories), `ident.py` (pairs of run directories, spin000 +
  main, JSON report) and `launch4.sh` (job list in parallel, `run_case.py` of this repo by default) — previously
  only in `~/tem-dev/wet-work/scripts`.
- `tests/regression/reburn/variants.json` + `mk_variants.py`: the `model_settings` overlay of every variant
  (`s4_f1` from `leg`; `F1c_*`, `F1cdbg*`, `F1p_*`, `E1on*` from the smoke bases), so no case is made by hand.
- All scripts take `--cases/--runs/--regcases/--figdir/--work` (or `TEM_WORK`, default `~/tem-dev`) instead of
  hard-coded `~/tem-dev` and `/tmp` paths; `split_run.py` imports this repo's `tests/regression/run_case.py`.
- `run_all.sh`: builds `19af750f`/`06c41688`/`f94fceaf` in `git worktree`s (with `scripts/get_external.sh`),
  makes the cases, runs the identity and balance sets and compares with the committed summary JSONs
  (`compare_summaries.py`). Expected: 88/88 stage folders identical and PASS for every summary; about
  1.5 core-h. Details in `tests/regression/reburn/README.md`.
- New: `docs/figures/fire_balance/identity_summary.json` (`ident.py` on the existing runs: fb2 vs fb3 42/42,
  fb3 vs fb4 46/46, 7266 files) and `test_summary.png` (`plot_test_summary.py`).

Validation (≈ 0.11 core-h, all in a new scratch directory `~/tem-dev/wet-work/rbval`; existing case and run
directories only read):

| check | result |
|---|---|
| `mk_variants.py --verify-against`: rebuilt s4_f1 (8 cases) and all 16 smoke variants vs the hand-made ones | **all identical**: parsed and byte-level config.js (paths normalised) and every other file by md5 (39–41 files per case) |
| `mk_smoke_cases.py` + `mk_explicit_case.py` rebuilt F0g/F1c/F1p/E1 | **identical** (configs, `fri-fire.nc` and explicit-fire NetCDF edits, `cmt_firepar.txt`, output specs) |
| analysis scripts with the new options on the existing runs (outputs to scratch) | `compare_summaries.py` **PASS** for all 5 summary JSONs; regenerated `restart_continuity.png`, `n_leak_gaps.png`, `c_balance_gap.png` byte-identical to the committed ones |
| `run_all.sh --dry-run` (all sets) | plan of 3 builds, 67 identity + 23 balance jobs (16 new, the rest reused), 3 splits |
| `run_all.sh --bins ... --sets cases,identity,check` on leg/TVC_Warming_nfeed + F0g (6 runs) | **8/8** stage folders identical (fb2 vs fb3, fb3 vs fb4); the fresh runs from the regenerated cases are also identical to the existing runs (6/6, including F0g from the pre-fix Reburn binary vs fb2) |

The build step (worktrees, `get_external.sh`, `make`) was checked in dry-run only, so as not to add worktrees to
the shared repository here.

![test summary](figures/fire_balance/test_summary.png)

`test_summary.png`: bit identity per comparison; N balance gap sums with the fixes off and on; spin-up → main
re-deposition ratio; restart continuity (max |restarted − continuous|) — all from the committed JSONs.
