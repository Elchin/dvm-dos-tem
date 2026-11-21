# WORK NOTES: snow_thermal_closure (and snow_enth_trace)

Branch `fix/snow-thermal-closure` off 3d282843 (fix/richards-mass-guard tip), worktree `~/tem-dev/snowleak-wt` (calib-elchin). Local commits only, not pushed.
Elchin (verbatim): "Approve the snowmelt fix code change as a local commit with tests, and the replay runs at about 0.31 core-h".

## 1. The leak, localised
`snow_enth_trace` (logging only) on the CMT06 1945 d147 day (W0_X1c09_alt_AB, alive start, fix A+B, snowmelt_heat_cap on) printed the snow cells at four phases of the daily thermal step. Before the solve (after `enthApplyMeltDirect`):

| cell | T (C) | liq | ice | H of the layer (J/m2) | H the solver assigns at that T (J/m2) |
|---|---|---|---|---|---|
| 0 snow | 0 | 16.55 | 0 | 0 | 0 |
| 1 snow | -0.0059 | 0 | 49.64 | -16.566e6 | -9.781e6 |

The difference (6.785e6 J/m2) is the whole daily residual. Cell 1 is off its freezing curve: it sits at T = -0.0059 C, inside the free-water band (-eps, 0) with eps = 0.01 C, but holds no liquid. The solver's H(T) on that band assumes ~41% of the water is liquid. The solve then sees +6.785e6 J that the layer never had, and the scatter (`liq = liq(cell, T)`) turns about 20 kg of ice into water without Lf being charged; the pack melts out early.

Cause (file/line, fix/richards-mass-guard code): `Soil_Env::enthApplyMeltDirect`, src/Soil_Env_enthalpy.cpp, cold-content branch (about line 109): `if (E < need) { l->tem += E / C; q.mcc += E; E = 0.0; break; }`. It warms the cell and leaves, skipping the "back onto the freezing curve at the same H" step at the end of the loop body (about lines 121-125). The comment there says T < 0 with liq = 0 is on the curve, but the free-water curve spans [-eps, 0]. It happens on a melt day when the energy left after melting out the top cell is smaller than the cold content of the next cell and ends that cell's warming inside the band: a 2-layer pack (top melts out, second warms to -eps..0). It is not in the hand-off or the scatter, and not in the solver or layer merge, so the Form B "solver or layer merge" fallback was not needed. All 20 known spike runs have snowmelt_heat_cap = 1.

## 2. The change
- `model_settings.snow_thermal_closure` (bool, default false; needs thermal_solver 1 and snowmelt_heat_cap true, FATAL otherwise). When on: after the direct melt and before the solve, every snow cell that is off its freezing curve (same `onCurve` test as step 3) is put back on it at the same layer enthalpy (`equilibrate`, H-conserving, mass unchanged). A day with no off-curve cell runs the same code as before. New CSV columns `snleak`, `snleak_n` (only when on): the enthalpy the legacy path would have created (H_solver(T) - H_layer) and the number of cells touched.
- `model_settings.snow_enth_trace` (int 0/1/2, default 0; thermal_solver 1 only): 1 prints the snow cells (T, liq, ice, layer H, solver H) after the direct melt, after the scatter and after the hand-off to stderr on days with |res| > 1e3 J/m2; 2 prints every day. Logging only.
- Files: include/ModelData.h, src/ModelData.cpp (parse, FATALs), include/Soil_Env.h, src/Cohort.cpp (pass-through, CSV), src/Soil_Env_enthalpy.cpp, config/config.js (docs), tests/unit/test_snow_thermal_closure.cpp.
- No change to snowmelt_heat_cap, the solver, fix A/B, or anything in parameters/.

## 3. Tests (tests/unit/test_snow_thermal_closure.cpp, 5 new cases)
1. Defaults off; parsing; FATAL without thermal_solver 1 / without snowmelt_heat_cap; trace range 0-2.
2. Reproducer column (top cell melts out, remainder warms cold second cell into the band): off |res| > 1e5 J/m2 (the leak), on |res| <= 1e-6 W/m2, `snleak` equals the legacy residual (rel 1e-6), mass closes to 1e-14, budget identity holds, pack keeps more ice.
3. Bit-identity: 3 pack shapes x 4 melt energies, switch on == off bitwise (res, src, gtop, mmelt, every layer T/liq/ice) on days with snleak_n == 0; with snleak_n > 0 the day closes.
4. Trace 0/1/2 give identical states; trace 2 prints.
5. Fuzz: 200 random 1-4 layer packs, random melt energy and air T: |res| gate holds and mass closes with the switch on.
Full suite: 476 cases, 476 passed (471 before + 5 new), assertions 3,172,994 passed and the 8 known failing assertions, status SUCCESS (build_snowleak_unittest.log).

## 4. Replays (calib-elchin, single-cell, 1 core, binary dvmdostem_snowleak_v1)
Reference runs are the j4 outputs from binary 43c28e8f. "Identical" = every variable of every output .nc equal (only run_status total_runtime differs) and every energy_daily.csv column equal.

| Case | Residual, switch off (new binary) | Residual, switch on | Non-leak days / runs |
|---|---|---|---|
| W0_X1c09_alt_AB tr+sc (1945 d147) | 78.53 W/m2 (6.785e6 J), identical to old run | 1.9e-7 J (max over run 2.7e-8 W/m2); snleak 6.785e6 J, 1 day, 1 cell | off == old, bit-identical; on == old up to row 16207 (first differing row = the leak day) |
| N0_X3_fl03 eq stage 1 (first 145 y) | old run: 78.8 and 19.5 W/m2 on rows 50144, 51239 | no spike day; max 1.1e-8 W/m2; snleak 6.81e6 J on 1 day | on == old up to row 50144 |
| N0_X1c09_fl03 eq stage 1 (first 185 y) | old run: 181.3 W/m2 row 66934 | no spike day; max 1.1e-8; snleak 1.57e7 J | on == old up to row 66934 |
| N0_X3_fl10 eq stage 1 (first 355 y) | old run: 58.4 W/m2 row 129349 | no spike day; max 1.1e-8; snleak 5.04e6 J | on == old up to row 129349 |
| R0_X3_fl03 own (tr+sc, R0 stage-1 restart) | no spike; off identical to R0 | identical to R0 and to off, snleak 0 | all days |
| R0_X3_fl03_alive, R0_X3_fl03_dead (tr+sc, frozen refs) | no spike; off identical | identical to old and off | all days |

So for the CMT06 X3 fl0.03 calibration the switch is a no-op: switch off is bit-identical to R0, and switch on gives the same output because none of its tr+sc runs has an off-curve snow cell. After a leak day the on-run trajectory differs from the old one (it keeps ~20 mm of snow a day longer), as expected. Not replayed: the j3 legacy-kit cases (S1_sd250 tr, S3_X1_sompr3 sc, S3_X3_micb0.2 sc) and P0_X3_K200 (stage-1 row 417693 would need a 1,140-year spin-up); the mechanism is the same in all 20 spike runs, which share the signature (melt day, nsnow 2, cap on).

Figure: `snow_thermal_closure_replays.png` (residual on the leak day off vs on; stage-1 max |res|; snleak vs old residual).

## 5. Cost
Replay stage: 0.200 core-h (cap 0.31): trace runs 0.03, W0 on/off 0.033, R0 identity set 0.093, stage-1 replays 0.042. Build and unit tests are not on the ledger. CMT06 ledger 3.1046 -> 3.3042 of 4.5 (stop 4.0).

## 6. Status and risks
- Default off; production use of the switch changes eq trajectories after the first leak day (about 0.7 per 1000 eq years), so it belongs in a re-spin, not as a silent change. It is a physics change (snow keeps its ice on those days) and needs Elchin's approval before production use.
- The fix is a re-equilibration of one cell at constant enthalpy; a more structural cure (re-equilibrating inside `enthApplyMeltDirect`'s cold-content branch for every user, switch-independent) would change legacy results and was not done.
- Only snowmelt_heat_cap = 1 is covered; with the legacy S[0] source no off-curve state was seen.
- Not pushed.
