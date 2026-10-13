# Mean-preserving monthly-to-daily forcing (branch `rts-light-m2d`)

## Problem
`Climate::monthly2daily` interpolates linearly between month midpoints. That is
not mean-preserving: near the convex July peak the daily curve averages below
the monthly mean. On Gydan ERA5 (1941-2023) it loses about 10.7 % of thaw
degree-days (835 vs 933 degC-d), which is about 5.5 % of Stefan ALT. Other
problems in the old code:
- `interpolation_range` used `data.at(curr_jan+13)` as "next Jan", which is
  actually next February.
- It always returned 366 days, with no leap-year handling.
- nirr and girr overshot below zero and were clipped, which biased their means.

## Method: Killworth (1996), chosen over Epstein (1991)
- Daily values are piecewise linear between month-midpoint pseudo-values `p_k`.
  The `p_k` are chosen so that the mean of the discrete, day-centred daily
  values in each month equals the monthly mean exactly.
- The constraint is tridiagonal (weights about 1/8, 3/4, 1/8). We solve it with
  the Thomas algorithm, O(n), over a 36-month window: previous year, target
  year, next year. Window-edge errors decay like 0.17^k, about 1e-9 after 12
  months, so consecutive years join continuously.
- The calendar weights are built once and cached in a static object. Each call
  then costs about 36x3 flops plus 365 evaluations, the same order as the old
  `temutil::resample` path.
- Why not Epstein: the harmonic method needs a 12-term Fourier fit plus a
  sin(x)/x correction for every day (more trig per call). It is periodic within
  a year, so it gives jumps at year boundaries in transient runs. It also still
  needs a non-negativity fix.
- Non-negativity (vapo, nirr, par, girr): if a month has negative daily values,
  they are set to 0 and the rest of that month is rescaled so the mean is
  unchanged. This applies around polar-night zeros. tair is not clipped.
- girr (12-value potential radiation) uses the same method, so the nirr/girr
  cloudiness ratio stays consistent. Precipitation is unchanged
  (`calculate_daily_prec`).

## Calendar / leap years
The model is noleap throughout: `DINY = 365`, `temutil::day_of_year`, DAYINM,
and the output time axes all assume it. The old 366th element was never read.
The new path returns exactly 365 values. Real leap years are not represented.
That would need calendar changes across Runner, EnvData and the outputs, so it
is out of scope.

## Switch
`model_settings.monthly2daily_mean_preserving`:
- `true` (default on this branch): the new path.
- `false`: the legacy code, bit-identical to the old output. This deliberately
  keeps the +13 quirk and the 366-day container. The flag is read into
  `ModelData::m2d_mean_preserving` and copied to `Climate` in the `Cohort`
  constructor.

## Known edge cases
- At the historic-to-projected switch, each file's series is padded with its
  own first or last year. The join is therefore only as smooth as the two
  files' Dec and Jan values allow (the legacy code behaves the same way).
- Equilibrium and pre-run stages use the periodic 12-month climatology.

## Tests (`make rts-test`, CI `.github/workflows/rts-tests.yml`)
- C++ unit tests (`tests/rts/test_rts_units.cpp`):
  - daily means equal monthly means within 1e-6 (relative) for tair, vapo,
    nirr and par, including exact-zero months;
  - non-negativity;
  - Dec 31 to Jan 1 continuity across years and in the climatology;
  - constant input gives constant output;
  - TDD is closer to the monthly-step TDD than the legacy result;
  - the legacy routine is bit-identical to a verbatim copy of the old code.
- Regression references: `tests/rts/reference/{baseline,rts,rts_ice}` are on
  the new path. `tests/rts/reference/legacy_m2d/` holds the old references.
  The flag-off runs must match them at zero tolerance.
- The synthetic ice profile top moved from 0.3 m to 0.4 m. The new ALD of
  about 0.25 m otherwise reached the ice without any RTS.

## Results (default params, n-factor not calibrated; RTS off; 1986-2020 mean ALT)
| Case | ALT old (m) | ALT new (m) | change | model TDD 1986-2015 old/new (degC-d) |
|---|---|---|---|---|
| Gydan CMT05 | 0.247 | 0.272 | +10.2 % | 823 / 903 |
| Gydan CMT70 | 0.193 | 0.212 | +9.7 % | 823 / 903 |
| Toolik synthetic | 0.216 | 0.239 | +10.6 % | 810 / 904 |

- Gydan CMT05, RTS on: total-C change at the end of the run is -3.56 kgC/m2
  with the new path, against -3.51 with the old.
- Runtime: the interpolation itself is 5.6 vs 13.8 us per call, so it is
  faster. Whole-run wall time changes because of the physics:
  - CMT05 stage 1: 75 to 79 s;
  - Toolik stage 1: 75 to 79 s;
  - CMT70 stage 1: 122 to 192 s (+57 %, reproduced in a serial rerun). The
    likely cause is more hourly freeze/thaw substeps in the solver. This has
    not been profiled.
- Scripts: /workspace/russia-inputs/m2d_compare.py; results in
  m2d_runs/results.json.
