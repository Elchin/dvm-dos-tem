# Stage B compute-cost benchmark (enthalpy solver vs legacy thermal step)

Nothing here is used by the model build.

- `bench_patch.py TREE` instruments a throw-away copy of the model: it times
  `soilenv.updateDailyGroundT` (the whole legacy thermal step) with `std::chrono::steady_clock`,
  counts Crank-Nicolson tridiagonal solves, and, with `DVM_THERMAL_BENCH=prefix`, dumps each day's
  start-of-day column (T, liq, ice, solid capacity, texture, organic fraction, unfrozen and frozen
  conductivity), the surface temperature and the legacy end-of-day T to `prefix.<pid>.csv`.
- `bench_enthalpy.cpp` replays those days through `include/EnthalpySolver.h` (NCZ with SFCC 0 and 1,
  Katzenelson with SFCC 0) using the legacy top conductance (cn = 100), zero bottom flux and dt = 1 d.
  It reports iterations, tridiagonal solves and µs per column-day (mean, p95, max) for all days,
  phase-change days and other days, next to the legacy numbers from the same days.
- `run_bench.sh TREE OUTDIR` builds -O0 and -O2 copies and runs the calib8 TVC and HPC NoWarming_nfeed
  cases with the regression protocol (`tests/regression/run_case.py`, 20-yr spin + 75-yr main),
  using at most 4 concurrent processes.
