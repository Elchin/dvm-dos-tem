# Reburn smoke / identity helpers (calib-elchin paths)

Used for docs/WORK_NOTES_reburn.md. The scripts copy the TVC_Warming_nfeed regression case into new directories, so the originals and repo `parameters/` are not touched. The paths are the calib-elchin ones (`~/tem-dev/...`); edit them for another machine.

- `mk_smoke_cases.py`: builds F0g / F1c / F1p. F0g is `fire 0` with the stage asking for fire; F1c is classic fire; F1p is process fire with a 21-line `cmt_firepar.txt`. All three use an FRI copy (15 yr, day 196, severity 3) and enable the BURN* outputs.
- `mk_explicit_case.py`: builds E1, the explicit-fire case for the classic-fire identity check (new binary with `fire 1` vs the pre-Reburn binary).
- `analyze_smoke.py`: finds fire events, computes the annual C/N balance residuals, and writes the figures in `docs/figures/reburn/`.

Each case is run with `tests/regression/run_case.py BINARY CASEDIR RUNDIR` (20-yr spin-up plus 75-yr main). Outputs are compared with `identdir.py` (every variable of every NetCDF file, `total_runtime` ignored).
