# WORK NOTES: `model_settings.eq_nunlim_years` (N-unlimited early equilibrium)

2026-11-20. Branch `fix/richards-mass-guard` (local only, **not pushed**), on top of 0c40fd48 (Richards switches incl. fix B).
Elchin (verbatim): "Approve the eq_nunlim_years switch, local commit with tests".
Plan: `~/cal-cmt06/PLAN_spinup_fixB_recal.md` section 1. Diagnosis: `~/cal-cmt06/NOTES_vegflip.md`.

## Why

CMT06 Toolik spin-up is bistable (NOTES_vegflip): during eq years 2-10 plant N either keeps up with C (alive, GPP/INGPP ~0.34)
or locks up (dead, ~0.05). With no N inputs (NINPUT = BNF = NLOST = 0) the low-N state is absorbing, so the outcome of a
calibration case depended on the solver path (fix B flipped X3 from alive to dead) rather than on the parameters.
Running the first K eq years N-unlimited lets C pools (leaf, litter, SOM) build up before the N cycle switches on.

## What

- `model_settings.eq_nunlim_years` (int, **default 0 = unchanged legacy path, bit-identical**), parsed with
  `bnf::valid_int(v, 0, 100000)`; anything else (negative, > 100000, non-integer, bool, string) is FATAL.
- `K > 0` (eq stage only): eq years `[0, K)` run with `nfeed` off, then `nfeed` is restored to `stage_settings.eq.nfeed`
  and eq years `[K, eq_years)` run as before (`TEM.cpp`, eq branch; split computed by `eq_nunlim_segments()` in `ModelData.h`).
- Refused (FATAL) unless `0 < K < eq years`, fire is off in eq (dsb), and `K % output_interval == 0` (keeps the netCDF year
  index of the two `run_years` segments on the output grid).
- pr, sp, tr and sc are untouched. No physics or parameter change; no change to restart files' layout.

## Tests

- `tests/unit/test_eq_nunlim.cpp`: parse (default 0; accepts 0, 1, 100, 200, 100000; FATAL on -1, 100001, 100.5, true, "100")
  and segment split (K = 0 single legacy segment; K in {1, 100, 200, 999} two segments; invalid K / fire on / output_interval
  0 / K not a multiple of output_interval give no split).
- Full suite: 471 test cases, 471 passed (the 8 failed assertions are the known pre-existing expected-failure assertions).
- Full-model identity (kit v4 P0, short stage 1 `-p 100 -e 30 -s 0`): see section "P0 results".

## Not changed / noted and deferred (Elchin: "Note the open model questions and proceed")

1. No N inputs at CMT06 Toolik (NINPUT = BNF = NLOST = 0): the N cycle is closed and the low-N state absorbing. Deferred.
2. The soil N safety net fires ~4 months every year (~4000 events per 1000-yr eq) in alive and dead runs alike. Deferred.
3. NETNMIN is exactly 0.000 in the alive eq while plant uptake is ~0.52 g/m2/yr. Deferred (would need its own approval).

These are model questions, not fixed by this switch; the switch only removes the dependence of the spin-up outcome on the
early-eq path.

## P0 results

(filled in after the kit v4 P0 pilot; see below)
