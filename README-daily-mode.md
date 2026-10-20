# dvm-dos-tem opt-in daily climate mode (SC stage)

Source: copy of /home/ejafarov/tvc_hpc_mip/dvm-dos-tem (master 419a4eb3, clean; identical src/include
to tvc_hpc_mip/dvm-dos-tem.zip and to the source of the monthly-run binary). Branch: daily-climate-mode.
Patch vs pristine: daily-mode.patch (git diff 419a4eb3 -- src include).

## Enable
config.js -> "model_settings": { ..., "climate_timestep": "daily" }   (absent/"monthly" = original behavior)
IO.proj_climate_file must be a daily file: 365 steps/yr (365_day calendar), precip in mm/day,
length a multiple of 365 and >= sc_yrs*365, no fill/NaN/negative precip/nirr/vapor_press. Fails loudly otherwise.

## Build (native, plain ssh)
./build-conda.sh        # activates conda env /home/ejafarov/miniconda3/envs/dvmdostem-build and runs make
(env: conda-forge gxx 13.4, libboost 1.83, libnetcdf 4.9.2, jsoncpp 1.9.8, openblas/lapacke 0.3.34; rpath baked in,
 so ./dvmdostem runs without activating the env)

## Run (SC only, from restart-tr.nc placed in output_dir)
./dvmdostem -p 0 -e 0 -s 0 -t 0 -n 150 -f config/config.js -l disabled --no-output-cleanup
(this tree has no --restart-run flag; --no-output-cleanup keeps output/restart-tr.nc, which SC loads
 when IO.restart_from is empty)
