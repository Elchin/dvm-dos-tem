{
 "general": {
  "output_global_attributes": {
   "run_name": "Demo dvmdostem run",
   "description": "A basic test run. Data is not scientifically sound."
  }
 },
 "IO": {
  "parameter_dir": "/home/ejafarov/tem-dev/wet-work/cmt31_nup/runs/N2_B1_mb015/params/",
  "hist_climate_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/historic-climate.nc",
  "proj_climate_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/projected-climate_ssp1_2_6_access_cm2.nc",
  "veg_class_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/vegetation.nc",
  "drainage_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/drainage.nc",
  "soil_texture_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/soil-texture.nc",
  "co2_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/co2.nc",
  "proj_co2_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/projected-co2_ssp1_2_6.nc",
  "ch4_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/ch4.nc",
  "proj_ch4_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/projected-ch4_ssp1.nc",
  "runmask_file": "/home/ejafarov/tem-dev/usbzb/runs/Q1/run-mask.nc",
  "topo_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/topo.nc",
  "fri_fire_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/fri-fire.nc",
  "hist_exp_fire_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/historic-explicit-fire.nc",
  "proj_exp_fire_file": "/home/ejafarov/tem-dev/usbzb/archive/zenodo_archive/workflows/ch4_manuscript_files/input_forcings/baseline/projected-explicit-fire.nc",
  "output_dir": "/home/ejafarov/tem-dev/wet-work/cmt31_nup/runs/N2_B1_mb015/output/",
  "restart_from": "",
  "output_spec_file": "/home/ejafarov/tem-dev/wet-work/cmt31_nup/specs/usbzb_spec_Nup_m.csv",
  "output_monthly": 1,
  "output_nc_eq": 1,
  "output_nc_sp": 0,
  "output_nc_tr": 1,
  "output_nc_sc": 1,
  "output_interval": 1
 },
 "calibration-IO": {
  "unique_pid_tag": "",
  "caldata_tree_loc": "/tmp"
 },
 "stage_settings": {
  "inter_stage_pause": false,
  "pr": {
   "env": true,
   "bgc": false,
   "nfeed": false,
   "avlnflg": false,
   "baseline": false,
   "dsb": false,
   "dsl": false,
   "dyn_lai": false
  },
  "eq": {
   "env": true,
   "bgc": true,
   "nfeed": true,
   "avlnflg": true,
   "baseline": true,
   "dsb": false,
   "dsl": true,
   "dyn_lai": true
  },
  "sp": {
   "env": true,
   "bgc": true,
   "nfeed": true,
   "avlnflg": true,
   "baseline": true,
   "dsb": false,
   "dsl": true,
   "dyn_lai": true
  },
  "tr": {
   "env": true,
   "bgc": true,
   "nfeed": true,
   "avlnflg": true,
   "baseline": true,
   "dsb": false,
   "dsl": true,
   "dyn_lai": true
  },
  "sc": {
   "env": true,
   "bgc": true,
   "nfeed": true,
   "avlnflg": true,
   "baseline": true,
   "dsb": false,
   "dsl": true,
   "dyn_lai": true
  }
 },
 "model_settings": {
  "cell_timelimit": 0,
  "dynamic_lai": 1,
  "n_rootless_fix": 2,
  "n_som_move_fix": 1,
  "deep_carbon": 0,
  "n_nonvasc_extract": 1,
  "daily_c_split": 0,
  "baseline_start": 1901,
  "baseline_end": 1931,
  "wetland": {
   "enabled": true,
   "ch4": true,
   "run_on": true,
   "groundwater": true,
   "cmts": [
    6,
    31,
    32,
    55,
    61,
    75,
    92,
    93
   ],
   "interception_scheme": 0,
   "diag": 1,
   "runon_energy_mode": 1,
   "ch4_conserve": 1,
   "ch4_daily": 0,
   "keep_expelled_water": 1
  },
  "energy_diag": 1,
  "eq_co2_mode": 1,
  "thermal_solver": 1,
  "thermal_substeps": 4,
  "richards_dt_fix": 1,
  "thermal_sfcc": 1,
  "snow_empty_cleanup": 1,
  "frozen_layer_liq_cap": 1,
  "nonvasc_trans_withdraw": 1,
  "n_safetynet_conserve": 1,
  "bnf_scheme": 2,
  "n_leach_fix": 1
 }
}