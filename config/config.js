{
  "general": {
    "output_global_attributes": {
        "run_name": "Demo dvmdostem run",
        "description": "A basic test run. Data is not scientifically sound."
    }
  },

  "IO": {
    "parameter_dir":      "parameters/",
    "hist_climate_file":  "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/historic-climate.nc",
    "proj_climate_file":  "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/projected-climate.nc",
    "veg_class_file":     "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/vegetation.nc",
    "drainage_file":      "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/drainage.nc",
    "soil_texture_file":  "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/soil-texture.nc",
    "co2_file":           "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/co2.nc",
    "proj_co2_file":      "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/projected-co2.nc",
    "ch4_file":           "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/ch4.nc",
    "proj_ch4_file":      "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/projected-ch4.nc",
    "runmask_file":       "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/run-mask.nc",
    "topo_file":          "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/topo.nc",
    "fri_fire_file":      "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/fri-fire.nc",
    "hist_exp_fire_file": "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/historic-explicit-fire.nc",
    "proj_exp_fire_file": "demo-data/cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10/projected-explicit-fire.nc",
    "fire_fuel_model_file": "external/Fireweed/Input_Files/Standard_Fuel_Models_Table.csv",
    "fire_fosberg_a_file": "external/Fireweed/Input_Files/NWGC_Fosberg_Table_A.txt",
    "fire_fosberg_b_file": "external/Fireweed/Input_Files/NWGC_Fosberg_Table_B.txt",
    "fire_fosberg_c_file": "external/Fireweed/Input_Files/NWGC_Fosberg_Table_C.txt",
    "fire_fosberg_d_file": "external/Fireweed/Input_Files/NWGC_Fosberg_Table_D.txt",
    "output_dir":         "output/",
    "restart_from":       "", // default
    "output_spec_file":   "config/output_spec.csv",
    "output_monthly":     1, //JSON specific
    "output_nc_eq":       0,
    "output_nc_sp":       0,
    "output_nc_tr":       1,
    "output_nc_sc":       1, 
    "output_interval":    1 
  },

  // Define storage locations for json files generated and used
  // during calibration. The calibration-viewer.py program will read
  // this config file and these settings to determine where to look for
  // the json files. dvmdostem will create a directory tree that looks
  // like this:
  //  dvmdostem/
  //  └── calibration
  //      ├── daily
  //      │   ├── year_00000_daily_drivers.text
  //      │   └── year_00001_daily_drivers.text
  //      ├── monthly
  //      │   ├── 0000000.json
  //      │   ├── ...
  //      │   └── 0000011.json
  //      └── yearly
  //          └── 00000.json
  "calibration-IO": {
    "unique_pid_tag": "",
    "caldata_tree_loc": "/tmp"

    // NOTE: It is generally reccomended that the files be kept in the /tmp
    // directory so that the operating system will clean up the files, as the
    // output can be voluminous, especially when generating monthly or daily
    // files. The main reason to have this location configurable is so that
    // we will be able to run dvmdostem on Atlas under the control of PEST and
    // can use the compute node-specific /scratch directories and keep different
    // running instances of dvmdostem from overwriting eachothers json files.
  },
  "stage_settings": {
    "inter_stage_pause": false,
    //Per-stage module settings
    "pr": {
      "env": true,
      "bgc": false,
      "nfeed": false,
      "avlnflg": false,
      "baseline": false,
      "dsb": false,
      "dsl": false,
      //FW_NOTES: Fire will never be on in the pre-run so we don't need a switch?
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
      "fire": false,
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
      "fire": false,
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
      "fire": false,
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
      "fire": false,
      "dyn_lai": true
    }

    // maybe less confusing if these settings are only available from cmd line?
    //"tr_yrs": 109,
    //"sc_yrs": 100

    // ??
    //"restartfile_dir": "DATA/Toolik_10x10_30yrs/" // location for restart-XX.nc file

  },

  "module_settings": {
    "dsb": {
      "fire": {
        // Wildfire ignitions modes: [Note: In development.  Numbers may change!]
        // 0: Default/old: FRI for [PR/]EQ/SP, explicit for TR/SC.
        // 1: Use fire return interval for this stage if fire is on.
        // 2+: Reserved for future use.
        // Fire is alway off in the PR stage and can only be FRI for EQ/SP so no switches are
        // provided for these stages.
        "ignition_tr": 0,
        "ignition_sc": 0,

        // More future controls per stage can go here...
        //"severity_pr":
        //"severity_eq":
        // ...

        //Turn on the process based wildfire model:
        //The fire settings that follow only pertain if this switch is set to true.
        "fire_process_wildfire": true,

        //The maximum snow depth (cm) where a fire can start, i.e. the extinguishment limit.  Set to
        //a negative value to ignore snow presence:
        //FW_NOTE: This could be made CMT specific like vsmburn.
        "fire_max_snow": 1.0,
        //Calculate dynamic fuel moisture for dynamic fuel models:
        "fire_dynamic_fuel": true,
        //Should moss be treated as a dead fuel (true) or live fuel (false):
        "fire_moss_as_dead_fuel": true,
        //Calculate the fuel bed depth, AKA delta, (true) or use a fixed depth (false):
        "fire_calculate_delta": false,
        //Ground fire model paramters and settings:
        //The fraction of aboveground fire heat (HPA) that enters the soil:
        "fire_heat_frac_to_soil": 0.25,
        //The fraction of heat lost from the soil:
        "fire_gf_heat_loss_factor": 0.50,
        //The maximum heat penetration depth for the soil surface heat transfer approximation (cm):
        "fire_gf_surface_pd": 1.0,
        //The maximum heat penetration depth for the smoldering heat transfer approximation (cm):
        "fire_gf_smolder_pd": 3.0,
        //The soil layer thickness for the simulation (cm):
        "fire_gf_layer_thickness": 0.1,
        //A temporary parameter to force the fuel model choice:
        //Default value of -1 indicates to ignore the setting.
        "fire_temp_fm": -1,
        //The following parameters to provide alternate meterology at the time of the fire
        //(possibly temporary).
        //Fire wind speed (m/min), -1 to ingnore:
        "fire_windspeed": -1.0,
        //Air temperature at the time of the fire (C), -300 to ignore:
        "fire_tempair": -300.0,
        //Percent relative humidity at the time of fire (%), -1 to ingnore:
        "fire_rh_pct": -1.0,
        //The fraction of dead woody debris that is consumed with fire.  1.0 matches old model behavior.
        "fire_dwd_combust_frac": 1.0
      }
      // Future disturbances go here...
      //"insect": { ... },
      //"thermokarst":{ ... },
    }
  },

  "model_settings": {
    "cell_timelimit": 0, //Run time limit per cell in seconds. 0 for no limit.
    "dynamic_lai": 1,                   // from model (1) or from input (0)
    "n_rootless_fix": 2,                // N1 fix (v5o default 2; 0 = legacy: rootless-layer net mineralisation N lost)
    "n_som_move_fix": 1,                // N2 fix (v5o default 1; 0 = legacy: SOM N moved shallow->deep overwritten)
    "deep_carbon": 1,                   // deep soil C burial (deep_carbon_v2 scheme; default 1; 0 = legacy v5q burial, bit-identical)
    "n_nonvasc_extract": 1,             // N3 fix (v5o default 1; 0 = legacy: moss N uptake not extracted from soil)
    "daily_c_split": 0,                 // v5p: statistical daily split of the monthly C fluxes, output only (0 off, 1 on; the *_DS outputs need 1)
    "baseline_start": 1901,  //start year for baseline EQ climate
    "baseline_end": 1931,    //end year for baseline EQ climate
    // Wetland bundle (regional_wetland_deepc merge, phase 3). OFF by default; enabled = false is
    // bit-identical to the pre-merge fork. The old flat keys "ch4", "run-on", "groundwater"
    // are deprecated aliases (used only if this block is absent).
    "wetland": {
      "enabled": false,                       // master switch
      "ch4": true,                            // CH4 production/oxidation/transport (only if enabled)
      "run_on": false,                        // run-on hydrology (only if enabled)
      "groundwater": false,                   // calibrated groundwater recharge (only if enabled)
      "cmts": [6, 31, 32, 55, 61, 75, 92, 93],// CMTs treated as wetland
      "interception_scheme": 0,               // 0 = legacy/v5q interception (all cohorts); 1 = upstream CLM5-style formula
      "runon_energy_mode": 0,                 // 0 = legacy (recycle own runoff, puddle has no heat). 1 = no self-recycle; run-on/puddle carry energy at max(Tair, 0 C) only when runon_status is true; puddle ice is in the water budget and the restart
      "ch4_conserve": 0,                      // 0 = legacy CH4 bookkeeping. 1 = full-day dt, FV-CN diffusion, positivity limiter, fluxes scaled
      "ch4_daily": 0,                         // 0 = CH4Flux in the monthly loop (end-of-month state every day). 1 = daily env loop
      "reco_ch4ox_output": 0,                 // output only. 0 = legacy (CH4 oxidation CO2 not in RH/RECO/NEE). 1 = add oxidation CO2. 2 = also add CH4 efflux
      "keep_expelled_water": 0                // 0 = legacy (freeze-expansion / organic-shrink water leaves as qdrain). 1 = poorly drained wetland-CMT cells keep it (soil, else puddle, else qover)
    },
    // thermal_solver defaults to 0 (legacy Stefan/CN). With thermal_solver = 1,
    // thermal_substeps defaults to 4 (Elchin, 2026-11-10); allowed 1, 2, 4, 8, 12, 24 or 48.
    // Do not set thermal_substeps with thermal_solver = 0 (FATAL unless it is 1).
    // richards_dt_fix defaults to 0 (legacy, bit-identical) in the code and in this
    // production template. For new-solver work and the Stage E baseline set it to 1
    // explicitly (Elchin, 2026-11-10); production legacy stays 0 until recalibration.
    // New-solver example (do not uncomment in this production template):
    // "thermal_solver": 1,
    // "thermal_sfcc": 1,
    // "thermal_substeps": 4,
    // "richards_dt_fix": 1
    // frozen_layer_liq_cap (default 0 = legacy): 1 = liquid in partly frozen layers limited to the
    // pore space left by ice (validity check + Richards groundwater input)
    // snow_empty_cleanup defaults to 0 (legacy, bit-identical): 1 removes emptied snow layers
    // (no ice/liquid, zero thickness) at the end of each day's snow update, both solvers.
//    //"dynamic_climate": 0,
//    //"varied_co2": 0,
//    //"fire_severity_as_input": 0,    // fire sev. as input or ??
//    //"output_starting_year": -9999
  }

//  "output_switches": {
//    "daily_output": 0,
//    "monthly_output": 0,
//    "yearly_output": 1,
//    "summarized_output": 0,
//    "soil_climate_output": 0
//  }
}
