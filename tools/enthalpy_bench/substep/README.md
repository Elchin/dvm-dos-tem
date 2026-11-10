Sub-step anomaly study (design doc section 15d). Run dir on calib-elchin: ~/tem-dev/ssf-work.
- mkvar.sh NAME N: case-set variant (thermal_substeps N) from cs/HPCW_s24src (md5-checked copy of Stage 0 cs/ts1s24/HPC_Warming_nfeed, S23 C1 inputs)
- runv.sh BIN CS RUN: one run_case.py run (20-yr spin + 75-yr main); tab.py RUN...: ALD / closure / dump-day table; same.py A B: value identity of all outputs
- dbg_patch.py: diagnostic-only patch (ENTH_DBG_FILE dump of the soil-water window), not part of the model
- reg.sh: thermal_solver 0 bit-identity (dc0_calib8, dc1_prod, prod_t125_s23) and ts1 s1 calib8 REF vs NEW
