#!/bin/bash
# Stage D driver (runs after c2_vm.sh). <=4 concurrent runs, nice 5. Usage: stageD_vm.sh [phases: mk eval sub gd]
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
W=$HOME/tem-dev/esc-work; H=$HOME/tem-dev/esc/tests/regression; CS=$HOME/tem-dev/regress/cs; V=$W/cs; T=$W/tools; D=$W/stageD
P=" $* "; has(){ [[ "$P" == *" $1 "* ]]; }
while pgrep -f "[c]2_vm.sh" > /dev/null; do sleep 15; done
echo "$(date) stageD start [$*]"; uptime; mkdir -p $D
mkv(){ s=$1; n=$2; shift 2; [ -d $V/$n ] || python3 $T/mk_variant.py $s $V/$n "$@"; }
jw(){ while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 2; done; }
runset(){ bin=$1; cs=$2; out=$3; mkdir -p $out; for c in $(ls -d $cs/*/); do n=$(basename $c)
    [ -f $out/$n/meta.json ] && continue
    nice -n 5 python3 $H/run_case.py $bin $c $out/$n > $out/$n.out 2>&1 & jw; done; }
closeall(){ for out in "$@"; do for c in $(ls -d $out/*/); do n=$(basename $c)
    python3 $H/closure.py $c > $out/$n.closure.json 2>/dev/null; echo "$(basename $out) $n closure_rc=$?"; done; done; }
WC=TVC_Warming_nfeed,HPC_Warming_nfeed
if has mk; then
  mkv $CS/ed_dc1_calib8 ts1s0_dc1 ALL '{"thermal_solver":1,"thermal_sfcc":0}' '{"thermal_solver":1,"thermal_sfcc":0}'
  for k in 1 2 4 8 24; do mkv $CS/ed_dc0_calib8 dsub$k $WC "{\"thermal_solver\":1,\"thermal_substeps\":$k}" "{\"thermal_solver\":1,\"thermal_substeps\":$k}"; done
  mkv $CS/ed_dc0_calib8 dleg_o2 $WC '{}' '{}'
  mkv $CS/ed_dc0_calib8 gd_leg_c8 HPC_Warming_nfeed '{"soil_giveup_diag":1}' '{"soil_giveup_diag":1}'
  mkv $CS/ed_dc0_calib8 gd_ts1_c8 HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_diag":1}' '{"thermal_solver":1,"soil_giveup_diag":1}'
  mkv $CS/ed_dc1_prod gd_leg_prod HPC_Warming_nfeed '{"soil_giveup_diag":1}' '{"soil_giveup_diag":1}'
  mkv $CS/ed_dc1_prod gd_ts1_prod HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_diag":1}' '{"thermal_solver":1,"soil_giveup_diag":1}'
  mkv $CS/ed_dc0_calib8 gr_ts1_c8 HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_retry":1}' '{"thermal_solver":1,"soil_giveup_retry":1}'
  mkv $CS/ed_dc1_prod gr_ts1_prod HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_retry":1}' '{"thermal_solver":1,"soil_giveup_retry":1}'
  mkv $CS/ed_dc0_calib8 gd_ts1s0_c8 HPC_Warming_nfeed '{"thermal_solver":1,"thermal_sfcc":0,"soil_giveup_diag":1}' '{"thermal_solver":1,"thermal_sfcc":0,"soil_giveup_diag":1}'
  mkv $CS/ed_dc1_calib8 gd_leg_c8d1 HPC_Warming_nfeed '{"soil_giveup_diag":1}' '{"soil_giveup_diag":1}'
  mkv $CS/ed_dc1_calib8 gd_ts1_c8d1 HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_diag":1}' '{"thermal_solver":1,"soil_giveup_diag":1}'
  mkv $CS/ed_dc1_calib8 gr_ts1_c8d1 HPC_Warming_nfeed '{"thermal_solver":1,"soil_giveup_retry":1,"soil_giveup_diag":1}' '{"thermal_solver":1,"soil_giveup_retry":1,"soil_giveup_diag":1}'
fi
if has eval; then runset $W/c2_O0/dvmdostem $V/ts1s0_dc1 $D/ts1s0_dc1; fi
if has gd; then for s in gd_leg_c8 gd_ts1_c8 gd_leg_prod gd_ts1_prod gr_ts1_c8 gr_ts1_prod gd_ts1s0_c8 gd_leg_c8d1 gd_ts1_c8d1 gr_ts1_c8d1; do runset $W/c2_O0/dvmdostem $V/$s $D/$s; done; fi
wait
# substep timing: O2 build, at most 2 concurrent so timings are comparable
if has sub; then for s in dleg_o2 dsub1 dsub2 dsub4 dsub8 dsub24; do mkdir -p $D/$s
  for c in TVC_Warming_nfeed HPC_Warming_nfeed; do [ -f $D/$s/$c/meta.json ] && continue
    nice -n 5 python3 $H/run_case.py $W/c2_O2/dvmdostem $V/$s/$c $D/$s/$c > $D/$s/$c.out 2>&1 &
    while [ $(jobs -r | wc -l) -ge 2 ]; do sleep 2; done; done; done; wait; fi
closeall $(ls -d $D/*/)
echo "$(date) STAGEDDONE"
# analysis (Stage D report): obs eval, tables, substeps, give-up diagnostics
if has ana; then cd $W
  python3 $T/stageD_eval.py $D/obs_stageD.json $D/eval_c2b.json leg_dc0=leg_dc0 leg_dc1=leg_dc1 legprod_dc1=c2reg0_dc1_prod/A ts1_dc0=c2_ts1_dc0 ts1_dc1=c2_ts1_dc1 s0_dc0=c2_ts1s0_dc0 s0_dc1=$D/ts1s0_dc1 ts1prod_dc1=c2_ts1_dc1_prod
  python3 $T/dtab.py $D/eval_c2b.json NoWarming; python3 $T/subtab.py
  (cd $D && python3 $T/gd.py gd_leg_c8 gd_ts1_c8 gd_ts1s0_c8 gd_leg_c8d1 gd_ts1_c8d1 gr_ts1_c8d1 gd_leg_prod gd_ts1_prod gr_ts1_c8 gr_ts1_prod && python3 $T/gd2.py gd_ts1_c8 gd_leg_c8d1 gd_ts1_c8d1 gd_leg_prod gd_ts1_prod)
fi
