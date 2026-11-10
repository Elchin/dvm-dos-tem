#!/bin/bash
# Stage C extra verification (after esc_vm.sh). <=4 concurrent runs.
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
W=$HOME/tem-dev/esc-work; H=$HOME/tem-dev/esc/tests/regression; CS=$HOME/tem-dev/regress/cs; V=$W/cs; T=$W/tools
while pgrep -f "^bash esc_vm.sh" > /dev/null; do sleep 10; done
echo "$(date) extra start"; uptime
mkdir -p $V
S1='{"thermal_solver":1}'
mkv(){ s=$1; n=$2; shift 2; [ -d $V/$n ] || python3 $T/mk_variant.py $s $V/$n "$@"; }
mkv $CS/ed_dc0_calib8 ts1s0_dc0 ALL '{"thermal_solver":1,"thermal_sfcc":0}' '{"thermal_solver":1,"thermal_sfcc":0}'
mkv $CS/ed_dc1_prod   ts1_dc1_prod ALL "$S1" "$S1"
mkv $CS/ed_dc0_calib8 rt_l2e TVC_Warming_nfeed,HPC_Warming_nfeed '{}' "$S1"
mkv $CS/ed_dc0_calib8 rt_e2l TVC_Warming_nfeed,HPC_Warming_nfeed "$S1" '{"thermal_solver":0}'
mkv $CS/ed_dc0_calib8 geo1 HPC_Warming_nfeed '{"thermal_solver":1,"geothermal_flux_mode":1,"geothermal_flux":0.06}' '{"thermal_solver":1,"geothermal_flux_mode":1,"geothermal_flux":0.06}'
mkv $CS/ed_dc0_calib8 geo2 HPC_Warming_nfeed '{"thermal_solver":1,"geothermal_flux_mode":2}' '{"thermal_solver":1,"geothermal_flux_mode":2}' --geomap 0.06
mkv $CS/ed_dc0_calib8 face1 HPC_Warming_nfeed '{"thermal_solver":1,"thermal_face_mode":1}' '{"thermal_solver":1,"thermal_face_mode":1}'
mkv $CS/ed_dc0_calib8 sub4 HPC_Warming_nfeed '{"thermal_solver":1,"thermal_substeps":4}' '{"thermal_solver":1,"thermal_substeps":4}'
mkv $CS/ed_dc0_calib8 katz HPC_Warming_nfeed '{"thermal_solver":1,"thermal_sfcc":0,"thermal_nl_solver":1}' '{"thermal_solver":1,"thermal_sfcc":0,"thermal_nl_solver":1}'
mkv $CS/ed_dc0_calib8 ts1rf HPC_Warming_nfeed '{"thermal_solver":1,"stefan_front_resist_fix":1}' '{"thermal_solver":1,"stefan_front_resist_fix":1}'
jobs_wait(){ while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 2; done; }
for s in ts1s0_dc0 ts1_dc1_prod rt_l2e rt_e2l geo1 geo2 face1 sub4 katz ts1rf; do
  for c in $(python3 -c "import json;print(' '.join(json.load(open('$V/$s/caseset.json'))['cases']))"); do
    [ -f $W/x_$s/$c/meta.json ] && continue
    mkdir -p $W/x_$s; nice -n 5 python3 $H/run_case.py $W/O0/dvmdostem $V/$s/$c $W/x_$s/$c > $W/x_$s/$c.out 2>&1 &
    jobs_wait
  done
done; wait
for d in $W/x_*/; do for c in $(ls -d $d*/); do n=$(basename $c)
  python3 $H/closure.py $c > $d/$n.closure.json 2> $d/$n.closure.err; echo "$(basename $d) $n closure_rc=$?"; done; done
echo "$(date) EXTRADONE"
