#!/bin/bash
# Stage C verification on calib-elchin (<=4 cores, scratch dirs only). usage: esc_vm.sh STAGEC_COMMIT REF_COMMIT [phase...]
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
C=$1; REF=$2; shift 2; PH="${*:-build ut leg ts1 cmp time}"
W=$HOME/tem-dev/esc-work; mkdir -p $W; R=$HOME/tem-dev/dvm-dos-tem; H=$HOME/tem-dev/esc/tests/regression; CS=$HOME/tem-dev/regress/cs
LIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
mk(){ d=$1; shift; (cd $d && mkdir -p obj && nice -n 5 make -j4 CC="$CXX $*" SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" SITE_SPECIFIC_LIBS="$LIBS" GIT_SHA=stageC > build.log 2>&1); echo "build $d [$*] rc=$? errors=$(grep -c 'error:' $d/build.log) warnings=$(grep -c 'warning:' $d/build.log)"; }
has(){ [[ " $PH " == *" $1 "* ]]; }
runset(){ bin=$1; cs=$2; out=$3; mkdir -p $out; for c in $(ls -d $CS/$cs/*/); do n=$(basename $c)
    python3 $H/run_case.py $bin $c $out/$n > $out/$n.out 2>&1 &
    while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 2; done; done; wait; }
if has build; then
  echo "$(date) builds"; uptime
  for d in O0 O2 ref; do rm -rf $W/$d; mkdir -p $W/$d; done
  git -C $R archive $C | tar -x -C $W/O0; git -C $R archive $C | tar -x -C $W/O2; git -C $R archive $REF | tar -x -C $W/ref
  mk $W/ref; mk $W/O0; mk $W/O2 -O2
  md5sum $W/ref/dvmdostem $W/O0/dvmdostem $W/O2/dvmdostem
fi
if has ut; then
  echo "$(date) unit tests"
  (cd $W/O0 && nice -n 5 make -j4 CC="$CXX" SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" SITE_SPECIFIC_LIBS="$LIBS" GIT_SHA=stageC unittest > ut.log 2>&1); echo "unittest rc=$?"
  grep "test cases:\|assertions:\|Status:" $W/O0/ut.log
  (cd $W/O0/tests/unit && LD_LIBRARY_PATH=$CONDA_PREFIX/lib ./dvmdostem_unittest -ts=enthalpy_integration 2>&1 | tail -3)
  python3 $H/test_closure.py 2>&1 | tail -3
fi
if has leg; then
  for cs in dc0_calib8 dc1_calib8 dc0_prod dc1_prod; do
    echo "$(date) solver0 ref-vs-stageC $cs"; bash $H/run_regression.sh -j 4 $W/ref/dvmdostem $W/O0/dvmdostem $CS/$cs $W/reg0_$cs > $W/reg0_$cs.out 2>&1; echo "rc=$?"; grep "RESULT\|bit-identical" $W/reg0_$cs.out
  done
fi
if has ts1; then
  for s in dc0 dc1; do
    [ -d $CS/ts1_${s}_calib8 ] || python3 $HOME/tem-dev/regress/make_cs.py $CS/ed_${s}_calib8 $CS/ts1_${s}_calib8 --ms-json $HOME/tem-dev/regress/ms_ts1.json --name ts1_${s}_calib8
    echo "$(date) solver1 $s"; runset $W/O0/dvmdostem ts1_${s}_calib8 $W/ts1_${s}
    echo "$(date) legacy ed $s"; runset $W/O0/dvmdostem ed_${s}_calib8 $W/leg_${s}
    for c in $(ls -d $W/ts1_${s}/*/); do n=$(basename $c)
      python3 $H/closure.py $c > $W/ts1_${s}/$n.closure.json; echo "$s $n closure_rc=$?"; done
  done
fi
if has cmp; then
  for s in dc0 dc1; do python3 $HOME/tem-dev/compare_ts.py $W/leg_${s} $W/ts1_${s}; done
  python3 - <<P
import json,glob
for f in sorted(glob.glob("$W/ts1_*/*.closure.json")):
    r=json.load(open(f)); k=lambda x: r.get(x)
    print(f.split("/")[-2], r["case"], "PASS" if r["closure_pass"] else "FAIL", r["closure_fail"],
      "W%P=%.3f Cmax=%.4f N=%s res_max=%.2e resday_max=%.2e unattr_max=%.2e fail=%s solves=%s init=%.3g clip=%.3g bsnow=%.3g bwater=%.3g blayer=%.3g bfire=%.3g" % (
      k("W_resid_pctP"),k("C_resid_max_abs_annual"),k("N_resid_mean_gm2yr"),k("Enth_res_max_abs_Wm2"),k("Enth_resday_max_abs_Wm2"),k("Enth_unattr_max_abs_Wm2"),
      k("Enth_fallback_failures"),k("Enth_solves"),k("Enth_init_mean_Wm2"),k("Enth_clip_mean_Wm2"),k("Enth_bsnow_mean_Wm2"),k("Enth_bwater_mean_Wm2"),k("Enth_blayer_mean_Wm2"),k("Enth_bfire_mean_Wm2")))
P
fi
if has time; then
  [ -d $CS/ts1t_dc0_calib8 ] || python3 $HOME/tem-dev/regress/make_cs.py $CS/dc0_calib8 $CS/ts1t_dc0_calib8 --ms-json $HOME/tem-dev/regress/ms_ts1.json --name ts1t_dc0_calib8
  echo "$(date) timing (sequential, 1 process)"; uptime
  for b in O0 O2; do for cs in dc0_calib8 ts1t_dc0_calib8; do
    for c in $(ls -d $CS/$cs/*/); do n=$(basename $c); t0=$(date +%s.%N)
      python3 $H/run_case.py $W/$b/dvmdostem $c $W/time_${b}_${cs}/$n > /dev/null 2>&1; rc=$?; t1=$(date +%s.%N)
      echo "TIME $b $cs $n rc=$rc wall=$(echo "$t1 - $t0" | bc) load=$(cut -d' ' -f1 /proc/loadavg)"
    done; done; done
fi
echo "$(date) ESCDONE"
