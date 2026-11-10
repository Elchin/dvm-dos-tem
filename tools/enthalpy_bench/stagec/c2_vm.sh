#!/bin/bash
# Stage C2 verification. usage: c2_vm.sh COMMIT [phase...]   phases: build ut leg ts1 rst
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
C=$1; shift; PH="${*:-build ut leg ts1 rst}"
W=$HOME/tem-dev/esc-work; R=$HOME/tem-dev/esc; H=$R/tests/regression; CS=$HOME/tem-dev/regress/cs; V=$W/cs
LIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
mk(){ d=$1; shift; (cd $d && mkdir -p obj && nice -n 5 make -j4 CC="$CXX $*" SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" SITE_SPECIFIC_LIBS="$LIBS" GIT_SHA=stageC > build.log 2>&1); echo "build $d [$*] rc=$? errors=$(grep -c 'error' $d/build.log)"; }
has(){ [[ " $PH " == *" $1 "* ]]; }
runset(){ bin=$1; cs=$2; out=$3; mkdir -p $out; for c in $(ls -d $cs/*/); do n=$(basename $c)
    nice -n 5 python3 $H/run_case.py $bin $c $out/$n > $out/$n.out 2>&1 &
    while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 2; done; done; wait
  for c in $(ls -d $out/*/); do n=$(basename $c); python3 $H/closure.py $c > $out/$n.closure.json 2>/dev/null; echo "$(basename $out) $n closure_rc=$?"; done; }
if has build; then
  for d in c2_O0 c2_O2; do rm -rf $W/$d; mkdir -p $W/$d; git -C $R archive $C | tar -x -C $W/$d; done
  mk $W/c2_O0; mk $W/c2_O2 -O2; md5sum $W/c2_O0/dvmdostem $W/c2_O2/dvmdostem
fi
if has ut; then
  (cd $W/c2_O0 && nice -n 5 make -j4 CC="$CXX" SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" SITE_SPECIFIC_LIBS="$LIBS" GIT_SHA=stageC unittest > ut.log 2>&1); echo "unittest rc=$?"
  grep "test cases:\|assertions:\|Status:" $W/c2_O0/ut.log; grep -B2 -A8 "ERROR" $W/c2_O0/ut.log | grep -v KNOWN | head -40
fi
if has leg; then
  for cs in dc0_calib8 dc1_calib8 dc0_prod dc1_prod; do
    rm -rf $W/c2reg0_$cs; bash $H/run_regression.sh -j 4 $W/ref/dvmdostem $W/c2_O0/dvmdostem $CS/$cs $W/c2reg0_$cs > $W/c2reg0_$cs.out 2>&1
    echo "$(date) solver0 $cs rc=$? $(grep 'RESULT\|bit-identical' $W/c2reg0_$cs.out | tr '\n' ' ')"
  done
fi
if has ts1; then
  for s in ts1_dc0_calib8:ts1_dc0 ts1_dc1_calib8:ts1_dc1; do runset $W/c2_O0/dvmdostem $CS/${s%%:*} $W/c2_${s##*:}; done
  runset $W/c2_O0/dvmdostem $V/ts1_dc1_prod $W/c2_ts1_dc1_prod
  runset $W/c2_O0/dvmdostem $V/ts1s0_dc0 $W/c2_ts1s0_dc0
fi
if has rst; then
  for f in $W/c2_ts1_dc0/HPC_Warming_nfeed/main/output/restart-sc.nc $W/c2reg0_dc0_calib8/B/HPC_Warming_nfeed/main/output/restart-sc.nc; do
    echo "$f: $(ncdump -h $f | grep -c ':thermal_solver = 1')"; done
  runset $W/c2_O0/dvmdostem $V/rt_l2e $W/c2_rt_l2e; runset $W/c2_O0/dvmdostem $V/rt_e2l $W/c2_rt_e2l
  grep -h "WARNING: restart" $W/c2_rt_*/*/main/run.log $W/c2_rt_*/*/*.out 2>/dev/null | sort | uniq -c
fi
echo "$(date) C2DONE"
