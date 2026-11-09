#!/bin/bash
# Stage B compute-cost evaluation on calib-elchin (<= 4 cores). Never touches the branch tree's
# sources: every model build is a throw-away `git archive` copy instrumented by bench_patch.py.
# usage: run_bench.sh TREE OUTDIR [CASESET_DIR]
#   TREE    worktree of enthalpy/stageB (committed HEAD is used)
#   OUTDIR  scratch output (model copies, dumps, run dirs, results)
set -e
TREE=$(readlink -f $1); OUT=$(readlink -f $2); CS=${3:-$HOME/tem-dev/regress/cs/calib8}
CASES="TVC_NoWarming_nfeed HPC_NoWarming_nfeed"
mkdir -p $OUT
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
B=$TREE/tools/enthalpy_bench
for opt in O0 O2; do
  M=$OUT/model_$opt; rm -rf $M; mkdir -p $M/obj
  git -C $TREE archive HEAD | tar -x -C $M
  python3 $B/bench_patch.py $M
  extra=""; [ $opt = O2 ] && extra="-O2 -Wno-error=aggressive-loop-optimizations"  # pre-existing EnvData.cpp:952/997 z[il+1] read past MAX_SOI_LAY trips -Werror at -O2
  (cd $M && nice -n 5 make -j4 CC="$CXX $extra" SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" \
     SITE_SPECIFIC_LIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib" GIT_SHA=bench_$opt > build.log 2>&1)
  (cd $B && $CXX $extra -g -std=c++11 -I../../include bench_enthalpy.cpp ../../src/EnthalpySolver.cpp -o $OUT/bench_enthalpy_$opt)
done
H=$TREE/tests/regression
jobs=()
for opt in O0 O2; do for c in $CASES; do
  # dump run (first 3650 days of each of spin000 and main) and timing-only run (no dump)
  REG_EXTRA_ENV="DVM_THERMAL_BENCH=$OUT/dump_${opt}_${c} DVM_THERMAL_BENCH_DAYS=3650" \
    python3 $H/run_case.py $OUT/model_$opt/dvmdostem $CS/$c $OUT/rd_${opt}_${c} > $OUT/rd_${opt}_${c}.log 2>&1 &
  REG_EXTRA_ENV="DVM_THERMAL_BENCH=$OUT/time_${opt}_${c} DVM_THERMAL_BENCH_DAYS=0" \
    python3 $H/run_case.py $OUT/model_$opt/dvmdostem $CS/$c $OUT/rt_${opt}_${c} > $OUT/rt_${opt}_${c}.log 2>&1 &
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 5; done
done; done
wait
for opt in O0 O2; do for c in $CASES; do
  for r in rt rd; do for st in spin000 main; do
    echo "$opt $c $r $st $(grep -h THERMAL_BENCH $OUT/${r}_${opt}_${c}/$st/run.log) $(grep -h EXIT_CODE $OUT/${r}_${opt}_${c}/$st/run.log)"
  done; done
done; done > $OUT/inmodel_timing.txt
for opt in O0 O2; do
  $OUT/bench_enthalpy_$opt -r 3 $OUT/dump_${opt}_*.csv > $OUT/bench_$opt.txt
  for c in $CASES; do $OUT/bench_enthalpy_$opt -r 3 $OUT/dump_${opt}_${c}.*.csv > $OUT/bench_${opt}_${c}.txt; done
done
echo DONE
