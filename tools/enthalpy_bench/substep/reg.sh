#!/bin/bash
# waits for the s48 jobs, then thermal_solver=0 regression (ref stageD-eval O0 vs fixAB O0) and ts1 s1 set
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib
W=~/tem-dev/ssf-work; H=~/tem-dev/ssf/tests/regression; CS=~/tem-dev/regress/cs; REF=~/tem-dev/esc-work/c2_O0/dvmdostem; NEW=$W/fixAB/dvmdostem
while [ $(grep -c "_s48 rc=0" $W/runs/DONE.log) -lt 3 ]; do sleep 10; done
for cs in dc0_calib8 dc1_prod prod_t125_s23; do
  rm -rf $W/reg0_$cs; bash $H/run_regression.sh -j 4 --no-repeat $REF $NEW $CS/$cs $W/reg0_$cs > $W/reg0_$cs.out 2>&1
  echo "$(date +%T) ts0 $cs rc=$? $(grep -i 'RESULT\|bit-identical\|A==B' $W/reg0_$cs.out | tr '\n' ' ')" >> $W/reg.log
done
for b in REF NEW; do bin=${!b}; out=$W/ts1s1_$b; rm -rf $out; mkdir -p $out
  for c in $(ls -d $CS/ts1_dc0_calib8/*/); do n=$(basename $c)
    nice -n 5 python3 $H/run_case.py $bin $c $out/$n > $out/$n.out 2>&1 &
    while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 2; done; done; wait
  for c in $(ls -d $out/*/); do n=$(basename $c); python3 $H/closure.py $c > $out/$n.closure.json 2>/dev/null; echo "$b $n closure_rc=$?" >> $W/reg.log; done
done
echo "ALLDONE $(date +%T)" >> $W/reg.log
