#!/bin/bash
# runv.sh BIN CSNAME RUNNAME [ENV=..]
W=~/tem-dev/ssf-work; source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $HOME/miniconda3/envs/dvmdostem-build
export REG_LD_LIBRARY_PATH=$CONDA_PREFIX/lib; export REG_EXTRA_ENV="$4"
rm -rf $W/runs/$3; mkdir -p $W/runs
nice -n 5 python3 $W/run_case.py $1 $W/cs/$2 $W/runs/$3 > $W/runs/$3.out 2>&1; echo "DONE $3 rc=$? $(date +%T)" >> $W/runs/DONE.log
