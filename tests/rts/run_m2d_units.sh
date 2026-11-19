#!/bin/bash
# Build + run tests/rts/test_m2d_units.cpp against the worktree obj/*.o (build the model first).
set -u
cd "$(dirname "$0")/../.."
OUT=${1:-/tmp/test_m2d_units}
OBJS=$(ls obj/*.o | grep -vx "obj/TEM[.]o")
${CXX:-c++} -std=c++20 -g -DBOOST_ALL_DYN_LINK -Iinclude ${CONDA_PREFIX:+-I$CONDA_PREFIX/include} -I external/Fireweed/C++ -I external/Burnup/C++ \
  -o "$OUT" tests/rts/test_m2d_units.cpp $OBJS ${CONDA_PREFIX:+-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib} \
  -lnetcdf -lboost_system -lboost_filesystem -lboost_program_options -lboost_thread -lboost_log -ljsoncpp -lpthread -lreadline -llapacke || exit 3
"$OUT"
