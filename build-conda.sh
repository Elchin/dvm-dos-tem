#!/bin/bash
# Native build of dvmdostem using the conda env "dvmdostem-build"
# (conda-forge gxx 13, boost 1.83, netcdf 4.9.2, jsoncpp, openblas/lapacke).
# Usage: ./build-conda.sh [make args]   (run from the source dir)
set -e
source /home/ejafarov/miniconda3/etc/profile.d/conda.sh
conda activate /home/ejafarov/miniconda3/envs/dvmdostem-build
make -j16 CC="$CXX" \
  SITE_SPECIFIC_INCLUDES="-I$CONDA_PREFIX/include" \
  SITE_SPECIFIC_LIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib" "$@"
