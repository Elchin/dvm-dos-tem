# Stage B compute-cost benchmark (standalone; not part of the model build).
#   make OPT=-O0   or   make OPT=-O2
CC  ?= g++
OPT ?= -O0
bench_enthalpy$(SUFFIX): bench_enthalpy.cpp ../../src/EnthalpySolver.cpp ../../include/EnthalpySolver.h
	$(CC) $(OPT) -g -std=c++11 -I../../include bench_enthalpy.cpp ../../src/EnthalpySolver.cpp -o $@
clean:
	rm -f bench_enthalpy bench_enthalpy_O0 bench_enthalpy_O2
