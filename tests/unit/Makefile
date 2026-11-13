# Unit tests for dvm-dos-tem (doctest, header-only, vendored in this directory).
# Normally invoked from the top-level Makefile:  make unittest
# (which first (re)builds the model objects in ../../obj and passes CC/INCLUDES/LIBDIR/LIBS).
# Standalone (objects must already exist in ../../obj):
#   make -C tests/unit CC=g++ INCLUDES=... LIBDIR=... run
CC      ?= g++
ROOT    := ../..
TFLAGS  := -g -std=c++20 -fPIC -DBOOST_ALL_DYN_LINK -Wno-deprecated-declarations
LIBS    ?= -lnetcdf -lboost_system -lboost_filesystem -lboost_program_options \
           -lboost_thread -lboost_log -ljsoncpp -lpthread -lreadline -llapacke
MODEL_OBJECTS ?= $(filter-out $(ROOT)/obj/TEM.o,$(wildcard $(ROOT)/obj/*.o))
TEST_SRCS := $(wildcard test_*.cpp)
TEST_OBJS := $(patsubst %.cpp,obj/%.o,$(TEST_SRCS))
BIN     := dvmdostem_unittest

all: $(BIN)

obj/%.o: %.cpp test_access.h test_helpers.h nc_fixtures.h doctest.h $(wildcard $(ROOT)/include/*.h)
	@mkdir -p obj
	$(CC) -c $(TFLAGS) $(INCLUDES) -I$(ROOT)/include -I$(ROOT)/external/Fireweed/C++ -I$(ROOT)/external/Burnup/C++ $< -o $@

$(BIN): $(TEST_OBJS) $(MODEL_OBJECTS)
	$(CC) -o $@ $(TEST_OBJS) $(MODEL_OBJECTS) $(LIBDIR) $(LIBS)

run: $(BIN)
	./$(BIN) $(DOCTEST_ARGS)

clean:
	rm -rf obj $(BIN)

.PHONY: all run clean
