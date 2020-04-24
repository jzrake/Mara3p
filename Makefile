# =====================================================================
# Mara build system
# =====================================================================
#
#
# External library dependencies: HDF5, MPI
#
#
# Notes
# -----
#
# - A useful resource for techniques to process Makefile dependencies:
# www.microhowto.info/howto/automatically_generate_makefile_dependencies.html
#
# - Using -O0 rather than -O3 during development may reduce compilation time
# significantly.


# Build configuration (these may be over-ridden in the Makefile.in)
# =====================================================================
CXX      = mpicxx
CXXFLAGS = -std=c++17 -Isrc -Wall -O0 -MMD -MP
LDFLAGS  = -lhdf5
ALL_TARGETS = mara \
	examples/euler1d \
	examples/euler1d_moving_mesh \
	examples/computable \
	examples/euler1d_moving_mesh_plm \
	examples/hardware_check \
	examples/mhd3d \
	examples/serial_perf \
	problems/sedov \
	problems/sedov2d \
	problems/minidisk
TARGETS = $(ALL_TARGETS)

# Conditional compilation macros
MARA_COMPILE_EULER1D = 1
MARA_COMPILE_EULER2D = 1
MARA_COMPILE_LOCALLY_ISOTHERMAL = 1

# If a Makefile.in exists in this directory, then use it
-include Makefile.in


# Build macros
# =====================================================================
SRC := $(wildcard src/*.cpp) $(wildcard examples/*.cpp) $(wildcard problems/*.cpp)
OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)

OBJ_APP     := $(filter src/app_%.o,     $(OBJ))
OBJ_CORE    := $(filter src/core_%.o,    $(OBJ))
OBJ_MESH    := $(filter src/mesh_%.o,    $(OBJ))
OBJ_MODULE  := $(filter src/module_%.o,  $(OBJ))
OBJ_PROBLEM := $(filter src/problem_%.o, $(OBJ))
OBJ_TEST    := $(filter src/test_%.o,    $(OBJ))


# Build config
# =====================================================================
GIT_COMMIT := $(shell git rev-parse --short=7 HEAD)
MARA_H_TMP := $(shell mktemp -u make.XXXXXX)


# Build rules
# =====================================================================
default: $(TARGETS)

all: $(ALL_TARGETS)

mara: src/mara.o src/parallel_computable.o $(OBJ_APP) $(OBJ_CORE) $(OBJ_MESH) $(OBJ_MODULE) $(OBJ_PROBLEM) $(OBJ_TEST)
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/euler1d: examples/euler1d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/euler1d_moving_mesh: examples/euler1d_moving_mesh.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/euler1d_moving_mesh_plm: examples/euler1d_moving_mesh_plm.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/hardware_check: examples/hardware_check.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/mhd3d: examples/mhd3d.o src/scheme_mhd_v1.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/serial_perf: examples/serial_perf.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/computable: examples/computable.o src/parallel_computable.o
	$(CXX) -o $@ $^ $(LDFLAGS)

problems/sedov: problems/sedov.o
	$(CXX) -o $@ $^ $(LDFLAGS)

problems/sedov2d: problems/sedov2d.o src/scheme_sedov2d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

problems/minidisk: problems/minidisk.o problems/minidisk_io.o problems/minidisk_scheme.o src/parallel_computable.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(DEP) $(ALL_TARGETS)

Makefile.in:
	@touch $@

src/mara.hpp: Makefile.in
	@$(RM) $(MARA_H_TMP)
	@echo "#pragma once"                                                               >> $(MARA_H_TMP)
	@echo "#define MARA_GIT_COMMIT \"$(GIT_COMMIT)\""                                  >> $(MARA_H_TMP)
	@echo "#define MARA_INSTALL_PATH \"$(PWD)\""                                       >> $(MARA_H_TMP)
	@echo "#define MARA_COMPILE_EULER1D $(MARA_COMPILE_EULER1D)"                       >> $(MARA_H_TMP)
	@echo "#define MARA_COMPILE_EULER2D $(MARA_COMPILE_EULER2D)"                       >> $(MARA_H_TMP)
	@echo "#define MARA_COMPILE_LOCALLY_ISOTHERMAL $(MARA_COMPILE_LOCALLY_ISOTHERMAL)" >> $(MARA_H_TMP)
	@cmp -s $(MARA_H_TMP) $@ || (echo "[mara.hpp updated]"; cat $(MARA_H_TMP) > $@)
	@$(RM) $(MARA_H_TMP)

$(OBJ): src/mara.hpp

-include $(DEP)
