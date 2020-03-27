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
	examples/task_parallel \
	problems/sedov \
	problems/sedov2d
TARGETS = $(ALL_TARGETS)


# If a Makefile.in exists in this directory, then use it
-include Makefile.in


# Build macros
# =====================================================================
SRC := $(wildcard src/*.cpp) $(wildcard examples/*.cpp) $(wildcard problems/*.cpp)
OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)


# Build config
# =====================================================================
GIT_COMMIT := $(shell git rev-parse --short=7 HEAD)
MARA_H_TMP := $(shell mktemp -u make.XXXXXX)


# Build rules
# =====================================================================
default: $(TARGETS)
all: $(ALL_TARGETS)

mara: src/main.o src/core_unit_test.o src/test_app.o src/test_core.o src/test_mesh.o src/test_model.o src/test_physics.o
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

src/mara.hpp:
	@$(RM) $(MARA_H_TMP)
	@echo "#define MARA_GIT_COMMIT \"$(GIT_COMMIT)\""   >> $(MARA_H_TMP)
	@echo "#define MARA_INSTALL_PATH \"$(PWD)\""  >> $(MARA_H_TMP)
	@cmp -s $(MARA_H_TMP) $@ || (echo "[mara.hpp updated]"; cat $(MARA_H_TMP)>$@)
	@$(RM) $(MARA_H_TMP)

$(OBJ): src/mara.hpp

-include $(DEP)
