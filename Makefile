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


# Build configuration
# =====================================================================


# Default build macros (these may be over-ridden by entries in the Makefile.in)
CXX      = mpicxx
CXXFLAGS = -std=c++17 -Isrc -Wall -O0 -MMD -MP
LDFLAGS  = -lhdf5
ALL_TARGETS = mara \
	examples/euler1d \
	examples/euler1d_moving_mesh \
	examples/euler1d_moving_mesh_plm \
	examples/mhd3d \
	examples/mhd3d_tbp \
	examples/task_parallel \
	problems/sedov
TARGETS = $(ALL_TARGETS)


# If a Makefile.in exists in this directory, then use it
-include Makefile.in


# Build macros
# =====================================================================
SRC         := $(wildcard src/*.cpp) $(wildcard examples/*.cpp) $(wildcard problems/*.cpp)
OBJ         := $(SRC:%.cpp=%.o)
DEP         := $(SRC:%.cpp=%.d)


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

examples/mhd3d: examples/mhd3d.o src/scheme_mhd_v1.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/mhd3d_tbp: examples/mhd3d_tbp.o src/scheme_mhd_v2.o src/app_ui.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/task_parallel: examples/task_parallel.o
	$(CXX) -o $@ $^ $(LDFLAGS)

problems/sedov: problems/sedov.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(DEP) $(ALL_TARGETS)

-include $(DEP)
