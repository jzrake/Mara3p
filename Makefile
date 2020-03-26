# =====================================================================
# Copyright 2020, Jonathan Zrake
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# =====================================================================




# Default build macros (these may be over-ridden in the Makefile.in)
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
	problems/sedov2d
TARGETS = $(ALL_TARGETS)


# If a Makefile.in exists in this directory, then use it
-include Makefile.in


# Mara source variables
# =====================================================================
SRC := $(wildcard src/*.cpp) $(wildcard examples/*.cpp) $(wildcard problems/*.cpp)
OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)


# Lua source variables
# =====================================================================
LUA_SRC := $(filter-out src/lua/core/lua.c src/lua/core/luac.c, $(wildcard src/lua/core/*.c))
LUA_MOD := $(wildcard src/lua/modules/*.cpp)
LUA_OBJ := $(LUA_SRC:%.c=%.o) $(LUA_MOD:%.cpp=%.o)
LUA_DEP := $(LUA_SRC:%.c=%.d) $(LUA_MOD:%.cpp=%.d)


# Mara source variables
# =====================================================================
CORE_OBJ     = $(filter src/core_%.o,     $(OBJ))
PARALLEL_OBJ = $(filter src/parallel_%.o, $(OBJ))
PHYSICS_OBJ  = $(filter src/physics_%.o,  $(OBJ))
TEST_OBJ     = $(filter src/test_%.o,     $(OBJ))


# Build config
# =====================================================================
GIT_COMMIT := $(shell git rev-parse --short=7 HEAD)
MARA_H_TMP := $(shell mktemp -u make.XXXXXX)


# Build rules
# =====================================================================
default: $(TARGETS)

show:
	@echo LUA_OBJ: $(LUA_OBJ)
	@echo CORE_OBJ: $(CORE_OBJ)
	@echo TEST_OBJ: $(TEST_OBJ)
	@echo PHYSICS_OBJ: $(PHYSICS_OBJ)
	@echo PARALLEL_OBJ: $(PARALLEL_OBJ)

all: $(ALL_TARGETS)

mara: src/mara.o $(CORE_OBJ) $(PARALLEL_OBJ) $(PHYSICS_OBJ) $(TEST_OBJ) $(LUA_OBJ)
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

examples/mhd3d_tbp: examples/mhd3d_tbp.o src/scheme_mhd_v2.o src/scheme_mhd_rules.o src/app_ui.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/serial_perf: examples/serial_perf.o
	$(CXX) -o $@ $^ $(LDFLAGS)

examples/task_parallel: examples/task_parallel.o
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
	$(RM) $(OBJ) $(DEP) $(ALL_TARGETS) $(LUA_OBJ) $(LUA_DEP)

src/mara.hpp: .FORCE
	@$(RM) $(MARA_H_TMP)
	@echo "#define MARA_GIT_COMMIT \"$(GIT_COMMIT)\""        >> $(MARA_H_TMP)
	@echo "#define MARA_INSTALL_PATH \"$(PWD)\""             >> $(MARA_H_TMP)
	@cmp -s $(MARA_H_TMP) $@ || (echo "[mara.hpp updated]"; cat $(MARA_H_TMP) > $@)
	@$(RM) $(MARA_H_TMP)

.FORCE:

-include $(DEP)
