CXX := g++
CC  := gcc

CXXFLAGS := -O2 -std=c++17 -Iinclude -Isrc/io/xtc_c/reader/include -pthread
CFLAGS   := -O2 -Isrc/io/xtc_c/reader/include -pthread
LDFLAGS  := -pthread
LDLIBS   :=

# GCC < 9 needs explicit filesystem library
CXX_MAJOR := $(shell $(CXX) -dumpversion | cut -d. -f1)
ifeq ($(shell [ $(CXX_MAJOR) -lt 9 ] 2>/dev/null && echo 1),1)
LDLIBS += -lstdc++fs
endif

CPP_SRCS := \
  apps/xtc_all_properties/main.cpp \
	src/runtime/run_config.cpp \
	src/analysis/intrinsics/x_grid_cache.cpp \
	src/analysis/intrinsics/z_grid_cache.cpp \
	src/analysis/intrinsics/channel_roi.cpp \
	src/analysis/intrinsics/in_channel_mask.cpp \
  apps/xtc_all_properties/all_props/axial_profile_x_export.cpp \
  apps/xtc_all_properties/all_props/channel_count_xz.cpp \
  apps/xtc_all_properties/all_props/density_x.cpp \
  apps/xtc_all_properties/all_props/water_atom_density_x.cpp \
  apps/xtc_all_properties/all_props/density_z_in_x_channel.cpp \
  apps/xtc_all_properties/all_props/dipole_x.cpp \
  apps/xtc_all_properties/all_props/dipole_z_in_x_channel.cpp \
  apps/xtc_all_properties/all_props/coord_x.cpp \
  apps/xtc_all_properties/all_props/nacl_cluster_x.cpp \
  apps/xtc_all_properties/all_props/gating_flux.cpp \
  apps/xtc_all_properties/all_props/jump_msd.cpp

C_SRCS := \
  src/io/xtc_c/reader/src/xtc.c \
  src/io/xtc_c/reader/src/stub_trr.c

CPP_OBJS := $(addprefix build/,$(CPP_SRCS:.cpp=.o))
C_OBJS   := $(addprefix build/,$(C_SRCS:.c=.o))
OBJS     := $(CPP_OBJS) $(C_OBJS)

OUT := bin/simio_xtc_all_properties
MINI_TOPO_SRC := apps/mini_gro_topology/main.cpp
MINI_TOPO_BIN := bin/simio_mini_gro_topology
SCAN_XTC_SRC := apps/scan/scan_xtc.cpp
SCAN_XTC_OBJ := build/apps/scan/scan_xtc.o
SCAN_XTC_BIN := bin/simio_scan_xtc

.PHONY: all clean scan_xtc
all: $(OUT) $(MINI_TOPO_BIN) $(SCAN_XTC_BIN)

scan_xtc: $(SCAN_XTC_BIN)

$(OUT): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(MINI_TOPO_BIN): $(MINI_TOPO_SRC)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

$(SCAN_XTC_BIN): $(SCAN_XTC_OBJ) $(C_OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(SCAN_XTC_OBJ) $(C_OBJS) $(LDLIBS)

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build bin

# ---- unit tests (minimal) ----
bin/test_cache_store: tests/unit/test_cache_store.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Iinclude -o $@ $<

test_cache_store: bin/test_cache_store
	./bin/test_cache_store

bin/test_cache_store_strict: tests/unit/test_cache_store_strict.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Iinclude -o $@ $<

test_cache_store_strict: bin/test_cache_store_strict
	./bin/test_cache_store_strict
