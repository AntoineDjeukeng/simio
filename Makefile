CXX := g++
CC  := gcc

CXXFLAGS := -O2 -std=c++17 -Iinclude -Ireader/include -pthread
CFLAGS   := -O2 -Ireader/include -pthread
LDFLAGS  := -pthread
LDLIBS   :=

# GCC < 9 needs explicit filesystem library
CXX_MAJOR := $(shell $(CXX) -dumpversion | cut -d. -f1)
ifeq ($(shell [ $(CXX_MAJOR) -lt 9 ] 2>/dev/null && echo 1),1)
LDLIBS += -lstdc++fs
endif

CPP_SRCS := \
  examples/xtc_all_properties.cpp \
  examples/all_props/channel_count_xz.cpp \
  examples/all_props/density_x.cpp \
  examples/all_props/density_z_in_x_channel.cpp \
  examples/all_props/dipole_x.cpp \
  examples/all_props/dipole_z_in_x_channel.cpp \
  examples/all_props/coord_x.cpp \
  examples/all_props/gating_flux.cpp \
  examples/all_props/jump_msd.cpp

C_SRCS := \
  reader/src/xtc.c \
  reader/src/stub_trr.c

CPP_OBJS := $(addprefix build/,$(CPP_SRCS:.cpp=.o))
C_OBJS   := $(addprefix build/,$(C_SRCS:.c=.o))
OBJS     := $(CPP_OBJS) $(C_OBJS)

OUT := bin/simio_xtc_all_properties

.PHONY: all clean

all: $(OUT)

$(OUT): $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf build bin
