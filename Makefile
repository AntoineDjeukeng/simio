# Simple Makefile for simio skeleton
CC      ?= cc
AR      ?= ar
CFLAGS  ?= -std=c11 -O2 -Wall -Wextra -Werror -Iinclude
LDFLAGS ?=
LIB     := build/libsimio.a

# Collect sources
SRC := $(shell find src -name '*.c')
OBJ := $(patsubst src/%.c,build/obj/%.o,$(SRC))

# Tests -> bin/<name>
TEST_SRCS := $(shell find tests -maxdepth 1 -name '*.c')
TEST_BINS := $(patsubst tests/%.c,bin/%,$(TEST_SRCS))

# Examples -> bin/examples/<name>
EX_SRCS := $(shell find examples -maxdepth 1 -name '*.c' 2>/dev/null)
EX_BINS := $(patsubst examples/%.c,bin/examples/%,$(EX_SRCS))

.PHONY: all clean tests examples

all: $(LIB) $(TEST_BINS) $(EX_BINS)

$(LIB): $(OBJ)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

bin/%: tests/%.c $(LIB)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< $(LIB) -lpthread -lm -o $@

bin/examples/%: examples/%.c $(LIB)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< $(LIB) -lpthread -lm -o $@

# Convenience: build & run all tests
tests: $(TEST_BINS)
	@set -e; \
	for t in $(TEST_BINS); do \
	  echo $$t; $$t; \
	done

examples: $(EX_BINS)

clean:
	rm -rf build bin
