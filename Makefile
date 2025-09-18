# Simple Makefile for simio skeleton
CC      ?= cc
AR      ?= ar
CFLAGS  ?= -std=c11 -O2 -Wall -Wextra -Werror -Iinclude
LDFLAGS ?=
LIB     := build/libsimio.a

# Collect sources
SRC := $(shell find src -name '*.c')
OBJ := $(patsubst src/%.c,build/obj/%.o,$(SRC))

# All tests become binaries bin/<name>
TEST_SRCS := $(shell find tests -maxdepth 1 -name '*.c')
TEST_BINS := $(patsubst tests/%.c,bin/%,$(TEST_SRCS))

all: $(LIB) $(TEST_BINS)

$(LIB): $(OBJ)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

bin/%: tests/%.c $(LIB)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< $(LIB) -lpthread -lm -o $@

clean:
	rm -rf build bin
