# Simple Makefile for simio skeleton
CC      ?= cc
AR      ?= ar
CFLAGS  ?= -std=c11 -O2 -Wall -Wextra -Werror -Iinclude
LDFLAGS ?= 
LIB     := build/libsimio.a

# Collect sources
SRC := $(shell find src -name '*.c')
OBJ := $(patsubst src/%.c,build/obj/%.o,$(SRC))

TESTS := tests/test_build.c
TEST_BINS := bin/test_build

all: $(LIB) $(TEST_BINS)

$(LIB): $(OBJ)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

bin/test_build: $(LIB) $(TESTS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(TESTS) $(LIB) -lpthread -lm -o $@

clean:
	rm -rf build bin

.PHONY: all clean