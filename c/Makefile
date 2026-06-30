# Makefile for the Vektorregelung C library, demo and tests.
#
# Targets:
#   make            build the static library (build/libvektorregelung.a)
#   make demo       build and run the closed-loop simulation demo
#   make test       build and run the self-contained test suite
#   make clean      remove build artifacts
#
# The demo and tests are built with VR_REAL_DOUBLE so their results match the
# Python reference. The library object files use the default `float` type,
# matching a typical embedded single-precision FPU.

CC      ?= cc
CSTD    ?= -std=c99
CFLAGS  ?= $(CSTD) -O2 -Wall -Wextra -Wpedantic -Iinclude
LDLIBS  := -lm

BUILD   := build
SRC     := $(wildcard src/*.c)
OBJ     := $(patsubst src/%.c,$(BUILD)/%.o,$(SRC))
LIB     := $(BUILD)/libvektorregelung.a

# Host (double-precision) objects for demo/test so they match Python.
HOSTOBJ := $(patsubst src/%.c,$(BUILD)/host_%.o,$(SRC))

.PHONY: all demo test clean
all: $(LIB)

$(LIB): $(OBJ)
	@mkdir -p $(BUILD)
	$(AR) rcs $@ $^

$(BUILD)/%.o: src/%.c
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/host_%.o: src/%.c
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -DVR_REAL_DOUBLE -c $< -o $@

demo: $(HOSTOBJ) examples/sim_demo.c
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -DVR_REAL_DOUBLE examples/sim_demo.c $(HOSTOBJ) \
		-o $(BUILD)/sim_demo $(LDLIBS)
	./$(BUILD)/sim_demo

test: $(HOSTOBJ) tests/test_vektorregelung.c
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -DVR_REAL_DOUBLE tests/test_vektorregelung.c $(HOSTOBJ) \
		-o $(BUILD)/test_vektorregelung $(LDLIBS)
	./$(BUILD)/test_vektorregelung

clean:
	rm -rf $(BUILD)
