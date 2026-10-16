CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -pthread -Isrc -Isrc/robin-map/include
LDFLAGS  := -lpthread -lz
PREFIX   := $(CURDIR)

# parasail is a git submodule (src/parasail, pinned to v2.6.2). It is built on
# demand with CMake and exposed at src/parasail/lib/libparasail.so for linking
# and for the T3 rescue SIMD SW calls used by validate.
PARASAIL_SRC   := src/parasail
PARASAIL_BUILD := $(PARASAIL_SRC)/build
PARASAIL_INC   := -I$(PARASAIL_SRC)
PARASAIL_LIB   := $(PARASAIL_SRC)/lib/libparasail.so
PARASAIL_RPATH := -Wl,-rpath,$(CURDIR)/$(PARASAIL_SRC)/lib

# pigz is a git submodule (src/pigz, pinned to v2.8). Its binary src/pigz/pigz is
# what config.py points at for (de)compression.
PIGZ := src/pigz/pigz

.PHONY: all install clean parasail pigz

all: mertable validate assemble pigz

parasail: $(PARASAIL_LIB)

pigz: $(PIGZ)

$(PIGZ):
	$(MAKE) -C src/pigz pigz

$(PARASAIL_LIB): $(PARASAIL_SRC)/parasail.h
	cmake -S $(PARASAIL_SRC) -B $(PARASAIL_BUILD) -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
	cmake --build $(PARASAIL_BUILD) -j
	mkdir -p $(PARASAIL_SRC)/lib
	cp -a $$(find $(PARASAIL_BUILD) -name 'libparasail.so*' -type f -o -name 'libparasail.so*' -type l) $(PARASAIL_SRC)/lib/

mertable: src/01.mertable.cpp src/hip/mertable_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

validate: src/02.validate.cpp src/hip/validator_core.hpp src/hip/panel_core.hpp $(PARASAIL_LIB)
	$(CXX) $(CXXFLAGS) $(PARASAIL_INC) $< $(LDFLAGS) $(PARASAIL_LIB) $(PARASAIL_RPATH) -lm -o $@

assemble: src/03.assembler.cpp src/hip/validator_core.hpp src/hip/debruijn_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

install: all
	mkdir -p $(PREFIX)/bin
	mv mertable validate assemble $(PREFIX)/bin/

clean:
	rm -rf mertable validate assemble bin/mertable bin/validate bin/assemble
	-$(MAKE) -C src/pigz clean
