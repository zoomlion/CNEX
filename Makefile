CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -pthread -Isrc -Isrc/robin-map/include
LDFLAGS  := -lpthread -lz
PREFIX   := $(CURDIR)

# Local parasail vendor (2.6.1 from the pip package; identical API surface to
# upstream 2.6.2 for the SIMD SW calls used by the in-validate T3 rescue).
PARASAIL_INC := -Isrc/parasail/include
PARASAIL_LIB := src/parasail/lib/libparasail.so
PARASAIL_RPATH := -Wl,-rpath,$(CURDIR)/src/parasail/lib

.PHONY: all install clean

all: mertable validate assemble

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
