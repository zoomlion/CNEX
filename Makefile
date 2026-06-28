CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -pthread -Isrc -Isrc/robin-map/include
LDFLAGS  := -lpthread -lz
PREFIX   := $(CURDIR)

.PHONY: all install clean

all: validate assemble

validate: src/02.validate.cpp src/hip/validator_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

assemble: src/03.assembler.cpp src/hip/validator_core.hpp src/hip/debruijn_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

install: all
	mkdir -p $(PREFIX)/bin
	mv validate assemble $(PREFIX)/bin/

clean:
	rm -rf validate assemble bin/validate bin/assemble
