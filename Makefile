CXX      := g++
CXXFLAGS := -std=c++17 -O3 -march=native -pthread -Isrc -Isrc/robin-map/include
LDFLAGS  := -lpthread -lz
PREFIX   := $(CURDIR)

.PHONY: all install clean

all: mertable validate assemble

mertable: src/01.mertable.cpp src/hip/mertable_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

validate: src/02.validate.cpp src/hip/validator_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

assemble: src/03.assembler.cpp src/hip/validator_core.hpp src/hip/debruijn_core.hpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

install: all
	mkdir -p $(PREFIX)/bin
	mv mertable validate assemble $(PREFIX)/bin/

clean:
	rm -rf mertable validate assemble bin/mertable bin/validate bin/assemble
