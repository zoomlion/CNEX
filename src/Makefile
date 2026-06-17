CXX = g++
BASE_FLAGS = -std=c++17 -O3 -march=native -pthread -I. -Irobin-map/include
LDFLAGS = -lpthread -lz

TARGETS = 02.assembler
PROFILE_DIR = profile

all: 02.assembler

02.assembler: 02.assembler.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) $< $(LDFLAGS) -o $@

# ─── PGO instrumented build ───

pgo-gen: 02.assembler.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) -fprofile-generate $< $(LDFLAGS) -o 02.assembler.pgo

# ─── PGO optimized rebuild ───

pgo-use: 02.assembler.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) -flto -fprofile-use $< $(LDFLAGS) -o 02.assembler.pgo
	strip 02.assembler.pgo

clean:
	rm -f $(TARGETS) 02.assembler.pgo
	rm -rf $(PROFILE_DIR)

.PHONY: all clean pgo pgo-gen pgo-use
