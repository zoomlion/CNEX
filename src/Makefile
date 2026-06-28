CXX = g++
BASE_FLAGS = -std=c++17 -O3 -march=native -pthread -I. -Irobin-map/include
LDFLAGS = -lpthread -lz

TARGETS = 02.validate 03.assembler
PROFILE_DIR = profile

all: 02.validate 03.assembler

02.validate: 02.validate.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) $< $(LDFLAGS) -o $@

03.assembler: 03.assembler.cpp hip/validator_core.hpp hip/debruijn_core.hpp
	$(CXX) $(BASE_FLAGS) $< $(LDFLAGS) -o $@

# ─── PGO instrumented build ───

pgo-gen-02: 02.validate.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) -fprofile-generate $< $(LDFLAGS) -o 02.validate.pgo

pgo-gen-03: 03.assembler.cpp hip/validator_core.hpp hip/debruijn_core.hpp
	$(CXX) $(BASE_FLAGS) -fprofile-generate $< $(LDFLAGS) -o 03.assembler.pgo

# ─── PGO optimized rebuild ───

pgo-use-02: 02.validate.cpp hip/validator_core.hpp
	$(CXX) $(BASE_FLAGS) -flto -fprofile-use $< $(LDFLAGS) -o 02.validate.pgo
	strip 02.validate.pgo

pgo-use-03: 03.assembler.cpp hip/validator_core.hpp hip/debruijn_core.hpp
	$(CXX) $(BASE_FLAGS) -flto -fprofile-use $< $(LDFLAGS) -o 03.assembler.pgo
	strip 03.assembler.pgo

clean:
	rm -f $(TARGETS) *.pgo
	rm -rf $(PROFILE_DIR)

.PHONY: all clean pgo-gen-02 pgo-gen-03 pgo-use-02 pgo-use-03
