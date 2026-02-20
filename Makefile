# ─────────────────────────────────────────────────────────────────────────────
# Minimum-Cost Genome Planner — release Makefile
# ─────────────────────────────────────────────────────────────────────────────
# Requires: SDSL-lite v2.  Quick start (first time):
#   make install_sdsl   # clone + build SDSL v2.1.1 into ./deps/sdsl/
#   make                # compile all four binaries into ./bin/
#
# Using an existing SDSL install:
#   make SDSL_PREFIX=/path/to/your/sdsl/install
#
# Targets:
#   all           – build all four binaries into ./bin/
#   install_sdsl  – clone + build SDSL v2.1.1 into SDSL_PREFIX
#   clean         – remove ./bin/
# ─────────────────────────────────────────────────────────────────────────────

# Auto-detect SDSL:
#   1. ./deps/sdsl           (local, created by  make install_sdsl)
#   2. ../fm_index_benchmarks/local  (shared install in the sibling dir)
# Override any time with:  make SDSL_PREFIX=/custom/path
ifeq ($(origin SDSL_PREFIX), undefined)
  ifneq ($(wildcard $(abspath ./deps/sdsl)/include/sdsl/csa_wt.hpp),)
    SDSL_PREFIX := $(abspath ./deps/sdsl)
  else
    SDSL_PREFIX := $(abspath ../fm_index_benchmarks/local)
  endif
endif

# Force g++ to avoid NVHPC ABI mismatches (see README)
ifneq ($(shell command -v g++ 2>/dev/null),)
  CXX := g++
else
  CXX := c++
endif

# Detect OpenMP support
OMP_TEST := $(shell echo 'int main(){}' | $(CXX) -fopenmp -x c++ - -o /dev/null 2>&1)
ifneq ($(OMP_TEST),)
  OMP_FLAG :=
else
  OMP_FLAG := -fopenmp
endif

CXXFLAGS := -O3 -std=c++17 -I$(SDSL_PREFIX)/include
LDFLAGS  := -L$(SDSL_PREFIX)/lib -Wl,-rpath,$(SDSL_PREFIX)/lib \
            -lsdsl -ldivsufsort -ldivsufsort64 -pthread

BINDIR := bin
BINS   := $(BINDIR)/create_index \
          $(BINDIR)/genome_planner_flex \
          $(BINDIR)/greedy_planner_clean \
          $(BINDIR)/max_block_greedy_clean

.PHONY: all clean install_sdsl check_sdsl

all: check_sdsl $(BINS)

# ── dependency check ─────────────────────────────────────────────────────────
check_sdsl:
	@if [ ! -f "$(SDSL_PREFIX)/include/sdsl/csa_wt.hpp" ]; then \
	  echo "[ERROR] SDSL not found at $(SDSL_PREFIX)"; \
	  echo "        Run  make install_sdsl  to build SDSL automatically."; \
	  echo "        Or:  make SDSL_PREFIX=/path/to/existing/sdsl/install"; \
	  exit 1; \
	fi

# ── bin/ directory ────────────────────────────────────────────────────────────
$(BINDIR):
	mkdir -p $@

# ── binaries ─────────────────────────────────────────────────────────────────
$(BINDIR)/create_index: create_index.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(BINDIR)/genome_planner_flex: genome_planner_flex.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OMP_FLAG) $< -o $@ $(LDFLAGS)

$(BINDIR)/greedy_planner_clean: greedy_planner_clean.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OMP_FLAG) $< -o $@ $(LDFLAGS)

$(BINDIR)/max_block_greedy_clean: max_block_greedy_clean.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OMP_FLAG) $< -o $@ $(LDFLAGS)

# ── install SDSL ──────────────────────────────────────────────────────────────
install_sdsl:
	bash install_sdsl.sh $(SDSL_PREFIX)

# ── clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(BINDIR)
