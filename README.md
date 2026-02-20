# Minimum-Cost Genome Planner — release

This directory contains the four C++ tools that implement the minimum-cost
genome construction framework described in the paper.

---

## Tools

| Binary | Description |
|---|---|
| `create_index` | Build an FM-index over a source genome FASTA |
| `genome_planner_flex` | **Optimal DP planner** (linear or nonlinear synthesis cost) |
| `greedy_planner_clean` | Greedy baseline — Replication-First heuristic |
| `max_block_greedy_clean` | Greedy baseline — Max-Block heuristic |

All binaries accept `--help` for full parameter descriptions.

---

## Dependencies — SDSL-lite

All four tools require **[sdsl-lite v2](https://github.com/simongog/sdsl-lite)**
for FM-index construction and querying.

> **Known issue on clusters with NVHPC (nvc/nvc++).**  
> If your environment sets `CXX=nvc++`, mixing it with an SDSL build done
> under GCC will produce ABI errors at runtime (illegal instruction / abort).  
> The Makefile explicitly forces `g++`; if you need a different compiler set
> `CXX=g++ make` or `make FORCE_CXX=g++`.

> **SDSL version compatibility (`louds_tree.hpp`).**  
> Some older SDSL-lite checkouts (pre-v2.1.1) contain a bug in `louds_tree.hpp`
> where member names `m_select1`/`m_select0` were later renamed to
> `m_bv_select1`/`m_bv_select0`, causing a compile-time error when our code
> instantiates the FM-index templates.  
> `make install_sdsl` pins **v2.1.1** which includes the fix.  If you bring
> your own SDSL build, make sure it is v2.1.1 or a later commit from the
> `simongog/sdsl-lite` repository.

If SDSL is not yet installed, build it with the bundled script:

```bash
make install_sdsl        # clones sdsl-lite v2.1.1, installs into ./deps/sdsl/
```

Or point to an existing install:

```bash
make SDSL_PREFIX=/path/to/your/sdsl/install
```

---

## Building

```bash
cd release/
make            # compiles all four binaries into bin/
make clean      # removes bin/
```

The Makefile auto-detects `g++` and OpenMP availability.

---

## Quick usage

### Step 1 — Build an FM-index over the source genome

```bash
./bin/create_index source.fasta source.fm
```

### Step 2 — Run the planner of choice

```bash
# DP (optimal)
./bin/genome_planner_flex  <W> <target.fasta> <pcr> <join> <synth_linear> [synth_quad] source.fm

# Greedy replication-first
./bin/greedy_planner_clean <W> <target.fasta> <pcr> <join> <synth_linear> source.fm

# Greedy max-block
./bin/max_block_greedy_clean <W> <target.fasta> <pcr> <join> <synth_linear> source.fm
```

Output is CSV:  `filename, chromosome, length_bp, total_cost`  
followed by a `STATS_TOTAL` summary line.

---

## Parameter reference

| Parameter | Description |
|---|---|
| `W` | Maximum block length (bp). Reflects PCR amplification / synthesis length limits. |
| `target.fasta` | FASTA file of the genome to construct. |
| `pcr` (`c_reuse`) | Fixed cost per reused (PCR-amplified) block, independent of length. |
| `join` (`c_join`) | Fixed cost per junction between adjacent blocks. Not charged before the first block. |
| `synth_linear` (`c_s`) | Per-base synthesis cost coefficient (linear term). Synthesis cost = `c_s × L`. |
| `synth_quad` (`c_s2`) | *(DP only, optional)* Quadratic coefficient. Synthesis cost = `c_s × L + c_s2 × L²`. Omit or set to 0 for purely linear synthesis. |
| `source.fm` | FM-index file produced by `create_index`. |

---

## Examples

### Linear synthesis cost (viral experiments)

```bash
# W=500, pcr=5, join=1.5, synth=0.2 (linear only)
./bin/create_index source.fasta source.fm
./bin/genome_planner_flex 500 target.fasta 5 1.5 0.2 source.fm
./bin/greedy_planner_clean 500 target.fasta 5 1.5 0.2 source.fm
./bin/max_block_greedy_clean 500 target.fasta 5 1.5 0.2 source.fm
```

### Nonlinear synthesis cost — bacterial sweep

Parameters used in the paper (target: *E. coli* K-12 MG1655):

```bash
# W=1000, pcr=5, join=1.5, synth_linear=0.2, synth_quad=1e-4
./bin/create_index ecoli_source.fasta ecoli_source.fm
./bin/genome_planner_flex 1000 ecoli_target.fasta 5 1.5 0.2 1e-4 ecoli_source.fm
```

### Nonlinear synthesis cost — eukaryotic sweep

Parameters used in the paper (target: *S. cerevisiae* S288C):

```bash
# W=800, pcr=5, join=1.5, synth_linear=0.2, synth_quad=1e-4
./bin/create_index yeast_source.fasta yeast_source.fm
./bin/genome_planner_flex 800 scerevisiae_target.fasta 5 1.5 0.2 1e-4 yeast_source.fm
```

---

## Cost model

A block of length *L* is either **reused** (PCR-amplified) or **synthesized**:

- **Reuse** — the block occurs as an exact substring of the source genome  
  (tested via FM-index query).  Cost = `c_reuse` (fixed, independent of *L*).
- **Synthesis** — Cost = `c_s × L` (linear) or `c_s × L + c_s2 × L²` (nonlinear).
- **Join** — Cost = `c_join` per boundary between adjacent blocks  
  (not charged before the first block).

The DP planner finds the globally optimal partition of the target into blocks  
of length ≤ W that minimises total cost.  The greedy heuristics are provided  
as fast, approximate baselines.

---

## Notes

- **SDSL FM-index construction is single-threaded.** Parallelise by running  
  multiple `create_index` invocations simultaneously (e.g. via a Slurm job  
  array), not by giving more CPUs to a single build.
- The planner binaries use OpenMP internally to parallelise across chromosomes  
  when available.
- Non-ACGT characters in the FASTA are stripped before indexing; the planner  
  operates on the cleaned sequence.
