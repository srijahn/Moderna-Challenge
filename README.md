# WISER × Moderna Quantum Challenge — RNA Folding

Exploring RNA secondary structure prediction (minimum free energy folding) using
classical benchmarking (ViennaRNA) and quantum/quantum-inspired methods
(QUBO formulation + QAOA/CVaR-VQE in PennyLane).

## Requirements

### Software
- **Python 3.9–3.11** (ViennaRNA's Python bindings can be picky about newer
  versions — 3.11 or below is safest)
- pip (or conda)

### Python packages

| Package      | Used for                                  |
|--------------|--------------------------------------------|
| `ViennaRNA`  | Classical MFE structure/energy benchmark (`import RNA`); real thermodynamic evaluation of candidate structures (`structure_metrics.py`) |
| `pennylane`  | QAOA / CVaR-VQE quantum circuit simulation |
| `numpy`      | QUBO matrices, numerical work              |
| `pandas`     | Reading/writing CSV data and results       |
| `matplotlib` | Scaling and noise plots                    |

Install everything with:

```bash
pip install ViennaRNA pennylane numpy pandas matplotlib
```

```bash
pip install -r requirements.txt
```

> **Windows note:** `ViennaRNA` is a C-extension package and doesn't always
> have a prebuilt wheel for Windows. If `pip install ViennaRNA` fails, try:
> - `conda install -c bioconda viennarna`, or
> - running the project inside WSL (Windows Subsystem for Linux)

---

## Folder structure

The repo is already organized this way — every script reads/writes these
paths directly, nothing needs to be moved before running:

```
Moderna/
├── data/
│   └── sequences.csv            (legacy — see note below)
├── results/
│   ├── noise_analysis.csv
│   ├── final_results.csv
│   ├── batch_accuracy.csv
│   ├── scaling_analysis_real.csv
│   ├── cvar_vqe_scaling_analysis.csv
│   ├── statistical_benchmark.csv
│   ├── statistical_benchmark_summary.csv
│   ├── noise_plot.png
│   ├── scaling_plot_real.png
│   └── cvar_vqe_scaling_plot.png
├── *.py                     (all scripts)
├── Rna_basics.docx
├── Moderna_RNA_Folding_Quantum_Report.docx   (Task 7 final report)
└── Moderna - WISER Quantum Challenge [SHARED].pdf
```

> **Note:** every script writes only into `results/`

---

## How to run

Sanity-check ViennaRNA first:
```bash
python test_vienna.py
```

Then:

1. **Pairing / QUBO basics** (print output only, no files written)
   ```bash
   python pair_finder.py
   python conflict_detector.py
   python structure_generator.py
   python simple_qubo.py
   python qubo_energy.py
   python random_solver.py
   ```

2. **Quantum circuit demos**
   ```bash
   python quantum_test.py
   ```

   **Real RNA QUBO solvers (two independent methods):**
   ```bash
   python qaoa_rna_solver.py
   python cvar_vqe_rna_solver.py
   ```

   Both now run on curated sequences (`test_sequences.py`) confirmed to
   actually fold under ViennaRNA -- see `select_test_sequences.py` :
   for additional/different validation sequences.

3. **Analysis & plots**
   ```bash
   python compare_to_vienna.py
   python batch_accuracy.py
   python generate_final_results.py
   python final_summary.py
   python noise_simulation.py
   python plot_noise.py
   python scaling_analysis_real.py
   python plot_scaling_real.py
   python cvar_vqe_scaling_analysis.py
   python plot_cvar_vqe_scaling.py
   python runtime_analysis.py
   python compare_qubo_encodings.py       # optional: second (one-hot) QUBO encoding vs. current one
   python plot_qubo_encoding_comparison.py
   ```

4. **Broader benchmark (12 sequences, multiple trials, mean ± std)**
   ```bash
   python benchmark_sequences.py     # sanity-check the 12 curated sequences fold
   python statistical_benchmark.py   # full sweep -- can take a while
   ```
   QAOA gets slow past ~11 qubits, so
   running in chunks, e.g.:

   ```bash
   python statistical_benchmark.py --methods vqe            # cheap, all 12 sequences
   python statistical_benchmark.py --methods qaoa --append  # only sequences <= QAOA_MAX_QUBITS
   ```

---

## Script reference

| Script | Purpose |
|---|---|
| `generate_sequences.py` | *(legacy, unused by the current pipeline)* Generates 50 random 20-nt RNA sequences → `data/sequences.csv` |
| `test_vienna.py` | Single-sequence ViennaRNA sanity check |
| `structure_generator.py` | Prints dot-bracket structure + MFE for one example sequence |
| `pair_finder.py` / `rna_to_qubo.py` | Finds valid Watson-Crick candidate base pairs |
| `conflict_detector.py` | Flags candidate pairs that share a nucleotide (can't coexist) |
| `simple_qubo.py` / `qubo_energy.py` | Toy 5-variable QUBO matrix + energy function |
| `random_solver.py` | Brute-force/random search over the toy QUBO |
| `quantum_test.py` | Minimal 1-qubit PennyLane circuit |
| `rna_to_qubo_full.py` | Full QUBO: wobble pairs, min loop size, overlap + crossing penalties, brute-force solver |
| `select_test_sequences.py` | Reproducible search for short sequences that actually fold under ViennaRNA and fit a brute-forceable qubit budget |
| `test_sequences.py` | Curated 10 nt / 12 nt sequences (found via the search above) used consistently across the whole pipeline as the standard validation cases |
| `qaoa_rna_solver.py` | QAOA wired to the real RNA QUBO (method 1 of 2) — loops over all 12 `benchmark_sequences.py` sequences |
| `cvar_vqe_rna_solver.py` | CVaR-VQE (two-local ansatz) wired to the real RNA QUBO (method 2 of 2) — loops over all 12 `benchmark_sequences.py` sequences, decodes to dot-bracket and compares to ViennaRNA MFE structure |
| `structure_metrics.py` | Real ViennaRNA energy evaluation (`eval_structure`, same units as MFE) + base-pair precision/recall/F1/distance and Hamming distance between two structures |
| `benchmark_sequences.py` | 12 curated sequences (8–14 nt, varied GC content) confirmed to fold and sized for local quantum simulation — the standard sequence set used across the whole pipeline (superseded `test_sequences.py`'s 10 nt / 12 nt pair, though `test_sequences.py` is still kept for `select_test_sequences.py`-style ad hoc searches) |
| `statistical_benchmark.py` | Runs QAOA + CVaR-VQE across all 12 benchmark sequences with multiple independent trials each → `results/statistical_benchmark.csv` (every trial) and `results/statistical_benchmark_summary.csv` (mean ± std per sequence/method) |
| `scaling_analysis_real.py` / `plot_scaling_real.py` | Measured (not estimated) QAOA resource scaling: qubits, circuit depth, gates, forward runtime |
| `cvar_vqe_scaling_analysis.py` / `plot_cvar_vqe_scaling.py` | Measured CVaR-VQE resource scaling, plus energy gap + success rate up to `MAX_OPT_QUBITS` |
| `noise_simulation.py` | Real PennyLane depolarizing noise model (`default.mixed`) — QAOA solution quality vs. noise level, looped over all 12 benchmark sequences → `results/noise_analysis.csv` |
| `plot_noise.py` | Plots `results/noise_analysis.csv` (one line per sequence) → `results/noise_plot.png` |
| `plot_scaling.py` | Plots the superseded formula-based scaling numbers → `results/scaling_plot.png` (removed from this repo as stale/superseded; rerun only if you specifically want to reference the old estimate) |
| `runtime_analysis.py` | Times ViennaRNA folding for a few sequences |
| `compare_to_vienna.py` | Runs QAOA + CVaR-VQE on all 12 curated benchmark sequences and compares to ViennaRNA using real thermodynamic energy and base-pair metrics |
| `batch_accuracy.py` | Same comparison as above, looped over all 12 curated `benchmark_sequences.py` cases, one row per sequence → `results/batch_accuracy.csv` |
| `generate_final_results.py` | Runs both methods on all 12 curated benchmark sequences and writes a full real-metrics summary table, one row per sequence → `results/final_results.csv` |
| `final_summary.py` | Prints a human-readable per-sequence + aggregate project summary from `results/final_results.csv` and the scaling tables |
| `official_example_benchmark.py` | Classical-only ViennaRNA MFE benchmark on the 44 nt example sequence given directly in the challenge brief (Task 2) — not run through QAOA/CVaR-VQE, since it needs ~313 qubits, far past this project's measured feasibility ceiling |
| `rna_to_qubo_onehot.py` | Second QUBO encoding (optional advanced task): one-hot-per-position pairing variables (2 directed qubits per candidate pair + one-hot/consistency penalties) instead of `rna_to_qubo_full.py`'s pair-indicator variables |
| `compare_qubo_encodings.py` | Measures qubit count and penalty-term count for both encodings across all 14 curated test/benchmark sequences, and brute-force-validates that both encodings agree with each other and with ViennaRNA MFE wherever feasible → `results/qubo_encoding_comparison.csv` |
| `plot_qubo_encoding_comparison.py` | Plots `results/qubo_encoding_comparison.csv` (grouped bar chart, qubits per sequence per encoding) → `results/qubo_encoding_comparison_plot.png` |

