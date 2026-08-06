# WISER × Moderna Quantum Challenge — RNA Folding

Exploring RNA secondary structure prediction (minimum free energy folding) using
classical benchmarking (ViennaRNA) and quantum/quantum-inspired methods
(QUBO formulation + QAOA/CVaR-VQE in PennyLane).

Full write-up: [`Moderna_RNA_Folding_Quantum_Report.docx`](./Moderna_RNA_Folding_Quantum_Report.docx). This README summarizes the same project for anyone browsing the repository directly.

---

## The Challenge

Moderna's mRNA medicines depend on predictable secondary structure — it affects molecular stability, translation efficiency, and manufacturability. Classical tools like ViennaRNA find the exact Minimum Free Energy (MFE) structure via dynamic programming, but that approach only works for the pseudoknot-free case; richer structural constraints break the recursion and make exact folding NP-hard in general. The WISER × Moderna brief asks participants to formulate RNA folding as a quantum or quantum-inspired optimization problem, reproduce known MFE benchmark structures for small sequences, and analyze how quantum resource requirements scale with sequence length.

## Our Approach

Each candidate base pair becomes a binary decision variable in a QUBO (Quadratic Unconstrained Binary Optimization) objective: valid Watson-Crick and wobble (G-U) pairs are rewarded, a minimum hairpin-loop constraint is enforced at candidate-generation time, and overlapping or crossing (pseudoknotted) pairs are penalized so the optimizer never selects them. Two independent quantum methods — QAOA and CVaR-VQE, both implemented in PennyLane — solve this QUBO and are checked against a classical ViennaRNA ground truth, not just against each other. Running two methods with different objectives (plain expectation value vs. CVaR) against the same QUBO gives a built-in cross-check: agreement is evidence a solution reflects the underlying problem rather than one circuit's quirk.

**Alternatives considered:** the brief also lists quantum annealing, Grover-style search, and tensor-network-inspired methods as options. Annealing was set aside because the project's toolchain and available hardware access centered on gate-model simulators (the QUBO itself is annealer-agnostic, so it would carry over if hardware access becomes available). Grover-style search doesn't naturally fit a weighted quadratic cost landscape — it's built for locating a single marked item via a known oracle, not minimizing among many low-lying near-optimal states. Tensor-network methods are suited to compactly representing highly entangled states, but this problem's core difficulty is combinatorial (which base pairs to select), not a state-representation problem — adapting a tensor-network optimizer to this QUBO was judged out of scope for the timeline.

## Methods & Tools

- **Classical baseline:** ViennaRNA (`RNA.fold`, `RNA.fold_compound().eval_structure()`) for MFE structures and real thermodynamic energy evaluation
- **QUBO formulation:** pair-indicator encoding (one qubit per candidate base pair), brute-force validated against ViennaRNA; a second one-hot encoding was also implemented and cross-validated as an optional advanced task
- **Quantum solvers:** QAOA and CVaR-VQE (two-local ansatz), both in PennyLane, both optimized with derivative-free COBYLA
- **Evaluation:** base-pair precision/recall/F1 and real thermodynamic energy gap against ViennaRNA (not just internal QUBO score)
- **Robustness testing:** multi-trial statistical benchmarking (mean ± std), a real PennyLane depolarizing-noise model, and measured (not estimated) qubit/circuit-depth/runtime scaling up to 50 nt

## Results

Across a curated 12-sequence benchmark (8–14 nt, 33%–100% GC content):

| Metric | QAOA | CVaR-VQE |
|---|---|---|
| Mean success rate (multi-trial) | ≈48% | 100% |
| Exact ViennaRNA MFE match | 4 of 12 best-method wins | 8 of 12 best-method wins |
| Practical qubit ceiling | ~11–13 qubits | ~13–14 qubits |
| Mean optimization runtime | ~9.8 s | ~2.9 s |

CVaR-VQE reproduces its own QUBO's exact optimum on every sequence and trial tested, and that optimum matches ViennaRNA's real MFE structure exactly (F1 = 1.0, energy gap = 0.0 kcal/mol) on 11 of 12 sequences. Under depolarizing noise, QAOA's mean success probability falls from 100% (ideal) to 74% at moderate noise (p = 0.05) and 55% at high noise (p = 0.10). Resource scaling was measured (not formula-estimated) up to 50 nt (378 qubits); full quality-checked optimization is only run up to ~13–14 qubits, where statevector simulation is still tractable.

## Limitations & Next Steps

- **QUBO objective is imperfect:** a simple pair-count reward diverges from real ViennaRNA thermodynamics on a minority of sequences (2 of 12) — a limitation of the objective, not the optimizers. *Next step: a length- and stack-dependent reward closer to real nearest-neighbor thermodynamics.*
- **QAOA reliability gap:** ≈48% vs. 100% success rate under this project's optimizer/step budget. *Next step: sweep circuit depth (p), step budget, and restart count to see if the gap closes.*
- **Pseudoknots excluded by design**, matching ViennaRNA's own DP assumption but limiting applicability to the pseudoknot-free subset of real structures. *Next step: relax or reweight the crossing-pair penalty to allow scoring pseudoknotted structures.*
- **Simulator-only:** no results here speak to physical NISQ hardware behavior. *Next step: run the smallest validated cases (8–11 qubits) on real hardware.*
- **Qubit ceiling:** full optimization tops out around 13–14 qubits. *Next step: windowed/hierarchical QUBO decomposition to push past this without waiting on hardware qubit counts.*

## Team Contributions

- **Sri Jahnavi Chinthalapudi** — QUBO formulation and both encodings (pair-indicator and one-hot), both quantum solvers (QAOA and CVaR-VQE), the statistical, scaling, and noise-robustness analyses, structural/energy evaluation metrics, and final code clean-up across the repository.
- **Yagna Priya Gummadi** — Classical ViennaRNA benchmark implementation, and the final report and accompanying presentation.
- **Sreeneha Narayanam** — Curated benchmark sequence selection and the background review (`Rna_basics.docx`).

---

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
