# WISER × Moderna Quantum Challenge — RNA Folding

Exploring RNA secondary structure prediction (minimum free energy folding) using
classical benchmarking (ViennaRNA) and quantum/quantum-inspired methods
(QUBO formulation + QAOA/CVaR-VQE in PennyLane).

Submission deadline: **August 7, 2026**. Challenge brief:
`Moderna - WISER Quantum Challenge [SHARED].pdf`.

---

## Project status

- ✅ Classical benchmarking against ViennaRNA is real and working, wired
  through `structure_metrics.py` into `compare_to_vienna.py`,
  `batch_accuracy.py`, `generate_final_results.py`, and
  `statistical_benchmark.py`.
- ✅ Full QUBO formulation (`rna_to_qubo_full.py`) implements wobble (G-U)
  pairs, minimum hairpin loop size, and overlap/crossing penalties, and is
  brute-force validated against ViennaRNA MFE on short sequences.
- ✅ Two independent quantum methods run end-to-end on the real QUBO:
  QAOA (`qaoa_rna_solver.py`) and CVaR-VQE (`cvar_vqe_rna_solver.py`), the
  Week 2 Coder A/B split.
- ✅ **Real thermodynamic comparison, not just internal QUBO energy.**
  `structure_metrics.py` evaluates candidate structures with
  `RNA.fold_compound(sequence).eval_structure()` (same units, same model as
  the MFE -- Task 3 of the challenge) and reports base-pair-level
  precision/recall/F1, base-pair distance, and Hamming distance against the
  ViennaRNA reference structure. Wired into `compare_to_vienna.py`,
  `batch_accuracy.py`, `generate_final_results.py`, and `final_summary.py`.
  Previously these scripts compared the *internal QUBO energy* (an
  arbitrary pair-count reward, not real physics) directly against
  ViennaRNA's MFE (real kcal/mol) -- two different scales, not a real
  benchmark result.
- ✅ **Fixed a bug this uncovered: the original validation sequences didn't
  fold.** Both hardcoded test sequences ("GGUGCCGAAC", 10 nt, and a
  seed=15 random 15-mer) have a fully-unpaired real ViennaRNA MFE structure
  (0.00 kcal/mol) -- ViennaRNA says they don't fold at all. Every
  "validated against ViennaRNA" claim built on them was therefore only
  checking the quantum solvers against their own QUBO's brute-force
  optimum, never against a real, non-trivial MFE structure. Replaced with
  curated sequences (`test_sequences.py`, found via
  `select_test_sequences.py`) that ViennaRNA confirms fold. With real
  sequences, CVaR-VQE finds the **exact** ViennaRNA structure on the 10 nt
  case (F1 = 1.0, energy gap = 0.0 kcal/mol); on the 12 nt case it finds
  its own QUBO's exact optimum, but that optimum still doesn't match
  ViennaRNA's real MFE (F1 = 0.33) -- a genuine, reportable limitation of
  the simple pair-count QUBO objective, not a bug.
- ✅ **Benchmark breadth + statistical rigor.** `benchmark_sequences.py`
  provides 8 curated sequences (8–14 nt, 36%–100% GC content) confirmed to
  fold and sized for local simulation. `statistical_benchmark.py` runs
  both methods across all 8 with multiple independent trials each (5 for
  CVaR-VQE, 3 for QAOA) and reports mean ± std per sequence
  (`results/statistical_benchmark.csv` for every trial,
  `results/statistical_benchmark_summary.csv` for the per-sequence
  mean ± std) instead of a single-run anecdote. Headline numbers from the
  current run: CVaR-VQE hits its own
  QUBO optimum on 7/8 sequences (success rate 1.0) with mean F1 ≈ 0.6–1.0
  depending on sequence; QAOA (reduced step budget for this batch sweep,
  see the script's docstring) is markedly less reliable, with success
  rates as low as 0.0–0.67 and higher variance. QAOA is also only run up
  to 11 qubits in this sweep -- see the runtime finding below.
- ⚠️ **QAOA scaling/runtime constraint (new finding).** Full QAOA (150
  steps × 2 restarts, PennyLane backprop differentiation of a full
  statevector) took ~90s locally at 11 qubits but ~17 minutes at 17
  qubits -- backprop-based gradient descent scales badly with qubit count
  in a way CVaR-VQE (derivative-free COBYLA) doesn't. This is why the
  curated test/benchmark sequences are capped around 11–13 qubits, and
  it's worth a line in the final report/scaling discussion.
- ✅ Resource scaling (qubits, circuit depth, gates, runtime) is measured
  (not estimated) for both methods up to 50 nt via `scaling_analysis_real.py`
  and `cvar_vqe_scaling_analysis.py`. Solve-quality metrics (energy gap,
  success rate) are only measured up to their respective `MAX_OPT_QUBITS`
  cutoffs -- see each script's docstring for why (statevector simulation
  cost).
- ✅ **Repo cleanup (2026-07-29).** `results/` is the single source of
  truth for every generated file: `batch_accuracy.csv` and the
  `statistical_benchmark*.csv` files now live there alongside everything
  else, and the superseded `scaling_plot.png` / formula-based
  `scaling_analysis.py` stay dropped. The old `results/vienna_results.csv`
  was also removed -- it was output from a `benchmark.py` script that no
  longer exists in the repo and nothing else reads it; classical ViennaRNA
  comparison now happens inline (via `structure_metrics.py`) inside
  `compare_to_vienna.py`, `batch_accuracy.py`, `generate_final_results.py`,
  and `statistical_benchmark.py` instead of a separate bulk-benchmark step.
  `generate_sequences.py` / `data/sequences.csv` (the original 50 random
  20 nt sequences) are legacy from that same dropped step and aren't used
  by anything currently in the pipeline -- kept only as a reference for
  how the original random sequences were produced.
- ✅ Background review is `Rna_basics.docx`, covering all three parts Task 1
  asks for: the biological problem and MFE thermodynamics, the
  computational challenge (combinatorial growth of the fold space, why
  classical DP is fast but breaks down once pseudoknots are allowed), and
  the proposed quantum approach (QUBO encoding + QAOA/CVaR-VQE).
- ❌ Final report and presentation deck (Task 7 submission package) not
  started.

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
└── Moderna - WISER Quantum Challenge [SHARED].pdf
```

> **Note:** every script writes only into `results/` — it's the single
> source of truth for generated output, nothing is duplicated at the
> project root. `scaling_plot.png` (and its underlying formula-based
> `scaling_analysis.py` estimate) is dropped entirely rather than kept
> alongside `scaling_analysis_real.py` / `scaling_plot_real.png` — don't
> regenerate or cite it. `data/sequences.csv` (50 random 20 nt sequences)
> and the script that made it, `generate_sequences.py`, are no longer
> read by anything else in the pipeline — the original `benchmark.py`
> that consumed them is gone, and classical ViennaRNA comparison now
> happens per-curated-sequence via `structure_metrics.py` instead. Kept
> only for reference.

---

## How to run

Sanity-check ViennaRNA first:
```bash
python test_vienna.py
```
If it prints a sequence, dot-bracket structure, and MFE energy, you're set.

Then, in order:

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
   actually fold under ViennaRNA -- see `select_test_sequences.py` if you
   want to search for additional/different validation sequences.

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
   ```

4. **Broader benchmark (8 sequences, multiple trials, mean ± std)**
   ```bash
   python benchmark_sequences.py     # sanity-check the 8 curated sequences fold
   python statistical_benchmark.py   # full sweep -- can take a while
   ```
   QAOA gets slow past ~11 qubits, so
   running in chunks, e.g.:

   ```bash
   python statistical_benchmark.py --methods vqe            # cheap, all 8 sequences
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
| `qaoa_rna_solver.py` | QAOA wired to the real RNA QUBO (method 1 of 2) — loops over all 8 `benchmark_sequences.py` sequences |
| `cvar_vqe_rna_solver.py` | CVaR-VQE (two-local ansatz) wired to the real RNA QUBO (method 2 of 2) — loops over all 8 `benchmark_sequences.py` sequences, decodes to dot-bracket and compares to ViennaRNA MFE structure |
| `structure_metrics.py` | Real ViennaRNA energy evaluation (`eval_structure`, same units as MFE) + base-pair precision/recall/F1/distance and Hamming distance between two structures |
| `benchmark_sequences.py` | 8 curated sequences (8–14 nt, varied GC content) confirmed to fold and sized for local quantum simulation — the standard sequence set used across the whole pipeline (superseded `test_sequences.py`'s 10 nt / 12 nt pair, though `test_sequences.py` is still kept for `select_test_sequences.py`-style ad hoc searches) |
| `statistical_benchmark.py` | Runs QAOA + CVaR-VQE across all 8 benchmark sequences with multiple independent trials each → `results/statistical_benchmark.csv` (every trial) and `results/statistical_benchmark_summary.csv` (mean ± std per sequence/method) |
| `scaling_analysis_real.py` / `plot_scaling_real.py` | Measured (not estimated) QAOA resource scaling: qubits, circuit depth, gates, forward runtime |
| `cvar_vqe_scaling_analysis.py` / `plot_cvar_vqe_scaling.py` | Measured CVaR-VQE resource scaling, plus energy gap + success rate up to `MAX_OPT_QUBITS` |
| `noise_simulation.py` | Real PennyLane depolarizing noise model (`default.mixed`) — QAOA solution quality vs. noise level, looped over all 8 benchmark sequences → `results/noise_analysis.csv` |
| `plot_noise.py` | Plots `results/noise_analysis.csv` (one line per sequence) → `results/noise_plot.png` |
| `plot_scaling.py` | Plots the superseded formula-based scaling numbers → `results/scaling_plot.png` (removed from this repo as stale/superseded; rerun only if you specifically want to reference the old estimate) |
| `runtime_analysis.py` | Times ViennaRNA folding for a few sequences |
| `compare_to_vienna.py` | Runs QAOA + CVaR-VQE on all 8 curated benchmark sequences and compares to ViennaRNA using real thermodynamic energy and base-pair metrics |
| `batch_accuracy.py` | Same comparison as above, looped over all 8 curated `benchmark_sequences.py` cases, one row per sequence → `results/batch_accuracy.csv` |
| `generate_final_results.py` | Runs both methods on all 8 curated benchmark sequences and writes a full real-metrics summary table, one row per sequence → `results/final_results.csv` |
| `final_summary.py` | Prints a human-readable per-sequence + aggregate project summary from `results/final_results.csv` and the scaling tables |
| `official_example_benchmark.py` | Classical-only ViennaRNA MFE benchmark on the 44 nt example sequence given directly in the challenge brief (Task 2) — not run through QAOA/CVaR-VQE, since it needs ~313 qubits, far past this project's measured feasibility ceiling |

---

## Known issues / TODO

- [ ] Try a second QUBO encoding (e.g. one-hot per-position pairing variables instead of pair-indicator
      variables) and compare qubit count / constraint-enforcement tradeoffs against the current encoding
      (optional advanced task)
- [ ] Explicitly document the pseudoknot-exclusion rationale (crossing-pair penalty in
      `rna_to_qubo_full.py`) as a stated modeling assumption in the final report, not just in code
- [ ] Write final report and presentation deck
