# WISER × Moderna Quantum Challenge — RNA Folding

Exploring RNA secondary structure prediction (minimum free energy folding) using
classical benchmarking (ViennaRNA) and quantum/quantum-inspired methods
(QUBO formulation + QAOA in PennyLane).

Submission deadline: **August 7, 2026**. Full execution plan:
`WISER_Moderna_4Week_Plan.docx`.

---

## Project status

- ✅ Classical benchmark pipeline (ViennaRNA) is real and working.
- ✅ Full QUBO formulation (`rna_to_qubo_full.py`) implements wobble (G-U)
  pairs, minimum hairpin loop size, and overlap/crossing penalties, and is
  brute-force validated against ViennaRNA MFE on short sequences.
- ⚠️ `qaoa_demo.py` / `qaoa_optimizer.py` are still small 1–2 qubit toy
  circuits, not connected to the real RNA QUBO -- kept only as an earlier
  step-by-step reference. The real pipeline lives in `qaoa_rna_solver.py`.
- ✅ Two independent quantum methods run end-to-end on the real QUBO:
  QAOA (`qaoa_rna_solver.py`) and CVaR-VQE (`cvar_vqe_rna_solver.py`), the
  Week 2 Coder A/B split. Both validated against brute force on 10 nt and
  15 nt test sequences.
- ✅ `compare_to_vienna.py`, `batch_accuracy.py`, `generate_final_results.py`,
  `final_summary.py`, and `noise_simulation.py` are now wired to real
  QAOA/CVaR-VQE output and a real PennyLane depolarizing noise model --
  no more hardcoded placeholder numbers.
- ⚠️ `resource_estimator.py` / `scaling_analysis.py` still contain the old
  formula-based estimates internally. They now print a SUPERSEDED notice
  pointing to the real measured versions (`scaling_analysis_real.py` /
  `cvar_vqe_scaling_analysis.py`) and are kept only for reference -- don't
  cite their numbers.
- ✅ Resource scaling (qubits, circuit depth, gates, runtime) is measured
  (not estimated) for both methods up to 50 nt via `scaling_analysis_real.py`
  and `cvar_vqe_scaling_analysis.py`.
- ❌ Final report and presentation deck (Week 4 deliverables) not started.

---

## Requirements

### Software
- **Python 3.9–3.11** (ViennaRNA's Python bindings can be picky about newer
  versions — 3.11 or below is safest)
- pip (or conda)

### Python packages

| Package      | Used for                                  |
|--------------|--------------------------------------------|
| `ViennaRNA`  | Classical MFE structure/energy benchmark (`import RNA`) |
| `pennylane`  | QAOA / quantum circuit simulation          |
| `numpy`      | QUBO matrices, numerical work              |
| `pandas`     | Reading/writing CSV data and results       |
| `matplotlib` | Scaling and noise plots                    |

Install everything with:

```bash
pip install ViennaRNA pennylane numpy pandas matplotlib
```

or use the included `requirements.txt`:

```bash
pip install -r requirements.txt
```

> **Windows note:** `ViennaRNA` is a C-extension package and doesn't always
> have a prebuilt wheel for Windows. If `pip install ViennaRNA` fails, try:
> - `conda install -c bioconda viennarna`, or
> - running the project inside WSL (Windows Subsystem for Linux)

---

## Folder structure

The scripts expect this layout:

```
Moderna/
├── data/
│   └── sequences.csv
├── results/
│   ├── vienna_results.csv
│   ├── noise_analysis.csv
│   ├── final_results.csv
│   ├── noise_plot.png
│   └── scaling_plot.png
├── *.py                     (all scripts)
├── RNA_Basics.md
├── findings_week1.md
└── WISER_Moderna_4Week_Plan.docx
```

Create the folders and move the CSVs before running anything:

**macOS/Linux:**
```bash
mkdir data results
mv sequences.csv data/
mv vienna_results.csv noise_analysis.csv final_results.csv results/
```

**Windows (PowerShell):**
```powershell
mkdir data
mkdir results
Move-Item sequences.csv data\
Move-Item vienna_results.csv results\
Move-Item noise_analysis.csv results\
Move-Item final_results.csv results\
```

---

## How to run

Sanity-check ViennaRNA first:
```bash
python test_vienna.py
```
If it prints a sequence, dot-bracket structure, and MFE energy, you're set.

Then, in order:

1. **Generate data** (optional — `sequences.csv` already provided)
   ```bash
   python generate_sequences.py
   ```

2. **Classical benchmark** (writes `results/vienna_results.csv`)
   ```bash
   python benchmark.py
   ```

3. **Pairing / QUBO basics** (print output only, no files written)
   ```bash
   python pair_finder.py
   python conflict_detector.py
   python structure_generator.py
   python simple_qubo.py
   python qubo_energy.py
   python random_solver.py
   ```

4. **Quantum circuit demos**
   ```bash
   python quantum_test.py
   python qaoa_demo.py
   python qaoa_optimizer.py
   ```

   **Real RNA QUBO solvers (two independent methods, Week 2):**
   ```bash
   python qaoa_rna_solver.py
   python cvar_vqe_rna_solver.py
   ```

5. **Analysis & plots** (run after step 2, since these read `vienna_results.csv`)
   ```bash
   python batch_accuracy.py
   python noise_simulation.py
   python plot_noise.py
   python plot_scaling.py
   python scaling_analysis.py
   python resource_estimator.py
   python runtime_analysis.py
   python compare_to_vienna.py
   python generate_final_results.py
   python final_summary.py
   ```

---

## Script reference

| Script | Purpose |
|---|---|
| `generate_sequences.py` | Generates 50 random 20-nt RNA sequences → `data/sequences.csv` |
| `benchmark.py` | Runs ViennaRNA MFE fold on all sequences → `results/vienna_results.csv` |
| `test_vienna.py` | Single-sequence ViennaRNA sanity check |
| `structure_generator.py` | Prints dot-bracket structure + MFE for one example sequence |
| `pair_finder.py` / `rna_to_qubo.py` | Finds valid Watson-Crick candidate base pairs |
| `conflict_detector.py` | Flags candidate pairs that share a nucleotide (can't coexist) |
| `simple_qubo.py` / `qubo_energy.py` | Toy 5-variable QUBO matrix + energy function |
| `random_solver.py` | Brute-force/random search over the toy QUBO |
| `quantum_test.py` | Minimal 1-qubit PennyLane circuit |
| `qaoa_demo.py` / `qaoa_optimizer.py` | 2-qubit toy QAOA circuit in PennyLane |
| `rna_to_qubo_full.py` | Full QUBO: wobble pairs, min loop size, overlap + crossing penalties, brute-force solver |
| `qaoa_rna_solver.py` | QAOA wired to the real RNA QUBO (method 1 of 2) |
| `cvar_vqe_rna_solver.py` | CVaR-VQE (two-local ansatz) wired to the real RNA QUBO (method 2 of 2) — decodes to dot-bracket and compares to ViennaRNA MFE |
| `scaling_analysis_real.py` / `plot_scaling_real.py` | Measured (not estimated) QAOA resource scaling: qubits, circuit depth, gates, forward runtime |
| `cvar_vqe_scaling_analysis.py` / `plot_cvar_vqe_scaling.py` | Measured CVaR-VQE resource scaling, plus energy gap + success rate up to `MAX_OPT_QUBITS` |
| `noise_simulation.py` | Illustrative noise-level vs. success-probability table (placeholder math) |
| `plot_noise.py` | Plots `results/noise_analysis.csv` → `results/noise_plot.png` |
| `resource_estimator.py` / `scaling_analysis.py` | Estimated qubit/variable counts vs. sequence length (formula-based, not measured) |
| `plot_scaling.py` | Plots hardcoded scaling numbers → `results/scaling_plot.png` |
| `runtime_analysis.py` | Times ViennaRNA folding for a few sequences |
| `compare_to_vienna.py` / `batch_accuracy.py` | Compares ViennaRNA MFE to a placeholder "quantum energy" |
| `generate_final_results.py` | Writes a placeholder summary table → `results/final_results.csv` |
| `final_summary.py` | Prints a placeholder project summary to the terminal |

---

## Known issues / TODO

- [x] Add wobble (G-U) pairing and minimum hairpin loop constraint (≥3 nt) to the QUBO (`rna_to_qubo_full.py`)
- [x] Add crossing-pair and overlapping-pair penalty terms to the QUBO objective (`rna_to_qubo_full.py`)
- [x] Brute-force validate the QUBO on sequences up to 20 candidate pairs (`brute_force_solve()`,
      used as ground truth by both `qaoa_rna_solver.py` and `cvar_vqe_rna_solver.py`)
- [x] Wire QAOA to the real RNA QUBO instead of the 2-qubit toy problem (`qaoa_rna_solver.py`)
- [x] Implement a second independent method (`cvar_vqe_rna_solver.py`, CVaR-VQE with a two-local ansatz)
- [x] Record qubit count / circuit depth / runtime for CVaR-VQE across the full sequence-length ladder
      (`cvar_vqe_scaling_analysis.py` + `plot_cvar_vqe_scaling.py`; energy gap and success rate are also
      measured up to `MAX_OPT_QUBITS`, beyond which only resource counts are reported -- see the script's
      docstring for why)
- [x] Replace placeholder "quantum energy" (`-3.0`) with actual QAOA/CVaR-VQE output in
      `compare_to_vienna.py`, `batch_accuracy.py`, `generate_final_results.py`, `final_summary.py`
- [x] Replace the arbitrary noise multipliers in `noise_simulation.py` with a real
      PennyLane depolarizing noise model (`default.mixed` device)
- [x] Replace formula-based estimates in `resource_estimator.py` / `scaling_analysis.py`
      with numbers measured from actual circuit builds (`scaling_analysis_real.py` /
      `cvar_vqe_scaling_analysis.py`). The old scripts are kept for reference but now
      print a SUPERSEDED notice on run -- don't cite their numbers in the report.
- [ ] Write final report and presentation deck (Week 4)