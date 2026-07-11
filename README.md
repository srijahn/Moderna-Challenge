# WISER × Moderna Quantum Challenge — RNA Folding

Exploring RNA secondary structure prediction (minimum free energy folding) using
classical benchmarking (ViennaRNA) and quantum/quantum-inspired methods
(QUBO formulation + QAOA in PennyLane).

Submission deadline: **August 7, 2026**. Full execution plan:
`WISER_Moderna_4Week_Plan.docx`.

---

## Project status

- ✅ Classical benchmark pipeline (ViennaRNA) is real and working.
- ⚠️ QUBO formulation currently only detects valid base pairs (Watson-Crick);
  wobble pairs, crossing/overlap penalties, and loop-length constraints are
  **not yet implemented**.
- ⚠️ QAOA scripts (`qaoa_demo.py`, `qaoa_optimizer.py`) are small 1–2 qubit
  toy circuits, not yet connected to the real RNA QUBO.
- ⚠️ `noise_simulation.py`, `resource_estimator.py`, `scaling_analysis.py`,
  `generate_final_results.py`, `final_summary.py`, `compare_to_vienna.py`,
  and `batch_accuracy.py` currently use **hardcoded/placeholder numbers**
  (e.g. quantum energy fixed at `-3.0`, noise levels are arbitrary
  multipliers) rather than results computed from an actual quantum run.
  Treat any numbers/plots from these scripts as illustrative only until
  they're wired up to real QAOA output.
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

- [ ] Add wobble (G-U) pairing and minimum hairpin loop constraint (≥3 nt) to the QUBO
- [ ] Add crossing-pair and overlapping-pair penalty terms to the QUBO objective
- [ ] Brute-force validate the QUBO on the 10-nt sequence against ViennaRNA MFE
- [ ] Wire `qaoa_optimizer.py` to the real RNA QUBO instead of the 2-qubit toy problem
- [ ] Implement a second independent method (D-Wave Ocean SDK or CVaR-VQE)
- [ ] Replace placeholder "quantum energy" (`-3.0`) with actual QAOA output in
      `compare_to_vienna.py`, `batch_accuracy.py`, `generate_final_results.py`, `final_summary.py`
- [ ] Replace the arbitrary noise multipliers in `noise_simulation.py` with a real
      Qiskit Aer / PennyLane noise model
- [ ] Replace formula-based estimates in `resource_estimator.py` / `scaling_analysis.py`
      with numbers measured from actual circuit builds
- [ ] Write final report and presentation deck (Week 4)