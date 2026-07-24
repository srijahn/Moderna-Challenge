import pandas as pd
import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe

# --- Original approach ------------------------------------------------------
# Hand-written example values, not computed from any actual run:
# data = {
#     "Metric": [
#         "ViennaRNA MFE", "Quantum Candidate Energy", "Energy Gap",
#         "Approx Accuracy (%)", "QAOA Success Probability",
#         "Best Gamma", "Best Beta",
#     ],
#     "Value": [-4.0, -3.0, 1.0, 75.0, 0.4796, 0.7, 0.3],
# }

# Hardcoded short test sequence (the same 10 nt sequence used as the Week 2
# validation case in qaoa_rna_solver.py / cvar_vqe_rna_solver.py) so this
# script can actually run both real quantum solvers locally.
sequence = "GGUGCCGAAC"

structure, mfe = RNA.fold(sequence)

candidates = get_candidate_pairs(sequence)
Q = build_qubo(candidates)

print("Running QAOA...")
qaoa_bits, qaoa_prob, _ = run_qaoa(Q, n_layers=3, steps=150, n_restarts=2, top_k=15)
qaoa_energy = energy(qaoa_bits, Q)

print("Running CVaR-VQE...")
vqe_bits, vqe_prob, depth = run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200)
vqe_energy = energy(vqe_bits, Q)

best_x, best_e = brute_force_solve(Q)

# Best (lowest, i.e. most negative) QUBO energy of the two methods
quantum_energy = min(qaoa_energy, vqe_energy)

# NOTE: QUBO energy and ViennaRNA MFE (kcal/mol) are different scales --
# see the caveat in compare_to_vienna.py.
energy_gap = abs(mfe - quantum_energy)
accuracy = (abs(quantum_energy) / abs(mfe)) * 100 if mfe != 0 else 0.0

data = {
    "Metric": [
        "Sequence",
        "ViennaRNA MFE",
        "Brute-force QUBO Optimum",
        "QAOA QUBO Energy",
        "QAOA Success Probability",
        "CVaR-VQE QUBO Energy",
        "CVaR-VQE Success Probability",
        "CVaR-VQE Circuit Depth",
        "Best Quantum QUBO Energy",
        "Energy Gap (vs ViennaRNA MFE)",
        "Approx Accuracy (%)",
    ],
    "Value": [
        sequence,
        mfe,
        best_e,
        qaoa_energy,
        round(qaoa_prob, 4),
        vqe_energy,
        round(vqe_prob, 4),
        depth,
        quantum_energy,
        round(energy_gap, 3),
        round(accuracy, 2),
    ],
}

df = pd.DataFrame(data)

df.to_csv("results/final_results.csv", index=False)

print("\n========== FINAL RESULTS ==========\n")
print(df)

print("\nSaved:")
print("results/final_results.csv")
