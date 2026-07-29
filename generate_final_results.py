import pandas as pd
import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket
from structure_metrics import real_energy, base_pair_metrics
from test_sequences import TEST_SEQUENCE_10NT

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

# --- Original 10 nt validation sequence ("GGUGCCGAAC") -----------------
# Turned out to have a fully-unpaired real ViennaRNA MFE structure (i.e.
# it doesn't actually fold) -- see test_sequences.py for why this was
# replaced with a curated sequence that ViennaRNA confirms folds.
sequence = TEST_SEQUENCE_10NT

vienna_structure, mfe = RNA.fold(sequence)

candidates = get_candidate_pairs(sequence)
Q = build_qubo(candidates)

print("Running QAOA...")
qaoa_bits, qaoa_prob, _ = run_qaoa(Q, n_layers=3, steps=150, n_restarts=2, top_k=15)
qaoa_energy = energy(qaoa_bits, Q)
qaoa_pairs = [candidates[k] for k, b in enumerate(qaoa_bits) if b == 1]
qaoa_structure = pairs_to_dot_bracket(len(sequence), qaoa_pairs)

print("Running CVaR-VQE...")
vqe_bits, vqe_prob, depth = run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200)
vqe_energy = energy(vqe_bits, Q)
vqe_pairs = [candidates[k] for k, b in enumerate(vqe_bits) if b == 1]
vqe_structure = pairs_to_dot_bracket(len(sequence), vqe_pairs)

best_x, best_e = brute_force_solve(Q)
best_pairs = [candidates[k] for k, b in enumerate(best_x) if b == 1]
best_structure = pairs_to_dot_bracket(len(sequence), best_pairs)

# Best (lowest, i.e. most negative) *QUBO* energy of the two methods --
# internal solver comparison only, used to pick which structure to report
# as "the quantum candidate" below.
if qaoa_energy <= vqe_energy:
    quantum_energy, quantum_structure, quantum_method = qaoa_energy, qaoa_structure, "QAOA"
else:
    quantum_energy, quantum_structure, quantum_method = vqe_energy, vqe_structure, "CVaR-VQE"

# Real ViennaRNA-based comparison (Task 3): same units, same thermodynamic
# model as the MFE, unlike the QUBO energy above which is an arbitrary
# reward/penalty scheme used only to steer the quantum optimizers.
quantum_real_energy = real_energy(sequence, quantum_structure)
real_energy_gap = quantum_real_energy - mfe
metrics = base_pair_metrics(quantum_structure, vienna_structure)

data = {
    "Metric": [
        "Sequence",
        "ViennaRNA MFE Structure",
        "ViennaRNA MFE Energy (kcal/mol)",
        "Brute-force QUBO Optimum (internal)",
        "QAOA QUBO Energy (internal)",
        "QAOA Success Probability",
        "CVaR-VQE QUBO Energy (internal)",
        "CVaR-VQE Success Probability",
        "CVaR-VQE Circuit Depth",
        "Best Quantum Method",
        "Best Quantum Structure",
        "Best Quantum Real Energy (kcal/mol)",
        "Real Energy Gap vs. MFE (kcal/mol)",
        "Base-Pair Precision",
        "Base-Pair Recall",
        "Base-Pair F1",
        "Base-Pair Distance",
        "Hamming Distance",
    ],
    "Value": [
        sequence,
        vienna_structure,
        mfe,
        best_e,
        qaoa_energy,
        round(qaoa_prob, 4),
        vqe_energy,
        round(vqe_prob, 4),
        depth,
        quantum_method,
        quantum_structure,
        round(quantum_real_energy, 2),
        round(real_energy_gap, 2),
        round(metrics["precision"], 3),
        round(metrics["recall"], 3),
        round(metrics["f1"], 3),
        metrics["base_pair_distance"],
        metrics["hamming_distance"],
    ],
}

df = pd.DataFrame(data)

df.to_csv("results/final_results.csv", index=False)

print("\n========== FINAL RESULTS ==========\n")
print(df.to_string(index=False))

print("\nSaved:")
print("results/final_results.csv")
