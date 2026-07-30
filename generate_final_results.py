import pandas as pd
import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket
from structure_metrics import real_energy, base_pair_metrics
from benchmark_sequences import BENCHMARK_SEQUENCES

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

# --- Superseded: single fixed 10 nt sequence, one Metric/Value column pair -
# Was TEST_SEQUENCE_10NT only, written as a single Metric/Value table. Now
# loops over the full 8-sequence curated benchmark set (benchmark_sequences.py,
# 8-14 nt, 35.7%-100% GC content, all confirmed by ViennaRNA to fold) and
# writes one row per sequence instead, so results/final_results.csv covers
# the whole benchmark, not one anecdote. final_summary.py reads this same
# wide format.

rows = []

for label, sequence, _, _, _ in BENCHMARK_SEQUENCES:
    print(f"\n=== {label}: {sequence} ===")

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

    # Best (lowest, i.e. most negative) *QUBO* energy of the two methods --
    # internal solver comparison only, used to pick which structure to report
    # as "the quantum candidate" below.
    if qaoa_energy <= vqe_energy:
        quantum_energy, quantum_structure, quantum_method = qaoa_energy, qaoa_structure, "QAOA"
    else:
        quantum_energy, quantum_structure, quantum_method = vqe_energy, vqe_structure, "CVaR-VQE"

    # Real ViennaRNA-based comparison (Task 3): same units, same
    # thermodynamic model as the MFE, unlike the QUBO energy above which is
    # an arbitrary reward/penalty scheme used only to steer the quantum
    # optimizers.
    quantum_real_energy = real_energy(sequence, quantum_structure)
    real_energy_gap = quantum_real_energy - mfe
    metrics = base_pair_metrics(quantum_structure, vienna_structure)

    rows.append({
        "Label": label,
        "Sequence": sequence,
        "ViennaRNA MFE Structure": vienna_structure,
        "ViennaRNA MFE Energy (kcal/mol)": mfe,
        "Brute-force QUBO Optimum (internal)": best_e,
        "QAOA QUBO Energy (internal)": qaoa_energy,
        "QAOA Success Probability": round(qaoa_prob, 4),
        "CVaR-VQE QUBO Energy (internal)": vqe_energy,
        "CVaR-VQE Success Probability": round(vqe_prob, 4),
        "CVaR-VQE Circuit Depth": depth,
        "Best Quantum Method": quantum_method,
        "Best Quantum Structure": quantum_structure,
        "Best Quantum Real Energy (kcal/mol)": round(quantum_real_energy, 2),
        "Real Energy Gap vs. MFE (kcal/mol)": round(real_energy_gap, 2),
        "Base-Pair Precision": round(metrics["precision"], 3),
        "Base-Pair Recall": round(metrics["recall"], 3),
        "Base-Pair F1": round(metrics["f1"], 3),
        "Base-Pair Distance": metrics["base_pair_distance"],
        "Hamming Distance": metrics["hamming_distance"],
    })

df = pd.DataFrame(rows)

df.to_csv("results/final_results.csv", index=False)

print("\n========== FINAL RESULTS (all 8 benchmark sequences) ==========\n")
print(df.to_string(index=False))

print("\nSaved:")
print("results/final_results.csv")
