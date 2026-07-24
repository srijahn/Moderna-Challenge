import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket

# --- Original placeholder sequence -----------------------------------------
# This 20 nt sequence produces 50+ candidate base pairs (see
# findings_week1.md), i.e. 50+ qubits -- far past what local statevector
# simulation of QAOA/CVaR-VQE can handle. Left here for reference only.
# sequence = "GGUGCCGAACAGUAGCACUC"

# Hardcoded short test sequence (the same 10 nt sequence used as the Week 2
# validation case in qaoa_rna_solver.py / cvar_vqe_rna_solver.py) so this
# script can actually run both real quantum solvers locally.
sequence = "GGUGCCGAAC"

# ViennaRNA reference
structure, mfe = RNA.fold(sequence)

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

# Best (lowest, i.e. most negative) QUBO energy of the two independent
# quantum methods
if qaoa_energy <= vqe_energy:
    quantum_energy, quantum_structure, quantum_method = qaoa_energy, qaoa_structure, "QAOA"
else:
    quantum_energy, quantum_structure, quantum_method = vqe_energy, vqe_structure, "CVaR-VQE"

# NOTE: the QUBO energy (from build_qubo()'s pair-count reward/penalty
# scheme) and ViennaRNA's MFE (kcal/mol) are on different scales, so this
# gap/accuracy is a rough proxy for "did the quantum solver find a
# good structure", not a true energy comparison. See RNA_Basics.md.
energy_gap = abs(mfe - quantum_energy)
accuracy = (abs(quantum_energy) / abs(mfe)) * 100 if mfe != 0 else 0.0

print()
print("========== RNA Folding Comparison ==========")

print("\nSequence:")
print(sequence)

print("\nViennaRNA Structure:")
print(structure)

print("\nViennaRNA Energy:")
print(mfe)

print(f"\nQAOA QUBO Energy (p={qaoa_prob:.3f}):")
print(qaoa_energy, " | structure:", qaoa_structure)

print(f"\nCVaR-VQE QUBO Energy (p={vqe_prob:.3f}, circuit depth={depth}):")
print(vqe_energy, " | structure:", vqe_structure)

print(f"\nBest Quantum Candidate Energy ({quantum_method}):")
print(quantum_energy, " | structure:", quantum_structure)

print("\nEnergy Gap (vs ViennaRNA MFE -- different scales, rough proxy only):")
print(round(energy_gap, 3))

print("\nApprox Accuracy:")
print(round(accuracy, 2), "%")

print("\n============================================")
