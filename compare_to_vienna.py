import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket
from structure_metrics import real_energy, base_pair_metrics
from test_sequences import TEST_SEQUENCE_10NT

# --- Original placeholder sequence -----------------------------------------
# This 20 nt sequence produces 50+ candidate base pairs (see
# findings_week1.md), i.e. 50+ qubits -- far past what local statevector
# simulation of QAOA/CVaR-VQE can handle. Left here for reference only.
# sequence = "GGUGCCGAACAGUAGCACUC"

# --- Original 10 nt validation sequence ("GGUGCCGAAC") -----------------
# Turned out to have a fully-unpaired real ViennaRNA MFE structure (i.e.
# it doesn't actually fold) -- see test_sequences.py for why this was
# replaced with a curated sequence that ViennaRNA confirms folds.
sequence = TEST_SEQUENCE_10NT

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

# Best (lowest, i.e. most negative) *QUBO* energy of the two independent
# quantum methods. QUBO energy is only used internally to pick a winner
# between the two solvers -- it is NOT compared to ViennaRNA's MFE anymore
# (see the real-energy section below for that).
if qaoa_energy <= vqe_energy:
    quantum_energy, quantum_structure, quantum_method = qaoa_energy, qaoa_structure, "QAOA"
else:
    quantum_energy, quantum_structure, quantum_method = vqe_energy, vqe_structure, "CVaR-VQE"

# ---------------------------------------------------------------------------
# Real ViennaRNA-based comparison (Task 3 of the challenge PDF): evaluate
# each candidate structure's actual thermodynamic free energy with
# RNA.fold_compound(sequence).eval_structure(), so it's directly comparable
# (same units, same model) to the MFE -- unlike the internal QUBO energy
# above, which is an arbitrary reward/penalty scheme used only to steer the
# quantum optimizers and is not on the same scale as kcal/mol.
# ---------------------------------------------------------------------------
qaoa_real_energy = real_energy(sequence, qaoa_structure)
vqe_real_energy = real_energy(sequence, vqe_structure)
quantum_real_energy = real_energy(sequence, quantum_structure)

qaoa_metrics = base_pair_metrics(qaoa_structure, structure)
vqe_metrics = base_pair_metrics(vqe_structure, structure)
quantum_metrics = base_pair_metrics(quantum_structure, structure)

real_energy_gap = quantum_real_energy - mfe  # signed: 0 = exact MFE match, >0 = worse (less stable)

print()
print("========== RNA Folding Comparison ==========")

print("\nSequence:")
print(sequence)

print("\nViennaRNA MFE Structure:")
print(structure)

print("\nViennaRNA MFE Energy (kcal/mol):")
print(mfe)

print(f"\n--- QAOA (internal QUBO energy={qaoa_energy}, p={qaoa_prob:.3f}) ---")
print("Structure:            ", qaoa_structure)
print(f"Real energy (kcal/mol):", round(qaoa_real_energy, 2))
print("Base-pair metrics:     ", qaoa_metrics)

print(f"\n--- CVaR-VQE (internal QUBO energy={vqe_energy}, p={vqe_prob:.3f}, depth={depth}) ---")
print("Structure:            ", vqe_structure)
print(f"Real energy (kcal/mol):", round(vqe_real_energy, 2))
print("Base-pair metrics:     ", vqe_metrics)

print(f"\n--- Best quantum candidate: {quantum_method} ---")
print("Structure:                     ", quantum_structure)
print(f"Real energy (kcal/mol):        ", round(quantum_real_energy, 2))
print(f"Real energy gap vs. MFE (kcal/mol):", round(real_energy_gap, 2))
print("Base-pair precision / recall / F1:",
      round(quantum_metrics["precision"], 3),
      round(quantum_metrics["recall"], 3),
      round(quantum_metrics["f1"], 3))
print("Base-pair distance (0 = exact structural match):", quantum_metrics["base_pair_distance"])
print("Hamming distance (dot-bracket, 0 = exact string match):", quantum_metrics["hamming_distance"])

print("\n============================================")
