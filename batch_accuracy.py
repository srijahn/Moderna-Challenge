import pandas as pd
import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket
from structure_metrics import real_energy, base_pair_metrics
from benchmark_sequences import BENCHMARK_SEQUENCES

print("\n========== BATCH ACCURACY ANALYSIS ==========\n")

rows = []

for index, (label, sequence, _, _, _) in enumerate(BENCHMARK_SEQUENCES):
    vienna_structure, mfe = RNA.fold(sequence)

    candidates = get_candidate_pairs(sequence)
    Q = build_qubo(candidates)

    qaoa_bits, _, _ = run_qaoa(Q, n_layers=3, steps=150, n_restarts=2, top_k=15)
    qaoa_pairs = [candidates[k] for k, b in enumerate(qaoa_bits) if b == 1]
    qaoa_structure = pairs_to_dot_bracket(len(sequence), qaoa_pairs)
    qaoa_qubo_energy = energy(qaoa_bits, Q)

    vqe_bits, _, _ = run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200)
    vqe_pairs = [candidates[k] for k, b in enumerate(vqe_bits) if b == 1]
    vqe_structure = pairs_to_dot_bracket(len(sequence), vqe_pairs)
    vqe_qubo_energy = energy(vqe_bits, Q)

    # select the one with lower QUBO energy (internal solver comparison)
    if qaoa_qubo_energy <= vqe_qubo_energy:
        quantum_structure, quantum_method = qaoa_structure, "QAOA"
    else:
        quantum_structure, quantum_method = vqe_structure, "CVaR-VQE"

    # Real ViennaRNA-based comparison: same units, same thermodynamic model as the MFE, unlike the QUBO energy above.
    quantum_real_energy = real_energy(sequence, quantum_structure)
    real_energy_gap = quantum_real_energy - mfe
    metrics = base_pair_metrics(quantum_structure, vienna_structure)

    print(f"Sequence {index + 1}/{len(BENCHMARK_SEQUENCES)}: {label} ({sequence}, {len(sequence)} nt)")
    print(f"ViennaRNA MFE structure: {vienna_structure}  ({mfe:.2f} kcal/mol)")
    print(f"Best quantum method: {quantum_method}")
    print(f"Quantum structure:       {quantum_structure}  ({quantum_real_energy:.2f} kcal/mol)")
    print(f"Real energy gap vs. MFE (kcal/mol): {real_energy_gap:.2f}")
    print(f"Base-pair precision/recall/F1: {metrics['precision']:.3f} / "
          f"{metrics['recall']:.3f} / {metrics['f1']:.3f}")
    print(f"Base-pair distance: {metrics['base_pair_distance']}  "
          f"(0 = exact structural match)")
    print("--------------------------")

    rows.append({
        "label": label,
        "sequence": sequence,
        "length": len(sequence),
        "vienna_mfe_structure": vienna_structure,
        "vienna_mfe_energy": mfe,
        "quantum_method": quantum_method,
        "quantum_structure": quantum_structure,
        "quantum_real_energy": quantum_real_energy,
        "real_energy_gap": real_energy_gap,
        "precision": metrics["precision"],
        "recall": metrics["recall"],
        "f1": metrics["f1"],
        "base_pair_distance": metrics["base_pair_distance"],
        "hamming_distance": metrics["hamming_distance"],
    })

df = pd.DataFrame(rows)
df.to_csv("results/batch_accuracy.csv", index=False)

print("\n========== SUMMARY ==========")
print(f"Sequences tested: {len(BENCHMARK_SEQUENCES)}")
print(f"Mean real energy gap (kcal/mol): {df['real_energy_gap'].mean():.2f}")
print(f"Mean base-pair F1: {df['f1'].mean():.3f}")
print(f"Mean base-pair distance: {df['base_pair_distance'].mean():.2f}")
print("\nSaved: results/batch_accuracy.csv")
