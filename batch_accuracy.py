import pandas as pd
import RNA

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe

# --- Original approach ------------------------------------------------------
# Looped over all 50 generated 20 nt sequences in results/vienna_results.csv.
# Each 20 nt sequence has 50+ candidate base pairs (see findings_week1.md),
# i.e. 50+ qubits -- far past what local statevector simulation of
# QAOA/CVaR-VQE can handle. Left here for reference only.
# df = pd.read_csv("results/vienna_results.csv")

# Hardcoded short test sequences (the same ones validated end-to-end in
# qaoa_rna_solver.py / cvar_vqe_rna_solver.py) so both real quantum methods
# can actually run locally for this batch comparison.
# seq_15 below == generate_test_sequence(15, seed=15) from
# cvar_vqe_rna_solver.py's __main__ block; hardcoded here instead of
# regenerated so this script has no dependency on that helper.
sequences = [
    "GGUGCCGAAC",       # 10 nt
    "UAAUUAAUCUACGCC",  # 15 nt
]

print("\n========== BATCH ACCURACY ANALYSIS ==========\n")

rows = []
total_accuracy = 0

for index, sequence in enumerate(sequences):
    _, mfe_signed = RNA.fold(sequence)
    mfe = abs(mfe_signed)

    candidates = get_candidate_pairs(sequence)
    Q = build_qubo(candidates)

    qaoa_bits, _, _ = run_qaoa(Q, n_layers=3, steps=150, n_restarts=2, top_k=15)
    qaoa_energy_raw = energy(qaoa_bits, Q)

    vqe_bits, _, _ = run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200)
    vqe_energy_raw = energy(vqe_bits, Q)

    # Best (lowest, i.e. most negative) QUBO energy of the two methods
    quantum_energy = abs(min(qaoa_energy_raw, vqe_energy_raw))

    # NOTE: QUBO energy and ViennaRNA MFE (kcal/mol) are different scales --
    # see the caveat in compare_to_vienna.py. Same accuracy formula as the
    # original placeholder version, just fed real numbers now.
    energy_gap = abs(mfe - quantum_energy)

    if mfe != 0:
        accuracy = ((mfe - energy_gap) / mfe) * 100
    else:
        accuracy = 100

    total_accuracy += accuracy

    print(f"Sequence {index + 1} ({sequence}, {len(sequence)} nt)")
    print(f"MFE Energy: {-mfe}")
    print(f"QAOA QUBO Energy: {qaoa_energy_raw:.2f}")
    print(f"CVaR-VQE QUBO Energy: {vqe_energy_raw:.2f}")
    print(f"Best Quantum QUBO Energy: {-quantum_energy:.2f}")
    print(f"Accuracy: {accuracy:.2f}%")
    print("--------------------------")

    rows.append([sequence, -mfe, qaoa_energy_raw, vqe_energy_raw, -quantum_energy, accuracy])

avg_accuracy = total_accuracy / len(sequences)

print("\n========== FINAL RESULT ==========")
print(f"Total Sequences Tested: {len(sequences)}")
print(f"Average Accuracy: {avg_accuracy:.2f}%")
