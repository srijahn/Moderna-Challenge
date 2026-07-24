import pandas as pd

# ============================================================================
# SUPERSEDED: this script estimates qubit/variable counts with a made-up
# formula (variables = length^2 / 8), not the actual QUBO. Real, measured
# resource counts (from the actual candidate-pair QUBO + built QAOA/CVaR-VQE
# circuits) now live in:
#   - scaling_analysis_real.py       -> results/scaling_analysis_real.csv
#   - cvar_vqe_scaling_analysis.py   -> results/cvar_vqe_scaling_analysis.csv
# Kept here for reference / comparison only -- do not cite these numbers in
# the final report.
# ============================================================================
print("NOTE: this script uses a made-up estimation formula, not the real")
print("QUBO. See scaling_analysis_real.py / cvar_vqe_scaling_analysis.py")
print("for actual measured resource counts.\n")

print("\n========== QUANTUM RESOURCE ESTIMATION ==========\n")

rna_lengths = [20, 30, 40, 50]

for length in rna_lengths:

    variables = int((length * length) / 8)

    qubits = variables

    circuit_depth = variables * 2

    print(f"RNA Length: {length}")
    print(f"Optimization Variables: {variables}")
    print(f"Estimated Qubits: {qubits}")
    print(f"Estimated Circuit Depth: {circuit_depth}")
    print("-" * 40)

print("\nConclusion:")
print("Quantum resource requirements increase rapidly")
print("with RNA sequence length, making large-scale")
print("RNA folding challenging for current quantum devices.")