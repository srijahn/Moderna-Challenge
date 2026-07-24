import random

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

lengths = [20, 30, 40, 50]

print("\n========== Scaling Analysis ==========\n")

for length in lengths:

    estimated_variables = int((length * length) / 8)

    estimated_qubits = estimated_variables

    print(f"RNA Length: {length}")
    print(f"Variables: {estimated_variables}")
    print(f"Estimated Qubits: {estimated_qubits}")
    print()