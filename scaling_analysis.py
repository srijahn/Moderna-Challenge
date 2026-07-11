import random

lengths = [20, 30, 40, 50]

print("\n========== Scaling Analysis ==========\n")

for length in lengths:

    estimated_variables = int((length * length) / 8)

    estimated_qubits = estimated_variables

    print(f"RNA Length: {length}")
    print(f"Variables: {estimated_variables}")
    print(f"Estimated Qubits: {estimated_qubits}")
    print()