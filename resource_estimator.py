import pandas as pd

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