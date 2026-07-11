import RNA

sequence = "GGUGCCGAACAGUAGCACUC"

# ViennaRNA reference
structure, mfe = RNA.fold(sequence)

# Best energy from your quantum-inspired solver
quantum_energy = -3.0

energy_gap = abs(mfe - quantum_energy)

accuracy = (abs(quantum_energy) / abs(mfe)) * 100

print()
print("========== RNA Folding Comparison ==========")

print("\nSequence:")
print(sequence)

print("\nViennaRNA Structure:")
print(structure)

print("\nViennaRNA Energy:")
print(mfe)

print("\nQuantum Candidate Energy:")
print(quantum_energy)

print("\nEnergy Gap:")
print(round(energy_gap, 3))

print("\nApprox Accuracy:")
print(round(accuracy, 2), "%")

print("\n============================================")