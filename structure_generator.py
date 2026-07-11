import RNA

# Example RNA sequence
sequence = "GGUGCCGAACAGUAGCACUC"

# Generate Minimum Free Energy structure
structure, mfe = RNA.fold(sequence)

print("\n========== RNA STRUCTURE GENERATION ==========\n")

print("Sequence:")
print(sequence)

print("\nDot-Bracket Structure:")
print(structure)

print("\nMFE Energy:")
print(mfe)

# Count base pairs
pairs = structure.count("(")

print("\nNumber of Base Pairs:")
print(pairs)

print("\nStructure Length:")
print(len(structure))