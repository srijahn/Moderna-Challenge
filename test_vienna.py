import RNA

sequence = "AUGCGGAUAC"

structure, mfe = RNA.fold(sequence)

print("Sequence:", sequence)
print("Structure:", structure)
print("MFE:", mfe)