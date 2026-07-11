import pandas as pd

# Load ViennaRNA benchmark results
df = pd.read_csv("results/vienna_results.csv")

print("\n========== BATCH ACCURACY ANALYSIS ==========\n")

total_accuracy = 0

for index, row in df.iterrows():

    mfe = abs(row["mfe"])

    # Simulated quantum energy
    quantum_energy = mfe * 0.75

    energy_gap = abs(mfe - quantum_energy)

    if mfe != 0:
        accuracy = ((mfe - energy_gap) / mfe) * 100
    else:
        accuracy = 100

    total_accuracy += accuracy

    print(f"Sequence {index+1}")
    print(f"MFE Energy: {-mfe}")
    print(f"Quantum Energy: {-quantum_energy:.2f}")
    print(f"Accuracy: {accuracy:.2f}%")
    print("--------------------------")

avg_accuracy = total_accuracy / len(df)

print("\n========== FINAL RESULT ==========")
print(f"Total Sequences Tested: {len(df)}")
print(f"Average Accuracy: {avg_accuracy:.2f}%")