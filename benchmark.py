import RNA
import pandas as pd

df = pd.read_csv("data/sequences.csv")

results = []

for seq in df["sequence"]:
    structure, mfe = RNA.fold(seq)

    results.append({
        "sequence": seq,
        "structure": structure,
        "mfe": mfe
    })

results_df = pd.DataFrame(results)

results_df.to_csv(
    "results/vienna_results.csv",
    index=False
)

print(results_df.head())
print("\nBenchmark Complete!")