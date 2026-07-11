import pandas as pd

# Example values from your project

data = {
    "Metric": [
        "ViennaRNA MFE",
        "Quantum Candidate Energy",
        "Energy Gap",
        "Approx Accuracy (%)",
        "QAOA Success Probability",
        "Best Gamma",
        "Best Beta"
    ],

    "Value": [
        -4.0,
        -3.0,
        1.0,
        75.0,
        0.4796,
        0.7,
        0.3
    ]
}

df = pd.DataFrame(data)

df.to_csv("results/final_results.csv", index=False)

print("\n========== FINAL RESULTS ==========\n")
print(df)

print("\nSaved:")
print("results/final_results.csv")