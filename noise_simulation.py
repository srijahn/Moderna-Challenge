import pandas as pd

ideal_probability = 0.4796

noise_levels = {
    "No Noise": ideal_probability,
    "Low Noise": ideal_probability * 0.90,
    "Medium Noise": ideal_probability * 0.75,
    "High Noise": ideal_probability * 0.50
}

print("\n========== QAOA NOISE ANALYSIS ==========\n")

results = []

for level, prob in noise_levels.items():

    loss = ((ideal_probability - prob) / ideal_probability) * 100

    results.append([
        level,
        round(prob, 4),
        round(loss, 2)
    ])

    print(f"Noise Level: {level}")
    print(f"Success Probability: {prob:.4f}")
    print(f"Performance Loss: {loss:.2f}%")
    print("-" * 40)

df = pd.DataFrame(
    results,
    columns=[
        "Noise Level",
        "Success Probability",
        "Performance Loss (%)"
    ]
)

df.to_csv(
    "results/noise_analysis.csv",
    index=False
)

print("\nSaved:")
print("results/noise_analysis.csv")