import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/noise_analysis.csv")

plt.figure(figsize=(9, 6))

# NOTE: results/noise_analysis.csv now has one block of noise levels per
# benchmark sequence (see noise_simulation.py -- switched from a single
# fixed 10 nt sequence to looping over all 8 curated benchmark sequences).
# Plot one line per sequence instead of assuming a single global series.
if "Label" in df.columns:
    for label, group in df.groupby("Label", sort=False):
        plt.plot(group["Noise Level"], group["Success Probability"], marker="o", label=label)
    plt.legend(title="Sequence")
else:
    # Backward-compatible fallback for an older single-sequence CSV.
    plt.plot(df["Noise Level"], df["Success Probability"], marker="o")

plt.title("QAOA Performance Under Noise")
plt.xlabel("Noise Level")
plt.ylabel("Success Probability")
plt.grid(True)
plt.savefig("results/noise_plot.png")
plt.show()

print("Graph saved as results/noise_plot.png")
