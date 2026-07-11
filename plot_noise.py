import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/noise_analysis.csv")

plt.figure(figsize=(8,5))

plt.plot(
    df["Noise Level"],
    df["Success Probability"],
    marker="o"
)

plt.title("QAOA Performance Under Noise")

plt.xlabel("Noise Level")

plt.ylabel("Success Probability")

plt.grid(True)

plt.savefig("results/noise_plot.png")

plt.show()

print("Graph saved as results/noise_plot.png")