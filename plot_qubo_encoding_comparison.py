"""
Plots results/qubo_encoding_comparison.csv: qubit count per sequence for
the pair-indicator vs. one-hot encoding, side by side.
"""

import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("results/qubo_encoding_comparison.csv")
df = df.sort_values("n_candidate_pairs")

x = range(len(df))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar([i - width / 2 for i in x], df["pair_indicator_qubits"], width, label="pair-indicator")
ax.bar([i + width / 2 for i in x], df["onehot_qubits"], width, label="one-hot per-position")

ax.set_xticks(list(x))
ax.set_xticklabels(df["sequence_label"], rotation=45, ha="right")
ax.set_ylabel("Qubits")
ax.set_title("QUBO qubit count by encoding, per test sequence")
ax.legend()
fig.tight_layout()
fig.savefig("results/qubo_encoding_comparison_plot.png", dpi=150)
print("Wrote results/qubo_encoding_comparison_plot.png")
