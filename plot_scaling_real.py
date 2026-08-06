"""
Plots the REAL scaling data from scaling_analysis_real.py.
plot_scaling.py is plotted a hardcoded list of numbers.
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/scaling_analysis_real.csv")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# 1. Qubits required vs sequence length
axes[0, 0].plot(df["sequence_length"], df["n_qubits"], marker="o", color="tab:blue")
axes[0, 0].set_title("Qubits Required vs. RNA Length")
axes[0, 0].set_xlabel("RNA Length (nt)")
axes[0, 0].set_ylabel("Qubits (= candidate base pairs)")
axes[0, 0].grid(True)

# 2. Hamiltonian terms vs sequence length
axes[0, 1].plot(df["sequence_length"], df["n_hamiltonian_terms"], marker="o", color="tab:orange")
axes[0, 1].set_title("Cost Hamiltonian Terms vs. RNA Length")
axes[0, 1].set_xlabel("RNA Length (nt)")
axes[0, 1].set_ylabel("Number of Hamiltonian terms")
axes[0, 1].grid(True)

# 3. Circuit depth vs sequence length
axes[1, 0].plot(df["sequence_length"], df["circuit_depth"], marker="o", color="tab:green")
axes[1, 0].set_title("QAOA Circuit Depth vs. RNA Length")
axes[1, 0].set_xlabel("RNA Length (nt)")
axes[1, 0].set_ylabel("Circuit depth (decomposed gates)")
axes[1, 0].grid(True)

# 4. Forward-pass runtime vs sequence length (only where actually simulated)
sim_df = df[df["simulated"] == True]
not_sim_df = df[df["simulated"] == False]
axes[1, 1].plot(
    sim_df["sequence_length"], sim_df["forward_runtime_seconds"],
    marker="o", color="tab:red", label="measured"
)
if len(not_sim_df) > 0:
    ymax = sim_df["forward_runtime_seconds"].max() if len(sim_df) else 0
    axes[1, 1].scatter(
        not_sim_df["sequence_length"], [ymax] * len(not_sim_df),
        marker="x", color="gray",
        label="not simulated (exceeds qubit limit)"
    )
axes[1, 1].set_title("Forward-Pass Runtime vs. RNA Length")
axes[1, 1].set_xlabel("RNA Length (nt)")
axes[1, 1].set_ylabel("Runtime (s)")
axes[1, 1].legend()
axes[1, 1].grid(True)

plt.tight_layout()
plt.savefig("results/scaling_plot_real.png", dpi=150)
print("Saved results/scaling_plot_real.png")
