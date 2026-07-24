"""
Plots the CVaR-VQE scaling data from cvar_vqe_scaling_analysis.py.
Counterpart to plot_scaling_real.py, which plots the QAOA scaling data.
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/cvar_vqe_scaling_analysis.csv")
quality_df = df[df["quality_measured"] == True]

fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# 1. Qubits required vs sequence length
axes[0, 0].plot(df["sequence_length"], df["n_qubits"], marker="o", color="tab:blue")
axes[0, 0].set_title("Qubits Required vs. RNA Length")
axes[0, 0].set_xlabel("RNA Length (nt)")
axes[0, 0].set_ylabel("Qubits (= candidate base pairs)")
axes[0, 0].grid(True)

# 2. Circuit depth vs sequence length
axes[0, 1].plot(df["sequence_length"], df["circuit_depth"], marker="o", color="tab:green")
axes[0, 1].set_title("CVaR-VQE Circuit Depth vs. RNA Length")
axes[0, 1].set_xlabel("RNA Length (nt)")
axes[0, 1].set_ylabel("Circuit depth (decomposed gates)")
axes[0, 1].grid(True)

# 3. Forward-pass runtime vs sequence length (only where actually simulated)
sim_df = df[df["simulated"] == True]
not_sim_df = df[df["simulated"] == False]
axes[0, 2].plot(
    sim_df["sequence_length"], sim_df["forward_runtime_seconds"],
    marker="o", color="tab:red", label="measured"
)
if len(not_sim_df) > 0:
    ymax = sim_df["forward_runtime_seconds"].max() if len(sim_df) else 0
    axes[0, 2].scatter(
        not_sim_df["sequence_length"], [ymax] * len(not_sim_df),
        marker="x", color="gray", label="not simulated (exceeds qubit limit)"
    )
axes[0, 2].set_title("Forward-Pass Runtime vs. RNA Length")
axes[0, 2].set_xlabel("RNA Length (nt)")
axes[0, 2].set_ylabel("Runtime (s)")
axes[0, 2].legend()
axes[0, 2].grid(True)

# 4. Energy gap vs sequence length (only where quality was measured)
axes[1, 0].plot(
    quality_df["sequence_length"], quality_df["mean_energy_gap"],
    marker="o", color="tab:purple"
)
axes[1, 0].set_title("Mean Energy Gap vs. QUBO Optimum")
axes[1, 0].set_xlabel("RNA Length (nt)")
axes[1, 0].set_ylabel("Energy gap (QUBO units)")
axes[1, 0].grid(True)

# 5. Success rate vs sequence length
axes[1, 1].plot(
    quality_df["sequence_length"], quality_df["success_rate"] * 100,
    marker="o", color="tab:orange"
)
axes[1, 1].set_title("Success Rate vs. RNA Length")
axes[1, 1].set_xlabel("RNA Length (nt)")
axes[1, 1].set_ylabel("Success rate (%) -- exact QUBO optimum found")
axes[1, 1].set_ylim(-5, 105)
axes[1, 1].grid(True)

# 6. Full optimization runtime vs sequence length
axes[1, 2].plot(
    quality_df["sequence_length"], quality_df["mean_opt_runtime_seconds"],
    marker="o", color="tab:brown"
)
axes[1, 2].set_title("Mean Optimization Runtime vs. RNA Length")
axes[1, 2].set_xlabel("RNA Length (nt)")
axes[1, 2].set_ylabel("Runtime (s)")
axes[1, 2].grid(True)

note = (
    f"Quality metrics (energy gap, success rate, optimization runtime) are only measured up to "
    f"the qubit count where full CVaR-VQE optimization + brute-force ground truth stays tractable "
    f"(see MAX_OPT_QUBITS in cvar_vqe_scaling_analysis.py). Beyond that, only resource counts "
    f"(qubits/depth/gates) are shown -- reported honestly as simulator limits, not omitted."
)
fig.text(0.5, 0.01, note, ha="center", fontsize=8, wrap=True)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig("results/cvar_vqe_scaling_plot.png", dpi=150)
print("Saved results/cvar_vqe_scaling_plot.png")
