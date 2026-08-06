"""
QAOA noise analysis for the real RNA folding QUBO, using PennyLane's
built-in depolarizing noise channel (mixed-state simulation via the
`default.mixed` device).

Trains the real QAOA circuit (same formulation as qaoa_rna_solver.py)
noiselessly, then re-evaluates the *trained, fixed* circuit under
increasing single-qubit depolarizing noise applied after every QAOA layer.

Training uses derivative-free COBYLA (scipy.optimize.minimize).

Metric: instead of guessing at a "success probability", this tracks the
expectation value of the cost Hamiltonian (the actual quantity QAOA
optimizes) as noise increases, normalized against the maximally-mixed-state
expectation (which is 0, since the Ising Hamiltonian is a sum of traceless
Pauli terms). This gives a "solution quality" score that starts near 1.0 at
zero noise and decays toward 0.0 as noise pushes the state toward maximally
mixed -- a physically meaningful and honestly-labeled stand-in for the
"success probability" the original placeholder version guessed at.
"""

import numpy as np
import pandas as pd
import pennylane as qml
from scipy.optimize import minimize

from rna_to_qubo_full import get_candidate_pairs, build_qubo
from qaoa_rna_solver import qubo_to_ising, build_hamiltonians
from benchmark_sequences import BENCHMARK_SEQUENCES

# loops over the full 14-sequence curated benchmark set (benchmark_sequences.py,
# 8-14 nt, 35.7%-100% GC content, all confirmed by ViennaRNA to fold).

n_layers = 3

noise_levels = {
    "No Noise": 0.0,
    "Low Noise": 0.01,
    "Medium Noise": 0.05,
    "High Noise": 0.10,
}

all_rows = []

for label, sequence, _, _, _ in BENCHMARK_SEQUENCES:
    candidates = get_candidate_pairs(sequence)
    Q = build_qubo(candidates)
    n_qubits = Q.shape[0]

    h, J, offset = qubo_to_ising(Q)
    cost_h, mixer_h = build_hamiltonians(h, J, n_qubits)

    print(f"\n===== {label}: {sequence}  ({n_qubits} qubits) =====")

    def qaoa_layer(gamma, beta):
        qml.qaoa.cost_layer(gamma, cost_h)
        qml.qaoa.mixer_layer(beta, mixer_h)

    # --- 1. Train QAOA parameters noiselessly (ideal statevector device) ---
    dev_ideal = qml.device("default.qubit", wires=n_qubits)

    # No diff_method specified and plain numpy inputs below -- pure forward
    # simulation, no autodiff/backprop graph ever gets built.

    @qml.qnode(dev_ideal)
    def cost_function(params):
        gammas, betas = params[:n_layers], params[n_layers:]
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        qml.layer(qaoa_layer, n_layers, gammas, betas)
        return qml.expval(cost_h)

    def objective(flat_params):
        return float(cost_function(flat_params))

    print("Training QAOA (noiseless) to get fixed circuit parameters...")
    best_final_cost, best_params = None, None
    for r in range(2):
        rng = np.random.default_rng(r)
        x0 = np.concatenate([rng.uniform(0, 0.3, n_layers), rng.uniform(0, 0.3, n_layers)])
        result = minimize(
            objective, x0, method="COBYLA",
            options={"maxiter": 150, "rhobeg": 0.3},
        )
        final_cost = float(result.fun)
        if best_final_cost is None or final_cost < best_final_cost:
            best_final_cost, best_params = final_cost, result.x

    print(f"Trained. Noiseless expectation value: {best_final_cost:.4f}\n")

    # --- 2. Re-run the same fixed, trained circuit under depolarizing noise
    def noisy_expectation(depolarizing_p):
        dev_noisy = qml.device("default.mixed", wires=n_qubits)

        @qml.qnode(dev_noisy)
        def circuit(params):
            gammas, betas = params[:n_layers], params[n_layers:]
            for i in range(n_qubits):
                qml.Hadamard(wires=i)
            for layer in range(n_layers):
                qml.qaoa.cost_layer(gammas[layer], cost_h)
                qml.qaoa.mixer_layer(betas[layer], mixer_h)
                if depolarizing_p > 0:
                    for w in range(n_qubits):
                        qml.DepolarizingChannel(depolarizing_p, wires=w)
            return qml.expval(cost_h)

        return float(circuit(best_params))

    ideal_expectation = None

    for level, p in noise_levels.items():
        exp_val = noisy_expectation(p)

        if level == "No Noise":
            ideal_expectation = exp_val

        # Maximally-mixed-state expectation of a traceless Ising Hamiltonian
        # is 0, so quality = exp_val / ideal_expectation is ~1.0 at p=0 and
        # decays toward 0.0 as noise dominates.

        quality = exp_val / ideal_expectation if ideal_expectation else 1.0
        loss = (1 - quality) * 100

        all_rows.append([label, sequence, n_qubits, level, p, round(exp_val, 4), round(quality, 4), round(loss, 2)])

        print(f"Noise Level: {level} (depolarizing p={p})")
        print(f"Cost Hamiltonian Expectation: {exp_val:.4f}")
        print(f"Solution Quality (relative to noiseless): {quality:.4f}")
        print(f"Performance Loss: {loss:.2f}%")
        print("-" * 40)

# Column kept as "Success Probability". so plot_noise.py -- which reads that 
# exact column name -- keeps working unmodified. It's still a 0-1 "higher is 
# better" score, grounded in a real noisy-circuit simulation instead
# of assumed probability.
df = pd.DataFrame(
    all_rows,
    columns=[
        "Label",
        "Sequence",
        "Qubits",
        "Noise Level",
        "Depolarizing Probability",
        "Cost Hamiltonian Expectation",
        "Success Probability",
        "Performance Loss (%)",
    ],
)

df.to_csv("results/noise_analysis.csv", index=False)

print("\nSaved:")
print("results/noise_analysis.csv")
