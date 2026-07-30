"""
QAOA solver for the RNA folding QUBO (rna_to_qubo_full.py).

Replaces the toy 2-qubit circuit in qaoa_optimizer.py with a real QAOA run
over however many qubits the candidate-pair QUBO needs (one qubit per
candidate base pair).

Pipeline:
    sequence -> candidate pairs -> QUBO (Q matrix) -> Ising (h, J, offset)
    -> QAOA circuit -> trained parameters -> most likely bitstring
    -> decode back into selected base pairs -> compare to brute force
"""

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve
from benchmark_sequences import BENCHMARK_SEQUENCES


# ---------------------------------------------------------------------------
# 1. QUBO -> Ising conversion
# ---------------------------------------------------------------------------
def qubo_to_ising(Q):
    """
    Converts a QUBO matrix Q (x_i in {0,1}, energy = x^T Q x) into an
    equivalent Ising problem (z_i in {-1,+1}) via x_i = (1 - z_i) / 2.

    Returns:
        h       : array of linear (single-qubit PauliZ) coefficients
        J       : dict {(i,j): coeff} of quadratic (ZZ) coefficients, i<j
        offset  : constant term. QUBO_energy(x) = Ising_energy(z) + offset
    """
    n = Q.shape[0]
    Qs = (Q + Q.T) / 2  # symmetrize -- x^T Q x = x^T Qs x for any Q

    total_sum = Qs.sum()
    diag_sum = np.trace(Qs)
    offset = (total_sum + diag_sum) / 4

    row_sums = Qs.sum(axis=1)
    h = -row_sums / 2

    J = {}
    for i in range(n):
        for j in range(i + 1, n):
            if abs(Qs[i, j]) > 1e-12:
                J[(i, j)] = Qs[i, j] / 2

    return h, J, offset


# ---------------------------------------------------------------------------
# 2. Build QAOA cost + mixer Hamiltonians
# ---------------------------------------------------------------------------
def build_hamiltonians(h, J, n_qubits):
    coeffs = []
    ops = []

    for i, hi in enumerate(h):
        if abs(hi) > 1e-12:
            coeffs.append(hi)
            ops.append(qml.PauliZ(i))

    for (i, j), Jij in J.items():
        coeffs.append(Jij)
        ops.append(qml.PauliZ(i) @ qml.PauliZ(j))

    cost_h = qml.Hamiltonian(coeffs, ops)
    mixer_h = qml.Hamiltonian([1.0] * n_qubits, [qml.PauliX(i) for i in range(n_qubits)])
    return cost_h, mixer_h


def int_to_bits(idx, n_qubits):
    """int -> bitstring, wire 0 = most significant bit in PennyLane's probs
    ordering. Measurement outcome |0> has PauliZ eigenvalue +1, |1> has
    eigenvalue -1. Since x_i = (1 - z_i) / 2, x_i is just the measured bit."""
    return [(idx >> (n_qubits - 1 - k)) & 1 for k in range(n_qubits)]


# ---------------------------------------------------------------------------
# 3. QAOA circuit, with multi-restart + top-K classical post-selection
#
# Shallow QAOA (small n_layers) doesn't always peak sharply on the exact
# ground state -- that's expected/realistic NISQ-era behavior, not a bug.
# Standard practice is to (a) restart the optimizer from a few random seeds
# and keep the best, then (b) take the top-K most probable bitstrings from
# the final distribution and evaluate their *true* QUBO energy classically
# (cheap, since evaluating a candidate is trivial), returning the best one.
# ---------------------------------------------------------------------------
def run_qaoa(Q, n_layers=3, steps=150, step_size=0.03, n_restarts=2, top_k=15, seed_offset=0):
    # seed_offset lets callers run independent trials (e.g. for a statistical
    # benchmark across multiple runs): without it, every call restarts from
    # the exact same internal seeds (0, 1, ..., n_restarts-1) and would
    # silently return identical results every time it's called.
    n_qubits = Q.shape[0]
    h, J, offset = qubo_to_ising(Q)
    cost_h, mixer_h = build_hamiltonians(h, J, n_qubits)

    dev = qml.device("default.qubit", wires=n_qubits)

    def qaoa_layer(gamma, beta):
        qml.qaoa.cost_layer(gamma, cost_h)
        qml.qaoa.mixer_layer(beta, mixer_h)

    @qml.qnode(dev, diff_method="backprop")
    def cost_function(params):
        gammas, betas = params[0], params[1]
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        qml.layer(qaoa_layer, n_layers, gammas, betas)
        return qml.expval(cost_h)

    @qml.qnode(dev, diff_method="backprop")
    def probability_circuit(params):
        gammas, betas = params[0], params[1]
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        qml.layer(qaoa_layer, n_layers, gammas, betas)
        return qml.probs(wires=range(n_qubits))

    best_final_cost = None
    best_params = None

    for seed in range(n_restarts):
        np.random.seed(seed_offset + seed)
        params = pnp.array(
            [np.random.uniform(0, 0.3, n_layers), np.random.uniform(0, 0.3, n_layers)],
            requires_grad=True,
        )
        opt = qml.AdamOptimizer(stepsize=step_size)
        for step in range(steps):
            params = opt.step(cost_function, params)

        final_cost = float(cost_function(params))
        print(f"  restart {seed + 1}/{n_restarts}: final expectation = {final_cost:.4f}")

        if best_final_cost is None or final_cost < best_final_cost:
            best_final_cost = final_cost
            best_params = params

    probs = probability_circuit(best_params)
    top_indices = np.argsort(probs)[::-1][:top_k]

    # classically evaluate the true QUBO energy of each candidate and keep the best
    best_x_bits, best_qubo_energy, best_prob_seen = None, None, None
    for idx in top_indices:
        bits = int_to_bits(int(idx), n_qubits)
        e = energy(bits, Q)
        if best_qubo_energy is None or e < best_qubo_energy:
            best_qubo_energy = e
            best_x_bits = bits
            best_prob_seen = float(probs[idx])

    return best_x_bits, best_prob_seen, offset


# ---------------------------------------------------------------------------
# 4. Run + validate against brute force
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # NOTE: this used to run on a single fixed sequence (TEST_SEQUENCE_10NT).
    # Switched to loop over the full 8-sequence curated benchmark set
    # (benchmark_sequences.py, 8-14 nt, 35.7%-100% GC content, all confirmed
    # by ViennaRNA to fold and sized to stay under the ~90s/11-qubit runtime
    # observed here -- see README's QAOA scaling finding for why nothing
    # bigger is used at full step budget). Same set cvar_vqe_rna_solver.py
    # and batch_accuracy.py use, for direct comparability.
    summary_rows = []

    for label, sequence, _, _, _ in BENCHMARK_SEQUENCES:
        candidates = get_candidate_pairs(sequence)
        print(f"\n===== {label}: {sequence} =====")
        print(f"Candidate pairs (qubits needed): {len(candidates)}\n")

        Q = build_qubo(candidates)

        print("Running QAOA...")
        x_bits, best_prob, offset = run_qaoa(
            Q, n_layers=3, steps=150, n_restarts=2, top_k=15
        )

        qaoa_energy = energy(x_bits, Q)
        selected_qaoa = [candidates[k] for k, b in enumerate(x_bits) if b == 1]

        print(f"\nQAOA result (most probable bitstring, p={best_prob:.3f}):")
        print(f"  Energy: {qaoa_energy}")
        print("  Selected pairs:")
        for i, j, w in selected_qaoa:
            print(f"    ({i},{j}) {sequence[i]}-{sequence[j]}")

        # Ground truth from brute force (only feasible for small n_vars)
        best_x, best_e = brute_force_solve(Q)
        selected_bf = [candidates[k] for k, b in enumerate(best_x) if b == 1]

        print(f"\nBrute-force optimum (ground truth):")
        print(f"  Energy: {best_e}")
        print("  Selected pairs:")
        for i, j, w in selected_bf:
            print(f"    ({i},{j}) {sequence[i]}-{sequence[j]}")

        gap = abs(qaoa_energy - best_e)
        print(f"\nEnergy gap (QAOA vs. brute force): {gap}")
        if gap < 1e-6:
            print("QAOA found the exact optimum.")
        else:
            print("QAOA landed on a near-optimal (but not exact) solution.")
            print("This is normal for shallow QAOA circuits -- for the final")
            print("scaling analysis, report this gap honestly rather than hiding")
            print("it. To try closing it further: increase n_layers (circuit")
            print("depth), increase steps, or increase n_restarts/top_k above.")

        summary_rows.append((label, sequence, len(candidates), gap))

    print("\n===== Summary across all benchmark sequences =====")
    for label, sequence, n_qubits, gap in summary_rows:
        print(f"  {label:8} qubits={n_qubits:2d}  gap={gap:.4f}")
