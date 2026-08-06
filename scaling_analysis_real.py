"""
Real scaling analysis for the RNA-folding QAOA pipeline.


  1. Generates a reproducible test sequence for each length
  2. Runs it through the REAL pipeline: candidate pairs -> QUBO -> Ising
  3. Builds the REAL QAOA circuit and asks PennyLane for its actual
     resource counts (qubits, Hamiltonian terms, gate count, circuit depth)
     via qml.specs -- this is static resource analysis, not full statevector
     simulation, so it stays fast even at large qubit counts.
  4. For qubit counts small enough to actually simulate (<= MAX_SIM_QUBITS),
     also times a real forward-pass circuit execution on default.qubit, so
     we get genuine runtime numbers, not estimates.
"""

import random
import time

import numpy as np
import pandas as pd
import pennylane as qml
from pennylane import numpy as pnp

from rna_to_qubo_full import get_candidate_pairs, build_qubo
from qaoa_rna_solver import qubo_to_ising, build_hamiltonians

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
LENGTHS = [8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50]
N_LAYERS = 3                 # QAOA circuit depth (p) used for resource counting
MAX_SIM_QUBITS = 16          # statevector sim is 2^n -- cap real execution timing here


def generate_test_sequence(length, seed):
    """Deterministic random RNA sequence, so results are reproducible."""
    rng = random.Random(seed)
    return "".join(rng.choice("AUCG") for _ in range(length))


def build_circuit(n_qubits, cost_h, mixer_h, n_layers):
    dev = qml.device("default.qubit", wires=n_qubits)

    def qaoa_layer(gamma, beta):
        qml.qaoa.cost_layer(gamma, cost_h)
        qml.qaoa.mixer_layer(beta, mixer_h)

    @qml.qnode(dev)
    def circuit(params):
        for i in range(n_qubits):
            qml.Hadamard(wires=i)
        qml.layer(qaoa_layer, n_layers, params[0], params[1])
        return qml.expval(cost_h)

    return circuit


def measure_one_length(length, n_layers=N_LAYERS):
    seq = generate_test_sequence(length, seed=length)

    candidates = get_candidate_pairs(seq)
    n_qubits = len(candidates)  # REAL qubit requirement, not a formula

    Q = build_qubo(candidates)
    h, J, offset = qubo_to_ising(Q)
    n_ham_terms = int(np.sum(np.abs(h) > 1e-12)) + len(J)

    cost_h, mixer_h = build_hamiltonians(h, J, n_qubits)
    circuit = build_circuit(n_qubits, cost_h, mixer_h, n_layers)

    params = pnp.array([np.random.rand(n_layers), np.random.rand(n_layers)])

    # Static resource analysis -- fast even for large qubit counts, since it
    # doesn't build a full statevector.
    specs = qml.specs(circuit, level="device")(params)
    depth = specs.resources.depth
    num_gates = specs.resources.num_gates

    # Real forward-pass execution timing, only where simulation is feasible.
    if n_qubits <= MAX_SIM_QUBITS:
        t0 = time.time()
        circuit(params)
        forward_runtime = time.time() - t0
        simulated = True
    else:
        forward_runtime = None
        simulated = False

    return {
        "sequence_length": length,
        "n_qubits": n_qubits,
        "n_hamiltonian_terms": n_ham_terms,
        "circuit_depth": depth,
        "total_gates": num_gates,
        "forward_runtime_seconds": forward_runtime,
        "simulated": simulated,
    }


if __name__ == "__main__":
    print("\n========== REAL Scaling Analysis (measured, not estimated) ==========\n")

    results = []
    for length in LENGTHS:
        row = measure_one_length(length)
        results.append(row)
        sim_note = f"{row['forward_runtime_seconds']:.4f}s" if row["simulated"] else "not simulated (qubit count exceeds MAX_SIM_QUBITS)"
        print(
            f"length={length:>3}  qubits={row['n_qubits']:>4}  "
            f"H-terms={row['n_hamiltonian_terms']:>6}  depth={row['circuit_depth']:>5}  "
            f"gates={row['total_gates']:>6}  forward_runtime={sim_note}"
        )

    df = pd.DataFrame(results)
    df.to_csv("results/scaling_analysis_real.csv", index=False)
    print("\nSaved results/scaling_analysis_real.csv")
