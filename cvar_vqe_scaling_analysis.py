"""
Scaling analysis for the CVaR-VQE pipeline (cvar_vqe_rna_solver.py) --
the CVaR-VQE counterpart to scaling_analysis_real.py (which covers QAOA).

Together these two scripts give the Week 3 deliverable: run both quantum
methods across the sequence-length ladder and record, for each: qubit count,
circuit depth, number of QUBO variables, runtime, energy gap from the QUBO
optimum, and success rate.

Two tiers of measurement, same honesty principle as scaling_analysis_real.py
-- report the simulator limit rather than extrapolating past it:

  1. Resource counts (qubits, gate count, circuit depth) via qml.specs
     static analysis. Cheap, so measured across the FULL length ladder
     (up to 50 nt) even though nothing is actually simulated at that size.

  2. Solve quality (energy gap vs. the true QUBO optimum, success rate) via
     actually running CVaR-VQE end to end. This needs real statevector
     simulation plus a brute-force ground truth, so it's restricted to
     MAX_OPT_QUBITS, and run several independent trials per length to get
     a success rate rather than a single anecdote.

Note on "energy gap from MFE" (per the Week 3 plan wording): the ground
truth used here is the exact optimum of our own QUBO (via brute force),
not ViennaRNA's thermodynamic MFE directly -- those two only coincide if
the QUBO's pairing rules + weights are a faithful proxy for real folding
energetics, which is a modeling assumption, not a guarantee. We also print
ViennaRNA's MFE energy alongside for context, but for pseudo-random test
sequences at these lengths it is frequently 0.0 (no stable fold), which
would make a literal MFE energy gap meaningless as a quality signal.
"""

import random
import time

import numpy as np
import pandas as pd
import pennylane as qml

from rna_to_qubo_full import get_candidate_pairs, build_qubo
from cvar_vqe_rna_solver import (
    make_probability_circuit,
    all_bitstring_energies,
    run_cvar_vqe,
)

try:
    import RNA
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
LENGTHS = [8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50]  # same ladder as QAOA scaling
N_LAYERS = 3                # matches cvar_vqe_rna_solver.py default
MAX_SIM_QUBITS = 16         # forward-pass timing cutoff (statevector is 2^n)
MAX_OPT_QUBITS = 14         # full CVaR-VQE optimization + brute-force ground truth cutoff
                             # (kept below MAX_SIM_QUBITS since optimization runs hundreds
                             # of forward passes per trial, not just one)
N_TRIALS = 3                 # independent optimization runs per length, for success rate
TRIAL_MAXITER = 150          # COBYLA budget per trial (kept modest to cover the whole ladder)
ALPHA = 0.15                 # CVaR tail fraction, matches cvar_vqe_rna_solver.py default


def generate_test_sequence(length, seed):
    """Same deterministic generator as scaling_analysis_real.py, so QAOA and
    CVaR-VQE are measured on identical sequences at every length."""
    rng = random.Random(seed)
    return "".join(rng.choice("AUCG") for _ in range(length))


def measure_resources(n_qubits, Q):
    """Static resource counts -- fast even at large qubit counts, no
    simulation required. Mirrors scaling_analysis_real.py's approach."""
    circuit, _ = make_probability_circuit(n_qubits, N_LAYERS)
    params = np.random.uniform(0, 2 * np.pi, size=(N_LAYERS + 1, n_qubits))
    specs = qml.specs(circuit, level="device")(params)
    return specs.resources.depth, specs.resources.num_gates, circuit, params


def measure_forward_runtime(circuit, params):
    t0 = time.time()
    circuit(params)
    return time.time() - t0


def measure_solve_quality(Q, n_qubits, sequence):
    """Run CVaR-VQE end to end N_TRIALS times, compare each to the true QUBO
    optimum (via brute force / exhaustive enumeration), and report the
    average energy gap, success rate, and average wall-clock runtime."""
    energies, _ = all_bitstring_energies(Q)
    true_optimum = float(np.min(energies))

    gaps = []
    runtimes = []
    successes = 0

    for trial in range(N_TRIALS):
        t0 = time.time()
        x_bits, _, _ = run_cvar_vqe(
            Q, n_layers=N_LAYERS, alpha=ALPHA,
            n_restarts=1, maxiter=TRIAL_MAXITER, seed0=trial,
        )
        runtimes.append(time.time() - t0)

        bits_energy = float(np.array(x_bits) @ Q @ np.array(x_bits))
        gap = bits_energy - true_optimum
        gaps.append(gap)
        if gap < 1e-6:
            successes += 1

    vienna_mfe = None
    if HAVE_VIENNA:
        _, vienna_mfe = RNA.fold(sequence)

    return {
        "true_qubo_optimum": true_optimum,
        "mean_energy_gap": float(np.mean(gaps)),
        "success_rate": successes / N_TRIALS,
        "mean_opt_runtime_seconds": float(np.mean(runtimes)),
        "vienna_mfe_energy": vienna_mfe,
    }


def measure_one_length(length):
    seq = generate_test_sequence(length, seed=length)

    candidates = get_candidate_pairs(seq)
    n_qubits = len(candidates)  # REAL qubit requirement, not a formula

    Q = build_qubo(candidates)

    depth, num_gates, circuit, params = measure_resources(n_qubits, Q)

    if n_qubits <= MAX_SIM_QUBITS:
        forward_runtime = measure_forward_runtime(circuit, params)
        simulated = True
    else:
        forward_runtime = None
        simulated = False

    row = {
        "sequence_length": length,
        "n_qubits": n_qubits,
        "circuit_depth": depth,
        "total_gates": num_gates,
        "forward_runtime_seconds": forward_runtime,
        "simulated": simulated,
        "true_qubo_optimum": None,
        "mean_energy_gap": None,
        "success_rate": None,
        "mean_opt_runtime_seconds": None,
        "vienna_mfe_energy": None,
        "quality_measured": False,
    }

    if n_qubits <= MAX_OPT_QUBITS:
        quality = measure_solve_quality(Q, n_qubits, seq)
        row.update(quality)
        row["quality_measured"] = True

    return row


if __name__ == "__main__":
    print("\n===== CVaR-VQE Scaling Analysis (measured, not estimated) =====\n")

    results = []
    for length in LENGTHS:
        row = measure_one_length(length)
        results.append(row)

        sim_note = f"{row['forward_runtime_seconds']:.4f}s" if row["simulated"] else "not simulated"
        if row["quality_measured"]:
            quality_note = (
                f"gap={row['mean_energy_gap']:.3f}  "
                f"success={row['success_rate'] * 100:.0f}%  "
                f"opt_runtime={row['mean_opt_runtime_seconds']:.2f}s"
            )
        else:
            quality_note = "quality not measured (exceeds MAX_OPT_QUBITS)"

        print(
            f"length={length:>3}  qubits={row['n_qubits']:>3}  depth={row['circuit_depth']:>4}  "
            f"gates={row['total_gates']:>5}  forward={sim_note:>12}  {quality_note}"
        )

    df = pd.DataFrame(results)
    df.to_csv("results/cvar_vqe_scaling_analysis.csv", index=False)
    print("\nSaved results/cvar_vqe_scaling_analysis.csv")
