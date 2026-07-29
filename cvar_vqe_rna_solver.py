"""
CVaR-VQE solver for the RNA folding QUBO (rna_to_qubo_full.py).

This is the SECOND independent quantum method required by the Week 2 plan
(Coder A track) alongside qaoa_rna_solver.py (Coder B track / QAOA), so two
independent quantum results exist for every test sequence.

Method: Conditional Value-at-Risk VQE (Barkoutsos et al., 2020). Instead of
optimizing the full expectation value of the cost Hamiltonian (which is
dragged toward the *average* energy over all 2^n basis states -- most of
which are bad, since only a few bitstrings satisfy the pairing constraints),
CVaR-VQE only optimizes the average energy of the best alpha-fraction of the
measured distribution. This tends to converge faster on combinatorial
problems like this one, where the good solutions live in a small corner of
the space.

Ansatz: hardware-efficient two-local circuit (RY rotation layers + linear
CNOT entanglers), the same ansatz family referenced in Moderna's reference
demo.

Pipeline:
    sequence -> candidate pairs -> QUBO (Q matrix)
    -> two-local VQE ansatz -> CVaR-optimized parameters
    -> most likely low-energy bitstring -> decode to dot-bracket
    -> compare to ViennaRNA MFE (structure + energy) and to brute force
"""

import time

import numpy as np
import pennylane as qml
from scipy.optimize import minimize

from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve

try:
    import RNA
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False


# ---------------------------------------------------------------------------
# 1. Two-local ansatz (RY rotation layers + linear CNOT entanglers)
# ---------------------------------------------------------------------------
def two_local_circuit(params, n_qubits, n_layers):
    """params shape: (n_layers + 1, n_qubits) -- one RY rotation layer per
    entangling block, plus a final rotation layer with no entangler after it
    (standard TwoLocal convention)."""
    for layer in range(n_layers):
        for q in range(n_qubits):
            qml.RY(params[layer, q], wires=q)
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q + 1])
    # final rotation layer
    for q in range(n_qubits):
        qml.RY(params[n_layers, q], wires=q)


def make_probability_circuit(n_qubits, n_layers):
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def circuit(params):
        two_local_circuit(params, n_qubits, n_layers)
        return qml.probs(wires=range(n_qubits))

    return circuit, dev


# ---------------------------------------------------------------------------
# 2. Precompute the QUBO energy of every basis state (static -- depends only
#    on Q, not on the variational parameters). Feasible for the qubit counts
#    used here (<=~18), same limit already documented in qaoa_rna_solver.py /
#    scaling_analysis_real.py.
# ---------------------------------------------------------------------------
def all_bitstring_energies(Q):
    n = Q.shape[0]
    n_states = 2 ** n
    # bit matrix: row i = binary representation of i, wire 0 = MSB (matches
    # PennyLane's qml.probs ordering, same convention as int_to_bits() in
    # qaoa_rna_solver.py)
    X = ((np.arange(n_states)[:, None] >> np.arange(n - 1, -1, -1)) & 1).astype(float)
    energies = np.einsum("ij,jk,ik->i", X, Q, X)
    return energies, X


# ---------------------------------------------------------------------------
# 3. CVaR objective: weighted average energy of the alpha-fraction lowest-
#    energy probability mass in the measured distribution.
# ---------------------------------------------------------------------------
def cvar_from_probs(probs, sorted_order, sorted_energies, alpha):
    sorted_probs = probs[sorted_order]

    remaining = alpha
    total = 0.0
    weight_sum = 0.0
    for p, e in zip(sorted_probs, sorted_energies):
        w = min(p, remaining)
        if w <= 0:
            continue
        total += w * e
        weight_sum += w
        remaining -= w
        if remaining <= 1e-12:
            break

    if weight_sum <= 1e-12:
        return sorted_energies[0]
    return total / weight_sum


# ---------------------------------------------------------------------------
# 4. Run CVaR-VQE with multi-restart + top-K classical post-selection
#    (same post-selection idea as qaoa_rna_solver.py, for a fair comparison
#    between the two methods).
# ---------------------------------------------------------------------------
def run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200, top_k=15, seed0=0):
    n_qubits = Q.shape[0]
    energies, _ = all_bitstring_energies(Q)
    sorted_order = np.argsort(energies)
    sorted_energies = energies[sorted_order]

    circuit, dev = make_probability_circuit(n_qubits, n_layers)

    def objective(flat_params):
        params = flat_params.reshape(n_layers + 1, n_qubits)
        probs = np.array(circuit(params))
        return cvar_from_probs(probs, sorted_order, sorted_energies, alpha)

    best_cvar = None
    best_params = None

    for r in range(n_restarts):
        rng = np.random.default_rng(seed0 + r)
        x0 = rng.uniform(0, 2 * np.pi, size=(n_layers + 1) * n_qubits)

        result = minimize(
            objective, x0, method="COBYLA",
            options={"maxiter": maxiter, "rhobeg": 0.5},
        )
        print(f"  restart {r + 1}/{n_restarts}: final CVaR (alpha={alpha}) = {result.fun:.4f}")

        if best_cvar is None or result.fun < best_cvar:
            best_cvar = result.fun
            best_params = result.x

    final_params = best_params.reshape(n_layers + 1, n_qubits)
    probs = np.array(circuit(final_params))

    top_indices = np.argsort(probs)[::-1][:top_k]
    best_x_bits, best_qubo_energy, best_prob_seen = None, None, None
    for idx in top_indices:
        bits = [(int(idx) >> (n_qubits - 1 - k)) & 1 for k in range(n_qubits)]
        e = energy(bits, Q)
        if best_qubo_energy is None or e < best_qubo_energy:
            best_qubo_energy = e
            best_x_bits = bits
            best_prob_seen = float(probs[idx])

    depth = qml.specs(circuit, level="device")(final_params).resources.depth

    return best_x_bits, best_prob_seen, depth


# ---------------------------------------------------------------------------
# 5. Decode a set of selected (i, j, w) pairs into dot-bracket notation
# ---------------------------------------------------------------------------
def pairs_to_dot_bracket(length, selected_pairs):
    structure = ["."] * length
    for i, j, w in selected_pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)


# ---------------------------------------------------------------------------
# 6. Run on the Week 2 test sequences (10 nt and 15 nt), validate against
#    brute force, and compare to the ViennaRNA MFE structure/energy.
# ---------------------------------------------------------------------------
def evaluate_sequence(sequence, label):
    print(f"\n===== {label}: {sequence} (length {len(sequence)}) =====")

    candidates = get_candidate_pairs(sequence)
    n_qubits = len(candidates)
    print(f"Candidate pairs (qubits needed): {n_qubits}")

    Q = build_qubo(candidates)

    print("Running CVaR-VQE...")
    t0 = time.time()
    x_bits, best_prob, depth = run_cvar_vqe(Q, n_layers=3, alpha=0.15, n_restarts=2, maxiter=200)
    runtime = time.time() - t0

    vqe_energy = energy(x_bits, Q)
    selected_vqe = [candidates[k] for k, b in enumerate(x_bits) if b == 1]
    vqe_structure = pairs_to_dot_bracket(len(sequence), selected_vqe)

    print(f"\nCVaR-VQE result (p={best_prob:.3f}, circuit depth={depth}, runtime={runtime:.2f}s):")
    print(f"  QUBO energy: {vqe_energy}")
    print(f"  Dot-bracket: {vqe_structure}")

    best_x, best_e = brute_force_solve(Q)
    selected_bf = [candidates[k] for k, b in enumerate(best_x) if b == 1]
    bf_structure = pairs_to_dot_bracket(len(sequence), selected_bf)

    print(f"\nBrute-force optimum (ground truth for this QUBO):")
    print(f"  QUBO energy: {best_e}")
    print(f"  Dot-bracket: {bf_structure}")

    gap = abs(vqe_energy - best_e)
    print(f"\nQUBO energy gap (CVaR-VQE vs. brute force): {gap}")
    print("CVaR-VQE found the exact QUBO optimum." if gap < 1e-6
          else "CVaR-VQE landed on a near-optimal (but not exact) solution -- "
               "report this honestly in the Week 3 scaling analysis, same as qaoa_rna_solver.py.")

    if HAVE_VIENNA:
        vienna_structure, vienna_mfe = RNA.fold(sequence)
        match = "MATCH" if vienna_structure == vqe_structure else "differs"
        print(f"\nViennaRNA MFE structure: {vienna_structure}  (energy={vienna_mfe})")
        print(f"CVaR-VQE structure:      {vqe_structure}  ({match} vs. ViennaRNA)")
    else:
        print("\n(ViennaRNA not installed -- skipping direct MFE structure comparison.)")

    return {
        "sequence": sequence,
        "n_qubits": n_qubits,
        "vqe_energy": vqe_energy,
        "brute_force_energy": best_e,
        "energy_gap": gap,
        "circuit_depth": depth,
        "runtime_seconds": runtime,
    }


if __name__ == "__main__":
    from test_sequences import TEST_SEQUENCE_10NT, TEST_SEQUENCE_12NT

    # NOTE: this used to be a hardcoded "GGUGCCGAAC" (10 nt) plus a random
    # seed=15 15-mer. Both turned out to have a fully-unpaired real
    # ViennaRNA MFE structure -- i.e. neither actually folds -- so the
    # "compare to ViennaRNA MFE structure/energy" step below was always
    # comparing against an empty structure. Replaced with curated sequences
    # (test_sequences.py) that ViennaRNA confirms fold into a real hairpin.
    # Same 10 nt sequence used in qaoa_rna_solver.py, so the two independent
    # methods are directly comparable on identical input.
    seq_10 = TEST_SEQUENCE_10NT
    seq_12 = TEST_SEQUENCE_12NT

    results = [
        evaluate_sequence(seq_10, "10 nt sequence"),
        evaluate_sequence(seq_12, "12 nt sequence"),
    ]

    print("\n===== Summary =====")
    for r in results:
        print(
            f"len={len(r['sequence']):>2}  qubits={r['n_qubits']:>2}  "
            f"gap={r['energy_gap']:.4f}  depth={r['circuit_depth']:>3}  "
            f"runtime={r['runtime_seconds']:.2f}s"
        )
