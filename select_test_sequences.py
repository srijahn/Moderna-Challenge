"""

Need short sequences (deterministic,
seeded) for ones that:
  1. Have a non-trivial ViennaRNA MFE structure (at least one base pair).
  2. Stay within a qubit budget that brute_force_solve() can still validate
     exactly (<= ~20 candidate-pair variables), so quantum results can
     always be checked against a known-exact ground truth.

The winning sequences are hardcoded into test_sequences.py so every script
in the project (qaoa_rna_solver.py, cvar_vqe_rna_solver.py,
compare_to_vienna.py, batch_accuracy.py, generate_final_results.py,
noise_simulation.py) references the same, verified-folding validation set
instead of each picking its own arbitrary example.


"""

import random

import RNA

from rna_to_qubo_full import get_candidate_pairs


def search(length, n_tries=2000, max_qubits=20):
    """Try `n_tries` deterministic random sequences of the given length and
    return the ones that (a) actually fold under ViennaRNA and (b) fit
    within `max_qubits` candidate pairs, sorted by MFE (most stable first).
    """
    results = []
    for seed in range(n_tries):
        rng = random.Random(seed)
        seq = "".join(rng.choice("AUCG") for _ in range(length))

        structure, mfe = RNA.fold(seq)
        if structure.count("(") == 0:
            continue  # unfolded -- not a useful validation case

        candidates = get_candidate_pairs(seq)
        n_qubits = len(candidates)
        if n_qubits == 0 or n_qubits > max_qubits:
            continue

        results.append({
            "seed": seed,
            "sequence": seq,
            "structure": structure,
            "mfe": mfe,
            "n_qubits": n_qubits,
            "n_pairs": structure.count("("),
        })

    results.sort(key=lambda r: r["mfe"])  # most negative (most stable) first
    return results


if __name__ == "__main__":
    for length in (10, 12, 15, 16, 18):
        results = search(length)
        print(f"--- length {length}: {len(results)} folding candidates found ---")
        for r in results[:5]:
            print(
                f"  seed={r['seed']:5d}  qubits={r['n_qubits']:3d}  "
                f"pairs={r['n_pairs']}  mfe={r['mfe']:6.2f}  "
                f"{r['sequence']}  {r['structure']}"
            )
        print()
