"""
Compares the two QUBO encodings on the project's curated test sequences:

  1. Pair-indicator (rna_to_qubo_full.py)   -- 1 qubit per candidate pair
  2. One-hot per-position (rna_to_qubo_onehot.py) -- 2 qubits per candidate
     pair (one per direction) + explicit per-position one-hot penalties

For every sequence:
  - qubit count for each encoding (measured, not estimated)
  - number of quadratic penalty terms each encoding's constraint
    machinery adds (a rough proxy for "how much constraint enforcement
    costs" beyond the qubit count itself)
  - for sequences small enough to brute-force (<= BRUTE_FORCE_QUBIT_CAP
    qubits in the *larger* encoding): whether both encodings agree on
    the optimal structure, and whether that structure matches the real
    ViennaRNA MFE structure

Writes results/qubo_encoding_comparison.csv and prints a summary table.
"""

import numpy as np
import pandas as pd
import RNA

from rna_to_qubo_full import (
    MIN_LOOP_SIZE,
    build_qubo as build_qubo_pair_indicator,
    brute_force_solve as brute_force_pair_indicator,
    get_candidate_pairs,
)
from rna_to_qubo_onehot import (
    brute_force_solve_onehot,
    build_qubo_onehot,
    decode_solution,
)
from test_sequences import TEST_SEQUENCE_10NT, TEST_SEQUENCE_12NT
from benchmark_sequences import BENCHMARK_SEQUENCES

# Brute force is 2^n_vars -- cap on the *one-hot* qubit count (the larger
# of the two encodings), since that's the binding constraint for which
# sequences can get a real "did both encodings find the same optimum"
# check rather than just a qubit-count comparison.
BRUTE_FORCE_QUBIT_CAP = 22

SEQUENCES = [
    ("10nt (test_sequences)", TEST_SEQUENCE_10NT),
    ("12nt (test_sequences)", TEST_SEQUENCE_12NT),
] + [(label, seq) for label, seq, *_ in BENCHMARK_SEQUENCES]


def pairs_to_dot_bracket(pairs, length):
    db = ["."] * length
    for i, j in pairs:
        db[i] = "("
        db[j] = ")"
    return "".join(db)


def n_penalty_terms(Q):
    """Count off-diagonal nonzero entries above the diagonal -- i.e. how
    many distinct quadratic penalty/reward terms the QUBO actually
    encodes."""
    n = Q.shape[0]
    return int(np.count_nonzero(np.triu(Q, k=1)))


def run():
    rows = []

    for label, sequence in SEQUENCES:
        candidates = get_candidate_pairs(sequence, min_loop_size=MIN_LOOP_SIZE)
        n_candidates = len(candidates)

        Q_pair = build_qubo_pair_indicator(candidates)
        Q_onehot, var_index = build_qubo_onehot(candidates)

        pair_qubits = Q_pair.shape[0]
        onehot_qubits = Q_onehot.shape[0]

        row = {
            "sequence_label": label,
            "sequence": sequence,
            "length": len(sequence),
            "n_candidate_pairs": n_candidates,
            "pair_indicator_qubits": pair_qubits,
            "onehot_qubits": onehot_qubits,
            "qubit_ratio_onehot_over_pair": (
                onehot_qubits / pair_qubits if pair_qubits else float("nan")
            ),
            "pair_indicator_penalty_terms": n_penalty_terms(Q_pair),
            "onehot_penalty_terms": n_penalty_terms(Q_onehot),
        }

        vienna_structure, vienna_mfe = RNA.fold(sequence)
        row["vienna_mfe_structure"] = vienna_structure
        row["vienna_mfe_kcal_mol"] = vienna_mfe

        if onehot_qubits <= BRUTE_FORCE_QUBIT_CAP:
            best_x_pair, _ = brute_force_pair_indicator(Q_pair)
            pair_selected = {
                (i, j) for (i, j, w), bit in zip(candidates, best_x_pair) if bit == 1
            }
            pair_db = pairs_to_dot_bracket(pair_selected, len(sequence))

            best_x_onehot, _ = brute_force_solve_onehot(Q_onehot)
            onehot_selected = decode_solution(best_x_onehot, var_index)
            onehot_db = pairs_to_dot_bracket(onehot_selected, len(sequence))

            row["brute_forced"] = True
            row["pair_indicator_structure"] = pair_db
            row["onehot_structure"] = onehot_db
            row["encodings_agree"] = pair_selected == onehot_selected
            row["pair_indicator_matches_vienna"] = pair_db == vienna_structure
            row["onehot_matches_vienna"] = onehot_db == vienna_structure
        else:
            row["brute_forced"] = False
            row["pair_indicator_structure"] = None
            row["onehot_structure"] = None
            row["encodings_agree"] = None
            row["pair_indicator_matches_vienna"] = None
            row["onehot_matches_vienna"] = None

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("results/qubo_encoding_comparison.csv", index=False)

    print("===== QUBO encoding comparison: pair-indicator vs. one-hot per-position =====\n")
    print(
        f"{'sequence':24} {'len':>4} {'cands':>6} {'pair-qubits':>11} "
        f"{'onehot-qubits':>13} {'pair-terms':>10} {'onehot-terms':>12}  brute-forced"
    )
    for _, r in df.iterrows():
        print(
            f"{r['sequence_label']:24} {r['length']:4d} {r['n_candidate_pairs']:6d} "
            f"{r['pair_indicator_qubits']:11d} {r['onehot_qubits']:13d} "
            f"{r['pair_indicator_penalty_terms']:10d} {r['onehot_penalty_terms']:12d}  "
            f"{r['brute_forced']}"
        )

    brute_forced = df[df["brute_forced"]]
    if len(brute_forced):
        n_agree = brute_forced["encodings_agree"].sum()
        n_pair_ok = brute_forced["pair_indicator_matches_vienna"].sum()
        n_onehot_ok = brute_forced["onehot_matches_vienna"].sum()
        print(
            f"\nOf {len(brute_forced)} brute-forceable sequences: "
            f"encodings agree with each other on {n_agree}, "
            f"pair-indicator matches ViennaRNA MFE on {n_pair_ok}, "
            f"one-hot matches ViennaRNA MFE on {n_onehot_ok}."
        )

    avg_ratio = df["qubit_ratio_onehot_over_pair"].mean()
    print(f"\nMean qubit ratio (one-hot / pair-indicator): {avg_ratio:.2f}x")
    print("Wrote results/qubo_encoding_comparison.csv")


if __name__ == "__main__":
    run()
