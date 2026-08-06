"""
Statistical benchmark: run both quantum methods (QAOA, CVaR-VQE) across all
12 curated sequences in benchmark_sequences.py, with multiple independent
trials per sequence, and report mean +/- std instead of a single-run
anecdote.

For each (sequence, method) pair, a trial:
  1. Runs the solver (run_cvar_vqe / run_qaoa) with a different seed.
  2. Compares the returned bitstring's QUBO energy to the true QUBO optimum
     (via rna_to_qubo_full.brute_force_solve -- exact, since every sequence
     here was chosen to stay within brute-force range, <=20 qubits).
  3. Also decodes to dot-bracket and computes base-pair F1 against the real
     ViennaRNA MFE structure (structure_metrics.py), for context.
  "Success" = found the exact QUBO optimum (gap < 1e-6), same success
  criterion as cvar_vqe_scaling_analysis.py / scaling_analysis_real.py.


Usage:
    python statistical_benchmark.py                    # both methods, overwrite CSVs
    python statistical_benchmark.py --methods vqe       # CVaR-VQE only, all 8 sequences
    python statistical_benchmark.py --methods qaoa --append   # QAOA only, append to
                                                                # existing CSVs instead
                                                                # of overwriting

Output:
    results/statistical_benchmark.csv          -- one row per trial
    results/statistical_benchmark_summary.csv  -- mean +/- std per sequence/method
"""

import argparse
import os
import time

import numpy as np
import pandas as pd

from benchmark_sequences import BENCHMARK_SEQUENCES, gc_content
from rna_to_qubo_full import get_candidate_pairs, build_qubo, energy, brute_force_solve
from qaoa_rna_solver import run_qaoa
from cvar_vqe_rna_solver import run_cvar_vqe, pairs_to_dot_bracket
from structure_metrics import base_pair_metrics

try:
    import RNA
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
N_TRIALS_VQE = 5     # independent CVaR-VQE trials per sequence
N_TRIALS_QAOA = 3    # independent QAOA trials per sequence (slower -- fewer trials)

QAOA_MAX_QUBITS = 11        # sequences above this qubit count are skipped for QAOA
                             # in this sweep (13 nt / 13-qubit entry) 
QAOA_STEPS = 100              # reduced from qaoa_rna_solver.py's default of 150,
QAOA_RESTARTS = 2           # and 1 restart instead of 2 -- batch-sweep runtime tradeoff

VQE_MAXITER = 150            # reduced from cvar_vqe_rna_solver.py's default of 200
VQE_RESTARTS = 2             # and 1 restart instead of 2, same tradeoff as QAOA above
VQE_ALPHA = 0.15             # matches cvar_vqe_rna_solver.py default
VQE_LAYERS = 3               # matches cvar_vqe_rna_solver.py default
QAOA_LAYERS = 3               # matches qaoa_rna_solver.py default

RESULTS_DIR = "results"
TRIALS_CSV = os.path.join(RESULTS_DIR, "statistical_benchmark.csv")
SUMMARY_CSV = os.path.join(RESULTS_DIR, "statistical_benchmark_summary.csv")


# ---------------------------------------------------------------------------
# Per-trial runner
# ---------------------------------------------------------------------------
def run_vqe_trial(Q, trial_idx):
    t0 = time.time()
    x_bits, best_prob, depth = run_cvar_vqe(
        Q, n_layers=VQE_LAYERS, alpha=VQE_ALPHA,
        n_restarts=VQE_RESTARTS, maxiter=VQE_MAXITER, seed0=100 * trial_idx,
    )
    runtime = time.time() - t0
    return x_bits, runtime


def run_qaoa_trial(Q, trial_idx):
    t0 = time.time()
    x_bits, best_prob, _offset = run_qaoa(
        Q, n_layers=QAOA_LAYERS, steps=QAOA_STEPS,
        n_restarts=QAOA_RESTARTS, top_k=15, seed_offset=100 * trial_idx,
    )
    runtime = time.time() - t0
    return x_bits, runtime


TRIAL_RUNNERS = {"vqe": run_vqe_trial, "qaoa": run_qaoa_trial}
N_TRIALS = {"vqe": N_TRIALS_VQE, "qaoa": N_TRIALS_QAOA}


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
def run_sweep(methods):
    rows = []

    for label, sequence, expected_structure, expected_mfe, seed in BENCHMARK_SEQUENCES:
        candidates = get_candidate_pairs(sequence)
        n_qubits = len(candidates)
        Q = build_qubo(candidates)

        best_x, true_optimum = brute_force_solve(Q)
        vienna_structure, vienna_mfe = (RNA.fold(sequence) if HAVE_VIENNA
                                         else (None, None))

        for method in methods:
            if method == "qaoa" and n_qubits > QAOA_MAX_QUBITS:
                print(f"[{label}] skipping QAOA -- {n_qubits} qubits > "
                      f"QAOA_MAX_QUBITS ({QAOA_MAX_QUBITS})")
                continue

            n_trials = N_TRIALS[method]
            print(f"[{label}] {sequence} ({n_qubits} qubits) -- running "
                  f"{method.upper()} x {n_trials} trials")

            for trial in range(n_trials):
                x_bits, runtime = TRIAL_RUNNERS[method](Q, trial)

                trial_energy = energy(x_bits, Q)
                gap = trial_energy - true_optimum
                success = gap < 1e-6

                predicted_structure = pairs_to_dot_bracket(
                    len(sequence), [candidates[k] for k, b in enumerate(x_bits) if b == 1]
                )
                if HAVE_VIENNA:
                    metrics = base_pair_metrics(predicted_structure, vienna_structure)
                    f1 = metrics["f1"]
                else:
                    f1 = None

                rows.append({
                    "sequence_label": label,
                    "sequence": sequence,
                    "length": len(sequence),
                    "gc_percent": gc_content(sequence),
                    "n_qubits": n_qubits,
                    "method": method,
                    "trial": trial,
                    "qubo_energy": trial_energy,
                    "true_qubo_optimum": true_optimum,
                    "energy_gap": gap,
                    "success": success,
                    "runtime_seconds": runtime,
                    "f1_vs_vienna_mfe": f1,
                    "vienna_mfe_kcal_mol": vienna_mfe,
                })

    return pd.DataFrame(rows)


def summarize(df):
    grouped = df.groupby(["sequence_label", "length", "gc_percent", "n_qubits", "method"])
    summary = grouped.agg(
        n_trials=("trial", "count"),
        mean_energy_gap=("energy_gap", "mean"),
        std_energy_gap=("energy_gap", "std"),
        success_rate=("success", "mean"),
        mean_f1_vs_vienna_mfe=("f1_vs_vienna_mfe", "mean"),
        std_f1_vs_vienna_mfe=("f1_vs_vienna_mfe", "std"),
        mean_runtime_seconds=("runtime_seconds", "mean"),
    ).reset_index()
    # std of a single trial is NaN, not an error -- leave as-is, it's informative
    # (tells you n_trials==1 for that row) rather than silently zero-filling it.
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Multi-trial QAOA / CVaR-VQE benchmark across the 8 curated sequences."
    )
    parser.add_argument(
        "--methods", choices=["vqe", "qaoa", "both"], default="both",
        help="Which method(s) to run this invocation (default: both). "
             "QAOA is slower -- consider running vqe and qaoa separately, "
             "e.g. `--methods vqe` then `--methods qaoa --append`.",
    )
    parser.add_argument(
        "--append", action="store_true",
        help="Append to existing results/statistical_benchmark*.csv instead of "
             "overwriting (use when running --methods vqe and --methods qaoa "
             "as two separate invocations).",
    )
    args = parser.parse_args()

    methods = ["vqe", "qaoa"] if args.methods == "both" else [args.methods]

    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f"===== Statistical benchmark: methods={methods} =====\n")
    df = run_sweep(methods)

    if args.append and os.path.exists(TRIALS_CSV):
        existing = pd.read_csv(TRIALS_CSV)
        # drop any previous rows for the method(s) just re-run, so re-running
        # doesn't create duplicate trial rows
        existing = existing[~existing["method"].isin(methods)]
        df = pd.concat([existing, df], ignore_index=True)

    df.to_csv(TRIALS_CSV, index=False)
    print(f"\nWrote {len(df)} trial rows -> {TRIALS_CSV}")

    summary = summarize(df)
    summary.to_csv(SUMMARY_CSV, index=False)
    print(f"Wrote {len(summary)} summary rows -> {SUMMARY_CSV}\n")

    print(summary.to_string(index=False))
