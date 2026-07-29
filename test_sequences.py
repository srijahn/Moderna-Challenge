"""
Curated RNA test sequences for validating the quantum solvers.

These replace the project's original hardcoded 10 nt / 15 nt examples,
which (per structure_metrics.py's base-pair metrics) had a fully-unpaired
real ViennaRNA MFE structure -- i.e. they don't actually fold, so
"validating" the quantum solvers against them was really just validating
against the QUBO's own brute-force optimum, not against real thermodynamics.

Selected via select_test_sequences.py, which searches deterministic random
sequences for ones that (a) ViennaRNA actually folds into a non-trivial
structure, and (b) stay within a candidate-pair budget that brute-force can
still validate exactly. Both entries below are single-hairpin structures
(the simplest possible non-trivial fold), chosen with clean, unambiguous
ViennaRNA MFE structures so a "did the quantum solver find the right
structure" comparison is meaningful.

Every script that needs a short validation sequence should import from
here, so the whole project is validated against the same, known-folding
sequences instead of each script picking its own arbitrary example.
"""

# 10 nt: found by select_test_sequences.py, seed=439.
# ViennaRNA MFE: "(((....)))" at -2.60 kcal/mol -- a clean 3-bp stem with a
# 4 nt loop. 11 candidate pairs (qubits) -- brute-force validates in ~0.01s.
TEST_SEQUENCE_10NT = "GCCGCGCGGC"

# 12 nt: found by select_test_sequences.py, seed=1191.
# ViennaRNA MFE: "(((.....)))." at -4.10 kcal/mol -- a clean 3-bp stem with
# a 5 nt loop. 11 candidate pairs (qubits) -- brute-force validates in
# ~0.01s.
#
# Note: a 15 nt folding sequence found by the same search (seed=1662,
# "CCCGGAAAUAGCGGA") was tried first, but at 17 candidate-pair qubits, full
# QAOA (150 steps x 2 restarts, backprop diff on a 17-qubit statevector)
# took ~17 minutes locally -- too slow to be a practical validation case.
# This 12 nt / 11-qubit sequence gives a genuinely different length from
# the 10 nt case above while keeping a full QAOA run under ~90s.
TEST_SEQUENCE_12NT = "GCCAAAUGGGCG"


if __name__ == "__main__":
    import RNA

    for name, seq in [("10 nt", TEST_SEQUENCE_10NT), ("12 nt", TEST_SEQUENCE_12NT)]:
        structure, mfe = RNA.fold(seq)
        status = "folds" if structure.count("(") > 0 else "DOES NOT FOLD (!)"
        print(f"{name}: {seq}  ->  {structure}  ({mfe:.2f} kcal/mol)  [{status}]")
