"""
Classical-only benchmark on the official example sequence given in the
challenge brief itself (Moderna - WISER Quantum Challenge [SHARED].pdf,
Task 2: "Classical benchmark generation"):

    GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG  (44 nt)

This is NOT run through QAOA or CVaR-VQE. get_candidate_pairs() on this
sequence returns far more candidate pairs than this project's measured
QAOA feasibility ceiling (the README's runtime finding: ~17 min at 17
qubits on local statevector simulation, and that's before CVaR-VQE's
extra circuit-depth cost). Truncating the sequence to fit would mean this
is no longer actually "the official example," so instead this script only
does the classical half explicitly asked for in Task 2/3: ViennaRNA MFE
structure + energy, reported honestly as classical-only.

If you want a quantum-comparable run on real, official-example-scale RNA,
that requires either (a) a smaller/truncated fragment (explicitly labeled
as a fragment, not the full example) or (b) a QUBO decomposition /
window-based approach beyond this project's current scope -- see
findings_week1.md for the same tradeoff already hit with the original 20
nt placeholder sequence.
"""

import RNA

from rna_to_qubo_full import get_candidate_pairs

OFFICIAL_EXAMPLE_SEQUENCE = "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG"

if __name__ == "__main__":
    sequence = OFFICIAL_EXAMPLE_SEQUENCE
    structure, mfe = RNA.fold(sequence)
    candidates = get_candidate_pairs(sequence)

    print("========== Official Challenge-Brief Example Sequence ==========\n")
    print(f"Sequence ({len(sequence)} nt):")
    print(sequence)
    print("\nViennaRNA MFE Structure:")
    print(structure)
    print("\nViennaRNA MFE Energy (kcal/mol):")
    print(mfe)
    print(f"\nCandidate base pairs (qubits QAOA/CVaR-VQE would need): {len(candidates)}")
    print(
        "\nThis is well past this project's measured QAOA feasibility ceiling "
        "(see README scaling finding), so QAOA/CVaR-VQE are intentionally not "
        "run on this sequence here -- see this file's module docstring."
    )
