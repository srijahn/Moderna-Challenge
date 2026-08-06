"""
Classical-only benchmark on the official example sequence given in the
challenge brief itself (Moderna - WISER Quantum Challenge [SHARED].pdf,

    GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG  (44 nt)

This is NOT run through QAOA or CVaR-VQE. get_candidate_pairs() on this
sequence returns far more candidate pairs than the project's measured
QAOA feasibility ceiling. Truncating the sequence to fit would mean this
is no longer actually "the official example," so instead this script only
does the classical half: ViennaRNA MFE structure + energy.

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
