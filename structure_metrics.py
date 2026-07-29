"""
Shared structure-comparison metrics.

The old comparison scripts (compare_to_vienna.py, batch_accuracy.py,
generate_final_results.py) computed an "energy gap" and "accuracy" by
comparing the internal QUBO energy (an arbitrary pair-count reward scheme --
see rna_to_qubo_full.py's WC_WEIGHT / WOBBLE_WEIGHT / PENALTY) directly
against ViennaRNA's MFE, which is a real thermodynamic free energy in
kcal/mol. Those two numbers are not on the same scale, so the resulting
"accuracy %" was a rough proxy at best -- flagged honestly in code comments,
but not a real benchmark result.

This module provides the two things Task 3 of the challenge actually asks
for:

  1. real_energy(sequence, dot_bracket) -- evaluate a *candidate* structure
     with ViennaRNA's own thermodynamic model, so it's directly comparable
     (same units, same model) to RNA.fold()'s MFE.
  2. base_pair_metrics(predicted, reference) -- precision/recall/F1 over
     base pairs, base-pair distance, and Hamming distance between two
     dot-bracket strings, so "did the quantum solver find (approximately)
     the right structure" can be reported structurally, not just by energy.
"""

import RNA


def real_energy(sequence, dot_bracket):
    """Evaluate the real ViennaRNA free energy (kcal/mol) of a candidate
    dot-bracket structure for the given sequence. Directly comparable to
    RNA.fold(sequence)'s MFE, since both come from the same ViennaRNA
    thermodynamic model (nearest-neighbor energy parameters) -- unlike the
    internal QUBO energy, which is an arbitrary reward/penalty scheme used
    only to steer the quantum optimizers.
    """
    fc = RNA.fold_compound(sequence)
    return fc.eval_structure(dot_bracket)


def dot_bracket_to_pairs(dot_bracket):
    """Dot-bracket string -> set of (i, j) index pairs, i < j.

    Only handles the non-crossing '(' / ')' / '.' alphabet used by this
    project (pseudoknots are excluded by construction in
    rna_to_qubo_full.py, so '[', ']', etc. are not expected here).
    """
    stack = []
    pairs = set()
    for idx, ch in enumerate(dot_bracket):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise ValueError(f"Unbalanced dot-bracket string: {dot_bracket!r}")
            i = stack.pop()
            pairs.add((i, idx))
        elif ch != ".":
            raise ValueError(
                f"Unexpected character {ch!r} in dot-bracket string "
                f"(only '(', ')', '.' are supported)"
            )
    if stack:
        raise ValueError(f"Unbalanced dot-bracket string: {dot_bracket!r}")
    return pairs


def base_pair_metrics(predicted_db, reference_db):
    """Structural comparison between a predicted structure and a reference
    structure (typically the ViennaRNA MFE structure), both in dot-bracket
    notation over the same sequence.

    Returns a dict with:
      precision            -- fraction of predicted pairs also in reference
      recall               -- fraction of reference pairs also predicted
      f1                   -- harmonic mean of precision and recall
      base_pair_distance   -- |predicted \\ reference| + |reference \\ predicted|
                              (the standard RNA base-pair distance metric;
                              0 means an exact structural match)
      hamming_distance     -- number of dot-bracket character positions
                              that differ (paired vs. unpaired mismatches,
                              including pairs shifted to a different partner)
      n_predicted_pairs / n_reference_pairs -- for context
    """
    if len(predicted_db) != len(reference_db):
        raise ValueError("Structures must be the same length to compare.")

    pred_pairs = dot_bracket_to_pairs(predicted_db)
    ref_pairs = dot_bracket_to_pairs(reference_db)

    tp = len(pred_pairs & ref_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)

    if pred_pairs:
        precision = tp / len(pred_pairs)
    else:
        precision = 1.0 if not ref_pairs else 0.0

    if ref_pairs:
        recall = tp / len(ref_pairs)
    else:
        recall = 1.0 if not pred_pairs else 0.0

    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "base_pair_distance": fp + fn,
        "hamming_distance": sum(a != b for a, b in zip(predicted_db, reference_db)),
        "n_predicted_pairs": len(pred_pairs),
        "n_reference_pairs": len(ref_pairs),
    }


if __name__ == "__main__":
    # Small sanity check using the Task 3 example from the challenge PDF.
    sequence = "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG"
    mfe_structure, mfe = RNA.fold(sequence)

    candidate = ".................(((....))).................."[: len(sequence)]

    print(f"Sequence:          {sequence}")
    print(f"ViennaRNA MFE:     {mfe_structure}  ({mfe:.2f} kcal/mol)")
    print(f"Candidate:         {candidate}")
    print(f"Candidate energy:  {real_energy(sequence, candidate):.2f} kcal/mol")
    print(f"Metrics:           {base_pair_metrics(candidate, mfe_structure)}")
