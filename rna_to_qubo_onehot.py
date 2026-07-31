"""
Second RNA secondary structure -> QUBO encoding: one-hot per-position
pairing variables, instead of the one-variable-per-candidate-pair
("pair-indicator") encoding used by rna_to_qubo_full.py.

Optional advanced TODO item from the README:
    "Try a second QUBO encoding (e.g. one-hot per-position pairing
    variables instead of pair-indicator variables) and compare qubit
    count / constraint-enforcement tradeoffs against the current
    encoding."

-------------------------------------------------------------------------
Encoding comparison
-------------------------------------------------------------------------

Pair-indicator (rna_to_qubo_full.py):
    One binary variable x_k per candidate base pair k = (i, j).
    x_k = 1 means "positions i and j are paired".
    Structural validity (each position paired at most once, no
    pseudoknots) is enforced by a single flat penalty term between every
    *conflicting pair of candidate pairs* (O(n_candidates^2) potential
    penalty terms, decided once per (k, l) at build time).
    Qubit count = n_candidates.

One-hot per-position (this file):
    Two binary variables per candidate pair k = (i, j): y_{i->j} ("i's
    choice is to pair with j") and y_{j->i} ("j's choice is to pair with
    i"). Every position i owns one "choice group" {y_{i->j} : j a
    candidate partner of i}; "i stays unpaired" is the implicit all-zero
    state of that group, so a position is structurally valid whenever
    *at most one* variable in its group is 1 -- enforced by a one-hot
    ("at-most-one") penalty within each group.
    A second penalty term ties the two directions of the same physical
    pair together: y_{i->j} and y_{j->i} must agree (both 0 or both 1),
    enforced by a consistency penalty.
    Qubit count = 2 * n_candidates (every candidate pair now needs a
    variable on *each* side).

So the one-hot encoding roughly doubles the qubit count for the same
sequence, in exchange for expressing "each position pairs with at most
one partner" as an explicit, local one-hot constraint per position
(closer to how one-hot encodings are usually written for combinatorial
problems) rather than as pairwise penalties scattered across the whole
candidate-pair conflict graph. The trade only makes sense if that
locality is worth the extra qubits -- see compare_qubo_encodings.py for
measured numbers.
-------------------------------------------------------------------------
"""

import numpy as np

from rna_to_qubo_full import (
    MIN_LOOP_SIZE,
    PENALTY,
    WC_WEIGHT,
    WOBBLE_WEIGHT,
    get_candidate_pairs,
    pairs_cross,
)

SEQUENCE = "GGUGCCGAACAGUAGCACUC"

ONEHOT_PENALTY = PENALTY        # at-most-one-partner penalty per position
CONSISTENCY_PENALTY = PENALTY   # y_{i->j} == y_{j->i} penalty
CROSSING_PENALTY = PENALTY      # pseudoknot-exclusion penalty


# ---------------------------------------------------------------------------
# 1. Variable layout: two directed variables per candidate pair
# ---------------------------------------------------------------------------
def build_variable_index(candidates):
    """
    candidates: list of (i, j, w) with i < j, as returned by
    get_candidate_pairs().

    Returns:
      var_index: dict (from_pos, to_pos) -> qubit index, one entry for
                 each direction of each candidate pair (so 2 entries per
                 candidate).
      position_groups: dict position -> list of qubit indices that make
                 up that position's one-hot choice group.
      directed_pairs: list of (k, l, i, j, w) where k = var_index[i, j],
                 l = var_index[j, i] -- the two directed qubits standing
                 for the same physical candidate pair (i, j).
    """
    var_index = {}
    position_groups = {}
    directed_pairs = []

    for i, j, w in candidates:
        k = len(var_index)
        var_index[(i, j)] = k
        position_groups.setdefault(i, []).append(k)

        l = len(var_index)
        var_index[(j, i)] = l
        position_groups.setdefault(j, []).append(l)

        directed_pairs.append((k, l, i, j, w))

    return var_index, position_groups, directed_pairs


# ---------------------------------------------------------------------------
# 2. Build the QUBO matrix
# ---------------------------------------------------------------------------
def build_qubo_onehot(
    candidates,
    onehot_penalty=ONEHOT_PENALTY,
    consistency_penalty=CONSISTENCY_PENALTY,
    crossing_penalty=CROSSING_PENALTY,
):
    """
    candidates: list of (i, j, w), same format as rna_to_qubo_full's
                get_candidate_pairs() output.
    returns: (Q, var_index) where Q is the [n_vars, n_vars] QUBO matrix
             (energy(x) = x^T Q x) and var_index maps (from_pos, to_pos)
             -> qubit index, so a solution bitstring can be decoded back
             into base pairs.
    """
    var_index, position_groups, directed_pairs = build_variable_index(candidates)
    n_vars = len(var_index)
    Q = np.zeros((n_vars, n_vars))

    # --- reward term -------------------------------------------------
    # Split the pair's reward weight w across its two directed
    # variables, so that when consistency holds (both variables agree
    # at 1) the total reward for selecting pair (i, j) is still -w,
    # matching the pair-indicator encoding's reward for the same pair.
    for k, l, i, j, w in directed_pairs:
        Q[k, k] += -w / 2
        Q[l, l] += -w / 2

    # --- one-hot ("at most one partner") penalty per position --------
    # Within each position's choice group, penalize any two variables
    # both being 1 (i.e. position i "choosing" two different partners
    # at once).
    for i, group in position_groups.items():
        for a in range(len(group)):
            for b in range(a + 1, len(group)):
                k, l = group[a], group[b]
                Q[k, l] += onehot_penalty / 2
                Q[l, k] += onehot_penalty / 2

    # --- consistency penalty: y_{i->j} must equal y_{j->i} ------------
    # (y_k - y_l)^2 = y_k + y_l - 2*y_k*y_l  (using y^2 = y for binary y)
    for k, l, i, j, w in directed_pairs:
        Q[k, k] += consistency_penalty
        Q[l, l] += consistency_penalty
        Q[k, l] += -consistency_penalty
        Q[l, k] += -consistency_penalty

    # --- crossing (pseudoknot-exclusion) penalty ----------------------
    # Applied once per physical pair-of-pairs, on the canonical (i < j)
    # directed variable of each side -- consistency already ties that
    # variable to its mirror on the other side, so penalizing the
    # canonical direction is enough to discourage crossing pairs from
    # both being selected.
    for a in range(len(directed_pairs)):
        k_a, l_a, i_a, j_a, w_a = directed_pairs[a]
        for b in range(a + 1, len(directed_pairs)):
            k_b, l_b, i_b, j_b, w_b = directed_pairs[b]
            if pairs_cross((i_a, j_a, w_a), (i_b, j_b, w_b)):
                Q[k_a, k_b] += crossing_penalty / 2
                Q[k_b, k_a] += crossing_penalty / 2

    return Q, var_index


def energy(x, Q):
    x = np.array(x)
    return float(x @ Q @ x)


# ---------------------------------------------------------------------------
# 3. Decode a solution bitstring back into selected base pairs
# ---------------------------------------------------------------------------
def decode_solution(x, var_index):
    """
    Returns the set of (i, j) with i < j for which BOTH directed
    variables y_{i->j} and y_{j->i} are 1 (i.e. the consistency
    constraint actually held in this solution). Directed picks that
    disagree (only one side "chose" the pair) are not counted as a
    selected pair -- they indicate the consistency penalty wasn't
    satisfied, which brute_force_solve_onehot / a well-tuned annealer
    should avoid whenever it's cheaper to satisfy the constraint.
    """
    selected = set()
    for (from_pos, to_pos), idx in var_index.items():
        if from_pos < to_pos and x[idx] == 1:
            mirror_idx = var_index.get((to_pos, from_pos))
            if mirror_idx is not None and x[mirror_idx] == 1:
                selected.add((from_pos, to_pos))
    return selected


# ---------------------------------------------------------------------------
# 4. Brute-force solver (only feasible for small n_vars)
# ---------------------------------------------------------------------------
def brute_force_solve_onehot(Q, max_vars_for_bruteforce=24):
    n_vars = Q.shape[0]
    if n_vars > max_vars_for_bruteforce:
        raise ValueError(
            f"{n_vars} variables is too many for brute force "
            f"(2^{n_vars} combinations)."
        )

    best_energy = float("inf")
    best_x = None
    for bits in range(2 ** n_vars):
        x = [(bits >> k) & 1 for k in range(n_vars)]
        e = energy(x, Q)
        if e < best_energy:
            best_energy = e
            best_x = x

    return best_x, best_energy


# ---------------------------------------------------------------------------
# 5. Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    candidates = get_candidate_pairs(SEQUENCE, min_loop_size=MIN_LOOP_SIZE)
    print(f"Sequence: {SEQUENCE}")
    print(f"Candidate pairs: {len(candidates)}")

    Q, var_index = build_qubo_onehot(candidates)
    print(f"One-hot QUBO shape: {Q.shape}  "
          f"(pair-indicator encoding would use {len(candidates)} qubits)")

    if Q.shape[0] <= 24:
        best_x, best_e = brute_force_solve_onehot(Q)
        selected = decode_solution(best_x, var_index)
        print(f"\nBrute-force best energy: {best_e}")
        print("Selected pairs (consistency-satisfying only):")
        for i, j in sorted(selected):
            print(f"  ({i},{j}) {SEQUENCE[i]}-{SEQUENCE[j]}")
    else:
        print("\nToo many qubits for brute force on this sequence -- "
              "use compare_qubo_encodings.py's smaller test sequences instead.")
