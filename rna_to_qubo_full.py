"""
Full RNA secondary structure -> QUBO formulation.

Builds on rna_to_qubo.py / pair_finder.py / conflict_detector.py, but adds:
  1. Wobble (G-U) pairs, weighted lower than Watson-Crick pairs
  2. Minimum hairpin loop size (>= 3 unpaired nt between paired bases)
  3. Overlap penalty (two pairs sharing a nucleotide can't both be on)
  4. Crossing penalty (pseudoknots excluded -> non-crossing pairs only)

Objective (as a QUBO, to MINIMIZE):
    E(x) = -sum_k  w_k * x_k                         (reward valid pairs)
           + sum_{k<l in conflict} P * x_k * x_l      (penalize overlaps/crossings)

x_k in {0,1} = whether candidate pair k is selected in the structure.
"""

import numpy as np

# ---------------------------------------------------------------------------
# 1. Config
# ---------------------------------------------------------------------------
SEQUENCE = "GGUGCCGAACAGUAGCACUC"

MIN_LOOP_SIZE = 3      # minimum unpaired nt strictly between i and j
WC_WEIGHT = 2.0         # reward for a Watson-Crick pair (A-U, G-C)
WOBBLE_WEIGHT = 1.0     # reward for a wobble pair (G-U) -- weaker, so lower weight
PENALTY = 10.0          # conflict penalty; should be >> max possible reward
                         # so the optimizer never "profits" from breaking a constraint


# ---------------------------------------------------------------------------
# 2. Generate candidate pairs (Watson-Crick + wobble, respecting loop size)
# ---------------------------------------------------------------------------
def get_candidate_pairs(sequence, min_loop_size=MIN_LOOP_SIZE):
    """Returns a list of (i, j, weight) for every allowed candidate pair."""
    wc_pairs = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}
    wobble_pairs = {("G", "U"), ("U", "G")}

    candidates = []
    n = len(sequence)

    for i in range(n):
        for j in range(i + 1, n):
            # enforce minimum hairpin loop size
            if j - i <= min_loop_size:
                continue

            a, b = sequence[i], sequence[j]

            if (a, b) in wc_pairs:
                candidates.append((i, j, WC_WEIGHT))
            elif (a, b) in wobble_pairs:
                candidates.append((i, j, WOBBLE_WEIGHT))

    return candidates


# ---------------------------------------------------------------------------
# 3. Conflict detection: overlap + crossing
# ---------------------------------------------------------------------------
def pairs_overlap(p1, p2):
    """True if the two candidate pairs share a nucleotide index."""
    i1, j1, _ = p1
    i2, j2, _ = p2
    return len({i1, j1} & {i2, j2}) > 0


def pairs_cross(p1, p2):
    """True if the two candidate pairs form a pseudoknot (crossing pattern)."""
    i1, j1, _ = p1
    i2, j2, _ = p2
    # standard non-crossing condition check: pairs cross if exactly one endpoint
    # of p2 falls strictly inside the span of p1
    return (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)


# ---------------------------------------------------------------------------
# 4. Build the QUBO matrix
# ---------------------------------------------------------------------------
def build_qubo(candidates, penalty=PENALTY):
    """
    candidates: list of (i, j, weight)
    returns: Q (numpy array, shape [n_vars, n_vars])
             such that energy(x) = x^T Q x
    """
    n_vars = len(candidates)
    Q = np.zeros((n_vars, n_vars))

    # reward term: -weight on the diagonal (minimizing -weight*x rewards x=1)
    for k, (i, j, w) in enumerate(candidates):
        Q[k, k] = -w

    # penalty terms: overlap or crossing between pair k and pair l
    for k in range(n_vars):
        for l in range(k + 1, n_vars):
            if pairs_overlap(candidates[k], candidates[l]) or \
               pairs_cross(candidates[k], candidates[l]):
                # split penalty across Q[k,l] and Q[l,k] so x^T Q x adds it once
                Q[k, l] += penalty / 2
                Q[l, k] += penalty / 2

    return Q


def energy(x, Q):
    x = np.array(x)
    return float(x @ Q @ x)


# ---------------------------------------------------------------------------
# 5. Brute-force solver (only feasible for small n_vars, e.g. the 10 nt test
#    sequence from the plan -- use this to validate against ViennaRNA MFE)
# ---------------------------------------------------------------------------
def brute_force_solve(Q, max_vars_for_bruteforce=20):
    n_vars = Q.shape[0]
    if n_vars > max_vars_for_bruteforce:
        raise ValueError(
            f"{n_vars} variables is too many for brute force "
            f"(2^{n_vars} combinations). Use a QAOA/annealing solver instead."
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
# 6. Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    candidates = get_candidate_pairs(SEQUENCE)
    print(f"Sequence: {SEQUENCE}")
    print(f"Candidate pairs after loop-size filter: {len(candidates)}")
    for k, (i, j, w) in enumerate(candidates):
        kind = "WC" if w == WC_WEIGHT else "wobble"
        print(f"  x{k}: ({i},{j}) {SEQUENCE[i]}-{SEQUENCE[j]} [{kind}, w={w}]")

    Q = build_qubo(candidates)
    print(f"\nQUBO matrix shape: {Q.shape}")

    if len(candidates) <= 20:
        best_x, best_e = brute_force_solve(Q)
        selected = [candidates[k] for k, bit in enumerate(best_x) if bit == 1]
        print(f"\nBrute-force best energy: {best_e}")
        print("Selected pairs:")
        for i, j, w in selected:
            print(f"  ({i},{j}) {SEQUENCE[i]}-{SEQUENCE[j]}")
    else:
        print("\nToo many candidate pairs for brute force on this sequence "
              "-- try a shorter test sequence (~10 nt) to validate against "
              "ViennaRNA MFE, per the Week 1 plan.")