# Initial Findings

Sequence Length: 20

Candidate Base Pairs: 51

Conflicting Pair Selections: 225

Observation:

The number of possible pairings grows rapidly with sequence length.

Many candidate pairs cannot coexist because they share nucleotides.

This motivates formulating RNA folding as a constrained optimization problem suitable for QUBO and QAOA approaches.