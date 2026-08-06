
from rna_to_qubo_full import get_candidate_pairs

# Each entry: (label, sequence, expected_structure, expected_mfe, seed)
# seed = the select_test_sequences.py deterministic RNG seed that produced
# this sequence, for traceability / re-derivation.
BENCHMARK_SEQUENCES = [
    # label      sequence            structure          mfe     seed
    ("8nt",      "GCGAUAGC",         "((....))",        -0.10,  5524),
    ("9nt",      "GCAGGAGCA",        "((....)).",       -1.00,  681),
    ("10nt_a",   "GCCGCGCGGC",       "(((....)))",      -2.60,  439),
    ("10nt_b",   "CCAGAAAGGA",       "((.....)).",      -0.20,  2060),
    ("10nt_c",   "GCGCAUGCGC",       "(((....)))",      -1.30,  None),
    ("11nt",     "CGGAAGACCGA",      "(((....)))." ,    -2.30,  467),
    ("12nt",     "ACCACAGGUAAA",     "(((...)))...",    -0.30,  1949),
    ("13nt",     "GCGAAAAAUUCGA",    ".((((...)))).",   -1.10,  3219),
    ("14nt",     "AACCACAGGUAAAA",   ".(((...)))....",  -0.60,  1461),
    ("12nt_wobble", "CGCUAAAGGCGA",  "((((...))))." ,   -3.10,  3075),
    ("12nt_au",   "AAUGCUAAUGCA",    "..(((....)))",    -1.80,  5902),
    ("13nt_gc",   "GGCCAAGACGGCC",   "((((.....))))",   -5.80,  3402),
]


def gc_content(sequence):
    return 100.0 * sum(c in "GC" for c in sequence) / len(sequence)
 
 
if __name__ == "__main__":
    import RNA
 
    print("===== Benchmark sequence sanity check (vs. RNA.fold) =====\n")
    print(f"{'label':8} {'sequence':16} {'structure':16} {'mfe':>7}  "
          f"{'qubits':>6}  {'gc%':>5}  status")
 
    all_ok = True
    for label, seq, expected_structure, expected_mfe, seed in BENCHMARK_SEQUENCES:
        structure, mfe = RNA.fold(seq)
        n_qubits = len(get_candidate_pairs(seq))
        gc = gc_content(seq)
 
        folds = structure.count("(") > 0
        matches_expected = (structure == expected_structure and abs(mfe - expected_mfe) < 1e-6)
        status = "OK" if (folds and matches_expected) else "MISMATCH -- rerun select_test_sequences.py"
        if not (folds and matches_expected):
            all_ok = False
 
        print(f"{label:8} {seq:16} {structure:16} {mfe:7.2f}  "
              f"{n_qubits:6d}  {gc:5.1f}  {status}")
 
    print()
    if all_ok:
        print(f"All {len(BENCHMARK_SEQUENCES)} sequences confirmed folding and "
              f"match recorded structure/MFE. Safe to run statistical_benchmark.py.")
    else:
        print("Some sequences did not match the recorded structure/MFE -- this "
              "usually means a different ViennaRNA version/parameter set. "
              "Re-run select_test_sequences.py and update BENCHMARK_SEQUENCES "
              "before trusting statistical_benchmark.py's output.")