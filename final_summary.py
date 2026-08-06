
import pandas as pd

# reads one row per benchmark sequence written by generate_final_results.py
# and prints a per-sequence section plus an aggregate summary.

try:
    results = pd.read_csv("results/final_results.csv")
except FileNotFoundError:
    raise SystemExit(
        "results/final_results.csv not found -- run generate_final_results.py first."
    )

print("\n")
print("=" * 60)
print("WISER x MODERNA RNA FOLDING PROJECT")
print("=" * 60)

for _, row in results.iterrows():
    print(f"\n{'-' * 60}")
    print(f"{row['Label']}: {row['Sequence']}")
    print("-" * 60)

    print("\nCLASSICAL BENCHMARK")
    print("ViennaRNA MFE Structure:")
    print(row["ViennaRNA MFE Structure"])
    print("ViennaRNA MFE Energy (kcal/mol):")
    print(row["ViennaRNA MFE Energy (kcal/mol)"])

    print("\nQUANTUM RESULTS (internal QUBO energy -- NOT directly comparable to")
    print("the ViennaRNA kcal/mol values above; see COMPARISON section below for")
    print("the real, same-units comparison)")
    print("QAOA QUBO Energy:")
    print(f"{row['QAOA QUBO Energy (internal)']}  (success probability {row['QAOA Success Probability']})")
    print("CVaR-VQE QUBO Energy:")
    print(f"{row['CVaR-VQE QUBO Energy (internal)']}  (success probability {row['CVaR-VQE Success Probability']}, "
          f"circuit depth {row['CVaR-VQE Circuit Depth']})")
    print(f"Best quantum method (lower internal QUBO energy): {row['Best Quantum Method']}")

    print("\nCOMPARISON (real ViennaRNA thermodynamic model, via")
    print("RNA.fold_compound(sequence).eval_structure() -- same units as the MFE)")
    print("Best Quantum Structure:")
    print(row["Best Quantum Structure"])
    print(f"Best Quantum Real Energy = {row['Best Quantum Real Energy (kcal/mol)']} kcal/mol")
    print(f"Real Energy Gap vs. MFE = {row['Real Energy Gap vs. MFE (kcal/mol)']} kcal/mol  (0 = exact energy match)")
    print(f"Base-Pair Precision / Recall / F1 = {row['Base-Pair Precision']} / "
          f"{row['Base-Pair Recall']} / {row['Base-Pair F1']}")
    print(f"Base-Pair Distance = {row['Base-Pair Distance']}  (0 = exact structural match)")
    print(f"Hamming Distance = {row['Hamming Distance']}  (0 = exact dot-bracket string match)")

print("\n" + "=" * 60)
print(f"AGGREGATE (across all {len(results)} benchmark sequences)")
print("=" * 60)
print(f"Mean Real Energy Gap vs. MFE (kcal/mol): {results['Real Energy Gap vs. MFE (kcal/mol)'].mean():.2f}")
print(f"Mean Base-Pair F1: {results['Base-Pair F1'].mean():.3f}")
print(f"Mean Base-Pair Distance: {results['Base-Pair Distance'].mean():.2f}")
print(f"QAOA won on internal QUBO energy: {(results['Best Quantum Method'] == 'QAOA').sum()}/{len(results)} sequences")

print("\nSCALING ANALYSIS (measured -- see scaling_analysis_real.py / cvar_vqe_scaling_analysis.py)")
try:
    qaoa_scaling = pd.read_csv("results/scaling_analysis_real.csv")
    for length in [10, 20, 30, 40, 50]:
        row = qaoa_scaling[qaoa_scaling["sequence_length"] == length]
        if not row.empty:
            n_qubits = row.iloc[0]["n_qubits"]
            print(f"RNA Length {length} -> {int(n_qubits)} qubits (QAOA)")
except FileNotFoundError:
    print("(results/scaling_analysis_real.csv not found -- run scaling_analysis_real.py first)")

print("\nCONCLUSION")
print("QAOA and CVaR-VQE both ran end-to-end on the real RNA-folding QUBO")
print("(wobble pairs, loop-size, overlap and crossing constraints), not a")
print("simplified toy model, across all 8 curated benchmark sequences")
print("(benchmark_sequences.py, 8-14 nt, 35.7%-100% GC content) -- not just")
print("a single anecdote sequence -- and are compared to ViennaRNA using")
print("real thermodynamic energy (kcal/mol) and base-pair-level structural")
print("metrics, not just internal QUBO energy. Resource requirements")
print("increase rapidly with RNA sequence length -- see the scaling table above.")

print("\n" + "=" * 60)
