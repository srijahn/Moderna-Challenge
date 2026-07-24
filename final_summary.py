import pandas as pd

# --- Original approach ------------------------------------------------------
# Everything below was hardcoded print statements, not computed from any
# actual run (e.g. "Best Energy Found: -3.0", "RNA Length 20 -> 50
# Variables" from formula-based guesses). Left here for reference only:
#
# print("ViennaRNA Reference Energy:") ; print("-4.0")
# print("Best Energy Found:") ; print("-3.0")
# print("Gamma = 0.7") ; print("Beta  = 0.3")
# print("Success Probability") ; print("0.4796")
# print("Energy Gap = 1.0") ; print("Approx Accuracy = 75%")
# print("RNA Length 20 -> 50 Variables")
# print("RNA Length 30 -> 112 Variables")
# print("RNA Length 40 -> 200 Variables")
# print("RNA Length 50 -> 312 Variables")

# Reads the real numbers written by generate_final_results.py (run that
# script first) plus the measured scaling tables from
# scaling_analysis_real.py / cvar_vqe_scaling_analysis.py.
try:
    results = pd.read_csv("results/final_results.csv").set_index("Metric")["Value"]
except FileNotFoundError:
    raise SystemExit(
        "results/final_results.csv not found -- run generate_final_results.py first."
    )

print("\n")
print("=" * 60)
print("WISER x MODERNA RNA FOLDING PROJECT")
print("=" * 60)

print("\nCLASSICAL BENCHMARK")
print("Sequence:")
print(results["Sequence"])
print("ViennaRNA Reference Energy:")
print(results["ViennaRNA MFE"])

print("\nQUANTUM RESULTS")
print("QAOA QUBO Energy:")
print(f"{results['QAOA QUBO Energy']}  (success probability {results['QAOA Success Probability']})")
print("CVaR-VQE QUBO Energy:")
print(f"{results['CVaR-VQE QUBO Energy']}  (success probability {results['CVaR-VQE Success Probability']}, "
      f"circuit depth {results['CVaR-VQE Circuit Depth']})")
print("Best Quantum Candidate Energy (better of the two methods):")
print(results["Best Quantum QUBO Energy"])

print("\nCOMPARISON")
print(f"Energy Gap (vs ViennaRNA MFE) = {results['Energy Gap (vs ViennaRNA MFE)']}")
print(f"Approx Accuracy = {results['Approx Accuracy (%)']}%")

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
print("simplified toy model. Resource requirements increase rapidly with")
print("RNA sequence length -- see the scaling table above.")

print("\n" + "=" * 60)
