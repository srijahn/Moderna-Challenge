# Fundamentals of RNA Secondary Structure and MFE Prediction

## 1. RNA Composition & Base Pairing Rules
RNA is a single-stranded molecule made up of four primary chemical bases:
* **A** = Adenine
* **U** = Uracil
* **C** = Cytosine
* **G** = Guanine

Unlike DNA, single-stranded RNA folds back on itself to form 2D shapes (secondary structures) driven by hydrogen bonding between valid base pairs.

### Valid Base Pairs:
* **Canonical Watson-Crick Pairs:**
  * **A - U** (and **U - A**): Connected by **2 hydrogen bonds** (double bond).
  * **C - G** (and **G - C**): Connected by **3 hydrogen bonds** (triple bond, making C-G interactions significantly stronger and more thermodynamically stable).
* **Non-Canonical Wobble Pairs:**
  * **G - U** (and **U - G**): Connected by **2 hydrogen bonds** are weaker than C-G, but highly significant in thermodynamic structural stability.

---

## 2. Minimum Free Energy (MFE) & Structure Stability
RNA sequences spontaneously fold into specific secondary structures based on thermodynamic stability.

* **The Most Stable Structure:** The most stable structure adopted by an RNA sequence is called the **Minimum Free Energy (MFE)** structure.
* **Thermodynamic Basis:** Out of all possible folded configurations, the MFE structure minimizes the overall Gibbs Free Energy ($\Delta G$). The lower (more negative) the $\Delta G$, the more stable the fold.
* **ViennaRNA Benchmark:** **ViennaRNA** is the industry-standard software package that predicts MFE secondary structures using Turner's Nearest-Neighbor thermodynamic energy model.

---

## 3. Dot-Bracket Notation
Secondary structures are represented computationally using **Dot-Bracket Notation**, where string characters denote individual nucleotide states:

* `.` = **unpaired nucleotide**
* `(` = **opening pair** (nucleotide paired to a downstream position)
* `)` = **closing pair** (nucleotide paired to an upstream position)

### Standard Notation Examples:

1. **Unfolded Sequence (10 nt):**
   * **Sequence:**  `ACGUACGUAC`
   * **Structure:** `..........`
   * *Meaning:* All 10 nucleotides remain single-stranded and unpaired.

2. **Simple Stem-Loop / Hairpin (10 nt):**
   * **Sequence:**  `G C G C A U G C G C`
   * **Structure:** `( ( ( . . . . ) ) )`
   * *Meaning:*
     * Positions 1, 2, and 3 (`G C G`) pair with positions 10, 9, and 8 (`C G C`) via opening and closing brackets.
     * Positions 4, 5, 6, and 7 (`C A U G`) remain unpaired (`.`), forming a 4-nucleotide hairpin loop (note: hairpin loops require at least 3 unpaired bases due to physical steric constraints).

---

## 4. Why MFE Matters for mRNA Therapeutics
Predicting and controlling MFE is critical when designing mRNA vaccines and therapeutics:

1. **In Vivo Stability:** Unstructured single-stranded RNA regions (`.`) are rapidly degraded by cellular RNases. A stable MFE structure protects the mRNA and extends its half-life in target cells.
2. **Translation Efficiency:** If the 5' region is too tightly folded (extremely low $\Delta G$), cellular ribosomes cannot attach or initiate translation. Designing a balanced MFE ensures efficient protein expression.
3. **Manufacturability & Storage:** Unstable or incorrectly folded mRNA molecules tend to aggregate during manufacturing or lipid nanoparticle (LNP) encapsulation. Stable MFE designs ensure structural consistency across storage conditions.
