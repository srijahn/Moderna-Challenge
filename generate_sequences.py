import random
import pandas as pd

BASES = ['A', 'C', 'G', 'U']

# Curated, handpicked sequences targeting specific structural motifs
HANDPICKED_SEQUENCES = [
    {
        "sequence_id": "HP_01_Classic_Hairpin",
        "sequence": "GCGCAUGCGC",
        "length": 10,
        "type": "Handpicked - Hairpin"
    },
    {
        "sequence_id": "HP_02_Poly_CG_Stem",
        "sequence": "CGCGCGCGAAAAACGCGCGCG",
        "length": 21,
        "type": "Handpicked - High GC"
    },
    {
        "sequence_id": "HP_03_Wobble_Rich",
        "sequence": "GUGUGUGUAAAAAACACACACA",
        "length": 22,
        "type": "Handpicked - G-U Wobble"
    },
    {
        "sequence_id": "HP_04_Weak_AU_Loop",
        "sequence": "AUAUAUAUAAAAAAAUAUAUAU",
        "length": 22,
        "type": "Handpicked - Weak AU"
    },
    {
        "sequence_id": "HP_05_MicroRNA_let7",
        "sequence": "UGAGGUAGUAGGUUGUAUAGUU",
        "length": 22,
        "type": "Handpicked - Biological"
    },
    {
        "sequence_id": "HP_06_Variable_30nt",
        "sequence": "GGCCCCUUACCAAAGGGGCCAAAAAAAGGG",
        "length": 30,
        "type": "Handpicked - 30nt Multi-stem"
    }
]

def generate_random_sequence(min_len=15, max_len=35):
    """Generates an RNA sequence with random length between min_len and max_len."""
    length = random.randint(min_len, max_len)
    return ''.join(random.choice(BASES) for _ in range(length))

def build_dataset(num_synthetic=20):
    dataset = []
    
    # 1. Add handpicked curated sequences
    for item in HANDPICKED_SEQUENCES:
        dataset.append({
            'sequence_id': item['sequence_id'],
            'sequence': item['sequence'],
            'length': item['length'],
            'type': item['type']
        })
    
    # 2. Add variable-length synthetic sequences (15 to 35 nt)
    for i in range(1, num_synthetic + 1):
        seq = generate_random_sequence(min_len=15, max_len=35)
        dataset.append({
            'sequence_id': f"SYN_{i:02d}",
            'sequence': seq,
            'length': len(seq),
            'type': 'Synthetic'
        })
        
    df = pd.DataFrame(dataset)
    df.to_csv('data/sequences.csv', index=False)
    print(f"Dataset successfully created in 'data/sequences.csv' with {len(df)} total sequences.")

if __name__ == "__main__":
    build_dataset()
