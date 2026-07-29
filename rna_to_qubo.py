sequence = "GGUGCCGAACAGUAGCACUC"

pairs = []

for i in range(len(sequence)):
    for j in range(i + 1, len(sequence)):

        a = sequence[i]
        b = sequence[j]

        valid = (
            (a == 'A' and b == 'U') or
            (a == 'U' and b == 'A') or
            (a == 'C' and b == 'G') or
            (a == 'G' and b == 'C')
        )

        if valid:
            pairs.append((i, j))

print("Number of candidate pairs:", len(pairs))

for idx, pair in enumerate(pairs):
    print(f"x{idx} = {pair}")