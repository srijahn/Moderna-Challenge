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

conflicts = []

for p1 in range(len(pairs)):
    for p2 in range(p1 + 1, len(pairs)):

        a1, b1 = pairs[p1]
        a2, b2 = pairs[p2]

        overlap = (
            a1 == a2 or
            a1 == b2 or
            b1 == a2 or
            b1 == b2
        )

        if overlap:
            conflicts.append((p1, p2))

print("Number of conflicts:", len(conflicts))

print("\nFirst 20 conflicts:")

for c in conflicts[:20]:
    print(c)