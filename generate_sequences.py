import random
import csv

bases = ['A', 'U', 'C', 'G']

with open('data/sequences.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    writer.writerow(['sequence'])

    for i in range(50):
        seq = ''.join(random.choice(bases) for _ in range(20))
        writer.writerow([seq])

print("50 RNA sequences generated!")