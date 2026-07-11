import numpy as np
import random

Q = np.array([
    [-1, 5, 0, 0, 0],
    [5, -1, 0, 0, 0],
    [0, 0, -1, 5, 0],
    [0, 0, 5, -1, 0],
    [0, 0, 0, 0, -1]
])

def energy(x):
    x = np.array(x)
    return x.T @ Q @ x

best_energy = float("inf")
best_solution = None

for _ in range(1000):

    solution = [
        random.randint(0,1)
        for _ in range(5)
    ]

    e = energy(solution)

    if e < best_energy:
        best_energy = e
        best_solution = solution

print("Best Solution:")
print(best_solution)

print("Best Energy:")
print(best_energy)