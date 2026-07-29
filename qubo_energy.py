import numpy as np

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

solutions = [
    [0,0,0,0,0],
    [1,0,0,0,0],
    [1,1,0,0,0],
    [1,0,1,0,1],
    [0,0,0,0,1]
]

for s in solutions:
    print(s, "Energy =", energy(s))