import pennylane as qml
from pennylane import numpy as np

n_qubits = 2

dev = qml.device("default.qubit", wires=n_qubits)

@qml.qnode(dev)
def circuit(gamma, beta):

    for i in range(n_qubits):
        qml.Hadamard(wires=i)

    qml.CNOT(wires=[0,1])
    qml.RZ(2*gamma, wires=1)
    qml.CNOT(wires=[0,1])

    for i in range(n_qubits):
        qml.RX(2*beta, wires=i)

    return qml.probs(wires=[0,1])


gammas = [0.1,0.3,0.5,0.7,1.0]
betas = [0.1,0.3,0.5,0.7,1.0]

best_prob = 0
best_params = None

for gamma in gammas:
    for beta in betas:

        probs = circuit(gamma,beta)

        max_prob = np.max(probs)

        if max_prob > best_prob:
            best_prob = max_prob
            best_params = (gamma,beta)

print("Best Gamma, Beta:")
print(best_params)

print("Highest Probability:")
print(best_prob)