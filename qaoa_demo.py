import pennylane as qml
from pennylane import numpy as np

# 2 qubits

dev = qml.device("default.qubit", wires=2)

@qml.qnode(dev)
def circuit(gamma, beta):

    qml.Hadamard(wires=0)
    qml.Hadamard(wires=1)

    # Cost Hamiltonian
    qml.CNOT(wires=[0, 1])
    qml.RZ(2 * gamma, wires=1)
    qml.CNOT(wires=[0, 1])

    # Mixer Hamiltonian
    qml.RX(2 * beta, wires=0)
    qml.RX(2 * beta, wires=1)

    return qml.probs(wires=[0,1])


gamma = 0.5
beta = 0.3

result = circuit(gamma, beta)

print("QAOA Output:")
print(result)