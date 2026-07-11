import matplotlib.pyplot as plt

lengths = [20, 30, 40, 50]

variables = [50, 112, 200, 312]

plt.figure(figsize=(8,5))

plt.plot(lengths, variables, marker="o")

plt.title("RNA Length vs Optimization Variables")

plt.xlabel("RNA Length")

plt.ylabel("Number of Variables")

plt.grid(True)

plt.savefig("results/scaling_plot.png")

plt.show()

print("Graph saved as results/scaling_plot.png")