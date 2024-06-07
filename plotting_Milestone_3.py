import numpy as np
import matplotlib.pyplot as plt


save_fig = False

### DATA ###

x = ["", "", ""]
Theta0 = ["", "", ""]
Theta1 = ["", "", ""]
Theta2 = ["", "", ""]
Phi = ["", "", ""]
Psi = ["", "", ""]
delta_cdm = ["", "", ""]
delta_b = ["", "", ""]
v_cdm = ["", "", ""]
v_b = ["", "", ""]

k = [0.001, 0.01, 0.1]
colours = ["b", "orange", "green"]
linestyles = ["-", "--", "-."]

for i, filename in enumerate(["perturbations_k0.001.txt", "perturbations_k0.01.txt", "perturbations_k0.1.txt"]):
    x[i], Theta0[i], Theta1[i], Theta2[i], Phi[i], Psi[i], _, delta_cdm[i], delta_b[i], v_cdm[i], v_b[i], _, _, _, _ = np.loadtxt(filename, unpack=True)

### PLOTS ###

plt.title("delta-ene")
for i in range(3):
    plt.semilogy(x[i], 4 * abs(Theta0[i]), color=colours[i], linestyle=linestyles[0], label=f"delta gamma k = {k[i]}")
    plt.semilogy(x[i], abs(delta_b[i]), color=colours[i], linestyle=linestyles[1], label=f"delta b k = {k[i]}")
    plt.semilogy(x[i], abs(delta_cdm[i]), color=colours[i], linestyle=linestyles[2], label=f"delta cdm k = {k[i]}")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig()
plt.show()


plt.title("v-ene")
for i in range(3):
    plt.semilogy(x[i], 3 * abs(Theta1[i]), color=colours[i], linestyle=linestyles[0], label=f"v gamma, k = {k[i]}")
    plt.semilogy(x[i], abs(v_cdm[i]), color=colours[i], linestyle=linestyles[1], label=f"v cdm, k = {k[i]}")
    plt.semilogy(x[i], abs(v_b[i]), color=colours[i], linestyle=linestyles[2], label=f"v b, k = {k[i]}")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig()
plt.show()


plt.title("Theta2")
for i in range(3):
    plt.plot(x[i], Theta2[i], color=colours[i], label=f"Theta2, k = {k[i]}")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig()
plt.show()

plt.title("Psi")
for i in range(3):
    plt.plot(x[i], Phi[i], color=colours[i], label=f"Psi, k = {k[i]}")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig()
plt.show()

plt.title("Phi + Psi")
for i in range(3):
    plt.plot(x[i], Psi[i] + Phi[i], color=colours[i], label=f"Psi + Phi, k = {k[i]}")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig()
plt.show()
