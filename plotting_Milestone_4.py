import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const


# TODO: k_eq, compute


save_fig = False

### DATA ###

ell, C_ell = np.loadtxt("cells.txt", unpack=True)
ell_sw, C_ell_sw, = np.loadtxt("cells_sw.txt", unpack=True)
ell_isw, C_ell_isw, = np.loadtxt("cells_isw.txt", unpack=True)
ell_doppler, C_ell_doppler, = np.loadtxt("cells_doppler.txt", unpack=True)
ell_polarization, C_ell_polarization, = np.loadtxt("cells_polarization.txt", unpack=True)
ell_data, C_ell_data, _, _ = np.loadtxt("planck_cell_low.txt", unpack=True)

k, theta_k_2, theta_k_20, theta_k_120, matter_power_spectrum = np.loadtxt("functions_of_k.txt", unpack=True)
_, eta, _, _, _, _, _, _, _, _, _, _ = np.loadtxt("cosmology.txt", unpack=True)
eta_0 = eta[-1]

### Constants ###

Mpc = const.parsec * const.mega
h = const.h

### PLOTS ###

plt.title("C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell, C_ell)
plt.loglog(ell_sw, C_ell_sw, label="SW term, in source function")
plt.loglog(ell_isw, C_ell_isw, label="ISW term, in source function")
plt.loglog(ell_doppler, C_ell_doppler, label="Doppler term, in source function")
plt.loglog(ell_polarization, C_ell_polarization, label="Polarization term, in source function")
plt.loglog(ell_data, C_ell_data, ".", label="Planck 2018")
plt.xlabel("μK^2")
plt.ylabel("ℓ(ℓ+1)Cℓ/(2π)")
plt.legend()
if save_fig:
    plt.savefig(r"Plots\Milestone_4\C_ell.pdf")
plt.show()

plt.plot("Theta_k for ell values")
plt.plot(k * eta_0, theta_k_2, label="ell 2")
plt.plot(k * eta_0, theta_k_20, label="ell 20")
plt.plot(k * eta_0, theta_k_120, label="ell 120")
plt.xlabel("k * eta_0")
plt.legend()
if save_fig:
    plt.savefig(r"Plots\Milestone_4\Theta_k_ells.pdf")
plt.show()

plt.title("C_ell integrand")
plt.plot(k * eta_0, abs(theta_k_2**2) / (k * eta_0), label="ell 2")
plt.plot(k * eta_0, abs(theta_k_20**2) / (k * eta_0), label="ell 20")
plt.plot(k * eta_0, abs(theta_k_120**2) / (k * eta_0), label="ell 120")
plt.legend()
plt.xlabel("k * eta_0")
if save_fig:
    plt.savefig(r"Plots\Milestone_4\C_ell_integrad.pdf")
plt.show()

plt.title("Matter power spectrum")
plt.loglog(k * (h/Mpc), matter_power_spectrum * (Mpc**3/h**3))
plt.xlabel("h/Mpc")
plt.ylabel("Mpc^3/h^3")
if save_fig:
    plt.savefig(r"Plots\Milestone_4\matter_power_spectrum.pdf")
plt.show()