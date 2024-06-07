import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

# TODO: Ordne computing thingies

# Computed greier

# Present (e.g. make a table) the times (both x, redshift z and t) for when we have 
# 1) radiation-matter equality 
# 2) matter-dark energy equality
# 3) when the Universe starts to accelerate.

# Additionally present the age of the Universe today
# (t(x=0)) and the value for conformal time today (η(x=0)/c).

# Note that matter means baryons+CDM and radiation means photons+massless neutrinos.
# You can find the expressions for these times analytically and then use this to compute 
# the numerical values or
# solve numerically if you prefer this. (This)
# Mark these as vertical lines in the plots where relevant (if it's relevant for discussing the plot).
# These values will be useful also in future milestones to understand/explain results. (10p)


# TODO: Fortsett med plots

# Se Main.cpp for å gjøre det, kommentarer

# Plot of the luminosity distance versus redshift 
# with data from supernova observations (taken from Betoule et al. 2014)
# overplotted (see header line for file-format). (5p)

# Plots from MCMC fits to supernova data. 1) the 1σ constraints from MCMC 
# fits to supernova data in the ΩΛ vs ΩM plane 2) the posterior PDF 
# (the histogram of accepted samples) of the Hubble parameter H0. (10p)

# Remember to use sensible units in the plots
# (e.g. 100km/s/Mpc for H, Mpc/Gpc for distances and Gyr for times)
# Dette kan jeg gjøre nå i kveld

# Gi navn
# plot of luminosity distance versus redshift, gi navn

# two plots, 
# the 1σ constraints
# and the posterior PDF

x, eta, time, Hp, dHpdx, ddHpddx, OmegaB, OmegaCDM, OmegaLambda, OmegaR, OmegaNu, OmegaK = np.loadtxt("cosmology.txt", unpack=True)
# chi2, h_fitting, OmegaM_fitting, OmegaK_fitting = np.loadtxt("results_supernovafitting.txt", unpack=True, skiprows=1)
Mpc = const.parsec * const.mega
hecto = const.hecto
kilo = const.kilo
c = const.c

OmegaREL = OmegaR + OmegaNu
OmegaM = OmegaB + OmegaCDM

save_fig = False

# Plots that check that the code works properly

plt.title("(n * Hp) / c")
plt.plot(x, (eta * Hp) / c)
if save_fig:
    plt.savefig(r"Plots\Milestone_1\eta_Hp.pdf")
plt.show()

plt.title("1/Hp * dHpdx")
plt.plot(x, (1/Hp) * dHpdx)
if save_fig:
    plt.savefig(r"Plots\Milestone_1\dHp_dx.pdf")
plt.show()

plt.title("1/Hp * ddHpddx")
plt.plot(x, (1/Hp) * ddHpddx)
if save_fig:
    plt.savefig(r"Plots\Milestone_1\ddHpddx.pdf")
plt.show()

#########

plt.title("Hp")
plt.semilogy(x, Hp * (Mpc/(hecto * kilo)))
if save_fig:
    plt.savefig(r"Plots\Milestone_1\Hp.pdf")
plt.show()

plt.title("eta")
plt.semilogy(x, eta * Mpc / c)
if save_fig:
    plt.savefig(r"Plots\Milestone_1\eta.pdf")
plt.show()

plt.title("Cosmic time")
plt.semilogy(x, time)
if save_fig:
    plt.savefig(r"Plots\Milestone_1\time.pdf")
plt.show()

plt.title("Density parameters")
plt.plot(x, OmegaM, label = "Omega b + Omega CDM")
plt.plot(x, OmegaREL, label = "Omega R + Omega Nu")
plt.plot(x, OmegaLambda, label = "Omega Lambda")
plt.legend()
if save_fig:
    plt.savefig(r"Plots\Milestone_1\Omegas.pdf")
plt.show()