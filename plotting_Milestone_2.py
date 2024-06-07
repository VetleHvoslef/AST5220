import numpy as np
import matplotlib.pyplot as plt


# TODO: Ordne en computing python fil

# Computing, legg inn en python fil med alt, printer milestone 1 og gir alt 
# Compute the times (xdecoupling, zdecoupling, tdecoupling) for the last scattering surface surface and for when recombination happens.
# You can use the definition τ=1 or as the time when visibility function peaks for when last scattering happens.
# You can use the definition Xe=0.1 for when recombination happens.

# Compare these numbers with the prediction from the Saha equation. Compute the freeze-out abundance of free electrons today Xe(x=0)
# . If you include reionization compute these numbers with reionization turned off (15p).
# Compute the sound-horizon at decoupling rs (5p).

save_fig = False

x, Xe, Xe_Saha, ne, tau, dtaudx, ddtauddx, g_tilde, dgdx, ddgddx = np.loadtxt("recombination.txt", unpack=True)

plt.semilogy(x, tau, label="tau")
plt.semilogy(x, -dtaudx, label="-dtaudx")
plt.semilogy(x, ddtauddx, label="ddtauddx")
plt.legend()
plt.xlabel("x")
if save_fig:
    plt.savefig(r"Plots\Milestone_2\tau.pdf")
plt.show()

plt.plot(x, g_tilde/np.max(g_tilde), label="g_tilde")
plt.plot(x, dgdx/np.max(dgdx), label="dg_tildedx'")
plt.plot(x, ddgddx/np.max(ddgddx), label="ddg_tildeddx")
plt.xlim(-8, -6)
plt.xlabel("x")
plt.legend()
if save_fig:
    plt.savefig(r"Plots\Milestone_2\g_tilde_and_derivs.pdf")
plt.show()


plt.semilogy(x, Xe, label="Xe")
plt.semilogy(x, Xe_Saha, linestyle="--", label="Xe Saha")
plt.legend()
if save_fig:
    plt.savefig(r"Plots\Milestone_2\Xe.pdf")
plt.show()




