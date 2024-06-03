import numpy as np
import matplotlib.pyplot as plt

ell, C_ell = np.loadtxt("cells.txt", unpack=True)

ell_sw, C_ell_sw, = np.loadtxt("cells_sw.txt", unpack=True)
ell_isw, C_ell_isw, = np.loadtxt("cells_isw.txt", unpack=True)
ell_doppler, C_ell_doppler, = np.loadtxt("cells_doppler.txt", unpack=True)
ell_polarization, C_ell_polarization, = np.loadtxt("cells_polarization.txt", unpack=True)

plt.title("C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell, C_ell)
plt.savefig("C_ell.pdf")
plt.show()

plt.title("SW term, in source function, C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell_sw, C_ell_sw)
plt.savefig("SW_term.pdf")
plt.show()

plt.title("ISW term, in source function, C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell_isw, C_ell_isw)
plt.savefig("ISW_term.pdf")
plt.show()

plt.title("doppler term, in source function, C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell_doppler, C_ell_doppler)
plt.savefig("doppler_term.pdf")
plt.show()

plt.title("polarization term, in source function, C_ell * (ell (ell + 1)/2 pi)")
plt.loglog(ell_polarization, C_ell_polarization)
plt.savefig("polarization_term.pdf")
plt.show()