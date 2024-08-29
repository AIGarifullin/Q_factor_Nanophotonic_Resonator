import meep as mp
from meep import mpb
import math as mt
import numpy as np
from scipy.special import j1

resolution = 20  # pixels/a

lambda = 1.550e-6   # resonant wavelength
a = 0.515e-6        # units of um
r = 0.12e-6         # units of um
h = 0.22e-6         # units of um
w = 2.30e-6         # units of um

g = 2*np.pi/a       # reciprocal lattice vector

c = 3e8

ff = np.pi*r**2/(a*w)

nSi3N4 = 1.996
epsSi3N4 = nSi3N4**2
nAir = 1
epsAir = nAir**2

N = 10
ffmax = 0.5
ff = np.linspace(0, ffmax, N)

kappa_0 = ff/epsAir + (1 - ff)/epsSi3N4
kappa_1 = 2*ff*(1/epsAir - 1/epsSi3N4)*j1(g*r)/(g*r)

num_k = 20

###########################################################################

# The dependence of the resonant frequency from the radius of central hole:
for i in range(N - 1):
    R[i] = sqrt(ff[i]*a*w/np.pi)
    omega_res[i] = (1 - kappa_1[i]/(2*kappa_0[i]))*sqrt(kappa_0[i])*np.pi*c/a

# omega_res = 195.48 THz
 
# The dependence of gamma(j)
for j in range(N - 1):
    f[:,:,j] = np.genfromtxt("disprel_[j].dat", delimiter=",")
    omega[j] = f[num_k+1,1,j]*g*c
    delta[j] = 1 - omega[j]*a/(sqrt(kappa_0[j])*np.pi*c)
    gamma[j] = sqrt(kappa_1[j]**2/(4*kappa_0[j]**2) - delta[j]**2)


###########################################################################
    
fig, ax = plt.subplots()        

ax.plot(R, omega_res/(2*np.pi*1e12), color = 'red', label='TE-polarization')
ax.set_ylabel("Frequency, $\v$, THz", size=16)
ax.set_xlabel("Radius, R, nm", size=16)
ax.grid()
plt.legend()
plt.show()


