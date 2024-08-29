import meep as mp
from meep import mpb
import math as mt
import numpy as np
from scipy.special import j1
from subprocess import call
import matplotlib.pyplot as plt


resolution = 20  # pixels/a

wavelength = 1.550e-6   # resonant wavelength
a = 0.515e-6        # units of m
r = 0.316e-6         # units of m r = r_res = 0.316-6
h = 0.220e-6         # units of m
w = 2.300e-6         # units of m

g = 2*np.pi/a       # reciprocal lattice vector

c = 3e8

# fillin_factor = np.pi*r**2/(a*w)

nSi3N4 = 1.996
epsSi3N4 = nSi3N4**2
nAir = 1
epsAir = nAir**2

N = 100
ffmax = 0.3 # ffmax = ff_res = 0.265(at r_res = 0.316e-6)
ff = np.linspace(0.001, ffmax, N)

num_k = 20

R = np.zeros(N, float)
omega_res = np.zeros(N, float)
delta = np.zeros(N, float)
gamma = np.zeros(N, float)
kappa_0 = np.zeros(N, float)
kappa_1 = np.zeros(N, float)

###########################################################################

# The dependence of the resonant frequency from the radius of central hole:
for i in range(N):
    R = np.sqrt(ff[i]*a*w/np.pi)
    kappa_0[i] = ff[i]/epsAir + (1 - ff[i])/epsSi3N4
    kappa_1[i] = 2*ff[i]*(1/epsAir - 1/epsSi3N4)*j1(g*R)/(g*R)
    omega_res[i] = (1 - kappa_1[i]/(2*kappa_0[i]))*np.sqrt(kappa_0[i])*np.pi*c/a

R = np.sqrt(ff*a*w/np.pi)

# Checking the resonant frequency with specific value of the hole radius
r_res = 0.316e-6
ff_res = np.pi*r_res**2/(a*w)
kappa_0_res = ff_res/epsAir + (1 - ff_res)/epsSi3N4
kappa_1_res = 2*ff_res*(1/epsAir - 1/epsSi3N4)*j1(g*r_res)/(g*r_res)
Frequency_res = (1 - kappa_1_res/(2*kappa_0_res))*np.sqrt(kappa_0_res)*c/a/2e12
Lambda_res = c/(Frequency_res*1e12)

# omega_res = 195.48 THz
 
# The dependence of gamma(j)
# for j in range(N):
#     f[:,:,j] = np.genfromtxt("disprel_[j].dat", delimiter=",")
#     omega[j] = f[num_k+1,1,j]*g*c
#     kappa_0[j] = ff[j]/epsAir + (1 - ff[j])/epsSi3N4
#     kappa_1[j] = 2*ff[j]*(1/epsAir - 1/epsSi3N4)*j1(g*r)/(g*r)
#     delta[j] = 1 - omega[j]*a/(np.sqrt(kappa_0[j])*np.pi*c)
#     gamma[j] = np.sqrt(kappa_1[j]**2/(4*kappa_0[j]**2) - delta[j]**2)

###########################################################################
    
fig, ax = plt.subplots()        

ax.plot(R, omega_res/(2*np.pi*1e12), color = 'red', label='V(R)')
ax.plot(R, Frequency_res*np.ones(N, float), linestyle='dashed', color = 'blue', label='V_res = 195.46 THz')
ax.plot(r_res*np.ones(N, float), omega_res/(2*np.pi*1e12), linestyle='dashed', color = 'blue', label='R_res = 316 nm')
ax.set_ylabel("Frequency, v, THz", size=16)
ax.set_xlabel("Radius, R, nm", size=16)
ax.grid()
plt.legend()
plt.show()

fig, ax = plt.subplots()   

ax.plot(ff, omega_res/(2*np.pi*1e12), color = 'red', label='V(FF)')
ax.plot(ff, Frequency_res*np.ones(N, float), linestyle='dashed', color = 'blue', label='V_res = 195.46 THz')
ax.plot(ff_res*np.ones(N, float), omega_res/(2*np.pi*1e12), linestyle='dashed', color = 'blue', label='FF_res = 0.265')
ax.set_ylabel("Frequency, v, THz", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
plt.legend()
plt.show()

###########################################################################

# print("R =", R)
# print("Frequency(THz) =", omega_res/(2*np.pi*1e12))
print("ff_res =", ff_res, '\n' "Frequency_res(THZ) =", Frequency_res, '\n' "Lambda_res(m) =", Lambda_res)

# print("kappa_0=", kappa_0, '\n' "kappa_1=", kappa_1)