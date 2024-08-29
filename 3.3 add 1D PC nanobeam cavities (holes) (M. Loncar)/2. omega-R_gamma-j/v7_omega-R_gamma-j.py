import meep as mp
from meep import mpb
import math as mt
import numpy as np
from scipy.special import j1
from subprocess import call
import matplotlib.pyplot as plt


resolution = 20  # pixels/a

wavelength = 1.550e-6   # resonant wavelength
a = 0.600e-6        # units of m
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

N = 50
ffmax = 0.25 # ffmax = ff_res = 0.265(at r_res = 0.316e-6)
ff = np.linspace(0.01, ffmax, N)

# print("ff = ", ff, '\n' "ff[9] = ", ff[9])

num_k = 20

###########################################################################

R = np.zeros(N, float)
omega_res = np.zeros(N, float)
kappa_0 = np.zeros(N, float)
kappa_1 = np.zeros(N, float)

# omega_res = 195.48 THz
 
f = np.genfromtxt("modes_FF_30_6_21_k_05.dat", delimiter=",")

# The dependence of gamma(j)
omega_r = 2*np.pi*195.4725*1e12 
omega_air = f[:,2]*g*c
omega_diel = f[:,1]*g*c
omega_0 = omega_diel + (omega_air - omega_diel)/2
gamma = np.sqrt( (omega_air - omega_diel)**2/(omega_air + omega_diel)**2 - (omega_r - omega_0)**2/omega_0**2 + 0j)
    
# print("omega = ", omega, '\n' "delta = ", delta, '\n' "gamma = ", abs(gamma))
print("gamma = ", gamma.real)
print("ff = ", ff)

###########################################################################

fig, ax = plt.subplots()   

ax.scatter(ff, gamma.real, color = 'red', label='gamma(FF)')
ax.set_ylabel("gamma, a.u.", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
plt.xlim(0.095,0.18)
plt.legend()
plt.show()

fig, ax = plt.subplots()   

ax.plot(ff, omega_diel/(2*np.pi*1e12), color = 'red', label='omega_diel(FF)')
ax.plot(ff, omega_air/(2*np.pi*1e12), color = 'blue', label='omega_air(FF)')
ax.set_ylabel("Frequency, v, THz", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
plt.legend()
plt.show()

###########################################################################

# print("R_res[26] =", R_res[26])
# print("Frequency(THz) =", omega_res/(2*np.pi*1e12))
# print("ff_res =", ff_res, '\n' "Frequency_res(THZ) =", Frequency_res/1e12, '\n' "Lambda_res(m) =", Lambda_res)

###########################################################################

# v_THz = 0.323647*g*c/(2*np.pi*1e12)
# lambda_v = c/(v_THz*1e12)
# print("v_THz = ", v_THz, "lambda = ", lambda_v)