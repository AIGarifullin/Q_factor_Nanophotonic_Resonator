import meep as mp
from meep import mpb
import math as mt
import numpy as np
from scipy.special import j1
from subprocess import call
import matplotlib.pyplot as plt


resolution = 20  # pixels/a

wavelength = 1.550e-6   # resonant wavelength
a = 0.650e-6        # units of m
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

N = 30
ffmax = 0.3 # ffmax = ff_res = 0.265(at r_res = 0.316e-6)
ff = np.linspace(0.01, ffmax, N)

# print("ff = ", ff, '\n' "ff[9] = ", ff[9])

num_k = 20

###########################################################################

R = np.zeros(N, float)
omega_res = np.zeros(N, float)
kappa_0 = np.zeros(N, float)
kappa_1 = np.zeros(N, float)

# The dependence of the resonant frequency from the radius of central hole:
for i in range(N):
    R = np.sqrt(ff[i]*a*w/np.pi)
    kappa_0[i] = ff[i]/epsAir + (1 - ff[i])/epsSi3N4
    kappa_1[i] = 2*ff[i]*(1/epsAir - 1/epsSi3N4)*j1(g*R)/(g*R)
    omega_res[i] = (1 - kappa_1[i]/(2*kappa_0[i]))*np.sqrt(kappa_0[i])*np.pi*c/a
    
R_res = np.sqrt(ff*a*w/np.pi)

# Checking the resonant frequency with specific value of the hole radius
r_res = 0.316e-6
ff_res = np.pi*r_res**2/(a*w)
kappa_0_res = ff_res/epsAir + (1 - ff_res)/epsSi3N4
kappa_1_res = 2*ff_res*(1/epsAir - 1/epsSi3N4)*j1(g*r_res)/(g*r_res)
Frequency_res = (1 - kappa_1_res/(2*kappa_0_res))*np.sqrt(kappa_0_res)*c/a/2
Lambda_res = c/Frequency_res

# omega_res = 195.48 THz
 
# print(f[0,1])

f = np.genfromtxt("modes_FF_30_5_k_05.dat", delimiter=",")

# The dependence of gamma(j)
omega_r = 2*np.pi*Frequency_res 
omega_air = f[:,2]*g*c
omega_diel = f[:,1]*g*c
omega_0 = omega_diel + (omega_air - omega_diel)/2
gamma = np.sqrt( (omega_air - omega_diel)**2/(omega_air + omega_diel)**2 - (omega_r - omega_0)**2/omega_0**2 + 0j)
    
# print("kappa_0=", kappa_0, "kappa_1=", kappa_1)
# print("omega = ", omega, '\n' "delta = ", delta, '\n' "gamma = ", abs(gamma))
# print("gamma = ", abs(gamma))

###########################################################################

fig, ax = plt.subplots()

ax.plot(R_res, omega_res/(2*np.pi*1e12), color = 'red', label='V(R)')
ax.plot(R_res, Frequency_res/1e12*np.ones(N, float), linestyle='dashed', color = 'blue', label='V_res = 195.46 THz')
ax.plot(r_res*np.ones(N, float), omega_res/(2*np.pi*1e12), linestyle='dashed', color = 'blue', label='R_res = 316 nm')
ax.set_ylabel("Frequency, v, THz", size=16)
ax.set_xlabel("Radius, R, nm", size=16)
ax.grid()
plt.legend()
plt.show()

fig, ax = plt.subplots()   

ax.plot(ff, omega_res/(2*np.pi*1e12), color = 'red', label='V(FF)')
ax.plot(ff, Frequency_res/1e12*np.ones(N, float), linestyle='dashed', color = 'blue', label='V_res = 195.46 THz')
ax.plot(ff_res*np.ones(N, float), omega_res/(2*np.pi*1e12), linestyle='dashed', color = 'blue', label='FF_res = 0.265')
ax.set_ylabel("Frequency, v, THz", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
plt.legend()
plt.show()

fig, ax = plt.subplots()   

ax.plot(ff, kappa_1**2, color = 'red', label='kappa_1(FF)')
ax.plot(ff, kappa_0**2, color = 'blue', label='kappa_0(FF)')
ax.set_ylabel("kappa, a.u.", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
plt.legend()
plt.show()

fig, ax = plt.subplots()   

ax.plot(ff, abs(gamma), color = 'red', label='gamma(FF)')
ax.set_ylabel("gamma, a.u.", size=16)
ax.set_xlabel("Filling Factor, FF, a.u.", size=16)
ax.grid()
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

print("R_res[26] =", R_res[26])
# print("Frequency(THz) =", omega_res/(2*np.pi*1e12))
# print("ff_res =", ff_res, '\n' "Frequency_res(THZ) =", Frequency_res/1e12, '\n' "Lambda_res(m) =", Lambda_res)

# print("kappa_0=", kappa_0, '\n' "kappa_1=", kappa_1)
# print(ff)

###########################################################################



# fig, ax = plt.subplots()

# f = np.genfromtxt("modes_FF_9_k_05.dat", delimiter=",")
# ax.plot(f[:,0],f[:,1:],'b-', linewidth = 3)
# ax.plot(f[:,0],f[:,0],'g-', linewidth = 3, linestyle = '--', label='Light line, $sc_y$ = 8, $sc_z$ = 8')
# plt.xlabel("Wavevector, $k_x a/2\pi$, a.u.", size=16)
# plt.ylabel("Frequency, $\omega a/2\pi c$, a.u.", size=16)
# # plt.axis([0,0.5,0,0.4])
# # plt.xticks([t for t in np.arange(0,0.6,0.1)])
# # plt.yticks([t for t in np.arange(0,0.5,0.1)])
# ax.grid(True)
# plt.legend()
# plt.show()

v_THz = 0.323647*g*c/(2*np.pi*1e12)
lambda_v = c/(v_THz*1e12)
print("v_THz = ", v_THz, "lambda = ", lambda_v)