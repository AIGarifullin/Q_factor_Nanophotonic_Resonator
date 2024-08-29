import meep as mp
from meep import mpb
import math as mt
import numpy as np

resolution = 20  # pixels/a

wavelength = 1.550e-6   # resonant wavelength
a = 1.980e-6        # units of m
r = 0.316e-6         # units of m r = r_res = 0.316-6
h = 0.220e-6         # units of m
w = 2.300e-6         # units of m

g = 2*np.pi/a       # reciprocal lattice vector

c = 3e8

nSi3N4 = 1.996
epsSi3N4 = nSi3N4**2
nAir = 1
epsAir = nAir**2

N = 30
ffmax = 0.3 # ffmax = ff_res = 0.265(at r_res = 0.316e-6)
ff = np.linspace(0.01, ffmax, N)

Si3N4 = mp.Medium(index=nSi3N4)

sc_x = 1
sc_y = 6
sc_z = 6

geometry_lattice = mp.Lattice(size=mp.Vector3(sc_x,sc_y,sc_z))

# for j in range(N):
R = np.zeros(N, float)

R = np.sqrt(ff[19]*a*w/np.pi)/a # units of "a" ff[0], ... ff[N-1]
h = h/a          # units of "a"
w = w/a          # units of "a"

geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w,h), material=Si3N4),
           mp.Cylinder(center=mp.Vector3(), radius=R, height=mp.inf, material=mp.air) ]

num_k = 5

k_points = mp.interpolate(num_k, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)])

num_bands = 4 # We need only first dielectric band

ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
              geometry=geometry,
              k_points=k_points,
              resolution=resolution,
              num_bands=num_bands)

ms.run_yodd_zeven()
