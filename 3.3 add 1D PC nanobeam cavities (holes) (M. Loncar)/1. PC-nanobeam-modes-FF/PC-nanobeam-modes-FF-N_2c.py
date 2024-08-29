from __future__ import division 

import argparse
import meep as mp
from meep import mpb
import math as mt
import numpy as np

def str2bool(v):
    str1 = v.lower()
    return v.lower() in ('yes', 'true', 't', '1')

def main(args):
    resolution = 20  # pixels/a
    wavelength = 1.550e-6   # resonant wavelength

    a = 0.600e-6        # units of m
    h = 0.220e-6         # units of m
    w = 2.300e-6         # units of m

    g = 2*np.pi/a       # reciprocal lattice vector

    c = 3e8

    nSi3N4 = 1.996
    epsSi3N4 = nSi3N4**2
    nAir = 1
    epsAir = nAir**2

    N = args.N        # number of defected lattices on either side of the center (7)
    Nd = np.linspace(1, N, N) # number of tapering profiles of PC nanobeam
    ff = np.zeros(N)  # ff of PC nanobeam
    ffstart = 0.1765 # at f_res
    ffend = 0.1373 # at gamma_max

    for j in range(N):
        ff[j] = ffstart - (Nd[j] - 1)**2/N**2*(ffstart - ffend) # perturbed lattice constants for a defect
    
    Si3N4 = mp.Medium(index=nSi3N4)

    sc_x = 1
    sc_y = 6
    sc_z = 6

    geometry_lattice = mp.Lattice(size=mp.Vector3(sc_x,sc_y,sc_z))

    R = np.zeros(N, float)

    h_a = h/a          # units of "a"
    w_a = w/a          # units of "a"

    for j in range(N):
        print("j =", j)
        R = np.sqrt(ff[j]*a*w/np.pi)/a # units of "a"
        geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w_a,h_a), material=Si3N4),
                 mp.Cylinder(center=mp.Vector3(), radius=R, height=mp.inf, material=mp.air) ]
        num_k = 5
        k_points = mp.interpolate(num_k, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)])
        num_bands = 2 # We need only one (or two) dielectric band
        ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                        geometry=geometry,
                        k_points=k_points,
                        resolution=resolution,
                        num_bands=num_bands)
        ms.run_yodd_zeven()

    
    
if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('-N', type = int,  default = 20, help ='number of defected holes on either side of the center')
        args = parser.parse_args()
        main(args)