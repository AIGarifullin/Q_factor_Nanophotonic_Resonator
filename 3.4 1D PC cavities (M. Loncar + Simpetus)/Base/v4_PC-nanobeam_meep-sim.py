# Meep Tutorial: Hz-polarized transmission and reflection through a cavity
# formed by a periodic sequence of layers of a dielectric material,
# with a defect formed by 14 periods with a perturbed lattice spacing.

# This structure is based on one analyzed in:
#    T. Zhong, J. Rochman, J.M. Kindem, E. Miyazono, and A. Faraon
#    "High quality factor nanophotonic resonators in bulk rare-earth
#    doped crystals" Opt. Express, 24 (1), 536-544 (2016). 

# We should proceed investigations and achieve set goals!

from __future__ import division

import argparse
import meep as mp
from meep import mpb
import math as mt
import numpy as np


def str2bool(v):
    str1 = v.lower()
    return v.lower() in ('yes', 'true', 't', '1')

# Parameters for YVO crystal

def main(args):
    PixelAvg = str2bool(args.pixavg)
    
    resolution = 50  # pixels in um
    lambda_cen = 1.550   # resonant wavelength in um

    lambda_min = 1.46        # minimum source wavelength in um
    lambda_max = 1.66        # maximum source wavelength in um
    fmin = 1/lambda_max
    fmax = 1/lambda_min
    fcen = 0.5*(fmin + fmax) # pulse center frequency
    df = fmax - fmin
    nfreq = 500       # number of frequencies at which to compute flux
    
    a = 0.755        # units of um
    h = 0.142         # units of um
    w = 0.600         # units of um

    g = 2*np.pi/a       # reciprocal lattice vector

    nSi3N4 = 1.996
    epsSi3N4 = nSi3N4**2
    nAir = 1
    epsAir = nAir**2
    
    PixelAvg = str2bool(args.pixavg)
      
    N = args.N         # number of tapering profiles of holes on either side of the center (25)
    N_add = args.N_add # number of tapering profiles of holes on either side of the center (10)
    Nd = np.linspace(1, N, N) # number of tapering profile of PC nanobeam
    ff = np.zeros(N)  # ff of PC nanobeam
    # ffstart = 0.305862 # at f_res
    ffstart = 0.319310
    ffend = 0.144482 # at gamma_max
    R = np.zeros(N, float)

        # Calculate the filling fractions and radii of tapering profile of holes
    for j in range(N):
        ff[j] = ffstart - (Nd[j] - 1)**2/(N - 1)**2*(ffstart - ffend) # tapering profile of filling fraction
        R[j] = np.sqrt(ff[j]*a*w/np.pi) # units of um 
    
#############################################################################################    
    pad = 1           # padding between nanobeam and PML edge in um        
    dpml = 1          # PML thickness in um
        
    sx = 2*(pad + dpml + a*N + a*N_add) # size of cell in x direction with tapering profile and with additional holes
    sy = 2*(pad + dpml) + w       # size of cell in y direction (perpendicular to the 1D PC.)
    sz = 2*(pad + dpml) + h       # size of cell in z direction
    
    cell = mp.Vector3(sx, sy, sz)
    
#############################################################################################
 
    Si3N4 = mp.Medium(index=nSi3N4)
      
    geometry = [mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w,h), material=Si3N4)]
        
    # Create tapering profile of holes to the right side of the L = 0
    xc_r = 0.5*a # x-coordinate of first hole
    for j in range(N): # number of mirrors
        geometry.append(mp.Cylinder(center=mp.Vector3(xc_r,0,0), radius=R[j], height=mp.inf, material=mp.air))
        xc_r = xc_r + a
    
    # Create tapering profile of holes to the left side of the L = 0
    xc_l = -0.5*a # x-coordinate of first hole
    for j in range(N): # number of mirrors
        geometry.append(mp.Cylinder(center=mp.Vector3(xc_l,0,0), radius=R[j], height=mp.inf, material=mp.air))
        xc_l = xc_l - a
        
    # Create N_add additional mirror segments of holes with constant radius to right side of L=0
    Rend = np.sqrt(ffend*a*w/np.pi)
    # N_add = 10 # number of additional mirror segments of holes
    
    xc_r = xc_r # x-coordinate of the first additional hole
    for j in range(N_add): # number of mirrors
        geometry.append(mp.Cylinder(center=mp.Vector3(xc_r,0,0), radius=Rend, height=mp.inf, material=mp.air))
        xc_r = xc_r + a
        
    # Create N_add additional mirror segments of holes with constant radius to left side of L=0
    xc_l = xc_l # x-coordinate of the first additional hole
    for j in range(N_add): # number of mirrors
        geometry.append(mp.Cylinder(center=mp.Vector3(xc_l,0,0), radius=Rend, height=mp.inf, material=mp.air))
        xc_l = xc_l - a
        
#############################################################################################                
    
    sim = mp.Simulation(cell_size = cell,
                        eps_averaging = PixelAvg,
                        geometry = geometry,
                        sources = [],
                        boundary_layers = [mp.PML(dpml)],
                        resolution = resolution)

    if str2bool(args.resonant_modes):
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth = df),
                                     component = mp.Ey,
                                     center = mp.Vector3(0, 0, 0),
                                     size = mp.Vector3(0, 0, 0)))

        sim.run(mp.at_beginning(mp.output_epsilon),
                mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0, 0), fcen, df)),
                until_after_sources = 400)
        # sim.run(mp.at_every(1/fcen/20, mp.output_efield_y), until = 1/fcen)
    else:
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth = df),
                                     component = mp.Ey,
                                     center = mp.Vector3(-0.5*sx + dpml, 0, 0),
                                     size = mp.Vector3(0, 0, 0)))

        freg = mp.FluxRegion(center = mp.Vector3(0.5*sx - dpml, 0, 0),
                             size = mp.Vector3(0, w, h))

        # transmitted flux
        trans = sim.add_flux(fcen, df, nfreq, freg)

        vol = mp.Volume(mp.Vector3(), size = mp.Vector3(sx, sy))

        sim.run(mp.at_beginning(mp.output_epsilon),
                #mp.during_sources(mp.in_volume(vol, mp.to_appended("ey-slice", mp.at_every(0.4, mp.output_efield_y)))),
        mp.during_sources(mp.to_appended("ey-slice", mp.at_every(0.4, mp.output_efield_y))),
        until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(0.5*sx - dpml, 0, 0), 1e-3))

        sim.display_fluxes(trans)  # print out the flux spectrum 
    
         
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type = int,  default = 20, help ='number of defected holes on either side of the center')
    parser.add_argument('-N_add', type = int,  default = 10, help ='number of additional constant radius holes on either side of the center')
    parser.add_argument('-pixavg', choices=['true', 'false'], default = 'false', help ='Pixel averaging')
    parser.add_argument('-r', '--resonant_modes', choices=['true', 'false'], default = 'false', help = "Compute resonant modes. Default is transmission spectrum.")
    # parser.add_argument('-df', type = float, default = 0.2, help ='pulse frequency width')

    args = parser.parse_args()
    main(args)