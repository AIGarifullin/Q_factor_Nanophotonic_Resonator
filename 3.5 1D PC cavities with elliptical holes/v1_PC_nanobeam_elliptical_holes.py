from __future__ import division

import argparse
import math as mt
import numpy as np

import meep as mp
from meep import mpb

from add_geom_object import EllipticalCylinder

def str2bool(v):
    str1 = v.lower()
    return v.lower() in ('yes', 'true', 't', '1')

def main(args):
    # PixelAvg = str2bool(args.pixavg)
    resolution = 40    # pixels/um
    a = 0.515    # period, units of um
    w = 1.200    # nanobeam width, units of um
    h = 0.300    # nanobeam height, units of um
    
    H_d = 10    # height (down) of substrate, units of um
    H_u = 4    # upper-cladding thickness, units of um

    g = 2*np.pi/a       # reciprocal lattice vector
     
    N_W = args.N_W    # full number of waveguide central holes
    N_M = args.N_M    # number of mirror holes on either side of the center
    N_T = args.N_T    # number of tapering profiles of holes on either side of the center
    
    R_1_start = 0.500    # central (starting) major radius of elliptical hole
    R_2_start = 0.110    # central (starting) minor radius of elliptical hole
    ff_start = np.pi*R_1_start*R_2_start/(a*w)    # start filling fraction
    
    R_1_end = 0.200    # central (ending) major radius of elliptical hole
    R_2_end = 0.060    # central (ending) minor radius of elliptical hole
    ff_end = np.pi*R_1_end*R_2_end/(a*w)    # start filling fraction
    
    R_1 = np.zeros(N_T+2, float)    # major radii of tapering profile holes
    R_2 = np.zeros(N_T+2, float)    # minor radii of tapering profile holes
    ff = np.zeros(N_T+2, float)    # filling factor of tapering profile of PC nanobeam
    
    N_tap = np.linspace(0, N_T+1, N_T+2) # number of holes (+2 with max, min radii) in tapering profile of PC nanobeam
    
    R_2 = mp.interpolate(N_T, [R_2_start,R_2_end])

    # Calculation of the filling fractions and radii of tapering profile of holes
    for j in range(N_T+2):
        ff[j] = (ff_start - ff_end)*(N_tap[j] - N_tap[-1])**2/(N_T + 1)**2 + ff_end    # tapering profile of filling fraction
        R_1[j] = ff[j]*a*w/(np.pi*R_2[j])    # units of um
    
    
    R_1_left = R_1[-2:0:-1]
    R_2_left = R_2[-2:0:-1]
    
    R_1_right = R_1[1:-1]
    R_2_right = R_2[1:-1]
    
    dair = 1.00    # air padding
    dpml = 1.00    # PML thickness

    sx = (2*(N_T + N_M) + N_W)*a
    sy = dpml + dair + w + dair + dpml
    sz = dpml + dair + h + dair + dpml

    cell_size = mp.Vector3(sx, sy, sz)
    boundary_layers = [mp.PML(dpml)]

    nSi3N4 = 1.996
    Si3N4 = mp.Medium(index = nSi3N4)

    geometry = [mp.Block(material = Si3N4, center = mp.Vector3(), size = mp.Vector3(mp.inf, w, h))]

    for j in range(N_M):
        major_radius = R_1_end
        minor_radius = R_2_end
        x_0 = 0.5*sx - 0.5*a - j*a
        
        # x = np.linspace(-x_0, -x_0 + minor_radius, int(resolution*minor_radius + 1))
        # y = major_radius*np.sqrt(1 - x**2/(minor_radius)**2)
        
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_end,
                                           minor_radius = R_2_end,
                                           height = mp.inf,
                                           center = mp.Vector3(-x_0, 0, 0)))                                 
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_end,
                                           minor_radius = R_2_end,
                                           height = mp.inf,
                                           center = mp.Vector3(x_0, 0, 0)))
    
    dist_N_M = abs(-0.5*sx + 0.5*a + (N_M - 1)*a + 0.5*a)
    
    for j in range(N_T):
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_left[j],
                                           minor_radius = R_2_left[j],
                                           height = mp.inf,
                                           center = mp.Vector3(-dist_N_M + 0.5*a + j*a, 0, 0)))
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_right[j],
                                           minor_radius = R_2_right[j],
                                           height = mp.inf,
                                           center = mp.Vector3(dist_N_M - 0.5*a - j*a, 0, 0)))
    
    dist_N_M_T = abs(-dist_N_M + 0.5*a + (N_T - 1)*a + 0.5*a)
        
    for j in range(int(N_W/2)):
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_start,
                                           minor_radius = R_2_start,
                                           height = mp.inf,
                                           center = mp.Vector3(-dist_N_M_T + 0.5*a + j*a, 0, 0)))
        geometry.append(EllipticalCylinder(material = mp.air,
                                           major_radius = R_1_start,
                                           minor_radius = R_2_start,
                                           height = mp.inf,
                                           center = mp.Vector3(dist_N_M_T - 0.5*a - j*a, 0, 0)))
    
    lambda_cen = 1.55    # resonant wavelength in um
    lambda_min = 1.46    # minimum source wavelength
    lambda_max = 1.66    # maximum source wavelength
    
    fmin = 1/lambda_max
    fmax = 1/lambda_min
    fcen = 0.5*(fmin + fmax)    # pulse center frequency
    df = fmax - fmin
    
    nfreq = 500    # number of frequencies at which to compute flux
    

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth = df), component = mp.Ey, center = mp.Vector3())]

    symmetries = [mp.Mirror(mp.X, +1), mp.Mirror(mp.Y, -1), mp.Mirror(mp.Z, +1)]

    sim = mp.Simulation(resolution = resolution,
                        cell_size = cell_size,
                        boundary_layers = boundary_layers,
                        geometry = geometry,
                        sources = sources,
                        dimensions = 3,
                        symmetries = symmetries)

    sim.run(mp.in_volume(mp.Volume(center = mp.Vector3(), size = mp.Vector3(sx,sy,0)), mp.at_end(mp.output_epsilon, mp.output_efield_y)),
            mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(), fcen, df)),
            until_after_sources = 500)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('-pixavg', choices=['true', 'false'], default = 'false', help ='Pixel averaging')
    parser.add_argument('-N_M',
                        type=int,
                        default=60,
                        help='number of mirror holes on either side of the center (default: 60)')
    parser.add_argument('-N_T',
                        type=int,
                        default=80,
                        help='number of tapering profiles of holes on either side of the center (default: 80)')
    parser.add_argument('-N_W',
                        type=int,
                        default=160,
                        help='full number of waveguide central holes (default: 160)')
    args = parser.parse_args()
    main(args)