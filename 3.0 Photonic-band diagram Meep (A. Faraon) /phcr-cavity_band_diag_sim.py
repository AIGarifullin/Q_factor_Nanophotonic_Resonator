# Meep Tutorial: Hz-polarized transmission and reflection through a cavity
# formed by a periodic sequence of layers of a dielectric material,
# with a defect formed by 14 periods with a perturbed lattice spacing.

# This structure is based on one analyzed in:
#    T. Zhong, J. Rochman, J.M. Kindem, E. Miyazono, and A. Faraon
#    "High quality factor nanophotonic resonators in bulk rare-earth
#    doped crystals" Opt. Express, 24 (1), 536-544 (2016).

from __future__ import division

import argparse
import meep as mp
import math as mt
import numpy as np

def str2bool(v):
    str1 = v.lower()
    return v.lower() in ('yes', 'true', 't', '1')


# Parameters for YVO crystal

def main(args):
    
    defect = str2bool(args.defect)
    fcen = args.fcen  # pulse center frequency
    df = args.df      # pulse frequency width
    nfreq = 500       # number of frequencies at which to compute flux
    
    PixelAvg = str2bool(args.pixavg)

    resolution = 20   # pixels/um
    
    nYVO = 2.2        # refractive index of the 1D PC
    eps = nYVO**2     # dielectric constant of the 1D PC
    wl = 1/fcen       # wavelength of light (in micrometer) for 4_F_3/2 - 4_I_9/2 transition
  # wl = 1.064        # wavelength of light (in micrometer) for 4_F_3/2 - 4_I_11/2 transition
  #
    wB = 0.852*wl     # width of the 1D PC' layer - beam (in micrometer)
    a = 0.352*wl      # nonperturbed lattice constant of the PC (in micrometer)
    wG = 0.21*wl      # width of the groove between PC's layers
  #
    theta = 60        # nanobeam interior angle
    H = wB/2*mt.tan(mt.radians(theta)) # high of the 1D PC's beam
  #    
    h = 0.7*H         # groove depth, high of the PC's beam until the small intermediate layer
  #
    hg = 0.3*H        # high of the small intermediate layer - groove
    wBg = 0.3*wB     # width of the small intermediate layer - groove
    
    Ng = args.Ng      # number of grooves on either side of the center (20). Total number of periods and grooves is 40
    Nb = args.Nb      # number of PC's layers - beams on either side of the center (19). Total number of beams is 38 + 1 (at the 
                      # center)
    thb = a - wG      # unperturbed thickness of a beam
    
    N = args.N        # number of defected lattices on either side of the center (7)
    Nd = np.linspace(1, 7, N) # array
    
    xL = np.zeros(N)  # x coordinate of lattice
    if defect:
        xL = 0.05*a/49*(Nd - 1)**2 + 0.95*a # perturbed lattice constants for a defect
    else:
        xL = a*np.ones(N)
        
    thbd = np.zeros(N)
    thbd = xL - wG    # perturbed thicknesses of beams
    
    pad = 3*a         # padding between last groove and PML edge         
    dpml = 1          # PML thickness
    
    sxwd = 2*(pad + dpml + a*(Ng - N)) # size of cell in x direction without defect
    
    sxd = 0
    for i in range(N):
        sxd = xL[i] + sxd # size of cell in x direction with defect on one side of the center
    
    #sx = sxwd + 2*sxd - thbd[0]  # size of cell in x direction with defect without an overlapping central defect beam
    sx = a # size of cell in x direction with only one period cell
    #sz = 2*(pad + dpml) + H

    sy = args.sy      # size of cell in y direction (perpendicular to the 1D PC.)
    sz = args.sz      # size of cell in z direction
    
    cell = mp.Vector3(sx, sy, sz)

    pml_layersY = mp.PML(dpml, direction=mp.Y)
    pml_layersZ = mp.PML(dpml, direction=mp.Z)
    
    # Create grooves as a long triangle prism with epsilon = eps 
#     Cg = (sx - 2*(pad + dpml))/2  # one of coordinates of grooves
    Cg = 6 # one of coordinates of grooves
#     Cg = Cg.item()
    prs = mp.Prism(vertices = [mp.Vector3(-Cg, 0, 0), mp.Vector3(-Cg, wBg/2, hg), mp.Vector3(-Cg, -wBg/2, hg)], 
                             height=2*Cg,
                             axis=mp.Vector3(1, 0, 0),
                             center=mp.Vector3(0, 0, 2/3*hg),
                             material=mp.Medium(epsilon=eps))
    
    geometry = [prs]
    
    # Create central defect beam  with epsilon = eps
    Cb0 = thbd[0]/2  # one of coordinates of grooves
    Cb0 = Cb0.item()
    geometry.append(mp.Prism(vertices = [mp.Vector3(-Cb0, 0, 0), mp.Vector3(-Cb0, wB/2, H), mp.Vector3(-Cb0, -wB/2, H)], 
                            height=2*Cb0,
                            axis=mp.Vector3(1, 0, 0),
                            center=mp.Vector3(0, 0, 2/3*H),
                            material=mp.Medium(epsilon=eps)))  

    fcen_s = 1/0.880  # pulse center frequency
    df_s = 4.5     # pulse freq. width: large df = short impulse

    s = mp.Source(src=mp.GaussianSource(fcen_s, fwidth=df_s), 
                                        component=mp.Ey,
                                        center = mp.Vector3(0, 0, 2*H/3),
                                        size = mp.Vector3(0, 0, 0))
    
    sim = mp.Simulation(cell_size = cell,
                        eps_averaging=PixelAvg,
                        geometry = geometry,
                        sources = [s],
                        boundary_layers = [pml_layersY, pml_layersZ],
                        resolution = resolution)
    
#     kx = 0.4
#     sim.k_point = mp.Vector3(kx)

    sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0, 2*H/3), fcen_s, df_s)),
        until_after_sources=300)

#     sim.run(mp.at_every(1/fcen_s/20, mp.output_efield_y), until=1/fcen_s)
    sim.run(mp.at_every(1/fcen_s/20), until=1/fcen_s)
    
    k_interp = 19

    sim.run_k_points(300, mp.interpolate(k_interp, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)]))
    
#         sim.display_fluxes(trans)  # print out the flux spectrum

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--resonant_modes', choices=['true', 'false'], default = 'false',
                        help = "Compute resonant modes. Default is transmission spectrum.")
    parser.add_argument('-Ng', type = int, default = 20, help ='number of grooves on either side of the center')
    parser.add_argument('-Nb', type = int, default = 19, help ='number of beams on either side of the center')
    parser.add_argument('-N', type = int,  default = 7, help ='number of defected lattices on either side of the center')
    parser.add_argument('-sy', type = int, default = 6, help ='size of cell in y direction (perpendicular to 1D PC.)')
    parser.add_argument('-sz', type = int, default = 6, help ='size of cell in z direction (perpendicular to 1D PC.)')
    parser.add_argument('-fcen', type = float, default = 1/0.880, help ='pulse center frequency')
    parser.add_argument('-df', type = float, default = 0.2, help ='pulse frequency width')
    parser.add_argument('-pixavg', choices=['true', 'false'], default = 'false', help ='Pixel averaging')
    parser.add_argument('-defect', dest='defect', choices=['true', 'false'], help='Do we need defect?')
    
    args = parser.parse_args()
    main(args)
