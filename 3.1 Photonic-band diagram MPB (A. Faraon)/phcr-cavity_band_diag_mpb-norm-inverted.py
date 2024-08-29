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
from meep import mpb
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

    resolution = 32   # pixels/um
    
#     TE = str2bool(args.TE)
    
#     nSiO2 = 1.5
    nYVO = 2.2      # refractive index of the 1D PC
#     eps = nSiO2**2    # dielectric constant of the 1D PC 
    eps = nYVO**2   # dielectric constant of the 1D PC
    wl = 1/fcen       # wavelength of light (in micrometer) for 4_F_3/2 - 4_I_9/2 transition
  # wl = 1.064        # wavelength of light (in micrometer) for 4_F_3/2 - 4_I_11/2 transition
  #
    a = 0.352*wl      # nonperturbed lattice constant of the PC (in micrometer)
    wB = 0.852*wl     # width of the 1D PC' layer - beam (in micrometer)
    wG = 0.21*wl      # width of the groove between PC's layers
    theta = 60        # nanobeam interior angle
    H = wB/2*mt.tan(mt.radians(theta)) # high of the 1D PC's beam   
    h = 0.7*H         # groove depth, high of the PC's beam until the small intermediate layer
    hg = 0.3*H        # high of the small intermediate layer - groove
    wBg = 0.3*wB     # width of the small intermediate layer - groove
    
    Ng = args.Ng      # number of grooves on either side of the center (20). Total number of periods and grooves is 40
    Nb = args.Nb      # number of PC's layers - beams on either side of the center (19). Total number of beams is 38 + 1 (at the 
                      # center)    
    N = args.N        # number of defected lattices on either side of the center (7)
    Nd = np.linspace(1, 7, N) # array
    
    xL = np.zeros(N)  # x coordinate of lattice
    if defect:
        xL = 0.05*a/49*(Nd - 1)**2 + 0.95*a # perturbed lattice constants for a defect
    else:
        xL = a*np.ones(N)
        
    thbd = np.zeros(N)
    thbd = xL - wG    # perturbed thicknesses of beams
    
    # units of "a"
    wB = wB/a
    wG = wG/a   
    H = H/a   
    h = h/a     
    hg = hg/a
    wBg = wBg/a
    xL = xL/a
    thbd = thbd/a

    # multiplying by 10
#     wB = 2*wB
#     H = 2*H   
#     hg = 2*hg
#     wBg = 2*wBg
#     thbd = 2*thbd

# -----------------------------------------------------------------------------------------------------------------------    
#     Create inverted 1D PC:
    Cg = mp.inf # 2*a # one of coordinates of grooves
    
    w_c = 1.1 # width coefficient -- to the width of the Block of the inverted 1D PC
    h_c = 1.1 # height coefficient -- to the heigth of the Block of the inverted 1D PC
    
    prs = mp.Block(center=mp.Vector3(0, 0, 1/2*H), size=mp.Vector3(2*Cg, w_c*wB, h_c*H), material=mp.Medium(epsilon=eps))
    
    geometry = [prs]
    
# ------------------------------------------------------------------------------------------------------------------------    
#     Cg = Cg/a
#     Cg = thbd[0]/2
#     Cg = Cg.item()
#     geometry.append(mp.Prism(vertices = [mp.Vector3(-Cg, 0, 0), mp.Vector3(-Cg, wBg/2, hg), mp.Vector3(-Cg, -wBg/2, hg)], 
#                               height=mp.inf, #2*Cg,
#                               axis=mp.Vector3(1, 0, 0),
#                               center=mp.Vector3(0, 0, 2/3*hg),
#                               material=mp.Medium(epsilon=eps)))
    
# ------------------------------------------------------------------------------------------------------------------------    
    
    # Create central beam  with epsilon = eps
    Cb0 = thbd[0]/2  # one of coordinates of groove
    Cb0 = Cb0.item()
# ----------------------------------------------------------------------------------------------------------------------
#     geometry.append(mp.Prism(vertices = [mp.Vector3(-Cb0, 0, 0), mp.Vector3(-Cb0, wB/2, H), mp.Vector3(-Cb0, -wB/2, H)], 
#                              height= 2*Cb0, #thbd[0].item()
#                              axis=mp.Vector3(1, 0, 0),
#                              center=mp.Vector3(0, 0, 2/3*H),
#                              material=mp.Medium(epsilon=eps)))  
# ----------------------------------------------------------------------------------------------------------------------
    
    # Create a next air layer with epsilon = 1 to the right side from the central beam
    xlp = Cb0 #0 
    geometry.append(mp.Prism(vertices = [mp.Vector3(xlp, wBg/2, hg), mp.Vector3(xlp, wB/2, H), mp.Vector3(xlp, - wB/2, H), mp.Vector3(xlp, -wBg/2, hg)], 
                                 height=wG,
                                 axis=mp.Vector3(1, 0, 0),
                                 center=mp.Vector3(xlp + wG/2, 0, 2/3*H),
                                 material=mp.Medium(epsilon=1)))
    
    
    
    # Wavevector, number of bands
    num_bands = 4
    
    num_k = 10
    k_points = mp.interpolate(num_k, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)])
    
    sc_x = 1 #a 
    sc_y = 2
    sc_z = 2
    
    geometry_lattice = mp.Lattice(size=mp.Vector3(sc_x,sc_y,sc_z))
  
#     geometry_lattice = mp.Lattice(size=mp.Vector3(sc_x))
    
    ms = mpb.ModeSolver(num_bands = num_bands,
                        k_points = k_points,
                        geometry = geometry,
                        geometry_lattice = geometry_lattice,
                        resolution = resolution)
    
    
    ms.run_te()
    te_freqs = ms.all_freqs
    te_gaps = ms.gap_list
    
#     ms.run_tm()
#     tm_freqs = ms.all_freqs
#     tm_gaps = ms.gap_list
    
    
#     ms.run_yodd_zeven()
    
#     ms.init_params(p=mp.NO_PARITY,reset_fields=False)
#     ms.output_epsilon()
    
    
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
#     parser.add_argument('-TE', dest='TE', type=str, default='True', help='TE (True) or TM (False) mode (default: True)')
    
    args = parser.parse_args()
    main(args)