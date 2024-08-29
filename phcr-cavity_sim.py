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

    resolution = 30   # pixels/um
    
    nYVO = 2.2        # refractive index of the 1D PC
    eps = nYVO**2     # dielectric constant of the 1D PC
    wl = 1/fcen        # wavelength of light (in micrometer) for 4_F_3/2 - 4_I_9/2 transition
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
    
    sx = sxwd + 2*sxd - thbd[0]  # size of cell in x direction with defect without an overlapping central defect beam
    #sz = 2*(pad + dpml) + H

    sy = args.sy      # size of cell in y direction (perpendicular to the 1D PC.)
    sz = args.sz      # size of cell in z direction
    
    cell = mp.Vector3(sx, sy, sz)

    # Create grooves as a long triangle prism with epsilon = eps 
    Cg = (sx - 2*(pad + dpml))/2  # one of coordinates of grooves
    Cg = Cg.item()
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

    # Create defected beams by adding each beam with epsilon = eps to the right and to the left side from the central defect beam
    xlp = - thbd[0]/2
    xrp = 0 
    xrn = thbd[0]/2
    xln = 0
    for i in range(N - 1): # N - 1 because we have defined central defected beam earlier
        xlp =   (xlp + xL[i]).item()   # left coordinates of beams (positive direction)
#         xrp =   (xlp + thbd[i + 1]).item()   # right coordinates of beams (positive direction)
        xrn =   (xrn - xL[i]).item()   # right coordinates of beams (negative direction)
        xln =   (xrn - thbd[i + 1]).item()   # left coordinates of beams (negative direction)
        geometry.append(mp.Prism(vertices = [mp.Vector3(xlp, 0, 0), mp.Vector3(xlp, wB/2, H), mp.Vector3(xlp, - wB/2, H)], 
                                 height=thbd[i + 1].item(),
                                 axis=mp.Vector3(1, 0, 0),
                                 center=mp.Vector3((xlp + thbd[i + 1]/2).item(), 0, 2/3*H),
                                 material=mp.Medium(epsilon=eps)))
        geometry.append(mp.Prism(vertices = [mp.Vector3(xln, 0, 0), mp.Vector3(xln, wB/2, H), mp.Vector3(xln, - wB/2, H)], 
                                 height=thbd[i + 1].item(),
                                 axis=mp.Vector3(1, 0, 0),
                                 center=mp.Vector3((xrn - thbd[i + 1]/2).item(), 0, 2/3*H),
                                 material=mp.Medium(epsilon=eps)))

    # Create normal beams by adding each beam with epsilon = eps to the right and to the left side from the all defect lattices
    xlp = xlp + xL[N-1] # left coordinate of the first normal beam (positive direction)
    xrp = 0 # right coordinate of the first normal beam (negative direction)
    xrn = xrn - xL[N-1] # right coordinate of the first normal beam (negative direction)
    xln = 0 # left coordinate of the first normal beam (negative direction)
    T = 0
    
    for i in range(Ng - N): # Ng - N is total number of beams without defect
        xlp =   xlp + T   # left coordinates of unperturbed beams (positive direction)
#         xrp =   xlp + thb   # right coordinates of unperturbed beams (positive direction)
        xrn =   xrn - T   # right coordinates of unperturbed beams (negative direction)
        xln =   xrn - thb   # left coordinates of unperturbed beams (negative direction)
        T = a
        geometry.append(mp.Prism(vertices = [mp.Vector3(xlp, 0, 0), mp.Vector3(xlp, wB/2, H), mp.Vector3(xlp, - wB/2, H)], 
                                height=thb,
                                axis=mp.Vector3(1, 0, 0),
                                center=mp.Vector3(xlp + thb/2, 0, 2/3*H),
                                material=mp.Medium(epsilon=eps)))
        geometry.append(mp.Prism(vertices=[mp.Vector3(xln, 0, 0), mp.Vector3(xln, wB/2, H), mp.Vector3(xln, - wB/2, H)], 
                                 height=thb,
                                 axis=mp.Vector3(1, 0, 0),
                                 center=mp.Vector3(xrn - thb/2, 0, 2/3*H),
                                 material=mp.Medium(epsilon=eps)))

# Create pads at two sides of the 1D PC structure with epsilon = eps
    Cp = Cg + pad
    geometry.append(mp.Prism(vertices = [mp.Vector3(-Cp, 0, 0), mp.Vector3(-Cp, wB/2, H), mp.Vector3(-Cp, -wB/2, H)], 
                            height=pad,
                            axis=mp.Vector3(1, 0, 0),
                            center=mp.Vector3(-Cp + pad/2, 0, 2/3*H),
                            material=mp.Medium(epsilon=eps)))
    geometry.append(mp.Prism(vertices = [mp.Vector3(Cp, 0, 0), mp.Vector3(Cp, wB/2, H), mp.Vector3(Cp, -wB/2, H)], 
                            height=pad,
                            axis=mp.Vector3(1, 0, 0),
                            center=mp.Vector3(Cp - pad/2, 0, 2/3*H),
                            material=mp.Medium(epsilon=eps)))   

    sim = mp.Simulation(cell_size = cell,
                        eps_averaging=PixelAvg,
                        geometry = geometry,
                        sources = [],
                        boundary_layers = [mp.PML(dpml)],
                        resolution = resolution)

    if str2bool(args.resonant_modes):
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth = df),
                                     component = mp.Ey,
                                     center = mp.Vector3(0, 0, 2*H/3),
                                     size = mp.Vector3(0, 0, 0)))

        #sim.symmetries.append(mp.Mirror(mp.Y, phase = -1))
        #sim.symmetries.append(mp.Mirror(mp.X, phase = -1))

        sim.run(mp.at_beginning(mp.output_epsilon),
                mp.after_sources(mp.Harminv(mp.Ey, mp.Vector3(0, 0, 2*H/3), fcen, df)),
                until_after_sources = 400)

        #sim.run(mp.at_every(1/fcen/20, mp.output_efield_y), until = 1/fcen)  # simulation of field propatation
    else:
        sim.sources.append(mp.Source(mp.GaussianSource(fcen, fwidth = df),
                                     component = mp.Ey,
                                     center = mp.Vector3(-0.5*sx + dpml, 0, 2*H/3),
                                     size = mp.Vector3(0, 0, 0)))

        #sim.symmetries.append(mp.Mirror(mp.Y, phase = -1))

        freg = mp.FluxRegion(center = mp.Vector3(0.5*sx - dpml, 0, 2*H/3),
                             size = mp.Vector3(0, wB, H))

        # transmitted flux
        trans = sim.add_flux(fcen, df, nfreq, freg)

        vol = mp.Volume(mp.Vector3(), size = mp.Vector3(sx, sy))

        sim.run(mp.at_beginning(mp.output_epsilon),
                #mp.during_sources(mp.in_volume(vol, mp.to_appended("hz-slice", mp.at_every(0.4, mp.output_efield_y)))),
                mp.during_sources(mp.to_appended("ey-slice", mp.at_every(0.4, mp.output_efield_y))),
                until_after_sources = mp.stop_when_fields_decayed(50, mp.Ey, mp.Vector3(0.5*sx - dpml, 0, 2*H/3), 1e-3))

        sim.display_fluxes(trans)  # print out the flux spectrum

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