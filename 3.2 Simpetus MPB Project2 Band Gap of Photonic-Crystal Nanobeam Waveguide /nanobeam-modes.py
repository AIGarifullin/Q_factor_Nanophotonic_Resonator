import meep as mp
from meep import mpb

resolution = 20  # pixels/a

a = 0.43         # units of um
r = 0.12         # units of um
h = 0.22         # units of um
w = 0.50         # units of um

r = r/a          # units of "a"
h = h/a          # units of "a"
w = w/a          # units of "a"

nSi = 3.5
Si = mp.Medium(index=nSi)

sc_x = 1
sc_y = 8
sc_z = 8

geometry_lattice = mp.Lattice(size=mp.Vector3(sc_x,sc_y,sc_z))

geometry = [ mp.Block(center=mp.Vector3(), size=mp.Vector3(mp.inf,w,h), material=Si),
             mp.Cylinder(center=mp.Vector3(), radius=r, height=mp.inf, material=mp.air) ]

num_k = 20
k_points = mp.interpolate(num_k, [mp.Vector3(0,0,0), mp.Vector3(0.5,0,0)])

num_bands = 4

ms = mpb.ModeSolver(geometry_lattice=geometry_lattice,
                    geometry=geometry,
                    k_points=k_points,
                    resolution=resolution,
                    num_bands=num_bands)

ms.run_yodd_zeven()
# ms.run_yodd_zodd()
# ms.run_yeven_zodd()
# ms.run_yeven_zeven()

# ms.run_te()
# ms.run_tm()


# ms.init_params(p=mp.NO_PARITY,reset_fields=False)
# ms.output_epsilon()

