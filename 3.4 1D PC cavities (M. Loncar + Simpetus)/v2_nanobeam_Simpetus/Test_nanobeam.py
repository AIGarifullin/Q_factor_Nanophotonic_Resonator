import meep as mp

# Define the waveguide parameters
n = 2.0  # refractive index of the waveguide material
width = 2.3  # width of the waveguide (in micrometers)
height = 0.22  # height of the waveguide (in micrometers)
resonance_wavelength = 1.55  # resonance wavelength (in micrometers)
quality_factor_threshold = 100e6  # desired quality factor threshold

# Define the range of air hole radii, number of holes, periods, and cavity lengths to sweep
radii = [0.1, 0.2, 0.3, 0.4]  # in micrometers
num_holes = [10, 12, 14, 16]
periods = [1.0, 1.2, 1.4, 1.6]  # in micrometers
cavity_lengths = [0.1, 0.2, 0.3, 0.4]  # in micrometers

optimal_radius = None
optimal_num_holes = None
optimal_period = None
optimal_cavity_length = None
max_quality_factor = 0

# Perform parameter sweep to find the optimal radii, number of holes, period, and cavity length
for radius in radii:
    for num_hole in num_holes:
        for period in periods:
            for cavity_length in cavity_lengths:
                # Create the geometry of the photonic crystal waveguide with a cavity at the center
                geometry = [mp.Block(size=mp.Vector3(mp.inf, width, height),
                                    center=mp.Vector3(0, 0, 0),
                                    material=mp.Medium(index=n))]
                for i in range(num_hole):
                    hole_position = -width / 2 + (2*i+1) * (period / 2)
                    geometry.append(mp.Cylinder(radius=radius, material=mp.vacuum, center=mp.Vector3(hole_position, 0, 0)))
                geometry.append(mp.Block(size=mp.Vector3(cavity_length, width, height), center=mp.Vector3(0, 0, 0), material=mp.Medium(index=n)))

                # Setup the simulation environment
                resolution = 10  # number of pixels per distance unit
                cell_size = mp.Vector3(16, 8, 0)
                pml_layers = [mp.PML(1.0)]

                # Create the simulation object
                sim = mp.Simulation(cell_size=cell_size,
                                    boundary_layers=pml_layers,
                                    geometry=geometry,
                                    resolution=resolution)

                # Run the simulation
                sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, mp.Vector3(0, 0, 0), 1e-3))

                # Calculate the quality factor of the resonance mode
                freq_data = sim.get_eigenmode_coefficients(mp.Harminv(mp.Ez, mp.Vector3(), resonance_wavelength, 1/resonance_wavelength).harminv_cons[0])
                for f in freq_data:
                    if f[2] > max_quality_factor and f[1] > 1e-20:
                        max_quality_factor = f[2]
                        optimal_radius = radius
                        optimal_num_holes = num_hole
                        optimal_period = period
                        optimal_cavity_length = cavity_length

# Print the optimal parameters
print("Optimal radius:", optimal_radius, "micrometers")
print("Optimal number of holes:", optimal_num_holes)
print("Optimal period:", optimal_period, "micrometers")
print("Optimal cavity length:", optimal_cavity_length, "micrometers")
print("Max quality factor:", max_quality_factor)