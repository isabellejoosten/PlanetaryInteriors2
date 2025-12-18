R = 2575000.0                               # Total radius in meters
bulkDensity = 1881.4                        # kg/m^3
mu = 8978.1394e9                            # m^3/s^2
g = 1.352

# Layers thicknesses
thick_core = 2100000.0                      # m
thick_HPI = 4000.0                          # m
thick_ocean = 400000.0                      # m
thick_crust = 70000.0                       # m

# Layer radii
R_core = thick_core                         # m
R_HPI = R_core + thick_HPI                  # m
R_ocean = R_HPI + thick_ocean               # m
R_crust = R                                 # m

# Rheology
bulkModulus_core = 200.0e9                  # Pa
bulkModulus_HPI = 20.0e9                    # Pa
bulkModulus_ocean = 2.5e9                   # Pa
bulkModulus_crust = 10.0e9                  # Pa

shearModulus_core = 70.0e9                  # Pa
shearModulus_HPI = 10.0e9                   # Pa
shearModulus_ocean = 0.0                    # Pa
shearModulus_crust = 3.0e9                  # Pa

viscosity_core = 1.0e20                     # Pa*s
viscosity_HPI = 1.0e18                      # Pa*s
viscosity_ocean = 0.0                       # Pa*s
viscosity_crust = 1.0e15                    # Pa*s

