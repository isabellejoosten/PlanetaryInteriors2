R = 2575000                             # Total radius in meters
bulkDensity = 1881.4                    # kg/m^3
mu = 8978.1394e9                        # m^3/s^2
g = 1.352
MoI_factor = 0.3414                     # no

# Layers thicknesses
thick_core = 2101000                        # m
thick_HPI = 4000.0                          # m
thick_ocean = 400000.0                      # m
thick_crust = 70000.0                       # m

#Density
density_core = 2566.7                       #kg/m^3
density_HPI = 1340                      #kg/m^3
density_ocean = 1094                      #kg/m^3
density_crust = 925                     #kg/m^3

# Layer radii
R_core = thick_core                         # m
R_HPI = R_core + thick_HPI                  # m
R_ocean = R_HPI + thick_ocean               # m
R_crust = R  


# Rheology
bulkModulus_core = 200e9                # Pa
bulkModulus_HPI = 20e9                  # Pa
bulkModulus_ocean = 2.5e9               # Pa
bulkModulus_crust = 10e9                # Pa

shearModulus_core = 70e9                # Pa
shearModulus_HPI = 10e9                # Pa
shearModulus_ocean = 0                  # Pa
shearModulus_crust = 3e9              # Pa

viscosity_core = 1e20                  # Pa*s
viscosity_HPI = 1e18
viscosity_ocean = 0
viscosity_crust = 1e15

