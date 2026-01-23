import numpy as np
import matplotlib.pyplot as plt
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import functions
import sympy


# averages
bulkModulus_homogeneous = (p.bulkModulus_core*functions.sphereVolume(p.R_core) + p.bulkModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.bulkModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean) + p.bulkModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
shearModulus_homogeneous = (p.shearModulus_core*functions.sphereVolume(p.R_core) + p.shearModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.shearModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.shearModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
viscosity_homogeneous = (p.viscosity_core*functions.sphereVolume(p.R_core) + p.viscosity_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.viscosity_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.viscosity_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)

# Computing tides for homogeneous Titan, with various assumptions
#tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
#tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
#tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")

# Sensitivity study for the viscoelastic case
range_shearModulus = np.arange(0.1*shearModulus_homogeneous, 10.1*shearModulus_homogeneous, 0.001*shearModulus_homogeneous)
range_viscosity = np.arange(10.0e-5*viscosity_homogeneous, 10.0e5*viscosity_homogeneous, 100.0*viscosity_homogeneous)

from functions import (
    thick_core_sym, thick_HPI_sym, thick_ocean_sym, thick_crust_sym,
    density_core_sym, density_HPI_sym, density_ocean_sym, density_crust_sym
)

thick_core_range = np.linspace(0.75*p.thick_core, 1.25*p.thick_core, 100)
density_core_range = np.linspace(0.5*p.density_core,1.5*p.density_core, 100)

solutions = []

for tc in thick_core_range:
    for dc in density_core_range:

        try:
            to, dhpi, do = functions.FAST_MULTILAYER_SOLVER(
                tc,
                dc,
                p.thick_HPI,
                p.thick_crust,
                p.density_crust
            )
        except FloatingPointError:
            continue

        if not np.isfinite(to):
            continue
        if to <= 0:
            continue

        solutions.append((tc, dc, to))

solutions = np.array(solutions)

tc_arr   = solutions[:, 0] / 1e3   # km
dc_arr   = solutions[:, 1]
to_arr   = solutions[:, 2] / 1e3   # km

if len(solutions) == 0:
    raise RuntimeError("No valid solutions found")

tc_test = p.thick_core
dc_test = p.density_core

to_fast, dhpi_fast, do_fast = functions.FAST_MULTILAYER_SOLVER(
    tc_test, dc_test, p.thick_HPI, p.thick_crust, p.density_crust
)

print("Ocean thickness (fast):", to_fast)

plt.figure(figsize=(7, 6))
sc = plt.scatter(
    tc_arr,
    dc_arr,
    c=to_arr,
    cmap="viridis",
    s=15
)

plt.xlabel("Core thickness (km)")
plt.ylabel("Core density (kg m$^{-3}$)")
plt.colorbar(sc, label="Ocean thickness (km)")
plt.title("Titan interior solutions")
plt.tight_layout()
plt.show()



#tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
#tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
#tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")


#tide_multilayer = functions.multiLayerTitan(known)
#titan = delftide.TidalResponse.examples.titan()