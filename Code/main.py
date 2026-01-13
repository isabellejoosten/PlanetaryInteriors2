import numpy as np
import matplotlib.pyplot as plt 
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import functions
import sympy
from matplotlib.gridspec import GridSpec


# averages
bulkModulus_homogeneous = (p.bulkModulus_core*functions.sphereVolume(p.R_core) + p.bulkModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.bulkModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean) + p.bulkModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
shearModulus_homogeneous = (p.shearModulus_core*functions.sphereVolume(p.R_core) + p.shearModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.shearModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.shearModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
viscosity_homogeneous = (p.viscosity_core*functions.sphereVolume(p.R_core) + p.viscosity_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.viscosity_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.viscosity_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)

# Computing tides for homogeneous Titan, with various assumptions
fig = plt.figure(figsize=(14, 8))
gs = GridSpec(2, 4, figure=fig)
axes = [fig.add_subplot(gs[0, 0], sharex=None),
        fig.add_subplot(gs[0, 1], sharex=None),
        fig.add_subplot(gs[0, 2], sharex=None),
        fig.add_subplot(gs[1, 0], sharex=None),
        fig.add_subplot(gs[1, 1], sharex=None),
        fig.add_subplot(gs[1, 2], sharex=None)]
ax_big = fig.add_subplot(gs[:, 3])

tide_e = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e", plot=False)
tide_v = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v", plot=False)
tide_l = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l", plot=False)

tide_e.plot(axes=axes, ax_big=ax_big, label="Elastic")
tide_v.plot(axes=axes, ax_big=ax_big, label="Viscoelastic")
tide_l.plot(axes=axes, ax_big=ax_big, label="Liquid")

ax_big.legend()
plt.tight_layout()
plt.show()

# Sensitivity study for the viscoelastic case
range_shearModulus = np.arange(0.1*shearModulus_homogeneous, 10.1*shearModulus_homogeneous, 0.001*shearModulus_homogeneous)
range_viscosity = np.arange(10.0e-5*viscosity_homogeneous, 10.0e5*viscosity_homogeneous, 100.0*viscosity_homogeneous)

thick_core, thick_HPI, thick_ocean, thick_crust = sympy.symbols('thick_core thick_HPI thick_ocean thick_crust')
density_core, density_HPI, density_ocean, density_crust = sympy.symbols('density_core density_HPI density_ocean density_crust')

# currently only knows if exactly 1 thickness and 2 densities are commented out needs to be rewritten to make it work for any combination of 3
known = {
    thick_core: p.thick_core,
    thick_HPI: p.thick_HPI,
    #thick_ocean: p.thick_ocean,
    thick_crust: p.thick_crust,
    density_core: p.density_core,
    #density_HPI: p.density_HPI,
    #density_ocean: p.density_ocean,
    density_crust: p.density_crust
}
result = functions.MultiLayerSolver(known)
print(result)


#tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
#tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
#tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")


tide_multilayer = functions.multiLayerTitan(known)
titan = delftide.TidalResponse.examples.titan()
