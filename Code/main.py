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


R_core,density_core,density_ocean = functions.MultiLayerSolver()

print("Radius of the core is ", np.floor(R_core/10**3) , 'km')
print("Density of the core is ", np.floor(density_core) , 'kg/m^3')
print("Density of the ocean is ", np.floor(density_ocean) , 'kg/m^3')


'''
print("Averaged bulk modulus: ", bulkModulus_homogeneous/1000000000, " GPa")
print("Averaged shear modulus: ", shearModulus_homogeneous/1000000000, " GPa")
print("Averaged viscosity: ", viscosity_homogeneous, " Pa*s")
'''

#tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
#tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
#tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")


tide_multilayer = functions.multiLayerTitan()
titan = delftide.TidalResponse.examples.titan()