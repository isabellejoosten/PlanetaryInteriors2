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

R_core, density_core, density_ocean = sympy.symbols('R_core, density_core, density_ocean')
eq1 = sympy.Eq(functions.totalRadius(R_core, p.R_HPI, p.R_ocean, p.R_ocean),p.R)
eq2 = sympy.Eq(functions.bulkdensity(R_core, density_core, p.R_HPI, p.density_HPI, p.R_ocean, density_ocean, p.R_crust, p.density_crust),p.bulkDensity)
eq3 = sympy.Eq(functions.MomentOfInertiaPlanet(R_core, density_core, p.R_HPI, p.density_HPI, p.R_ocean, density_ocean, p.R_crust, p.density_crust),p.MoI_factor*p.R**2*functions.sphereVolume(p.R)*p.bulkDensity)

solution = sympy.solve((eq1,eq2,eq3),(R_core,density_core,density_ocean))
print(solution)

print("Averaged bulk modulus: ", bulkModulus_homogeneous/1000000000, " GPa")
print("Averaged shear modulus: ", shearModulus_homogeneous/1000000000, " GPa")
print("Averaged viscosity: ", viscosity_homogeneous, " Pa*s")

tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")
