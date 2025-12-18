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

# Single-layer model - viscoelastic
layers = [delftide.TidalLayer("Titan", thickness=p.R, density=p.bulkDensity, shear_modulus=shearModulus_homogeneous, viscosity=viscosity_homogeneous)]
omega = 4.56e-6
model = delftide.TidalInterior("Homogeneous Titan", layers)
tide = delftide.TidalResponse(model, omega)
tide.plot()
print(tide)

# Single-layer model - viscoelastic
layers = [delftide.TidalLayer("Titan", thickness=p.R, density=p.bulkDensity, shear_modulus=shearModulus_homogeneous, viscosity=viscosity_homogeneous)]
omega = 4.56e-6
model = delftide.TidalInterior("Homogeneous Titan", layers)
tide = delftide.TidalResponse(model, omega)
tide.plot()
print(tide)

k2_verification = functions.k2_analytical(shearModulus_homogeneous, p.bulkDensity, p.g, p.R, viscosity_homogeneous, omega)
print("k2 calculated analytically: ", k2_verification)
#print(float(tide.k2.real))
k2error_real = (float(k2_verification.real) - float(tide.k2.real)) / float(tide.k2.real) * 100
#k2error_im = (float(k2_verification.imag) - float(tide.k2.imag)) / float(tide.k2.imag) * 100
print("k2 error (real part): ", k2error_real, "%")
#print("k2 error (imaginary part): ", k2error_im, "%")

h2_verification = functions.h2_analytical(shearModulus_homogeneous, p.bulkDensity, p.g, p.R, viscosity_homogeneous, omega)
print("h2 calculated analytically: ", k2_verification)
h2error_real = (float(h2_verification.real) - float(tide.h2.real)) / float(tide.h2.real) * 100
#h2error_im = (float(h2_verification.imag) - float(tide.h2.imag)) / float(tide.h2.imag) * 100
print("h2 error (real part): ", h2error_real, "%")
#print("h2 error (imaginary part): ", h2error_im, "%")

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
