import numpy as np
import matplotlib.pyplot as plt
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import functions

# Single-layer model
bulkModulus_homogeneous = (p.bulkModulus_core*functions.sphereVolume(p.R_core) + p.bulkModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.bulkModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean) + p.bulkModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
shearModulus_homogeneous = (p.shearModulus_core*functions.sphereVolume(p.R_core) + p.shearModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.shearModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.shearModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
viscosity_homogeneous = (p.viscosity_core*functions.sphereVolume(p.R_core) + p.viscosity_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.viscosity_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.viscosity_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
layers = [delftide.TidalLayer("Titan", thickness=p.R, density=p.bulkDensity, shear_modulus=shearModulus_homogeneous, viscosity=viscosity_homogeneous)]
omega = 4.56e-6
model = delftide.TidalInterior("Homogeneous Titan", layers)
tide = delftide.TidalResponse(model, omega)
tide.plot()
print(tide)

k2_verification = functions.k2_analytical(shearModulus_homogeneous, p.bulkDensity, p.g, p.R, viscosity_homogeneous, omega)
print("k2 calculated analytically: ", k2_verification)
#k2_verification2 = functions.analytical_k_n(2,p.mu,p.bulkDensity,p.g,p.R)
#print("k2 calculated analytically with version 2", k2_verification2)
#print(float(tide.k2.real))
k2error_real = (float(k2_verification.real) - float(tide.k2.real)) / float(tide.k2.real) * 100
k2error_im = (float(k2_verification.imag) - float(tide.k2.imag)) / float(tide.k2.imag) * 100
print("k2 error (real part): ", k2error_real, "%")
print("k2 error (imaginary part): ", k2error_im, "%")

print("Averaged bulk modulus: ", bulkModulus_homogeneous/1000000000, " GPa")
print("Averaged shear modulus: ", shearModulus_homogeneous/1000000000, " GPa")
print("Averaged viscosity: ", viscosity_homogeneous, " Pa*s")

#functions.analytical_k_n(2,p.mu,p.bulkDensity,p.g,p.R)
