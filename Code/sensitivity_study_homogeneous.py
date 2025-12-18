import numpy as np
import matplotlib.pyplot as plt
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import functions

# averages
bulkModulus_homogeneous = (p.bulkModulus_core*functions.sphereVolume(p.R_core) + p.bulkModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.bulkModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean) + p.bulkModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
shearModulus_homogeneous = (p.shearModulus_core*functions.sphereVolume(p.R_core) + p.shearModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.shearModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.shearModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
viscosity_homogeneous = (p.viscosity_core*functions.sphereVolume(p.R_core) + p.viscosity_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.viscosity_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.viscosity_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)

# Setting up ranges to iterate over
range_shearModulus = np.arange(0.1*shearModulus_homogeneous, 10.1*shearModulus_homogeneous, 0.1*shearModulus_homogeneous)
range_viscosity = np.arange(10.0e-5*viscosity_homogeneous, 10.0e5*viscosity_homogeneous, 10000.0*viscosity_homogeneous)
print(len(range_shearModulus), len(range_viscosity))

shear_k2_real = []
shear_k2_im = []
shear_h2_real = []
shear_h2_im = []

for i in range_shearModulus:
    tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, i, viscosity_homogeneous, "v", plot=False)
    shear_k2_real.append(float(tide_viscoelastic.k2.real))
    shear_k2_im.append(abs(float(tide_viscoelastic.k2.imag)))
    shear_h2_real.append(float(tide_viscoelastic.h2.real))
    shear_h2_im.append(float(tide_viscoelastic.h2.imag))
shearlists = [shear_k2_real, shear_k2_im, shear_h2_real, shear_h2_im]

visc_k2_real = []
visc_k2_im = []
visc_h2_real = []
visc_h2_im = []

for i in range_viscosity:
    tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, i, "v", plot=False)
    visc_k2_real.append(float(tide_viscoelastic.k2.real))
    visc_k2_im.append(abs(float(tide_viscoelastic.k2.imag)))
    visc_h2_real.append(float(tide_viscoelastic.h2.real))
    visc_h2_im.append(float(tide_viscoelastic.h2.imag))
visclists = [visc_k2_real, visc_k2_im, visc_h2_real, visc_h2_im]

# plot for shear modulus
fig, axs = plt.subplots(2, 2, sharex=True)
fig.suptitle('Sensitivity of Love numbers to changes in shear modulus')
axs[0,0].set_title('k2, real part')
axs[1,0].set_xlabel('Shear modulus')
axs[1,1].set_xlabel('Shear modulus')
axs[0,0].set_ylabel('k2 (real)')
#axs[0,0].set_xscale('log')
#axs[0,0].set_yscale('log')
axs[0,0].plot(range_shearModulus, shear_k2_real, label='k2 real', color='blue')
#axs[0,1].set_xscale('log')
#axs[0,1].set_yscale('log')
axs[0,1].set_title('k2, imaginary part')
axs[0,1].set_ylabel('k2 (imaginary)')
axs[0,1].plot(range_shearModulus, shear_k2_im, label='k2 imaginary', color='red')

#axs[1,0].set_yscale('log')
#axs[1,0].set_xscale('log')
axs[1,0].set_title('h2, real part')
axs[1,0].set_ylabel('h2 (real)')
axs[1,0].plot(range_shearModulus, shear_h2_real, label='h2 real', color='green')
#axs[1,1].set_yscale('log')
axs[1,1].set_title('h2, imaginary part')
axs[1,1].set_ylabel('h2 (imaginary)')
#axs[1,1].set_xscale('log')
axs[1,1].plot(range_shearModulus, shear_h2_im, label='h2 imaginary', color='orange')

plt.show()

# plot for viscosity
fig, axs = plt.subplots(2, 2, sharex=True)
fig.suptitle('Sensitivity of Love numbers to changes in viscosity')
axs[0,0].set_title('k2, real part')
axs[1,0].set_xlabel('Viscosity')
axs[1,1].set_xlabel('Viscosity')
axs[0,0].set_ylabel('k2 (real)')
axs[0,0].set_xscale('log')
#axs[0,0].set_yscale('log')
axs[0,0].plot(range_viscosity, visc_k2_real, label='k2 real', color='blue')
axs[0,1].set_xscale('log')
#axs[0,1].set_yscale('log')
axs[0,1].set_title('k2, imaginary part')
axs[0,1].set_ylabel('k2 (imaginary)')
axs[0,1].plot(range_viscosity, visc_k2_im, label='k2 imaginary', color='red')

#axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title('h2, real part')
axs[1,0].set_ylabel('h2 (real)')
axs[1,0].plot(range_viscosity, visc_h2_real, label='h2 real', color='green')
#axs[1,1].set_yscale('log')
axs[1,1].set_title('h2, imaginary part')
axs[1,1].set_ylabel('h2 (imaginary)')
axs[1,1].set_xscale('log')
axs[1,1].plot(range_viscosity, visc_h2_im, label='h2 imaginary', color='orange')

plt.show()

# Derivatives
shearslopes = []
for array in shearlists:
    slope = functions.find_slope(array, range_shearModulus)
    shearslopes.append(slope)
viscslopes = []
for array in visclists:
    slope = functions.find_slope(array, range_viscosity)
    viscslopes.append(slope)

# plot for shear modulus slope
fig, axs = plt.subplots(2, 2, sharex=True)
fig.suptitle('Sensitivity of Love numbers to changes in shear modulus')
axs[0,0].set_title('k2, real part')
axs[1,0].set_xlabel('Shear modulus')
axs[1,1].set_xlabel('Shear modulus')
axs[0,0].set_ylabel('Slope (absolute value)')
axs[0,0].set_xscale('log')
#axs[0,0].set_yscale('log')
axs[0,0].plot(range_shearModulus, shearslopes[0], label='k2 real', color='blue')
axs[0,1].set_xscale('log')
#axs[0,1].set_yscale('log')
axs[0,1].set_title('k2, imaginary part')
axs[0,1].set_ylabel('Slope (absolute value)')
axs[0,1].plot(range_shearModulus, shearslopes[1], label='k2 imaginary', color='red')

#axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title('h2, real part')
axs[1,0].set_ylabel('Slope (absolute value)')
axs[1,0].plot(range_shearModulus, shearslopes[2], label='h2 real', color='green')
#axs[1,1].set_yscale('log')
axs[1,1].set_title('h2, imaginary part')
axs[1,1].set_ylabel('Slope (absolute value)')
axs[1,1].set_xscale('log')
axs[1,1].plot(range_shearModulus, shearslopes[3], label='h2 imaginary', color='orange')

plt.show()

# plot for viscosity slope
fig, axs = plt.subplots(2, 2, sharex=True)
fig.suptitle('Sensitivity of Love numbers to changes in viscosity')
axs[0,0].set_title('k2, real part')
axs[1,0].set_xlabel('Viscosity')
axs[1,1].set_xlabel('Viscosity')
axs[0,0].set_ylabel('Slope (absolute value)')
axs[0,0].set_xscale('log')
#axs[0,0].set_yscale('log')
axs[0,0].plot(range_viscosity, viscslopes[0], label='k2 real', color='blue')
axs[0,1].set_xscale('log')
#axs[0,1].set_yscale('log')
axs[0,1].set_title('k2, imaginary part')
axs[0,1].set_ylabel('Slope (absolute value)')
axs[0,1].plot(range_viscosity, viscslopes[1], label='k2 imaginary', color='red')

#axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title('h2, real part')
axs[1,0].set_ylabel('Slope (absolute value)')
axs[1,0].plot(range_viscosity, viscslopes[2], label='h2 real', color='green')
#axs[1,1].set_yscale('log')
axs[1,1].set_title('h2, imaginary part')
axs[1,1].set_ylabel('Slope (absolute value)')
axs[1,1].set_xscale('log')
axs[1,1].plot(range_viscosity, viscslopes[3], label='h2 imaginary', color='orange')

plt.show()