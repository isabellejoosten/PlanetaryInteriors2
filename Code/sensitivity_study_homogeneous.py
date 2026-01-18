import numpy as np
import matplotlib.pyplot as plt
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import functions
from matplotlib import ticker

# averages
bulkModulus_homogeneous = (p.bulkModulus_core*functions.sphereVolume(p.R_core) + p.bulkModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.bulkModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean) + p.bulkModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
shearModulus_homogeneous = (p.shearModulus_core*functions.sphereVolume(p.R_core) + p.shearModulus_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.shearModulus_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.shearModulus_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)
viscosity_homogeneous = (p.viscosity_core*functions.sphereVolume(p.R_core) + p.viscosity_HPI*functions.shellVolume(p.R_core, p.R_HPI) + p.viscosity_ocean*functions.shellVolume(p.R_HPI, p.R_ocean)  + p.viscosity_crust*functions.shellVolume(p.R_ocean, p.R)) / functions.sphereVolume(p.R)

# Setting up ranges to iterate over
range_shearModulus = np.arange(0.1*shearModulus_homogeneous, 10.1*shearModulus_homogeneous, 0.2*shearModulus_homogeneous)
range_viscosity = np.arange(10.0e-5*viscosity_homogeneous, 10.0e5*viscosity_homogeneous, 20000*viscosity_homogeneous)
print(len(range_shearModulus), len(range_viscosity))

all_k2_real = []
all_k2_im = []
all_h2_real = []
all_h2_im = []

for shearMod in range_shearModulus:
    k2_real = []
    k2_im = []
    h2_real = []
    h2_im = []
    for visc in range_viscosity:
        tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearMod, visc, "v", plot=False)
        k2_real.append(float(tide_viscoelastic.k2.real[0]))
        k2_im.append(abs(float(tide_viscoelastic.k2.imag[0])))
        h2_real.append(float(tide_viscoelastic.h2.real[0]))
        h2_im.append(float(tide_viscoelastic.h2.imag[0]))
    all_k2_real.append(k2_real)
    all_k2_im.append(k2_im)
    all_h2_real.append(h2_real)
    all_h2_im.append(h2_im)

print(len(range_shearModulus))
print(len(range_viscosity))

xlabels = []
for i in range_viscosity:
    xlabels.append(f"{i:.2e}")
ylabels = []
for i in range_shearModulus:
    ylabels.append(f"{i:.2e}")

# Plots for k2:
fig, axs1 = plt.subplots(1, 2, sharey=True)
im1 = axs1[0].imshow(all_k2_real, cmap='magma', aspect='equal', origin='lower', norm='log')
im2 = axs1[1].imshow(all_k2_im, cmap='viridis', aspect='equal', origin='lower', norm='log')
fig.colorbar(im1, ax=axs1[0], shrink=0.5)
fig.colorbar(im2, ax=axs1[1], shrink=0.5)

fig.suptitle("Sensitivity of k2 to variation in viscosity and shear modulus")

axs1[0].set_xticks(range(len(range_viscosity)), labels=xlabels, rotation=90)
axs1[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
axs1[0].set_yticks(range(len(range_shearModulus)), labels=ylabels)
axs1[0].yaxis.set_major_locator(ticker.MultipleLocator(5))
axs1[0].set_xlabel("Viscosity [Pa*s]")
axs1[0].set_ylabel("Shear Modulus [Pa]")
axs1[0].set_title("Real part")

axs1[1].set_xticks(range(len(range_viscosity)), labels=xlabels, rotation=90)
axs1[1].xaxis.set_major_locator(ticker.MultipleLocator(5))
axs1[1].set_yticks(range(len(range_shearModulus)), labels=ylabels)
axs1[1].yaxis.set_major_locator(ticker.MultipleLocator(5))
axs1[1].set_xlabel("Viscosity [Pa*s]")
axs1[1].set_ylabel("Shear Modulus [Pa]")
axs1[1].set_title("Imaginary part")

plt.show()


# plots for h2
fig, axs2 = plt.subplots(1, 2, sharey=True)
im1 = axs2[0].imshow(all_h2_real, cmap='magma', aspect='equal', origin='lower', norm='log')
im2 = axs2[1].imshow(all_h2_im, cmap='viridis', aspect='equal', origin='lower', norm='symlog')
fig.colorbar(im1, ax=axs2[0], shrink=0.5)
fig.colorbar(im2, ax=axs2[1], shrink=0.5)

fig.suptitle("Sensitivity of h2 to variation in viscosity and shear modulus")

axs2[0].set_xticks(range(len(range_viscosity)), labels=xlabels, rotation=90)
axs2[0].xaxis.set_major_locator(ticker.MultipleLocator(5))
axs2[0].set_yticks(range(len(range_shearModulus)), labels=ylabels)
axs2[0].yaxis.set_major_locator(ticker.MultipleLocator(5))
axs2[0].set_xlabel("Viscosity [Pa*s]")
axs2[0].set_ylabel("Shear Modulus [Pa]")
axs2[0].set_title("Real part")

axs2[1].set_xticks(range(len(range_viscosity)), labels=xlabels, rotation=90)
axs2[1].xaxis.set_major_locator(ticker.MultipleLocator(5))
axs2[1].set_yticks(range(len(range_shearModulus)), labels=ylabels)
axs2[1].yaxis.set_major_locator(ticker.MultipleLocator(5))
axs2[1].set_xlabel("Viscosity [Pa*s]")
axs2[1].set_ylabel("Shear Modulus [Pa]")
axs2[1].set_title("Imaginary part")


plt.show()