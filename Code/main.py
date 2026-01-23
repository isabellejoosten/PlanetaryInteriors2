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
# Define ranges

thick_core_range = np.linspace(0.75*p.thick_core, 1.25*p.thick_core, 500)  # km
density_ocean_range = np.linspace(0.5*p.density_ocean, 1.5*p.density_ocean, 500)  # kg/m³

solutions = []

# Solve for multi-layer parameters
for tc in thick_core_range:
    for do in density_ocean_range:
        try:
            to, dhpi, dc = functions.FAST_MULTILAYER_SOLVER(
                tc,
                do,
                p.thick_HPI,
                p.thick_crust,
                p.density_crust
            )
        except FloatingPointError:
            continue

        if to <= 0:
            continue

        solutions.append((tc, do, to))

solutions = np.array(solutions)
tc_arr = solutions[:, 0] / 1e3   # convert m → km
do_arr = solutions[:, 1]
to_arr = solutions[:, 2] / 1e3   # ocean thickness km

# Compute Love numbers
k2_arr = []

for tc, do, to in solutions:
    try:
        tide = functions.create_titan_model(tc, do)
    except FloatingPointError:
        k2_arr.append(np.nan)
        continue
    k2_arr.append(tide.k2.real)
k2_arr = np.array(k2_arr)
k2_arr = k2_arr.flatten()

# Define k2 range
k2_min_1 = 0.589 - 0.150
k2_max_1 = 0.589 + 0.150

k2_min_2 = 0.637 - 0.224
k2_max_2 = 0.637 + 0.224

k2_min_3 = 0.375 - 0.06
k2_max_3 = 0.375 + 0.06


# Filter solutions
mask = (k2_arr >= k2_min_1) & (k2_arr <= k2_max_1)
tc_filtered_1 = tc_arr[mask]
do_filtered_1 = do_arr[mask]
k2_filtered_1 = k2_arr[mask]

# Filter 2
mask = (k2_arr >= k2_min_2) & (k2_arr <= k2_max_2)

tc_filtered_2 = tc_arr[mask]
do_filtered_2 = do_arr[mask]
k2_filtered_2 = k2_arr[mask]

# Filter 3
mask = (k2_arr >= k2_min_3) & (k2_arr <= k2_max_3)

tc_filtered_3 = tc_arr[mask]
do_filtered_3 = do_arr[mask]
k2_filtered_3 = k2_arr[mask]


# Plot only filtered results of first observer k2
plt.figure(figsize=(8,6))
sc = plt.scatter(tc_filtered_1, do_filtered_1, c=k2_filtered_1, cmap="viridis", s=30)
plt.xlabel("Core thickness (km)")
plt.ylabel("Ocean density (kg/m³)")
#plt.title(f"Titan k2 = 0.589 +/- 0.150")
plt.colorbar(sc, label="k2 (real part)")
plt.tight_layout()
plt.show()

# Plot second observerd k2
plt.figure(figsize=(8,6))
sc = plt.scatter(tc_filtered_2, do_filtered_2, c=k2_filtered_2, cmap="viridis", s=30)
plt.xlabel("Core thickness (km)")
plt.ylabel("Ocean density (kg/m³)")
#plt.title(f"Titan k2 filtered in range 0.637 +/- 0.224")
plt.colorbar(sc, label="k2 (real part)")
plt.tight_layout()
plt.show()

# Plot third observerd k2
plt.figure(figsize=(8,6))
sc = plt.scatter(tc_filtered_3, do_filtered_3, c=k2_filtered_3, cmap="viridis", s=30)
plt.xlabel("Core thickness (km)")
plt.ylabel("Ocean density (kg/m³)")
#plt.title(f"Titan k2 filtered in range 0.375 +/- 0.06")
plt.colorbar(sc, label="k2 (real part)")
plt.tight_layout()
plt.show()







#tide_elastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "e")
#tide_viscoelastic = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "v")
#tide_liquid = functions.singleLayerTitan(p.R, p.bulkDensity, shearModulus_homogeneous, viscosity_homogeneous, "l")


#tide_multilayer = functions.multiLayerTitan(known)
#titan = delftide.TidalResponse.examples.titan()