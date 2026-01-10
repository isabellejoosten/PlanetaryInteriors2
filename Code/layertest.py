import numpy as np
import functions
import parameters as params
import matplotlib.pyplot as plt

# Setting up arrays
M, g, p, r, rho = functions.create_arrays()

# Setting the initial core and mantle boundaries. Can be changed by editing params.py
core_boundary = params.R_core
HPI_boundary = params.R_HPI
ocean_boundary = params.R_ocean

# Setting up a counter and a moment of inertia variable for the purpose of iterating over different core and mantle boundaries
simcount = 0
inertia = 0.0

# Start iteration
while abs((inertia-params.MoI_factor)/params.MoI_factor*100) > 1.0 or abs((M[-1]-params.M_observed)/params.M_observed*100) > 1.0: # continue iterating as long as moment of inertia or total mass deviate by more than 1% from observations
    simcount += 1
    print("Starting simulation ", simcount)
    # Set mantle, core, and crust (ocean) densities
    for i in range(len(r)):
        if r[i] <= core_boundary:
            rho[i] = params.density_core
        elif core_boundary < r[i] and r[i] <= HPI_boundary:
            rho[i] = params.density_HPI
        elif HPI_boundary < r[i] and r[i] <= ocean_boundary:
            rho[i] = params.density_ocean
        else:
            rho[i] = params.density_crust

    # Integrating
    for i in range(0, len(r)-1):
        M[i+1] = functions.Mass(M[i], r[i+1], rho[i+1])

    for i in range(len(r)):
        g[i] = functions.Gravity(params.gravConstant, M[i], r[i])

    for i in np.flip(range(1, len(r))):
        p[i-1] = functions.Pressure(p[i], rho[i], g[i])
    
    # Calculate moment of inertia and compare to observations
    inertia = functions.inertia(r, rho, params.delta_r, M)
    print('Residual moment of inertia: ', abs((inertia-params.MoI_factor)/params.MoI_factor*100), ' percent')
    print('Residual mass: ', abs((M[-1]-params.M_observed)/params.M_observed*100), ' percent')

    '''
    if inertia > inertia_observed:
        mantle_boundary += delta_r
    elif inertia < inertia_observed:
        mantle_boundary -= delta_r
    
    '''
    # Adjust boundaries based on the difference between the observed and calculated mass and moment of inertia.
    # Increase the core size if the mass is too small, decrease the core size if the mass is too large.
    # Increase the mantle thickness if the moment of inertia is too small, decrease the mantle thickness if the moment of inertia is too large.
    if M[-1] > params.M_observed:
        core_boundary -= params.delta_r
    elif M[-1] < params.M_observed:
        core_boundary += params.delta_r
    if inertia > params.inertia_observed:
        mantle_boundary -= params.delta_r
    elif inertia < params.inertia_observed:
        mantle_boundary += params.delta_r

# Plot pressure, mass, density, and gravity as a function of radius.   
fig, axs = plt.subplots(1, 4, sharey=True, layout='constrained')

ax = axs[0]
ax.plot(p/1000000000,r/1000)
ax.set_xlabel('Pressure [GPa]')
ax.set_ylabel('Radius [m]')

ax = axs[1]
ax.plot(M/10000000000000000000,r/1000)
ax.set_xlabel('Mass [1000 kg]')

ax = axs[2]
ax.plot(rho, r/1000)
ax.set_xlabel('Density [kg/m^3]')

ax = axs[3]
ax.plot(g,r/1000)
ax.set_xlabel('Gravity [m/s^2]')

plt.show()
plt.clf()

# Print final results
print('\n---SIMULATION COMPLETE---')
print("Total mass: ", M[-1], " kg")
print("Mass deviation: ", (M[-1]-params.M_observed)/params.M_observed*100, "percent of observed mass")
print("Center pressure: ", p[0]/1000000000, " GPa")
print("Gravitational acceleration at surface: ", g[-1], " m/s^2")

# I honestly don't remember what this is for
for i in range(len(r)):
    if r[i] <= core_boundary and r[i+1] > core_boundary: 
        print("Gravitational acceleration at core-mantle boundary: ", g[i], " m/s^2")
    elif r[i] <= mantle_boundary and r[i+1] > mantle_boundary:
        print("Gravitational acceleration at mantle-shell boundary: ", g[i], " m/s^2")

# More final results
print("Core radius: ", core_boundary/1000, " km")
print("Crust thickness: ", (params.rtotal-mantle_boundary)/1000, " km")
print("Moment of inertia: ", inertia)
print("Number of simulations ran: ", simcount)
