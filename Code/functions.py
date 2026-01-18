import math
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import sympy
import numpy as np
import parameters as params

def sphereVolume(radius):
    return 4/3*math.pi*radius**3

def shellVolume(innerRadius, outerRadius):
    return sphereVolume(outerRadius) - sphereVolume(innerRadius)

def k2_analytical(mu, rho, g, radius, viscosity, omega):
    a = mu/(rho*g*radius)
    b = 1/(1-(complex(0, mu))/(viscosity*omega))
    return 3/2*1/(1+(19/2)*a*b)

def h2_analytical(mu, rho, g, radius, viscosity, omega):
    a = mu/(rho*g*radius)
    b = 1/(1-(complex(0, mu))/(viscosity*omega))
    return 5/2*1/(1+(19/2)*a*b)

def computeError_complex(a, b):
    errorReal = (float(a.real) - float(b.real)) / float(b.real)
    errorIm = (float(a.imag) - float(b.imag)) / float(b.imag)

    return errorReal, errorIm

def singleLayerTitan(radius, bulkDensity, shearModulus, viscosity, type, plot=True):
    """Type = "e": elastic rheology
    Type = "v": viscoelastic rheology
    Type = "l": liquid rheology"""
    if type == "v":
        modelName = "Single-layer Titan, viscoelastic"
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, shear_modulus=shearModulus, viscosity=viscosity)]
        omega = 4.56e-6
        k2_verification = k2_analytical(shearModulus, bulkDensity, p.g, p.R, viscosity, omega)
        h2_verification = h2_analytical(shearModulus, bulkDensity, p.g, p.R, viscosity, omega)
    elif type == "e":
        modelName = "Single-layer Titan, elastic"
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, shear_modulus=shearModulus, viscosity=1e-20)]
        omega = 4.56e-6
        k2_verification = k2_analytical(shearModulus, bulkDensity, p.g, p.R, 1e-20, omega)
        h2_verification = h2_analytical(shearModulus, bulkDensity, p.g, p.R, 1e-20, omega)
    elif type == "l":
        modelName = "Single-layer Titan, liquid"
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, shear_modulus=shearModulus, viscosity=viscosity)]
        omega = 0.00000000000000000000001
        k2_verification = k2_analytical(shearModulus, bulkDensity, p.g, p.R, viscosity, omega)
        h2_verification = h2_analytical(shearModulus, bulkDensity, p.g, p.R, viscosity, omega)
    model = delftide.TidalInterior(modelName, layers)
    tide = delftide.TidalResponse(model, omega)
    if plot==True:
        tide.plot()

    print(modelName)
    print(tide.k2)
    print(tide.h2)
    k2errorReal, k2errorIm = computeError_complex(k2_verification, tide.k2)
    print("k2 error (real part): ", round(k2errorReal*100, 5), "%")
    print("k2 error (imaginary part): ", round(k2errorIm*100, 5), "%")

    h2errorReal, h2errorIm = computeError_complex(h2_verification, tide.h2)
    print("h2 error (real part): ", round(h2errorReal*100, 5), "%")
    print("h2 error (imaginary part): ", round(h2errorIm*100, 5), "%")
    
    return tide

def multiLayerTitan(known):
    thick_core, thick_HPI, thick_ocean, thick_crust, density_core, density_HPI, density_ocean, density_crust = MultiLayerSolver(known)
    layers = [
    delftide.TidalLayer("Core               ", thickness=thick_core , density=density_core, shear_modulus=p.shearModulus_core, viscosity=p.viscosity_core),
    delftide.TidalLayer("High Pressure Ice  ", thickness=thick_HPI, density=density_HPI, shear_modulus=p.shearModulus_HPI, viscosity=p.viscosity_HPI),
    delftide.TidalLayer("Subsurface Ocean   ", thickness=thick_ocean, density=density_ocean, rheology=rheologies.LiquidRheology()),
    delftide.TidalLayer("Crust              ", thickness=thick_crust, density=density_crust, shear_modulus=p.shearModulus_crust, viscosity=p.viscosity_crust),
    ]
    omega = 4.56e-6
    model = delftide.TidalInterior("Multi-layer Titan", layers)
    tide = delftide.TidalResponse(model, omega)
    tide.plot()
    print(tide)
    return model, tide

def MomentOfInertiaSphere(radius, density):
    return 2/5*sphereVolume(radius)*density*radius**2

def MomentOfInertiaShell(innerRadius, outerRadius, density):
    return 2/5*(shellVolume(innerRadius,outerRadius)*density)*((outerRadius**5-innerRadius**5)/(outerRadius**3-innerRadius**3))

def totalRadius(thick_core,thick_HPI,thick_ocean,thick_crust):
    return thick_core + thick_HPI + thick_ocean + thick_crust

def bulkdensity(thick_core, density_core, thick_HPI, density_HPI, thick_ocean, density_ocean, thick_crust, density_crust):
    R_core = thick_core
    R_HPI = R_core + thick_HPI
    R_ocean = R_HPI + thick_ocean
    R_crust = R_ocean + thick_crust
    return (sphereVolume(R_core)*density_core + shellVolume(R_core,R_HPI) * density_HPI + shellVolume(R_HPI,R_ocean) * density_ocean + shellVolume(R_ocean,R_crust) * density_crust) / sphereVolume(R_crust)

def MomentOfInertiaPlanet(thick_core, density_core, thick_HPI, density_HPI, thick_ocean, density_ocean, thick_crust, density_crust):
    R_core = thick_core
    R_HPI = R_core + thick_HPI
    R_ocean = R_HPI + thick_ocean
    R_crust = R_ocean + thick_crust
    return MomentOfInertiaSphere(R_core, density_core) + MomentOfInertiaShell(R_core,R_HPI,density_HPI) + MomentOfInertiaShell(R_HPI, R_ocean,density_ocean) + MomentOfInertiaShell(R_ocean,R_crust,density_crust)

def MultiLayerSolver(known_params: dict):
    #all_symbols = {thick_core, thick_HPI, thick_ocean, thick_crust, density_core, density_HPI, density_ocean, density_crust}
    thick_core, thick_HPI, thick_ocean, thick_crust = sympy.symbols('thick_core thick_HPI thick_ocean thick_crust')
    density_core, density_HPI, density_ocean, density_crust = sympy.symbols('density_core density_HPI density_ocean density_crust')
    eq_radius = sympy.Eq(totalRadius(thick_core, thick_HPI, thick_ocean, thick_crust),p.R)
    eq_density = sympy.Eq(bulkdensity(thick_core, density_core, thick_HPI, density_HPI, thick_ocean, density_ocean, thick_crust, density_crust),p.bulkDensity)
    eq_moi = sympy.Eq(MomentOfInertiaPlanet(thick_core, density_core, thick_HPI, density_HPI, thick_ocean, density_ocean, thick_crust, density_crust),p.MoI_factor*p.R**2*sphereVolume(p.R)*p.bulkDensity)
    equations = [eq_radius, eq_density, eq_moi]
    eqs = [eq.subs(known_params) for eq in equations]
    all_symbols = (thick_core, thick_HPI, thick_ocean, thick_crust, density_core, density_HPI, density_ocean, density_crust)
    missing = [s for s in all_symbols if s not in known_params]
    if len(missing) != len(eqs):
        raise ValueError(f"Need exactly {len(eqs)} missing variables, got {len(missing)}")

    solution = sympy.solve(eqs,missing,dict=True)
    if not solution:
        raise RuntimeError("No solution")
    sol = solution[0]
    
    def to_float(val):
        if val is None:
            return 0.0
        if isinstance(val, sympy.Basic):
            return float(val.evalf())
        return float(val)
    
    result = (
        to_float(known_params.get(thick_core, sol.get(thick_core, 0))),
        to_float(known_params.get(thick_HPI, sol.get(thick_HPI, 0))),
        to_float(known_params.get(thick_ocean, sol.get(thick_ocean, 0))),
        to_float(known_params.get(thick_crust, sol.get(thick_crust, 0))),
        to_float(known_params.get(density_core, sol.get(density_core, 0))),
        to_float(known_params.get(density_HPI, sol.get(density_HPI, 0))),
        to_float(known_params.get(density_ocean, sol.get(density_ocean, 0))),
        to_float(known_params.get(density_crust, sol.get(density_crust, 0))),
    )

    return result
'''
    R_core = solution[0][R_core]
    density_core = solution[0][density_core]
    density_ocean = solution[0][density_ocean]
    return float(R_core),float(density_core),float(density_ocean)
'''

def find_slope(array1, array2):
    '''Finds the absolute value of the slope of a curve for each data point.'''
    slope = []
    for i in range(len(array1)):
        if i == len(array1)-1:
            derivative = 0
        else:
            derivative = abs((array1[i+1]-array1[i])/(array2[i+1]-array2[i]))
        slope.append(derivative)
    return slope

# ---------------------- LAYER MODEL ADAPTED FROM ASSIGNMENT 1 ------------------------------
def create_arrays():
    '''Sets up all arrays needed to perform the integrations for the 1D model. Returns arrays of zeros for mass, pressure, and gravity. Returns an array of evenly spaced radii, and an array of densities based on the core and mantle boundaries specified in params.py.'''
    r = np.arange(0, params.R + params.delta_r, params.delta_r)
    rho = np.zeros(len(r))
    for i in range(len(r)):
        if r[i] <= params.R_core:
            rho[i] = params.density_core
        elif params.R_core < r[i] and r[i] <= params.R_HPI:
            rho[i] = params.density_HPI
        elif params.R_HPI < r[i] and r[i] <= params.R_ocean:
            rho[i] = params.density_ocean
        else:
            rho[i] = params.density_crust
    M = np.zeros(len(r))
    p = np.zeros(len(r))
    g = np.zeros(len(r))

    return M, g, p, r, rho

def Pressure(p, rho, g):
    '''Performs a simple numerical integration for the pressure at each radius increment.'''
    func = -rho*g
    p = p - func*params.delta_r
    return p

def Mass(M, r, rho):
    '''Performs a simple numerical integration for the planet mass.'''
    func = 4*np.pi*rho*r**2
    M += func*params.delta_r
    return M

def Gravity(G, M, r):
    '''Returns the gravitational acceleration at radius r.'''
    if r == 0:
        g = 0
    else:
        g = G*M/r**2
    return g

def inertia(r, rho, delta_r, M):
    '''Returns the mass moment of inertia of the planet from an array of densities at each radius increment.'''
    inertia = 0
    for i in range(len(r)):
        inertia += delta_r*rho[i]*r[i]**4

    return inertia*(8*np.pi)/(3*M[-1]*r[-1]*r[-1])

def precision_round(number, digits=3):
    power = "{:e}".format(number).split('e')[1]
    return round(number, -(int(power) - digits))