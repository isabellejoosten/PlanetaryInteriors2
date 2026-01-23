import math
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import sympy
import numpy as np
import parameters as params

# Thickness symbols
thick_core_sym, thick_HPI_sym, thick_ocean_sym, thick_crust_sym = sympy.symbols(
    'thick_core thick_HPI thick_ocean thick_crust', real=True
)
# Density symbols
density_core_sym, density_HPI_sym, density_ocean_sym, density_crust_sym = sympy.symbols(
    'density_core density_HPI density_ocean density_crust', real=True
)
ALL_SYMBOLS = (
    thick_core_sym, thick_HPI_sym, thick_ocean_sym, thick_crust_sym,
    density_core_sym, density_HPI_sym, density_ocean_sym, density_crust_sym
)



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
'''
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
'''



def MomentOfInertiaSphere(radius, density):
    return (8*np.pi/15) * density * radius**5

def MomentOfInertiaShell(innerRadius, outerRadius, density):
    return (8*np.pi/15) * density * (outerRadius**5 - innerRadius**5)

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

'''
def MultiLayerSolver(known_params: dict):
    # Thickness symbols
    thick_core_sym, thick_HPI_sym, thick_ocean_sym, thick_crust_sym = sympy.symbols(
        'thick_core thick_HPI thick_ocean thick_crust', real=True
    )
    # Density symbols
    density_core_sym, density_HPI_sym, density_ocean_sym, density_crust_sym = sympy.symbols(
        'density_core density_HPI density_ocean density_crust', real=True
    )
    ALL_SYMBOLS = (
        thick_core_sym, thick_HPI_sym, thick_ocean_sym, thick_crust_sym,
        density_core_sym, density_HPI_sym, density_ocean_sym, density_crust_sym
    )
    # ---------------- Equations ----------------
    eq_radius = sympy.Eq(
        totalRadius(
            thick_core_sym,
            thick_HPI_sym,
            thick_ocean_sym,
            thick_crust_sym
        ),
        p.R
    )

    eq_density = sympy.Eq(
        bulkdensity(
            thick_core_sym, density_core_sym,
            thick_HPI_sym, density_HPI_sym,
            thick_ocean_sym, density_ocean_sym,
            thick_crust_sym, density_crust_sym
        ),
        p.bulkDensity
    )

    eq_moi = sympy.Eq(
        MomentOfInertiaPlanet(
            thick_core_sym, density_core_sym,
            thick_HPI_sym, density_HPI_sym,
            thick_ocean_sym, density_ocean_sym,
            thick_crust_sym, density_crust_sym
        ),
        p.MoI_factor * p.R**2 * sphereVolume(p.R) * p.bulkDensity
    )

    equations = [eq_radius, eq_density, eq_moi]
    equations = [eq.subs(known_params) for eq in equations]
    missing = [s for s in ALL_SYMBOLS if s not in known_params]

    if len(missing) != 3:
        raise ValueError(
            f"MultiLayerSolver needs exactly 3 unknowns, got {len(missing)}"
        )

    solutions = sympy.solve(
        equations,
        missing,
        dict=True,
        simplify=True
    )

    if not solutions:
        raise RuntimeError("No solution found")
    sol = solutions[0]

    def val(sym):
        if sym in known_params:
            return float(known_params[sym])
        if sym in sol:
            return float(sol[sym].evalf())
        return 0.0

    return (
        val(thick_core_sym),
        val(thick_HPI_sym),
        val(thick_ocean_sym),
        val(thick_crust_sym),
        val(density_core_sym),
        val(density_HPI_sym),
        val(density_ocean_sym),
        val(density_crust_sym),
    )
'''

def _build_fast_solver():
    eqs = [
        totalRadius(
            thick_core_sym,
            thick_HPI_sym,
            thick_ocean_sym,
            thick_crust_sym
        ) - p.R,

        bulkdensity(
            thick_core_sym, density_core_sym,
            thick_HPI_sym, density_HPI_sym,
            thick_ocean_sym, density_ocean_sym,
            thick_crust_sym, density_crust_sym
        ) - p.bulkDensity,

        MomentOfInertiaPlanet(
            thick_core_sym, density_core_sym,
            thick_HPI_sym, density_HPI_sym,
            thick_ocean_sym, density_ocean_sym,
            thick_crust_sym, density_crust_sym
        ) - p.MoI_factor * p.R**2 * sphereVolume(p.R) * p.bulkDensity
    ]

    unknowns = (
        thick_ocean_sym,
        density_HPI_sym,
        density_core_sym
    )

    sol = sympy.solve(eqs, unknowns, dict=True)
    if not sol:
        raise RuntimeError("Symbolic multilayer solve failed")

    sol = sol[0]

    return sympy.lambdify(
        (
            thick_core_sym,
            density_ocean_sym,
            thick_HPI_sym,
            thick_crust_sym,
            density_crust_sym
        ),
        (
            sol[thick_ocean_sym],
            sol[density_HPI_sym],
            sol[density_core_sym]
        ),
        "numpy"
    )

FAST_MULTILAYER_SOLVER = _build_fast_solver()

def create_titan_model(tc, do, thpi=p.thick_HPI, tcr=p.thick_crust, dcr=p.density_crust):
    """Given core thickness & density and fixed HPI & crust, return delftide model."""
    
    # Solve for the remaining unknowns
    to, dhpi, dc = FAST_MULTILAYER_SOLVER(tc, do, thpi, tcr, dcr)

    layers = [
        delftide.TidalLayer("Core", thickness=tc, density=dc,
                            shear_modulus=p.shearModulus_core, viscosity=p.viscosity_core),
        delftide.TidalLayer("High Pressure Ice", thickness=dhpi, density=p.density_HPI,
                            shear_modulus=p.shearModulus_HPI, viscosity=p.viscosity_HPI),
        delftide.TidalLayer("Subsurface Ocean", thickness=to, density=do,
                            rheology=rheologies.LiquidRheology()),
        delftide.TidalLayer("Crust", thickness=tcr, density=dcr,
                            shear_modulus=p.shearModulus_crust, viscosity=p.viscosity_crust)
    ]
    
    omega = 4.56e-6  # Titan's orbital frequency
    model = delftide.TidalInterior("Titan multi-layer", layers)
    tide = delftide.TidalResponse(model, omega)
    
    return tide

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
'''
def precision_round(number, digits=3):
    power = "{:e}".format(number).split('e')[1]
    return round(number, -(int(power) - digits))
'''