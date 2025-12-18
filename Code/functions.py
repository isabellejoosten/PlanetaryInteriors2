import math
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p
import sympy

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

def singleLayerTitan(radius, bulkDensity, shearModulus, viscosity, type):
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
    print(tide)

def multiLayerTitan():
    R_core,density_core,density_ocean = MultiLayerSolver()
    layers = [
    delftide.TidalLayer("Core               ", thickness=R_core , density=density_core, shear_modulus=p.shearModulus_core, viscosity=p.viscosity_core),
    delftide.TidalLayer("High Pressure Ice  ", thickness=p.thick_HPI, density=p.density_HPI, shear_modulus=p.shearModulus_HPI, viscosity=p.viscosity_HPI),
    delftide.TidalLayer("Subsurface Ocean   ", thickness=p.thick_ocean, density=density_ocean, rheology=rheologies.LiquidRheology()),
    delftide.TidalLayer("Crust              ", thickness=p.thick_crust, density=p.density_crust, shear_modulus=p.shearModulus_crust, viscosity=p.viscosity_crust),
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

def totalRadius(R_core,R_HPI,R_ocean,R_crust):
    return R_core + R_HPI + R_ocean + R_crust

def bulkdensity(R_core, density_core, R_HPI, density_HPI, R_ocean, density_ocean, R_crust, density_crust):
    return (sphereVolume(R_core)*density_core + shellVolume(R_core,R_HPI) * density_HPI + shellVolume(R_HPI,R_ocean) * density_ocean + shellVolume(R_ocean,R_crust) * density_crust) / sphereVolume(R_crust)

def MomentOfInertiaPlanet(R_core, density_core, R_HPI, density_HPI, R_ocean, density_ocean, R_crust, density_crust):
    return MomentOfInertiaSphere(R_core, density_core) + MomentOfInertiaShell(R_core,R_HPI,density_HPI) + MomentOfInertiaShell(R_HPI, R_ocean,density_ocean) + MomentOfInertiaShell(R_ocean,R_crust,density_crust)

def MultiLayerSolver():
    R_core, density_core, density_ocean = sympy.symbols('R_core, density_core, density_ocean')
    eq1 = sympy.Eq(totalRadius(R_core, p.thick_HPI, p.thick_ocean, p.thick_crust),p.R)
    eq2 = sympy.Eq(bulkdensity(R_core, density_core, p.R_HPI, p.density_HPI, p.R_ocean, density_ocean, p.R_crust, p.density_crust),p.bulkDensity)
    eq3 = sympy.Eq(MomentOfInertiaPlanet(R_core, density_core, p.R_HPI, p.density_HPI, p.R_ocean, density_ocean, p.R_crust, p.density_crust),p.MoI_factor*p.R**2*sphereVolume(p.R)*p.bulkDensity)

    solution = sympy.solve((eq1,eq2,eq3),R_core,density_core,density_ocean,dict=True)
    R_core = solution[0][R_core]
    density_core = solution[0][density_core]
    density_ocean = solution[0][density_ocean]
    return float(R_core),float(density_core),float(density_ocean)

