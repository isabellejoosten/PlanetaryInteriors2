import math
import matrix_propagator_incompressible_machinery as delftide
import rheologies
import parameters as p

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
    elif type == "e":
        modelName = "Single-layer Titan, elastic"
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, rheology=rheologies.ElasticRheology()),]
        omega = 4.56e-6
        k2_verification = k2_analytical(1.0, bulkDensity, p.g, p.R, viscosity, omega)
    elif type == "l":
        modelName = "Single-layer Titan, liquid"
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, rheology=rheologies.LiquidRheology()),]
        omega = 1.0
        k2_verification = k2_analytical(shearModulus, bulkDensity, p.g, p.R, 1.0, omega)
    model = delftide.TidalInterior(modelName, layers)
    tide = delftide.TidalResponse(model, omega)
    tide.plot()
    print(tide)


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
