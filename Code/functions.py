import math
import matrix_propagator_incompressible_machinery as delftide
import rheologies

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

def singleLayerTitan(radius, bulkDensity, shearModulus, viscosity, omega, modelName, type="viscoelastic"):
    if type == "viscoeleastic":
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, shear_modulus=shearModulus, viscosity=viscosity)]
        omega = 4.56e-6
    elif type == "elastic":
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, rheology=rheologies.ElasticRheology())]
        omega = 4.56e-6
    elif type == "liquid":
        layers = [delftide.TidalLayer(modelName, thickness=radius, density=bulkDensity, rheology=rheologies.LiquidRheology())]
        omega = 0
    model = delftide.TidalInterior(modelName, layers)
    tide = delftide.TidalResponse(model, omega)
    tide.plot()
    print(tide)