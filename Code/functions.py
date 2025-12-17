import math

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
