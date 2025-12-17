import math

def sphereVolume(radius):
    return 4/3*math.pi*radius**3

def shellVolume(innerRadius, outerRadius):
    return sphereVolume(outerRadius) - sphereVolume(innerRadius)

def analytical_k_n(n,mu,bulkDensity,g,R):
    mu_eff_hat = mu / (bulkDensity * g * R)
    mu_n_hat = (2*n**2+4*n+3)/n * mu_eff_hat
    k_n = 2 / (2*(n-1)) * (1 / (1 + mu_n_hat))
    return k_n

def k2_analytical(mu, rho, g, radius, viscosity, omega):
    a = mu/(rho*g*radius)
    b = 1/(1-(complex(0, mu))/(viscosity*omega))

    return 3/2*1/(1+(19/2)*a*b)
