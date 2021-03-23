include("moduli-conversion.jl")


K = 247.0e9
mu = 126.0e9
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0/rhos
V0c = 1.0/rhoc
lambda = lame_lambda(K, mu)
theta0 = -0.067
ΔG0 = -6.952e3
outer_radius = 1.0

M = K*mu/(K+4mu/3)
N = (6K+5mu)/(9K+3mu)
L = mu/(3K+mu)
V0 = 0.5*(V0s+V0c)

t1 = M*V0*(1-N)*theta0^2
beta2 = M*V0*(1-L)*theta0^2
alpha2 = -(ΔG0 - t1)

alpha = sqrt(alpha2)
beta = sqrt(beta2)

m = beta2/(alpha2*outer_radius^2)
