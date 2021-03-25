using DifferentialEquations
include("moduli-conversion.jl")

function interface_derivative(R, p, t)
    G = p[1]
    gamma2 = p[2]
    return -1.0 + (R^2 - gamma2)/G
end

K = 247.0e9
mu = 126.0e9
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc
lambda = lame_lambda(K, mu)
theta0 = -0.067
ΔG0Jmol = -14351.0
molarmass = 0.147
ΔG0 = ΔG0Jmol / molarmass
outer_radius = 1.0

M = K * mu / (K + 4mu / 3)
N = (6K + 5mu) / (9K + 3mu)
L = mu / (3K + mu)
V0 = 0.5 * (V0s + V0c)

beta2 = M * V0 * (1 - L) * theta0^2
alpha2 = M * V0 * (1 - N) * theta0^2
gamma2 = alpha2/beta2

G = -ΔG0/beta2

Tmax = 1.2
tspan = (0.0, Tmax)
prob = ODEProblem(interface_derivative, 1.0, tspan, [G,gamma2])
sol = solve(prob)

t = 0:1e-2:Tmax
numericalsol = sol.(t)

using PyPlot
fig, ax = PyPlot.subplots()
ax.plot(t, numericalsol, color = "black")
ax.plot([0, 1], [1, 0], "--", color = "black")
ax.plot()
ax.set_xlabel(L"t/t_K")
ax.set_ylabel(L"R/b")
ax.grid()
fig
# fig.savefig("normalized-interface-position.png")
