using Test
using LinearAlgebra
using PyPlot
include("plane-strain-solver.jl")
include("moduli-conversion.jl")
PS = PlaneStrainSolver

function solver_potential_difference(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
    V0s,
)
    solver = PS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )

    cse = V0c * PS.core_strain_energy(solver)
    ccw = PS.core_compression_work(solver, V0c)

    sse = V0s * PS.shell_strain_energy(solver, inner_radius)
    scw = PS.shell_compression_work(solver, inner_radius, V0s)

    pd = (sse - scw) - (cse - ccw)
end

function direct_simple_potential_difference(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
    V0s,
)

    K = bulk_modulus(lambda,mu)
    M = K*mu/(K+4mu/3)
    N = (6K+5mu)/(9*(K+mu/3))
    L = mu/(3*(K+mu/3))

    V0 = 0.5*(V0c+V0s)

    pd = M*V0*(N-1)*theta0^2 + M*V0*(1-L)*theta0^2*(inner_radius/outer_radius)^2

    return pd
end

K = 197.0               # GPa
mu = 86.0              # GPa
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3           # Kg/m^3
V0s = 1.0/rhos
V0c = 1.0/rhoc
lambda = lame_lambda(K, mu)
theta0 = -0.067

ΔG0 = -6.952e3     # J/Kg

outer_radius = 1.0
dx = outer_radius/1e3
inner_radius = dx:dx:outer_radius


solverpd =
    solver_potential_difference.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
        V0s,
    )
directpd =
    direct_simple_potential_difference.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
        V0s,
    )

pderror = maximum(abs.(solverpd - directpd))

@test pderror < abs(theta0^3)


pd = (solverpd*1e9 .+ ΔG0)/(-ΔG0)
using PyPlot
fig,ax = PyPlot.subplots()
ax.plot(inner_radius,pd)
ax.grid()
ax.set_ylabel(L"[\Phi]/(\Delta G_0)")
ax.set_xlabel("R/b")
fig.tight_layout()
fig.savefig("normalized-potential-difference.png")
