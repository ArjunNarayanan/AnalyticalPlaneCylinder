using PyPlot
include("hydrostatic-solver.jl")
include("moduli-conversion.jl")
HS = HydrostaticSolver

function solver_potential_difference(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
    V0s,
)

    solver = HS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )

    cse = HS.core_strain_energy(solver)
    ccw = HS.core_compression_work(solver,V0c)

    sse = HS.shell_strain_energy(solver,inner_radius)
    scw = HS.shell_compression_work(solver,inner_radius,V0s)

    return (sse - scw) - (cse - ccw)
end

K = 197.0e9
mu = 86.0e9
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0/rhos
V0c = 1.0/rhoc
lambda = lame_lambda(K, mu)
theta0 = -0.067

ΔG0 = -6.952e3

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

pd = (solverpd .+ ΔG0)/abs(ΔG0)

using PyPlot

fig,ax = PyPlot.subplots()
ax.plot(inner_radius,pd)
ax.grid()
fig
