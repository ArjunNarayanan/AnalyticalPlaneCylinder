using Test
using LinearAlgebra
using PyPlot
include("cylindrical-solver.jl")

function solver_potential_difference(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
    V0s,
)
    solver = CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )

    cse = V0c * core_strain_energy(solver)
    ccw = core_compression_work(solver, V0c)

    sse = V0s * shell_strain_energy(solver, inner_radius)
    scw = shell_compression_work(solver, inner_radius, V0s)

    pd = (sse - scw) - (cse - ccw)
end

function direct_potential_difference(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
    V0s,
)

    K = bulk_modulus(lambda, mu)
    gamma = (inner_radius / outer_radius)^2

    P1 = -K * mu / (K + 4mu / 3) * (V0s - V0c)
    P2 =
        K * mu / (2 * (K + mu / 3) * (K + 4mu / 3)^2) *
        (V0s * (36K^2 + 51K * mu + 40mu^2) / 27 + V0c * K * mu)
    P3 = K * mu / (K + 4mu / 3) * (V0s - V0c)
    P4 =
        K * mu^2 / ((K + mu / 3) * (K + 4mu / 3)^2) *
        (2 / 3 * V0s * (K - 2mu / 3) - V0c * K)
    P5 = -(K * mu)^2 / (2 * (K + mu / 3) * (K + 4mu / 3)^2) * (V0s - V0c)

    return P1 * theta0 +
           P2 * theta0^2 +
           P3 * theta0 * gamma +
           P4 * theta0^2 * gamma +
           P5 * theta0^2 * gamma^2
end


K = 247.0
mu = 126.0
V0s = 0.8
V0c = 0.5
lambda = lame_lambda(K, mu)
theta0 = -0.067

outer_radius = 1.0
inner_radius = 1e-3:1e-3:outer_radius


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
    direct_potential_difference.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
        V0s,
    )

pderror = norm(solverpd - directpd)/length(inner_radius)

@test pderror < 10eps()
