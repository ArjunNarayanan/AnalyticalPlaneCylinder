using Test
using LinearAlgebra
using PyPlot
include("plane-strain-solver.jl")
include("moduli-conversion.jl")
PS = PlaneStrainSolver

function core_displacement_coefficient(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)
    A1c =
        K * theta0 * mu / (2 * (lambda + mu) * (lambda + 2mu)) *
        (1 - inner_radius^2 / outer_radius^2)
    return A1c
end

function shell_displacement_coefficient1(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)
    A1s =
        K * theta0 * mu / (2 * (lambda + mu) * (lambda + 2mu)) *
        ((lambda + 2mu) / mu - inner_radius^2 / outer_radius^2)
    return A1s
end

function shell_displacement_coefficient2(inner_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)
    A2s = -K * theta0 * inner_radius^2 / (2 * (lambda + 2mu))
    return A2s
end

function solver_core_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
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
    return PS.core_strain_energy(solver)
end

function direct_core_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)

    A1c = core_displacement_coefficient(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    return 2 * (lambda + mu) * (A1c)^2
end

function solver_shell_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
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
    return PS.shell_strain_energy(solver, inner_radius)
end

function direct_shell_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)

    A1s = shell_displacement_coefficient1(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    A2s = shell_displacement_coefficient2(inner_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)

    return 2 * (lambda + mu) * (A1s)^2 + 2 * mu * (A2s / inner_radius^2)^2 -
           2 * K * theta0 * A1s + 0.5 * K * theta0^2
end

function solver_core_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
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
    return PS.core_compression_work(solver, V0)
end

function direct_core_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
)

    A1c = core_displacement_coefficient(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    return V0 * (2 * (lambda + mu) * A1c + 4 * (lambda + mu) * A1c^2)
end

function solver_shell_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
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
    return PS.shell_compression_work(solver, inner_radius, V0)
end

function direct_shell_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
)

    A1s = shell_displacement_coefficient1(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    A2s = shell_displacement_coefficient2(inner_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)
    cw =
        V0 * (
            -(1 - theta0) * K * theta0 +
            ((1 - theta0) * 2(lambda + mu) - 2 * K * theta0) * A1s -
            (1 - theta0) * 2 * mu * A2s / inner_radius^2 -
            4mu * A1s * A2s / inner_radius^2 + 4 * (lambda + mu) * (A1s^2)
        )
    return cw
end

function solver_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
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
    se = V0 * PS.core_strain_energy(solver)
    cw = PS.core_compression_work(solver, V0)
    return se - cw
end

function direct_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
)
    A1c = core_displacement_coefficient(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    return -2 * V0 * (lambda + mu) * (A1c + A1c^2)
end

function solver_shell_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
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
    se = V0 * PS.shell_strain_energy(solver, inner_radius)
    cw = PS.shell_compression_work(solver, inner_radius, V0)
    return se - cw
end

function direct_shell_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0,
)
    K = bulk_modulus(lambda, mu)
    A1s = shell_displacement_coefficient1(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    A2s = shell_displacement_coefficient2(inner_radius, lambda, mu, theta0)

    p =
        V0 * (
            (1 - theta0 / 2) * K * theta0 -
            2 * (lambda + mu) * A1s * (1 - theta0 + A1s) +
            2 * mu * A2s / inner_radius^2 *
            (1 - theta0 + A2s / inner_radius^2) +
            4mu * A1s * A2s / inner_radius^2
        )
end

K = 247.0
mu = 126.0
V0s = 0.8
V0c = 0.5
lambda = lame_lambda(K, mu)
theta0 = -0.067

outer_radius = 1.0
inner_radius = 1e-3:1e-3:outer_radius

solvercorestrainenergy =
    solver_core_strain_energy.(inner_radius, outer_radius, lambda, mu, theta0)
directcorestrainenergy =
    direct_core_strain_energy.(inner_radius, outer_radius, lambda, mu, theta0)

corestrainenergyerror =
    norm(solvercorestrainenergy - directcorestrainenergy) / length(inner_radius)
@test corestrainenergyerror < 10eps()

solvershellstrainenergy =
    solver_shell_strain_energy.(inner_radius, outer_radius, lambda, mu, theta0)
directshellstrainenergy =
    direct_shell_strain_energy.(inner_radius, outer_radius, lambda, mu, theta0)

shellstrainenergyerror =
    norm(solvershellstrainenergy - directshellstrainenergy) /
    length(inner_radius)
@test shellstrainenergyerror < 10eps()

solvercorecompressionwork =
    solver_core_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )
directcorecompressionwork =
    direct_core_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )
corecompressionworkerror =
    norm(solvercorecompressionwork - directcorecompressionwork) /
    length(inner_radius)
@test corecompressionworkerror < 10eps()

solvershellcompressionwork =
    solver_shell_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0s,
    )
directshellcompressionwork =
    direct_shell_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0s,
    )
shellcompressionworkerror =
    norm(solvershellcompressionwork - directshellcompressionwork) /
    length(inner_radius)
@test shellcompressionworkerror < 10eps()

solvercorepotential =
    solver_core_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0c)
directcorepotential =
    direct_core_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0c)
corepotentialerror =
    norm(solvercorepotential - directcorepotential) / length(inner_radius)
@test corepotentialerror < 10eps()

solvershellpotential =
    solver_shell_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0s)
directshellpotential =
    direct_shell_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0s)
shellpotentialerror = norm(solvershellpotential - directshellpotential)/length(inner_radius)
@test shellpotentialerror < 10eps()
