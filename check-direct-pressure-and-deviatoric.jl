using LinearAlgebra
using PyPlot
include("cylindrical-solver.jl")

function direct_core_pressure(inner_radius, outer_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)
    C = -K^2 * theta0 * mu / ((lambda + 2mu) * (lambda + mu))
    return C * (1 - inner_radius^2 / outer_radius^2)
end

function solver_core_pressure(inner_radius, outer_radius, lambda, mu, theta0)
    solver = CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    stress = core_stress(solver)
    return -1 / 3 * sum(stress)
end

function solver_core_deviatoric_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
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
    stress = core_stress(solver)
    p = -1 / 3 * sum(stress)
    return stress .+ p
end

function direct_core_deviatoric_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)
    C = 1 / 3 * K * theta0 * mu^2 / ((lambda + 2mu) * (lambda + mu))

    srr = C * (1 - inner_radius^2 / outer_radius^2)
    stt = srr
    szz = -2 * srr

    return [srr, stt, szz]
end

function direct_shell_pressure(inner_radius, outer_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)
    C = 1 / 3 * K * theta0 * mu / (lambda + mu)
    D = K^2 * theta0 * mu / ((lambda + 2mu) * (lambda + mu))
    return C + D * inner_radius^2 / outer_radius^2
end

function solver_shell_pressure(inner_radius, outer_radius, lambda, mu, theta0)
    solver = CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    stress = shell_stress(solver, inner_radius)
    p = -1 / 3 * sum(stress)
    return p
end

function direct_shell_deviatoric_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)

    K = bulk_modulus(lambda, mu)
    C = K * theta0 * mu / ((lambda + mu) * (lambda + 2mu))

    srr = C / 3 * (4lambda + 5mu) - C / 3 * mu * inner_radius^2 / outer_radius^2
    stt = -C / 3 * (2lambda + mu) - C / 3 * mu * inner_radius^2 / outer_radius^2
    szz =
        -2 / 3 * C * (lambda + 2mu) +
        2 / 3 * C * mu * inner_radius^2 / outer_radius^2

    return [srr, stt, szz]
end

function solver_shell_deviatoric_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
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
    stress = shell_stress(solver, inner_radius)
    p = -1 / 3 * sum(stress)
    return stress .+ p
end

K = 247.0
mu = 126.0
lambda = lame_lambda(K, mu)
theta0 = -0.067

outer_radius = 1.0
inner_radius = 1e-3:1e-3:outer_radius

solvercorepressure =
    solver_core_pressure.(inner_radius, outer_radius, lambda, mu, theta0)
directcorepressure =
    direct_core_pressure.(inner_radius, outer_radius, lambda, mu, theta0)

directcoredevstress =
    direct_core_deviatoric_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
solvercoredevstress =
    solver_core_deviatoric_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

directshellpressure =
    direct_shell_pressure.(inner_radius, outer_radius, lambda, mu, theta0)
solvershellpressure =
    solver_shell_pressure.(inner_radius, outer_radius, lambda, mu, theta0)

directshelldevstress =
    direct_shell_deviatoric_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
solvershelldevstress =
    solver_shell_deviatoric_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )


using Test

errorcorepressure =
    norm(solvercorepressure - directcorepressure) / length(inner_radius)
@test errorcorepressure < 10eps()

errorcoredevstress =
    sum(norm.(directcoredevstress - solvercoredevstress)) / length(inner_radius)
@test errorcoredevstress < 10eps()

errorshellpressure =
    norm(directshellpressure - solvershellpressure) / length(inner_radius)
@test errorshellpressure < 10eps()

errorshelldevstress =
    sum(norm.(directshelldevstress - solvershelldevstress)) /
    length(inner_radius)
@test errorshelldevstress < 10eps()
