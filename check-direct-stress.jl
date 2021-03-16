using LinearAlgebra
using PyPlot
include("plane-strain-solver.jl")
include("moduli-conversion.jl")
PS = PlaneStrainSolver

function direct_core_in_plane_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)

    K = bulk_modulus(lambda, mu)
    return K * theta0 * mu / (lambda + 2mu) *
           (1 - inner_radius^2 / outer_radius^2)
end

function direct_core_out_of_plane_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)

    return K * theta0 * lambda * mu / ((lambda + mu) * (lambda + 2mu)) *
           (1 - inner_radius^2 / outer_radius^2)
end

function direct_shell_radial_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)
    return K * theta0 * mu / (lambda + 2mu) *
           (1 - inner_radius^2 / outer_radius^2)
end

function direct_shell_circumferential_stress(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)

    return -K * theta0 * mu / (lambda + 2mu) *
           (1 + inner_radius^2 / outer_radius^2)
end

function direct_shell_out_of_plane_stress(inner_radius,outer_radius,lambda,mu,theta0)
    K = bulk_modulus(lambda,mu)
    C3 = K*theta0*lambda*mu/((lambda+mu)*(lambda+2mu))
    C4 = (lambda+2mu)/lambda
    return -C3*(C4 + inner_radius^2/outer_radius^2)
end

function solver_core_stress(inner_radius, outer_radius, lambda, mu, theta0)
    solver = PS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    return PS.core_stress(solver)
end

function solver_shell_stress(inner_radius, outer_radius, lambda, mu, theta0)
    solver = PS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    return PS.shell_stress(solver, inner_radius)
end

K = 247.0
mu = 126.0
lambda = lame_lambda(K, mu)
theta0 = -0.067

outer_radius = 1.0
inner_radius = 1e-3:1e-3:outer_radius

solvercorestress =
    solver_core_stress.(inner_radius, outer_radius, lambda, mu, theta0)
solvershellstress =
    solver_shell_stress.(inner_radius, outer_radius, lambda, mu, theta0)

direct_core_srr =
    direct_core_in_plane_stress.(inner_radius, outer_radius, lambda, mu, theta0)
direct_core_szz =
    direct_core_out_of_plane_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

direct_shell_srr =
    direct_shell_radial_stress.(inner_radius, outer_radius, lambda, mu, theta0)
direct_shell_stt =
    direct_shell_circumferential_stress.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
direct_shell_szz = direct_shell_out_of_plane_stress.(inner_radius,outer_radius,lambda,mu,theta0)

solver_shell_srr = [s[1] for s in solvershellstress]
solver_shell_stt = [s[2] for s in solvershellstress]
solver_shell_szz = [s[3] for s in solvershellstress]

errorshellsrr = norm(solver_shell_srr - direct_shell_srr)
errorshellstt = norm(solver_shell_stt - direct_shell_stt)
errorshellszz = norm(solver_shell_szz - direct_shell_szz)

using Test
@test errorshellsrr < 1e3eps()
@test errorshellszz < 1e3eps()
@test errorshellstt < 1e3eps()

solver_core_srr = [s[1] for s in solvercorestress]
solver_core_szz = [s[3] for s in solvercorestress]

errorcoresrr = norm(direct_core_srr - solver_core_srr)
errorcoreszz = norm(direct_core_szz - solver_core_szz)

@test errorcoresrr < 1e3eps()
@test errorcoreszz < 1e3eps()
