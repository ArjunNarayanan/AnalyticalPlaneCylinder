include("hydrostatic-solver.jl")
include("utilities.jl")
include("moduli-conversion.jl")
HS = HydrostaticSolver

function direct_coefficients(inner_radius, outer_radius, lambda, mu, theta0)
    K = bulk_modulus(lambda, mu)
    C1 = (lambda - 2mu) / (lambda + 2mu)
    C2 = -K / (2 * (lambda + 2mu))

    R2 = inner_radius^2
    b2 = outer_radius^2

    A1c = -C1 / 6 * theta0 * (1 - R2 / b2)
    A1s = theta0 / 3 * (1 + C1 / 2 * R2 / b2)
    A2s = C2 * theta0 * R2
    B = theta0 / 3 * (1 - R2 / b2)

    return [A1c, A1s, A2s, B]
end

function solver_coefficients(inner_radius, outer_radius, lambda, mu, theta0)
    solver = HS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    return [solver.A1c, solver.A1s, solver.A2s, solver.B]
end

K = 192.0
mu = 87.0
lambda = lame_lambda(K, mu)
theta0 = -0.067

outer_radius = 1.0
dx = outer_radius/1e3
inner_radius = dx:dx:outer_radius

directcoeffs = direct_coefficients.(inner_radius,outer_radius,lambda,mu,theta0)
solvercoeffs = solver_coefficients.(inner_radius,outer_radius,lambda,mu,theta0)

difference = (directcoeffs .- solvercoeffs)
err = [maximum(abs.(d)) for d in difference]
maxerr = maximum(err)

using Test
@test maxerr < 10eps()
