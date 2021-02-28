using LinearAlgebra
using PyPlot
include("cylindrical-solver.jl")

function coefficient_matrix_determinant(lambda, mu, inner_radius)
    return -4.0 * (lambda + mu) * (lambda + 2mu) / inner_radius
end

function core_displacement_coefficient(inner_radius,outer_radius,lambda,mu,theta0)
    K = bulk_modulus(lambda,mu)
    A1c = K*theta0*mu/(2*(lambda+mu)*(lambda+2mu))*(1 - inner_radius^2/outer_radius^2)
    return A1c
end

function shell_displacement_coefficient1(inner_radius,outer_radius,lambda,mu,theta0)
    K = bulk_modulus(lambda,mu)
    A1s =
        K * theta0 * mu / (2 * (lambda + mu) * (lambda + 2mu)) *
        ((lambda + 2mu) / mu - inner_radius^2 / outer_radius^2)
    return A1s
end

function shell_displacement_coefficient2(inner_radius,lambda,mu,theta0)
    K = bulk_modulus(lambda,mu)
    A2s = -K*theta0*inner_radius^2/(2*(lambda+2mu))
    return A2s
end

function solver_coefficients(inner_radius, outer_radius, lambda, mu, theta0)
    solver = CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )

    return solver.A1c, solver.A1s, solver.A2s
end

K = 247.0
mu = 126.0
lambda = lame_lambda(K, mu)

theta0 = -0.067
outer_radius = 1.0


inner_radius = 1e-3:1e-3:outer_radius
# direct_A1c = core_displacement_coefficient.(inner_radius,outer_radius,lambda,mu,theta0)
direct_A1s = shell_displacement_coefficient1.(inner_radius,outer_radius,lambda,mu,theta0)
direct_A2s = shell_displacement_coefficient2.(inner_radius,lambda,mu,theta0)

solver_coeffs = solver_coefficients.(inner_radius,outer_radius,lambda,mu,theta0)
# solver_A1c = [s[1] for s in solver_coeffs]
solver_A1s = [s[2] for s in solver_coeffs]
solver_A2s = [s[3] for s in solver_coeffs]

# fig,ax = PyPlot.subplots()
# ax.plot(inner_radius,solver_A1c,label="solver")
# ax.plot(inner_radius,direct_A1c,label="direct")
# ax.legend()
# ax.grid()
# fig

# fig,ax = PyPlot.subplots()
# ax.plot(inner_radius,solver_A1s,label="solver")
# ax.plot(inner_radius,direct_A1s,label="direct")
# ax.legend()
# ax.grid()
# fig

fig,ax = PyPlot.subplots()
ax.plot(inner_radius,solver_A2s,label="solver")
ax.plot(inner_radius,direct_A2s,label="direct")
ax.legend()
ax.grid()
fig
