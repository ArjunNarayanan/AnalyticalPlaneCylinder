using PyPlot
include("hydrostatic-solver.jl")
include("utilities.jl")
include("moduli-conversion.jl")
HS = HydrostaticSolver

function deviatoric_strain_ratio(
    inner_radius,
    outer_radius,
    ls,
    ms,
    lc,
    mc,
    theta0,
)

    solver =
        HS.CylindricalSolver(inner_radius, outer_radius, ls, ms, lc, mc, theta0)

    coreinvariant = second_invariant(deviatoric_strain(HS.core_strain(solver)))
    shellinvariant = second_invariant(
        deviatoric_strain(HS.shell_strain(solver, inner_radius)),
    )
    return coreinvariant/shellinvariant
end

K1, K2 = 247.0e9, 192.0e9    # GPa
mu1, mu2 = 126.0e9, 87.0e9   # GPa
# K1,K2 = 192.0e9,192.0e9
# mu1,mu2 = 87.0e9,87.0e9
lambda1 = lame_lambda(K1, mu1)
lambda2 = lame_lambda(K2, mu2)
theta0 = -0.067

outer_radius = 1.0e-3
dx = outer_radius/1e3
inner_radius = dx:dx:outer_radius

devinvariantratio =
    deviatoric_strain_ratio.(
        inner_radius,
        outer_radius,
        lambda1,
        mu1,
        lambda2,
        mu2,
        theta0,
    )


fig,ax = PyPlot.subplots()
ax.plot(inner_radius,devinvariantratio)
ax.grid()
ax.set_xlabel("Interface position")
ax.set_ylabel("Deviatoric strain invariant ratio")
fig
