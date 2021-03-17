using PyPlot
include("hydrostatic-solver.jl")
include("plane-strain-solver.jl")
include("plane-stress-solver.jl")
include("moduli-conversion.jl")
include("utilities.jl")
PS = PlaneStressSolver
PE = PlaneStrainSolver
HS = HydrostaticSolver

function deviatoric_strain_ratio(
    inner_radius,
    outer_radius,
    ls,
    ms,
    lc,
    mc,
    theta0,
    mod,
)

    solver = mod.CylindricalSolver(
        inner_radius,
        outer_radius,
        ls,
        ms,
        lc,
        mc,
        theta0,
    )

    coreinvariant = second_invariant(deviatoric_strain(mod.core_strain(solver)))
    shellinvariant = second_invariant(
        deviatoric_strain(mod.shell_strain(solver, inner_radius)),
    )
    return coreinvariant / shellinvariant
end

K1, K2 = 247.0, 192.0    # GPa
mu1, mu2 = 126.0, 87.0   # GPa
# K1,K2 = 192.0e9,192.0e9
# mu1,mu2 = 87.0e9,87.0e9
lambda1 = lame_lambda(K1, mu1)
lambda2 = lame_lambda(K2, mu2)
theta0 = -0.067

outer_radius = 1.0
dx = outer_radius / 1e3
inner_radius = dx:dx:outer_radius

planestrainratio = [
    deviatoric_strain_ratio(
        r,
        outer_radius,
        lambda1,
        mu1,
        lambda2,
        mu2,
        theta0,
        PE,
    ) for r in inner_radius
]
planestressratio = [
    deviatoric_strain_ratio(
        r,
        outer_radius,
        lambda1,
        mu1,
        lambda2,
        mu2,
        theta0,
        PS,
    ) for r in inner_radius
]
hydrostaticratio = [
    deviatoric_strain_ratio(
        r,
        outer_radius,
        lambda1,
        mu1,
        lambda2,
        mu2,
        theta0,
        HS,
    ) for r in inner_radius
]

fig, ax = PyPlot.subplots()
ax.plot(
    inner_radius,
    planestrainratio,
    label = "plane-strain",
    color = "black",
    linestyle = "solid",
)
ax.plot(
    inner_radius,
    planestressratio,
    label = "plane-stress",
    color = "black",
    linestyle = "dashed",
)
ax.plot(
    inner_radius,
    hydrostaticratio,
    label = "hydrostatic",
    color = "black",
    linestyle = "dotted",
)
ax.legend()
ax.grid()
ax.set_xlabel("Interface position")
ax.set_ylabel("Deviatoric strain invariant ratio")
fig.savefig("invariant-ratio-comparison.png")

maxplanestrain = maximum(planestrainratio)
maxplanestress = maximum(planestressratio)
maxhydrostatic = maximum(hydrostaticratio)
