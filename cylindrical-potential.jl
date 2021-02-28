using PyPlot
include("cylindrical-solver.jl")

function normalized_potential_difference(interfaceradius)
    K1, K2 = 247.0, 192.0    # GPa
    mu1, mu2 = 126.0, 87.0   # GPa

    lambda1 = lame_lambda(K1, mu1)
    lambda2 = lame_lambda(K2, mu2)

    rho1 = 3.93e3           # Kg/m^3
    rho2 = 3.68e3           # Kg/m^3
    V01 = 1.0 / rho1
    V02 = 1.0 / rho2

    ΔG0 = -6.95e-3

    theta0 = -0.067

    outerradius = 1.0

    analyticalsolution = CylindricalSolver(
        interfaceradius,
        outerradius,
        lambda1,
        mu1,
        lambda2,
        mu2,
        theta0,
    )


    shellpressure = shell_pressure(analyticalsolution, interfaceradius)
    shelldevstress =
        shell_deviatoric_stress(analyticalsolution, interfaceradius)
    shellpotential = potential(shellpressure, shelldevstress, K1, mu1, V01)

    corepressure = core_pressure(analyticalsolution)
    coredevstress = core_deviatoric_stress(analyticalsolution)
    corepotential = potential(corepressure, coredevstress, K2, mu2, V02)

    potentialdifference = 1.0 + (shellpotential - corepotential)/ΔG0

    return potentialdifference
end

interfaceposition = range(0.0,stop=1.0,length=100)
potentialdifference = normalized_potential_difference.(interfaceposition)

fig,ax = PyPlot.subplots()
ax.plot(interfaceposition,potentialdifference)
ax.grid()
ax.set_ylim(0.9,1.1)
fig
