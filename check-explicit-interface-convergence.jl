using DifferentialEquations
include("moduli-conversion.jl")

function interface_derivative(R,p,t)
    a2,b2 = p[1],p[2]

    return -(a2 - b2*R^2)
end

function explicit_solution(t,a,b)
    gamma = b/a
    return  (1.0 - tanh(a*b*t)/gamma)/(1 - gamma*tanh(a*b*t))
end

function solve_and_compute_error(dt,Tmax,a2,b2)
    tspan = (0.0,Tmax)
    prob = ODEProblem(interface_derivative,1.0,tspan,[a2,b2])
    sol = solve(prob,RK4(),dt=dt,adaptive=false)

    a = sqrt(a2)
    b = sqrt(b2)

    exactsol = explicit_solution.(sol.t,a,b)

    maxerror = maximum(abs.(exactsol-sol.u))
    return maxerror
end

function convergence_rate(err,dx)
    return (diff(log.(err))) ./(diff(log.(dx)))
end

K = 247.0e9
mu = 126.0e9
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc
theta0 = -0.067

ΔG0Jmol = -14351.0
molarmass = 0.147
ΔG0 = ΔG0Jmol / molarmass

M = K*mu/(K+4mu/3)
V0 = 0.5*(V0c+V0s)
N = (6K + 5mu)/(9K+3mu)
L = mu / (3K + mu)
V0 = 0.5 * (V0s + V0c)

beta2 = M * V0 * (1 - L) * theta0^2
alpha2 = M * V0 * (1 - N) * theta0^2
gamma2 = alpha2/beta2

G = -ΔG0/beta2
a2 = 1 + gamma2/G
b2 = 1.0/G

a = sqrt(a2)
b = sqrt(b2)

outer_radius = 1.0

Tmax = 1.0/(2*a*b)*log((a+b)/(a-b))

powers = [4,5,6,7,8,9,10]
dt = 1.0 ./ (2 .^ powers)

err = solve_and_compute_error.(dt,Tmax,a2,b2)
rate = convergence_rate(err,dt)

# tspan = (0.0,Tmax)
# prob = ODEProblem(interface_derivative,1.0,tspan,[a2,b2])
# sol = solve(prob,RK4(),dt=0.1,adaptive=false)
#
# numericalsol = sol.u
# exactsol = explicit_solution.(sol.t,a,b)

# using PyPlot
# fig,ax = PyPlot.subplots()
# ax.plot(sol.t,exactsol,color="black",label="exact")
# ax.scatter(sol.t,numericalsol,color="black",label="numerical")
# ax.plot([0,1],[1,0],"--",color="black")
# ax.legend()
# ax.grid()
# fig.savefig("plane-strain-exact-vs-numerical.png")
