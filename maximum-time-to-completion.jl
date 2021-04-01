using PyPlot
include("moduli-conversion.jl")

function exact_solution(t,f,c)
    t = tanh(c^2*f*t)
    return (1.0 - t/f)/(1.0 - t*f)
end

function plot_interface(t,R)
    fig,ax = PyPlot.subplots()
    ax.plot(t,R,color="black")
    ax.plot([0,1],[1,0],linestyle="dashed",color="black")
    ax.grid()
    return fig
end

K = 247.0e9
# mu = 126.0e9
mu = 3/2*K
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc
lambda = lame_lambda(K, mu)
theta0 = -0.067

V0 = 0.5*(V0c+V0s)
ΔG0Jmol = -14351.0
molarmass = 0.147
ΔG0 = ΔG0Jmol/molarmass

M = K * mu / (K + 4mu / 3)
N = (6K + 5mu) / (9K + 3mu)
L = mu / (3K + mu)

beta2 = M * V0 * (1 - L) * theta0^2
alpha2 = M * V0 * (1 - N) * theta0^2
gamma2 = alpha2/beta2

G = -ΔG0/beta2

c2 = 1 + inv(G)*gamma2
g2 = inv(G)

c = sqrt(c2)
g = sqrt(g2)

f = g/c

Tmax = 2.0
trange = 0:1e-3:Tmax

R = exact_solution.(trange,f,c)

plot_interface(trange,R)
