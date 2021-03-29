include("moduli-conversion.jl")

function explicit_solution(t,a,b)
    gamma = b/a
    return  (1.0 - tanh(a*b*t)/gamma)/(1 - gamma*tanh(a*b*t))
end

K = 247.0e9
mu = 126.0e9
rhos = 3.93e3           # Kg/m^3
rhoc = 3.68e3
V0s = 1.0 / rhos
V0c = 1.0 / rhoc
theta0 = -0.067
V0 = 0.5*(V0c+V0s)

ΔG0Jmol = -14351.0
molarmass = 0.147
ΔG0 = ΔG0Jmol / molarmass

M = K*mu/(K+4mu/3)
N = (6K + 5mu)/(9K+3mu)
L = mu / (3K + mu)

beta2 = M * V0 * (1 - L) * theta0^2
alpha2 = M * V0 * (1 - N) * theta0^2
gamma2 = alpha2/beta2

G = -ΔG0/beta2
a2 = 1 + gamma2/G
b2 = 1.0/G

psa = sqrt(a2)
psb = sqrt(b2)

C = mu*K/(K+4mu/3)
beta2 = 4/3*C*V0*theta0^2
alpha2 = 2/3*C*V0*theta0^2
gamma2 = alpha2/beta2

G = -ΔG0/beta2
a2 = 1 + gamma2/G
b2 = 1.0/G

hsa = sqrt(a2)
hsb = sqrt(b2)

pTmax = 1.0/(2*psa*psb)*log((psa+psb)/(psa-psb))
hTmax = 1.0/(2*hsa*hsb)*log((hsa+hsb)/(hsa-hsb))

pt = 0:1e-2:pTmax
ht = 0:1e-2:hTmax
psol = explicit_solution.(pt,psa,psb)
hsol = explicit_solution.(ht,hsa,hsb)

using PyPlot
fig,ax = PyPlot.subplots()
ax.plot(pt,psol,color="black",linestyle="dotted",label="plane-strain")
ax.plot(ht,hsol,color="black",label="homogeneous")
ax.plot([0,1],[1,0],"--",color="black")
ax.legend()
ax.grid()
fig.savefig("comapre-plane-strain-homogeneous.png")
