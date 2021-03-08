import sympy
from sympy import symbols, Rational

theta0, lamda, mu, gamma, V01, V02 = symbols("theta0 lambda mu gamma V01 V02")

K = lamda + Rational(2,3)*mu
C = K*theta0*mu/((lamda+mu)*(lamda+2*mu))

corepressure = -C*K*(1-gamma)

coresrr = Rational(1,3)*C*mu*(1-gamma)
corestt = coresrr
coreszz = -2*coresrr

V1 = V01*(1 - corepressure/K)
coredevnorm = coresrr**2 + corestt**2 + coreszz**2

corephi = corepressure*V01 - (corepressure**2)*V01/(2*K) - V1*coresrr + V01/(4*mu)*coredevnorm


shellpressure = Rational(1,3)*C*(lamda+2*mu) + C*K*gamma

V2 = V02*(1 - shellpressure/K)

shellsrr = Rational(1,3)*C*(4*lamda+5*mu) - Rational(1,3)*C*mu*gamma
shellstt = -Rational(1,3)*C*(2*lamda +mu) - Rational(1,3)*C*mu*gamma
shellszz = -Rational(2,3)*C*(lamda+2*mu) + Rational(2,3)*C*mu*gamma

shelldevnorm = shellsrr**2 + shellstt**2 + shellszz**2

shellphi = shellpressure*V02 - (shellpressure**2)*V02/(2*K) - V2*shellsrr + V02/(4*mu)*shelldevnorm
