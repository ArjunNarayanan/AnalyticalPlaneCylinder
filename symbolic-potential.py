import sympy

C, l, m, K, g, V0 = sympy.symbols("C l m K g V0")

p = 1/3*C*(l+2*m) + C*K*g
srr = 1/3*C*(4*l + 5*m) - 1/3*C*m*g
stt = -1/3*C*(2*l+m) - 1/3*C*m*g
szz = -2/3*C*(l+2*m) + 2/3*C*m*g

V = V0*(1 - p/K)

devnorm = srr**2 + stt**2 + szz**2

phi = p*V0 - (p**2)*V0/(2*K) - V*srr + V0/(4*m)*devnorm
