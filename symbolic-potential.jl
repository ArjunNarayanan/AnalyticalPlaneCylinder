using Symbolics

@variables C l m K g V0

p = 1//3*C*(l+2m) + C*K*g
srr = 1//3*C*(4l+5m) - 1//3*C*m*g
stt = -1//3*C*(2l+m) - 1//3*C*m*g
szz = -2//3*C*(l+2m) + 2//3*C*m*g
V = simplify(V0*(1 - p/K))

devnorm = simplify(srr^2 + stt^2 + szz^2)
