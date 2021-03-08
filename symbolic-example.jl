using Symbolics

@variables t x y
D = Differential(t)

z = t+t^2
D(z)

B = simplify.([t^2+t+t^2    2t+4t
               x+y+y+2t     x^2 - x^2 + y^2])
