function bulk_modulus(lambda, mu)
    return lambda + 2mu / 3
end

function lame_lambda(K,mu)
    return K - 2mu/3
end

function complex_potential(lambda, mu)
    K = bulk_modulus(lambda, mu)
    C1 = K * mu / (2 * (lambda + mu) * (lambda + 2mu))
    C2 = (lambda + 2mu) / mu
    C3 = K / (2 * (lambda + 2mu))

    p =
        -K / 2 - 2 * (lambda + mu) * C1 * C2 * (C1 * C2 - 1) +
        2mu * C3 * (1 + C3) - 4mu * C1 * C2 * C3
end

function simple_potential(lambda,mu)
    K = bulk_modulus(lambda,mu)

    p = K*mu*(36*K^2+51K*mu+40*mu^2)/(54*(K+mu/3)*(K+4mu/3)^2)
    return p
end


K = 247.0
mu = 126.0
lambda = lame_lambda(K, mu)

p1 = complex_potential(lambda,mu)
p2 = simple_potential(lambda,mu)
