module HydrostaticSolver

function bulk_modulus(lambda, mu)
    return lambda + 2mu / 3
end

function coefficient_matrix(
    inner_radius,
    outer_radius,
    ls,
    ms,
    lc,
    mc,
)
    R = inner_radius
    b = outer_radius
    R2 = R^2
    b2 = b^2

    row1 = [R, -R, -1 / R, 0]
    row2 = [2(lc + mc), -2(ls + ms), 2ms / R2, lc - ls]
    row3 = [0, 2(ls + ms), -2ms / b2, ls]
    row4 = [
        2lc * R2,
        2ls * (b2 - R2),
        0,
        (lc + 2mc)*R2 + (ls + 2ms) * (b2 - R2),
    ]

    op = vcat(row1', row2', row3', row4')
end

function coefficient_rhs(inner_radius, outer_radius, ls, ms, theta0)
    K = bulk_modulus(ls, ms)
    v = [
        0,
        -K * theta0,
        K * theta0,
        K * theta0 * (outer_radius^2 - inner_radius^2),
    ]
    return v
end

struct CylindricalSolver
    inner_radius::Any
    outer_radius::Any
    ls::Any
    ms::Any
    lc::Any
    mc::Any
    theta0::Any
    A1c::Any
    A1s::Any
    A2s::Any
    B::Any
    function CylindricalSolver(
        inner_radius,
        outer_radius,
        ls,
        ms,
        lc,
        mc,
        theta0,
    )
        m = coefficient_matrix(
            inner_radius,
            outer_radius,
            ls,
            ms,
            lc,
            mc,
        )
        r = coefficient_rhs(inner_radius, outer_radius, ls, ms, theta0)
        A1c, A1s, A2s, B = m \ r
        new(
            inner_radius,
            outer_radius,
            ls,
            ms,
            lc,
            mc,
            theta0,
            A1c,
            A1s,
            A2s,
            B,
        )
    end
end

function core_strain(solver::CylindricalSolver)
    return [solver.A1c,solver.A1c,solver.B]
end

function shell_strain(solver::CylindricalSolver,r)
    err = solver.A1s - solver.A2s/r^2
    ett = solver.A1s + solver.A2s/r^2
    ezz = solver.B
    return [err,ett,ezz]
end

end
