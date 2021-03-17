module PlaneStressSolver
include("moduli-conversion.jl")

function coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
    R = inner_radius
    b = outer_radius
    R2 = R^2
    b2 = b^2

    row1 = [R, -R, -1 / R, 0, 0]
    row2 = [2 * (lc + mc), -2 * (ls + ms), 2ms / R2, lc, -ls]
    row3 = [0, 2 * (ls + ms), -2ms / b2, 0, ls]
    row4 = [2lc, 0, 0, (lc + 2mc), 0]
    row5 = [0, 2ls, 0, 0, ls + 2ms]

    return vcat(row1',row2', row3', row4', row5')
end

function coefficient_rhs(ls, ms, theta0)
    K = bulk_modulus(ls, ms)
    v = [0, -K * theta0, K * theta0, 0, K * theta0]
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
    B1c::Any
    B1s::Any
    function CylindricalSolver(
        inner_radius,
        outer_radius,
        ls,
        ms,
        lc,
        mc,
        theta0,
    )
        m = coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
        r = coefficient_rhs(ls, ms, theta0)
        A1c, A1s, A2s, B1c, B1s = m \ r
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
            B1c,
            B1s,
        )
    end
end

function core_strain(solver::CylindricalSolver)
    return [solver.A1c,solver.A1c,solver.B1c]
end

function shell_strain(solver::CylindricalSolver,r)
    err = solver.A1s - solver.A2s/r^2
    ett = solver.A1s + solver.A2s/r^2
    ezz = solver.B1s
    return [err,ett,ezz]
end

end
