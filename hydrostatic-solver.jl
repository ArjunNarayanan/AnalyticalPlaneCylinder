module HydrostaticSolver

function bulk_modulus(lambda, mu)
    return lambda + 2mu / 3
end

function coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
    R = inner_radius
    b = outer_radius
    R2 = R^2
    b2 = b^2

    row1 = [R, -R, -1 / R, 0]
    row2 = [2(lc + mc), -2(ls + ms), 2ms / R2, lc - ls]
    row3 = [0, 2(ls + ms), -2ms / b2, ls]
    row4 =
        [2lc * R2, 2ls * (b2 - R2), 0, (lc + 2mc) * R2 + (ls + 2ms) * (b2 - R2)]

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
        m = coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
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
    return [solver.A1c, solver.A1c, solver.B]
end

function shell_strain(solver::CylindricalSolver, r)
    err = solver.A1s - solver.A2s / r^2 - solver.theta0 / 3
    ett = solver.A1s + solver.A2s / r^2 - solver.theta0 / 3
    ezz = solver.B - solver.theta0 / 3
    return [err, ett, ezz]
end

function core_stress(solver::CylindricalSolver)
    srr = 2 * (solver.lc + solver.mc) * solver.A1c + solver.lc * solver.B
    stt = srr
    szz = 2 * solver.lc * solver.A1c + (solver.lc + 2 * solver.mc) * solver.B

    return [srr, stt, szz]
end

function shell_stress(solver::CylindricalSolver, r)
    KT = bulk_modulus(solver.ls, solver.ms) * solver.theta0

    srr =
        2 * (solver.ls + solver.ms) * solver.A1s -
        2 * solver.ms * solver.A2s / r^2 + solver.ls * solver.B - KT
    stt =
        2 * (solver.ls + solver.ms) * solver.A1s +
        2 * solver.ms * solver.A2s / r^2 +
        solver.ls * solver.B - KT
    szz =
        2 * solver.ls * solver.A1s + (solver.ls + 2 * solver.ms) * solver.B - KT

    return [srr, stt, szz]
end

function core_strain_energy(solver::CylindricalSolver)
    stress = core_stress(solver)
    strain = core_strain(solver)
    return 0.5 * sum(stress .* strain)
end

function shell_strain_energy(solver::CylindricalSolver, r)
    stress = shell_stress(solver, r)
    strain = shell_strain(solver, r)

    return 0.5 * sum(stress .* strain)
end

function core_compression_work(solver::CylindricalSolver, V0)
    V = V0 * (1.0 + sum(core_strain(solver)))
    srr = core_stress(solver)[1]
    return V * srr
end

function shell_compression_work(solver::CylindricalSolver, r, V0)
    V = V0 * (1.0 + sum(shell_strain(solver, r)))
    srr = shell_stress(solver, r)[1]
    return V * srr
end



end
