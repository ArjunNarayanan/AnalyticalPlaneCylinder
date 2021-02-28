function bulk_modulus(l, m)
    return l + 2m / 3
end

function lame_lambda(k, m)
    return k - 2m / 3
end

function analytical_coefficient_matrix(inradius, outradius, ls, ms, lc, mc)
    a = zeros(3, 3)
    a[1, 1] = inradius
    a[1, 2] = -inradius
    a[1, 3] = -1.0 / inradius
    a[2, 1] = 2 * (lc + mc)
    a[2, 2] = -2 * (ls + ms)
    a[2, 3] = 2ms / inradius^2
    a[3, 2] = 2(ls + ms)
    a[3, 3] = -2ms / outradius^2
    return a
end

function analytical_coefficient_rhs(ls, ms, theta0)
    r = zeros(3)
    Ks = bulk_modulus(ls, ms)
    r[2] = -Ks * theta0
    r[3] = Ks * theta0
    return r
end

struct CylindricalSolver
    inradius::Any
    outradius::Any
    A1c::Any
    A1s::Any
    A2s::Any
    ls::Any
    ms::Any
    lc::Any
    mc::Any
    theta0::Any
    function CylindricalSolver(
        inradius,
        outradius,
        ls,
        ms,
        lc,
        mc,
        theta0,
    )
        a = analytical_coefficient_matrix(inradius, outradius, ls, ms, lc, mc)
        r = analytical_coefficient_rhs(ls, ms, theta0)
        coeffs = a \ r
        new(
            inradius,
            outradius,
            coeffs[1],
            coeffs[2],
            coeffs[3],
            ls,
            ms,
            lc,
            mc,
            theta0,
        )
    end
end

function radial_displacement(A::CylindricalSolver, r)
    if r <= A.inradius
        return A.A1c * r
    else
        return A.A1s * r + A.A2s / r
    end
end

function shell_radial_stress(ls, ms, theta0, A1, A2, r)
    return (ls + 2ms) * (A1 - A2 / r^2) + ls * (A1 + A2 / r^2) -
           (ls + 2ms / 3) * theta0
end

function shell_circumferential_stress(ls, ms, theta0, A1, A2, r)
    return ls * (A1 - A2 / r^2) + (ls + 2ms) * (A1 + A2 / r^2) -
           (ls + 2ms / 3) * theta0
end

function shell_out_of_plane_stress(ls, ms, A1, theta0)
    return 2 * ls * A1 - (ls + 2ms / 3) * theta0
end

function core_in_plane_stress(lc, mc, A1)
    return (lc + 2mc) * A1 + lc * A1
end

function core_out_of_plane_stress(lc, A1)
    return 2 * lc * A1
end

function pressure(stress)
    return -1.0 / 3.0 * sum(stress)
end

function deviatoric_stress(stress)
    p = pressure(stress)
    return stress .+ p
end

function shell_stress(A::CylindricalSolver, r)

    srr = shell_radial_stress(A.ls, A.ms, A.theta0, A.A1s, A.A2s, r)
    stt = shell_circumferential_stress(A.ls, A.ms, A.theta0, A.A1s, A.A2s, r)
    s33 = shell_out_of_plane_stress(A.ls, A.ms, A.A1s, A.theta0)

    return [srr, stt, s33]
end

function shell_pressure(A::CylindricalSolver, r)
    stress = shell_stress(A, r)
    return pressure(stress)
end

function shell_deviatoric_stress(A::CylindricalSolver, r)
    stress = shell_stress(A, r)
    return deviatoric_stress(stress)
end

function core_in_plane_stress(A::CylindricalSolver)
    return core_in_plane_stress(A.lc, A.mc, A.A1c)
end

function core_stress(A::CylindricalSolver)
    s11 = core_in_plane_stress(A.lc, A.mc, A.A1c)
    s33 = core_out_of_plane_stress(A.lc, A.A1c)
    return [s11, s11, s33]
end

function core_pressure(A::CylindricalSolver)
    stress = core_stress(A)
    return pressure(stress)
end

function core_deviatoric_stress(A::CylindricalSolver)
    stress = core_stress(A)
    return deviatoric_stress(stress)
end

function potential(
    pressure,
    devstress,
    bulkmodulus,
    shearmodulus,
    specificvolume0,
)
    specificvolume = specificvolume0 * (1.0 - pressure / bulkmodulus)
    devstressnorm = sum(devstress .^ 2)

    p1 = pressure * specificvolume0
    p2 = -(pressure^2) * specificvolume0 / (2bulkmodulus)
    p3 = -specificvolume * devstress[1]
    p4 = specificvolume0 / (4shearmodulus) * devstressnorm

    return p1 + p2 + p3 + p4
end
