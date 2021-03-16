function pressure(stress)
    return -1.0 / 3.0 * sum(stress)
end

function deviatoric_stress(stress)
    p = pressure(stress)
    return stress .+ p
end
