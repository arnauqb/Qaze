export gravity, SS_temperature, NT_temperature, gravity_cgs, eddington_luminosity, bolometric_luminosity, mass_accretion_rate

"Gravitational force in units of c^2/R_g"
function gravity(r, z, bh::BlackHoleStruct)
    d = sqrt(r^2 + z^2)
    constant = -1. / d^2
    result = constant * [r / d, abs(z) / d]
    return result
end

"Accretion disc temperature profile following Shakura & Sunyaev model."
function SS_temperature(r, bh::BlackHoleStruct)
    a = 3 * G * bh.M * M_SUN * mass_accretion_rate(bh)
    b = 8 * pi * (r * bh.R_g)^3 * SIGMA_SB
    T = (a/b)^(0.25)
    return T
end

function NT_temperature(r, bh::BlackHoleStruct)
    T_SS = SS_temperature(r, bh)
    T_NT = T_SS * nt_rel_factors(r, bh.spin, bh.isco)
    return T_NT
end

"Gravitational force in cgs units."
function gravity_cgs(r, z, bh::BlackHoleStruct)
    return gravity(r,z,bh) *  C^2 / bh.R_g
end

function eddington_luminosity(bh::BlackHoleStruct)
    constant = 4 * pi * M_P * C^3 / SIGMA_T
    return constant * bh.R_g
end

function bolometric_luminosity(bh::BlackHoleStruct)
    bol_lumin = bh.mdot * eddington_luminosity(bh)
    return bol_lumin
end

function mass_accretion_rate(bh::BlackHoleStruct)
    bol_lumin = bolometric_luminosity(bh)
    Mdot = bol_lumin / (bh.eta * C^2)
    return Mdot
end