using Roots
export cak_density, cak_z0, cak_z0_aux

function cak_z0_aux(z, r, wind::WindStruct)
    gravity_z = gravity_cgs(r, z, wind.bh)[2]
    T_gas = wind.sed.disk_core_temperature(r) 
    cs2 = K_B * T_gas  / (wind.config["disk"]["mu"] * M_P) 
    gas_pressure = cs2 / (z * wind.bh.R_g)
    T_rad = SS_temperature(r, wind.bh)
    rad_pressure = (SIGMA_T * SIGMA_SB) / (M_P * C) * T_rad^4
    return gravity_z + gas_pressure + rad_pressure
end

function local_eddington_ratio(r, z, wind::WindStruct)
    T_rad = SS_temperature(r, wind.bh)
    rad_pressure = (SIGMA_T * SIGMA_SB) / (M_P * C) * T_rad^4
    T_gas = wind.sed.disk_core_temperature(r) 
    cs2 = K_B * T_gas  / (wind.config["disk"]["mu"] * M_P) 
    gas_pressure = cs2 / (z * wind.bh.R_g)
    return rad_pressure / gas_pressure
end

function cak_z0(r, wind::WindStruct)
    #z0 = find_zero(z->cak_z0_aux(z, r, wind), (1, 100.), Bisection())
    z0 = fzero(z->cak_z0_aux(z, r, wind), 50.) 
    return z0
end

function cak_density(r, wind::WindStruct)
    ALPHA = 0.6
    K = 0.03
    z_0 = cak_z0(r, wind)
    gamma = local_eddington_ratio(r, z_0, wind)
    T = SS_temperature(r, wind.bh)
    constant = 1 / (SIGMA_T / M_P * thermal_velocity(T)) * G * wind.bh.M * M_SUN * ALPHA 
    constant *= (1. - ALPHA)^((1-ALPHA) / ALPHA)
    a = z_0 / r^3 / wind.bh.R_g^2
    b = (K * gamma)^(1 / ALPHA)
    return constant * a * b
end