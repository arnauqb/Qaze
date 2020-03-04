using Roots
export cak_density, cak_z0, cak_z0_aux, local_eddington_ratio, radiation_pressure, gas_pressure, cak_mdot

function radiation_pressure(z, r, wind::WindStruct)
    #T_rad = SS_temperature(r, wind.bh)
    T_rad = wind.sed.disk_nt_temperature4(r)^(0.25)
    rad_pressure = (SIGMA_T * SIGMA_SB) / (M_P * C) * T_rad^4
    return rad_pressure
end

function gas_pressure(z, r, wind::WindStruct)
    T_gas = wind.sed.disk_core_temperature(r) 
    cs2 = K_B * T_gas  / (wind.config["disk"]["mu"] * M_P) 
    gas_pressure = cs2 / (z * wind.bh.R_g)
end

function cak_z0_aux(z, r, wind::WindStruct)
    gravity_z = gravity_cgs(r, z, wind.bh)[2]
    return gravity_z + radiation_pressure(z, r, wind) + gas_pressure(z,r,wind)
end

function local_eddington_ratio(r, z, wind::WindStruct)
    gravity_z = gravity_cgs(r, z, wind.bh)[2]
    return radiation_pressure(z, r, wind) / abs(gravity_z)
end

function cak_z0(r, wind::WindStruct)
    #z0 = find_zero(z->cak_z0_aux(z, r, wind), 50., ) 
    z0 = minimum(find_zeros(z->cak_z0_aux(z, r, wind), 0, 100))
    try
        @assert z0 > 0
    catch
        println("CAK z0 unphyiscal")
        println(r)
        println(z0)
        throw(DomainError)
    end
    return z0
end

"Mass loss rate from Pereyra 2015"
function cak_mdot(r, wind::WindStruct)
    ALPHA = 0.6
    #K = 0.003
    K = 0.03
    T = SS_temperature(r, wind.bh)
    vth = thermal_velocity(T) * C
    z0 = 3 * SIGMA_T / M_P * mass_accretion_rate(wind.bh) / (8 * π * C)
    factor1 = G * wind.bh.M * M_SUN * z0 / (SIGMA_T / M_P * vth * r^3 * wind.bh.R_g^3)
    fmfactor = ALPHA * (1-ALPHA)^((1.0-ALPHA)/ALPHA) * K^(1/ALPHA)
    mdot = factor1 * fmfactor
    return mdot
end

function cak_density(r, wind::WindStruct)
    mdot = cak_mdot(r, wind)
    v_th = thermal_velocity(SS_temperature(r, wind.bh)) * C
    n =   mdot / v_th / (wind.config["disk"]["mu"] * M_P) 
    try
        @assert n >= 0
    catch
        println(mdot)
        println(v_th)
        println(wind.config["disk"]["mu"])
        throw(DomainError)
    end
    return n
end

# Kure model

