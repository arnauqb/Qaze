using Roots
using Optim: optimize, Brent
const ALPHA = 0.6
export cak_surface_mloss, cak_density, cak_normalized_mdot, cak_characteristic_mloss

function compute_temperature(r_0, wind)
    r0_arg = get_index(wind.grids.disk_range, r_0)
    localmdot = wind.grids.mdot[r0_arg]
    if localmdot == 0
        println("WARNING")
        Throw(DomainError)
    end
    t4 = wind.sed.disk_nt_temperature4(r_0)  / wind.bh.mdot * localmdot
    T = t4^(1/4) #SS_temperature(r_0, wind.bh)
    return T
end

"Nozzle function defined in Pereyra et al. (2004) (Paper I)"
function cak_nozzle_function(z, r_0, wind::WindStruct)
    x = z[1] / r_0
    T = compute_temperature(r_0, wind) 
    b = thermal_velocity(25e3) * C
    M = wind.bh.M * M_SUN
    R0 = r_0 * wind.bh.R_g
    numerator = (1 + x^2)^(1 / ALPHA)
    denom_1 = x / (1+x^2)^(3/2)
    denom_2 = - SIGMA_E * SIGMA_SB * T^4 / (G * M * C) * R0^2
    denom_3 = - 4 * b^2 * R0 * x / (G * M)
    denom = (denom_1 + denom_2 + denom_3)^((1-ALPHA) / ALPHA)
    return numerator / denom
end

function cak_characteristic_mloss(r_0, wind::WindStruct, K=0.03)
    T = compute_temperature(r_0, wind) 
    b = thermal_velocity(25e3) * C
    M = wind.bh.M * M_SUN
    R0 = r_0 * wind.bh.R_g
    constant = ALPHA * (1-ALPHA)^((1-ALPHA)/ALPHA) / (b * SIGMA_E)
    term_1 = G * M / R0^2
    term_2 = (SIGMA_E * SIGMA_SB * T^4 * K * R0^2 / (C * G * M))^(1/ALPHA)
    return constant * term_1 * term_2
end

function cak_normalized_mdot(r_0, wind::WindStruct)
    f(z) = cak_nozzle_function(z, r_0, wind)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    return mdot
end

function cak_surface_mloss(r_0, wind::WindStruct, K=0.03)
    f(z) = cak_nozzle_function(z, r_0, wind)
    mdot = optimize(f, 0, 2 * r_0, Brent()).minimum
    Sigma = cak_characteristic_mloss(r_0, wind, K) * mdot
    return Sigma
end

function cak_density(r_0, wind::WindStruct, K=nothing)
    if K === nothing
        K = wind.config["radiation"]["cak_K"]
    end
    Sigma = cak_surface_mloss(r_0, wind, K)
    T = compute_temperature(r_0, wind) 
    b = thermal_velocity(T) * C
    return Sigma / b / M_P
end