include("constants.jl")
include("utils.jl")
include("structures.jl")
include("integrate.jl")

function thermal_velocity(T, mu = 1.)
    v = sqrt(K_B * T / (mu * M_P)) / C
    return v
end

function eddington_luminosity(bh::BlackHole)
    constant = 4 * pi * M_P * C^3 / SIGMA_T
    return constant * bh.R_g
end

function initialize_uv_fraction(wind::Wind)
    if wind.config["radiation"]["disk_uv_fraction"]
        wind.grids.disk_range[:], wind.grids.uv_fractions[:] = wind.sed.compute_uv_fractions(
            inner_radius=wind.bh.disk_r_min,
            outer_radius=wind.bh.disk_r_max,
            return_all=false,
            log_spaced=true
        )
    else
        wind.grids.disk_r_range = 10 .^(range(log10(wind.bh.disk_r_min), stop = log10(wind.bh.disk_r_max), length=wind.grids.n_disk))
        wind.grids.uv_fractions = wind.radiation.f_uv .* ones(Float64, wind.grids.n_disk)
    end
end

function nt_rel_factors(r, astar, isco)
    yms = sqrt(isco)
    y1 = 2 * cos((acos(astar) - pi) / 3)
    y2 = 2 * cos((acos(astar) + pi) / 3)
    y3 = -2 * cos(acos(astar) / 3)
    y = sqrt(r);
    C = 1 - 3 / r + 2 * astar / r^1.5
    B = 3 * (y1 - astar)^2 * log((y - y1) / (yms - y1)) / (y * y1 * (y1 - y2) * (y1 - y3))
    B += 3 * (y2 - astar)^2 * log((y - y2) / (yms - y2)) / (y * y2 * (y2 - y1) * (y2 - y3))
    B += 3 * (y3 - astar)^2 * log((y - y3) / (yms - y3)) / (y * y3 * (y3 - y1) * (y3 - y2))
    A = 1 - yms / y - 3 * astar * log(y / yms) / (2 * y)
    factor = (A - B) / C
    return factor
end

function opacity_x(xi)
    if xi <= 1e5
        return 100
    else
        return 1
    end
end

function compute_tau_x(r, z, wind::Wind)
    r_arg = get_index(wind.grids.r_range, r)
    z_arg = get_index(wind.grids.z_range, z)
    line_coords = drawline(1,1,r_arg,z_arg)
    tau = 0.
    tau_length = size(line_coords)[1]-1
    #for (k, row) in enumerate(eachrow(line_coords)[:aux])
    for k in 1:tau_length
        row = line_coords[k,:]
        i, j = row
        rp = wind.grids.r_range[i]
        zp = wind.grids.z_range[j]
        d = sqrt(rp^2 + zp^2) * wind.bh.R_g
        density = wind.grids.density[i,j]
        xi0 = wind.radiation.xray_luminosity / (density * d^2)
        xi = xi0
        for dummy in 1:3
            tau_prov = (tau + density * opacity_x(xi)) * d / k * SIGMA_T
            xi = xi0 * exp(-tau_prov)
        end
        tau = tau + (density * opacity_x(xi))
    end
    d = sqrt(r^2 + z^2) * wind.bh.R_g
    tau =  tau / tau_length *  SIGMA_T * d
    @assert tau >= 0
    return tau
end

function ionization_parameter(r, z, density, wind::Wind)
    d2 = (r^2 + z^2) * wind.bh.R_g^2
    tau_x = compute_tau_x(r,z,wind)
    xi = wind.radiation.xray_luminosity * exp(-tau_x) / (density * d2) + 1e-20
    @assert xi >= 0
    return xi
end

function compute_tau_eff(density, dv_dr, v_th)
    if dv_dr == 0
        return 1.
    end
    t = density * SIGMA_T * abs(v_th / dv_dr)
    return t
end

function force_multiplier_k(xi)
    k = 0.03 + 0.385 * exp(-1.4 * xi^0.6)
    return k
end

function force_multiplier_eta(xi)
    if (log10(xi) < 0.5)
        aux = 6.9 * exp(0.16 * xi^0.4)
        eta_max = 10^aux
    else
        aux = 9.1 * exp(-7.96e-3 * xi)
        eta_max = 10^aux
    end
end

function force_multiplier(t, xi)
    ALPHA = 0.6
    TAU_MAX_TOL = 1e-3
    k = force_multiplier_k(xi)
    eta = force_multiplier_eta(xi)
    tau_max = t * eta
    if tau_max < TAU_MAX_TOL
        aux = (1 - ALPHA) * (tau_max^ALPHA)
    else
        aux = ((1 + tau_max)^(1-ALPHA) - 1) / ((tau_max)^(1-ALPHA))
    end
    fm = k * t^(-ALPHA) * aux
    @assert fm >= 0
    return fm
end

function force_radiation(r, z, fm, wind::Wind ; include_tau_uv = false)
    if z < 1e-3
        return [0.,0.]
    end
    int_values = integrate(r, z, wind, include_tau_uv=include_tau_uv)
    force = wind.radiation.force_constant * (1. + fm) * int_values
    return force
end
