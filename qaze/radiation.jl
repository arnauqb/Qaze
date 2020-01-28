using Cubature
include("constants.jl")
include("utils.jl")
include("structures.jl")

function eddington_luminosity(bh::BlackHole)
    constant = 4 * pi * M_P * C^3 / SIGMA_T
    return constant * bh.R_g
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

function tau_uv_disk_blob(wind::Wind, r_d, phi_d, r, z)
    grids = wind.grids
    line_length = sqrt(r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    r_arg = searchsortedlast(grids.r_range, r)
    z_arg = searchsortedlast(grids.z_range, z)
    r_d_arg = searchsortedlast(grids.r_range, r_d)
    dr = abs(r_arg - r_d_arg)
    dz = abs(z_arg - 1)
    length = max(dr, dz) + 1
    line_coords = drawline(r_d_arg, 1, r_arg, z_arg)
    tau = 0.
    for row in eachrow(line_coords)
        i,j = row
        density = grids.density[i,j]
        fm = grids.fm[i,j]
        tau += density * (1 + fm)
    end
    tau = tau / length * line_length * wind.bh.R_g * SIGMA_T
    return tau
end

function integrate_z_kernel(x, r, z, wind::Wind)
    r_d, phi_d = x
    r_d_arg = searchsortedlast(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = exp(-tau_uv) / delta^2
    f_uv = wind.grids.uv_fraction[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end
function integrate_z(r, z, wind::Wind)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    (val, err) = hcubature(x->integrate_z_kernel(x, r, z, wind),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["adaptive_integral_epsrel"],
            abstol=0.,
            )
    val *= 2 * z^2
    return val
end

function integrate_r_kernel(x, r, z, wind::Wind)
    r_d, phi_d = x
    r_d_arg = searchsortedlast(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = (r - r_d * cos(phi_d)) * exp(-tau_uv) / delta^2
    f_uv = wind.grids.uv_fraction[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end

function integrate_r(r, z, wind::Wind)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    (val, err) = hcubature(x->integrate_r_kernel(x, r, z, wind),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["adaptive_integral_epsrel"],
            abstol=0.,
            )
    val *= 2 * z
    return val
end

function opacity_x(xi)
    if xi <= 1e5
        return 100
    else
        return 1
    end
end

function tau_x(r, z, wind::Wind)
    r_arg = searchsortedlast(wind.grids.r_range, r)
    z_arg = searchsortedlast(wind.grids.z_range, z)
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
    return tau
end







