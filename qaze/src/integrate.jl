using Cubature
include("utils.jl")

function tau_uv_disk_blob(wind::Wind, r_d, phi_d, r, z)
    grids = wind.grids
    line_length = sqrt(r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    r_arg = get_index(grids.r_range, r)
    z_arg = get_index(grids.z_range, z)
    r_d_arg = get_index(grids.r_range, r_d)
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
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = exp(-tau_uv) / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end

function integrate_z_notau_kernel(x, r, z, wind::Wind)
    r_d, phi_d = x
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = 1 / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end

function integrate_z(r, z, wind::Wind ; include_tau_uv=true)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    reltol = wind.config["radiation"]["integral_rtol"]
    if include_tau_uv
        (val, err) = hcubature(x->integrate_z_kernel(x, r, z, wind),
                xmin,
                xmax,
                reltol = reltol,
                abstol=0.,
                )
    else
        (val, err) = hcubature(x->integrate_z_notau_kernel(x, r, z, wind),
                xmin,
                xmax,
                reltol = reltol,
                abstol=0.,
                )
    end
    val *= 2 * z^2
    return val
end

function integrate_r_kernel(x, r, z, wind::Wind)
    r_d, phi_d = x
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = (r - r_d * cos(phi_d)) * exp(-tau_uv) / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end

function integrate_r_notau_kernel(x, r, z, wind::Wind)
    r_d, phi_d = x
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part = (r - r_d * cos(phi_d)) / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part * phi_part
end

function integrate_r(r, z, wind::Wind ; include_tau_uv=true)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    reltol = wind.config["radiation"]["integral_rtol"]
    if include_tau_uv
        (val, err) = hcubature(x->integrate_r_kernel(x, r, z, wind),
                xmin,
                xmax,
                reltol = reltol,
                abstol=0.,
                )
    else
        (val, err) = hcubature(x->integrate_r_notau_kernel(x, r, z, wind),
                xmin,
                xmax,
                reltol = reltol,
                abstol=0.,
                )
    end
    val *= 2 * z
    return val
end

function integrate_kernel(x, v, r, z, wind)
    r_d, phi_d = x
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    v[1] = (r - r_d * cos(phi_d)) * aux
    v[2] = aux
end

function integrate_notau_kernel(x, v, r, z, wind)
    r_d, phi_d = x
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  1. / delta^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    v[1] = (r - r_d * cos(phi_d)) * aux
    v[2] = aux
end


function integrate(r, z, wind::Wind; include_tau_uv=true)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    if include_tau_uv
        (val, err) = hcubature(2, 
                (x,v) ->integrate_kernel(x, v, r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
    else
        (val, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(x, v, r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
    end
    val .*= [2 * z, 2 * z^2]
    return val
end