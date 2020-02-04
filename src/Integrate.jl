export tau_uv_disk_blob, integrate_kernel, integrate_notau_kernel, integrate
using Cubature
#using HCubature

function tau_uv_disk_blob(wind::WindStruct, r_d, phi_d, r, z)
    grids = wind.grids
    line_length = sqrt(r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    r_arg = get_index(grids.r_range, r)
    z_arg = get_index(grids.z_range, z)
    r_d_arg = get_index(grids.r_range, r_d)
    dr = abs(r_arg - r_d_arg)
    dz = abs(z_arg - 1)
    line_arg_length = max(dr, dz) + 1
    line_coords = drawline(r_d_arg, 1, max(r_arg-1,1), max(z_arg-1,1))
    tau = 0.
    for row in eachrow(line_coords)
        i,j = row
        r1 = wind.grids.r_range[i]
        r2 = wind.grids.r_range[min(i+1, wind.grids.n_r)]
        z1 = wind.grids.z_range[j]
        z2 = wind.grids.z_range[min(j+1, wind.grids.n_z)]
        delta_d = sqrt((r2-r1)^2 + (z2-z1)^2) * wind.bh.R_g
        density = grids.density[i,j]
        fm = grids.fm[i,j]
        tau += density * (1 + fm) * delta_d
    end
    #tau = tau / length * line_length * wind.bh.R_g * SIGMA_T
    tau = tau * SIGMA_T
    return tau
end

function integrate_kernel(v, r_d, phi_d, r, z, wind)
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

function integrate_notau_kernel(v, r_d, phi_d, r, z, wind)
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


function integrate(r, z, wind::WindStruct; include_tau_uv=true)
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    if include_tau_uv
        (val, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
    else
        (val, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
    end
   # if include_tau_uv
   #     (val, err) = hcubature(x ->integrate_kernel(x[0], x[1], r, z, wind),
   #             xmin,
   #             xmax,
   #             rtol = wind.config["radiation"]["integral_rtol"],
   #             atol=0.,
   #             )
   # else
   #     (val, err) = hcubature(x ->integrate_notau_kernel(x[0], x[1], r, z, wind),
   #             xmin,
   #             xmax,
   #             rtol = wind.config["radiation"]["integral_rtol"],
   #             atol=0.,
   #             )
    val .*= [2 * z, 2 * z^2]
    return val
end