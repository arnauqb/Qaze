export tau_uv_disk_blob,
       tau_uv_disk_blob_fromstreamline,
       integrate_kernel, integrate_notau_kernel, integrate,
       integrate_fromstreamline_kernel, integrate_fromstreamline
using Cubature

"Computes the UV optical depth inside a leaf, taking into account the opacity
boost by the force multiplier."
function compute_tauuv_leaf(point, intersection, leaf, wind::WindStruct)
    deltad = distance2d(point, intersection) * wind.bh.R_g
    density = leaf.data[1]
    fm = leaf.data[2]
    tauuv = density * deltad * SIGMA_T * (1.0 + fm)
    return tauuv
end

"""
Compute the UV optical depth from a disc patch located at (r_d, phi_d),
until a gas element at (r,z). 
"""
function tau_uv_disk_blob(r_d, phi_d, r, z, wind::WindStruct)
    r_d > r ? backwards = true : backwards = false
    point1 = [r_d, wind.z_0]
    point1leaf = findleaf(wind.quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(wind.quadtree, point2)
    if point1leaf == point2leaf
        tauuv = compute_tauuv_leaf(point1, point2, point1leaf, wind)
        delta = sqrt(r^2 + z^2 + r_d^2 - 2 * r * r_d * cos(phi_d))
        tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
        return tauuv
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    tauuv = compute_tauuv_leaf(point1, intersection, point1leaf, wind)
    currentpoint = intersection
    backwards && (currentpoint[1] -= 1e-8)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf, wind)
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(wind.quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf, wind)
    delta = sqrt(r^2 + z^2 + r_d^2 - 2 * r * r_d * cos(phi_d))
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
    return tauuv
end

function integrate_kernel(v, r_d, phi_d, r, z, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, wind)
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
    v[1] = aux * (r - r_d * cos(phi_d)) 
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
    val .*= [2 * z, 2 * z^2]
    return val
end

function tau_uv_disk_blob_fromstreamline(p, psi, r, r_d, z, delta, wind::WindStruct)
    r_d > r ? backwards=true : backwards=false
    point1 = [r_d, wind.z_0]
    point1leaf = findleaf(wind.quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(wind.quadtree, point2)
    if point1leaf == point2leaf
        tauuv = compute_tauuv_leaf(point1, point2, point1leaf, wind)
        tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
        return tauuv
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    tauuv = compute_tauuv_leaf(point1, intersection, point1leaf, wind)
    currentpoint = intersection
    backwards && (currentpoint[1] -= 1e-8)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf, wind)
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(wind.quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf, wind)
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
    return tauuv
end

function integrate_fromstreamline_notau_kernel(v, p, psi, r, z, wind, r_max)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    common_part = nt * mdot * f_uv * p / (p^2 + z^2)^2 / r_d^3
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline_kernel(v, p, psi, r, z, wind, r_max)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    delta = sqrt(p^2 + z^2)
    tauuv = tau_uv_disk_blob_fromstreamline(p, psi, r, r_d, z, delta, wind)
    common_part = nt * mdot * f_uv * p / delta^4 / r_d^3 * exp(-tauuv)
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline(r, z, wind::WindStruct; include_tau_uv=true)
    r_max = wind.config["disk"]["outer_radius"]
    xmin = (0., 0.)
    xmax = (r + r_max, pi)
    if include_tau_uv
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_kernel(v, x[1], x[2], r, z, wind, r_max),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol=0.,
        )
    else
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_notau_kernel(v, x[1], x[2], r, z, wind, r_max),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol=0.,
        )
    end
    return [2*z, 2*z^2] .* val
end

function integrate_split(r, z, wind::WindStruct; include_tau_uv=true, eps=0)
    if include_tau_uv
        xmin = (wind.config["disk"]["inner_radius"], 0.)
        xmax = (r-eps, pi)
        (val1, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
        xmin = (r+eps, 0.)
        xmax = (wind.config["disk"]["outer_radius"], pi)
        (val2, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
        val = val1 + val2
    else
        xmin = (wind.config["disk"]["inner_radius"], 0.)
        xmax = (r-eps, pi)
        (val1, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
        xmin = (r+eps, 0.)
        xmax = (wind.config["disk"]["outer_radius"], pi)
        (val2, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
        val = val1 + val2
    end
    val .*= [2 * z, 2 * z^2]
    return val
end
