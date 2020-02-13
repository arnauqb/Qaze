export tau_uv_disk_blob, integrate_kernel, integrate_notau_kernel, integrate, compute_tauuv_tree
using Cubature
#using HCubature

function tau_uv_disk_blob_old(wind::WindStruct, r_d, phi_d, r, z, return_list=false)
    if (z<=0. || r < 0)
        return 0
    end
    rp_arg = get_index(wind.grids.r_range, r_d)
    rd_arg = get_index(wind.grids.r_range, r_d)
    zp_arg = get_index(wind.grids.z_range, 0.)
    r_arg = get_index(wind.grids.r_range, r)
    z_arg = get_index(wind.grids.z_range, z)
    compute_rp(zp) = r_d + (zp / z * (r - r_d))
    compute_zp(rp) = z / (r - r_d)  * (rp - r_d)
    deltad = 0.
    r1 = r_d
    z1 = 0.
    z2 = 0.
    r2 = r1
    tau = 0.
    #r_list = [convert(Float64,r1)]
    #z_list = [convert(Float64,z1)]
    step_r = convert(Int, sign(r - r_d))
    lambda_r = 0.
    lambda_z = 0.
    r2_candidate = 0.
    z2_candidate = 0.
    line_length = abs(r_arg - rd_arg) + abs(z_arg - zp_arg)
    delta_d_total = 0.
    counter = 1
    if rd_arg == r_arg
        vertical = true
    else
        vertical = false
    end
    while(true)
        if counter >= line_length
            break
        end
        r1 = r2
        z1 = z2
        density = wind.grids.density[rp_arg, zp_arg]
        fm = wind.grids.fm[rp_arg, zp_arg]
        if vertical
            zp_arg += 1
            z2 = wind.grids.z_range[zp_arg]
            deltad = (z2-z1) * wind.bh.R_g
            delta_d_total += deltad
            tau += density * SIGMA_T * (1 + fm) * deltad
            counter += 1
            #push!(r_list, r1)
            #push!(r_list, z1)
        else
            try
                r2_candidate = wind.grids.r_range[rp_arg + step_r]
                lambda_r = (r2_candidate - r1) / (r - r_d)
            catch
                lambda_r = Inf32
            end
            try
                z2_candidate = wind.grids.z_range[zp_arg + 1]
                lambda_z = (z2_candidate - z1) / z
            catch
                lambda_z = Inf32
            end
            if lambda_r < lambda_z
                r2 = r2_candidate 
                z2 = compute_zp(r2)
                rp_arg += step_r
            elseif lambda_r == lambda_z
                r2 = r2_candidate 
                z2 = z2_candidate 
                rp_arg += step_r
                zp_arg += 1
                counter += 1
            elseif lambda_r > lambda_z
                z2 = z2_candidate 
                r2 = compute_rp(z2)
                zp_arg += 1
            end
            counter +=1
            deltad = sqrt((r2-r1)^2 + (z2-z1)^2) * wind.bh.R_g
            delta_d_total += deltad
            #push!(r_list, r1)
            #push!(z_list, z1)
            tau += density * SIGMA_T * (1. + fm) * deltad 
        end
    end
    # add last bit
    density = wind.grids.density[r_arg, z_arg]
    fm = wind.grids.fm[r_arg, z_arg]
    if !vertical
        deltad = sqrt((r-r2)^2 + (z-z2)^2) * wind.bh.R_g
        delta_d_total += deltad
    else
        delta_d_total += (z-z2) * wind.bh.R_g
    end
    tau += density * SIGMA_T * deltad * (1 + fm) 
    # normalize
    #push!(r_list, r)
    #push!(z_list, z)
    delta = sqrt(r^2 + z^2 + r_d^2 - 2 * r * r_d * cos(phi_d))
    try
        if vertical 
            @assert isapprox(delta_d_total, z * wind.bh.R_g,  atol=0, rtol=1e-1)
        else
            asd = sqrt((r-r_d)^2 + z^2)
            @assert isapprox(delta_d_total, asd * wind.bh.R_g,  atol=0, rtol=1e-1)
        end
    catch
        println("distances do not match!")
        println("r_d: $r_d, phi_d: $phi_d r:$r, z:$z")
        asd = sqrt((r-r_d)^2 + z^2)
        println("z : $z")
        println(delta_d_total / wind.bh.R_g)
        println(asd)
    end
    tau = tau / sqrt((r-r_d)^2 + z^2) * delta
    return tau
    #if return_list
    #    return tau, r_list, z_list
    #else
    #    return tau
    #end
end

function compute_tauuv_leaf(point, intersection, leaf, wind::WindStruct)
    deltad = distance2d(point, intersection) * wind.bh.R_g
    density = leaf.data[1]
    tauuv = density * deltad * SIGMA_T
    return tauuv
end

function tau_uv_disk_blob(r_d, phi_d, r, z, wind::WindStruct)
    backwards = false
    if r_d > r
        backwards = true
    end
    point1 = [r_d, 0.]
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
    if backwards
        currentpoint[1] -= 1e-8
    end
    currentleaf = findleaf(wind.quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf, wind)
        currentpoint = intersection
        if backwards
            currentpoint[1] -= 1e-8
        end
        currentleaf = findleaf(wind.quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf, wind)
    delta = sqrt(r^2 + z^2 + r_d^2 - 2 * r * r_d * cos(phi_d))
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
    return tauuv
end

function integrate_kernel(v, r_d, phi_d, r, z, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    #tau_uv = tau_uv_disk_blob(wind, r_d, phi_d, r, z)
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