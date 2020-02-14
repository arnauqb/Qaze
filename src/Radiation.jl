export thermal_velocity, eddington_luminosity, initialize_uv_fraction, nt_rel_factors, opacity_x, compute_taux_grid, ionization_parameter,
       compute_tau_eff, force_multiplier, force_radiation, compute_taux_leaf, compute_tau_x

function thermal_velocity(T, mu = 1.)
    v = sqrt(K_B * T / (mu * M_P)) / C
    return v
end

function initialize_uv_fraction(wind::WindStruct)
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

function compute_tau_x_old(r, z, wind::WindStruct)
    r_arg = get_index(wind.grids.r_range, r)
    z_arg = get_index(wind.grids.z_range, z)
    line_coords = drawline(1,1, r_arg, z_arg)
    println("r: $r, z: $z")
    println(line_coords)
    tau = 0.
    tau_length = size(line_coords)[1]
    #println("r: $r, z: $z")
    #println(line_coords)
    for k in 1:(tau_length-1)
        row = line_coords[k,:]
        i, j = row
        rp1 = wind.grids.r_range[i]
        rp2 = wind.grids.r_range[min(i+1, wind.grids.n_r)]
        zp1 = wind.grids.z_range[j]
        zp2 = wind.grids.z_range[min(j+1, wind.grids.n_z)]
        println("rp1: $rp1 rp2: $rp2, zp1: $zp1, zp2: $zp2")
        d2 = sqrt(rp2^2 + zp2^2) * wind.bh.R_g
        delta_d = sqrt((rp2-rp1)^2 + (zp2-zp1)^2) * wind.bh.R_g
        density = wind.grids.density[i,j]
        #println("r: $rp z :$zp den: $density")
        xi0 = wind.radiation.xray_luminosity / (density * d2^2)
        xi = xi0
        for dummy in 1:2
            dtau = density * opacity_x(xi) * delta_d * SIGMA_T
            tau_prov = tau + dtau
            xi = xi0 * exp(-tau_prov)
        end
        println(delta_d / wind.bh.R_g)
        tau = tau + density * opacity_x(xi) * delta_d * SIGMA_T
    end
    #d = sqrt(r^2 + z^2) * wind.bh.R_g
    @assert tau >= 0
    return tau
end

function compute_taux_grid(r, z, wind::WindStruct)
    if r <= wind.grids.r_range[1] || z < 0
        return 0.
    end
    r1 = 0.
    z1 = 0.
    r_arg = get_index(wind.grids.r_range, r)
    z_arg = get_index(wind.grids.z_range, z)
    m = z / r 
    deltad = 0.
    r2 = wind.grids.r_range[1]
    z2 = m * r2
    tau = sqrt(r2^2 + z2^2) * wind.bh.R_g * SIGMA_T * wind.grids.n_vacuum
    rp_arg = 1
    zp_arg = get_index(wind.grids.z_range, z2)
    lambda_r = 0.
    lambda_z = 0.
    max_steps = max(abs(r_arg - 1), abs(z_arg - 1))
    delta_d_total = sqrt(r2^2 + z2^2) * wind.bh.R_g
    counter = 1
    #while ((rp_arg < r_arg) || (zp_arg < z_arg))
    while(true)
        if counter >= max_steps
            break
        end
        r1 = r2
        z1 = z2
        density = wind.grids.density[rp_arg, zp_arg]
        try
            r2_candidate = wind.grids.r_range[rp_arg + 1]
            lambda_r = (r2_candidate - r1) / r
        catch
            lambda_r = Inf32
            nothing
        end
        try
            z2_candidate = wind.grids.z_range[zp_arg + 1]
            lambda_z = (z2_candidate - z1) / z
        catch
            lambda_z = Inf32
        end
        if lambda_r < lambda_z
            rp_arg += 1
            r2 = wind.grids.r_range[rp_arg]
            z2 = m * r2 
        elseif lambda_r == lambda_z
            rp_arg += 1
            zp_arg += 1
            r2 = wind.grids.r_range[rp_arg]
            z2 = wind.grids.z_range[zp_arg]
            counter += 1
        elseif lambda_r > lambda_z
            zp_arg += 1
            z2 = wind.grids.z_range[zp_arg]
            r2 = z2 / m
        end
        deltad = sqrt((r2-r1)^2 + (z2-z1)^2) * wind.bh.R_g
        delta_d_total += deltad
        d2 = (r2^2 + z2^2) * wind.bh.R_g^2
        xi0 = wind.radiation.xray_luminosity / (density * d2)
        xi = xi0
        for dummy in 1:2
            dtau = density * opacity_x(xi) * deltad * SIGMA_T
            tau_prov = tau + dtau
            xi = xi0 * exp(-tau_prov)
        end
        tau += density * SIGMA_T * deltad * opacity_x(xi)
        counter += 1
    end
    # add last bit
    density = wind.grids.density[r_arg, z_arg]
    deltad = sqrt((r-r2)^2 + (z-z2)^2) * wind.bh.R_g
    delta_d_total += deltad
    d2 = (r^2 + z^2) * wind.bh.R_g^2
    try
        @assert isapprox(delta_d_total, sqrt(d2), atol=0, rtol=5e-1)
    catch
        println("distances do not match!")
        println("r : $r z : $z")
        println("deltad_total: $(delta_d_total/wind.bh.R_g) d2 : $(sqrt(d2)/wind.bh.R_g)")
        @assert isapprox(delta_d_total, sqrt(d2), atol=0, rtol=1e-2)
    end
    xi0 = wind.radiation.xray_luminosity / (density * d2)
    xi = xi0
    for dummy in 1:2
        dtau = density * opacity_x(xi) * deltad * SIGMA_T
        tau_prov = tau + dtau
        xi = xi0 * exp(-tau_prov)
    end
    tau += density * SIGMA_T * deltad * opacity_x(xi)
    return tau
end

function compute_taux_leaf(point, intersection, leaf, wind::WindStruct)
    deltad = distance2d(point, intersection) * wind.bh.R_g
    d = distance2d([0.,wind.config["wind"]["z_0"]], point) * wind.bh.R_g
    density = leaf.data[1]
    xi0 = wind.radiation.xray_luminosity / (density * d^2)
    taux = 0.
    xi = xi0
    for i = 1:2
        taux = density * opacity_x(xi) * deltad * SIGMA_T
        xi = xi0 * exp(-taux)
    end
    #println("point :$point, den : $density")
    taux = density * opacity_x(xi) * deltad * SIGMA_T
    return taux
end

function compute_tau_x(r, z, wind::WindStruct ; return_coords = false)
    z = max(z, wind.config["wind"]["z_0"])
    @assert z >= 0
    #r = max(r,0.)
    point1 = [0.0, wind.config["wind"]["z_0"]]
    coords_list = [point1]
    point1leaf = findleaf(wind.quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(wind.quadtree, point2)
    if point1leaf == point2leaf
        push!(coords_list, point2)
        taux = compute_taux_leaf(point1, point2, point1leaf, wind)
        return taux
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    taux = compute_taux_leaf(point1, intersection, point1leaf, wind)
    currentpoint = intersection
    push!(coords_list, intersection)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        push!(coords_list, intersection)
        taux += compute_taux_leaf(currentpoint, intersection, currentleaf, wind)
        #println("r : $r,z : $z")
        currentpoint = intersection
        currentleaf = findleaf(wind.quadtree, currentpoint)
    end
    taux += compute_taux_leaf(currentpoint, point2, currentleaf, wind)
    push!(coords_list, point2)
    if return_coords
        return taux, coords_list
    else
        return taux
    end
end

function ionization_parameter(r, z, density, tau_x, wind::WindStruct)
    d2 = (r^2 + z^2) * wind.bh.R_g^2
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

function force_radiation(r, z, fm, wind::WindStruct ; include_tau_uv = false)
    if (z <= 1e-3 || wind.config["wind"]["gravity_only"])
        return [0.,0.]
    end
    int_values = integrate(r, z, wind, include_tau_uv=include_tau_uv)
    if wind.config["wind"]["nofm"]
        fm = 0.
    end
    force = wind.radiation.force_constant * (1. + fm) * int_values
    return force
end