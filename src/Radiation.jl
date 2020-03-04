export thermal_velocity,
       eddington_luminosity,
       initialize_uv_fraction!,
       update_mdot_grid!,
       nt_rel_factors,
       opacity_x,
       compute_taux_leaf,
       compute_tau_x,
       ionization_parameter,
       compute_tau_eff,
       force_multiplier,
       force_radiation

function thermal_velocity(T, mu = 1.)
    v = sqrt(K_B * T / (mu * M_P)) / C
    return v
end

"""
Computes the radial fraction of UV flux. 
If the flag disk_uv_fraction is on, then the fraction is computed using the QSOSED model
in Kubota and Done 2018.
Otherwise, the fraction is constant throughout radius equal to f_uv
"""
function initialize_uv_fraction!(wind::WindStruct)
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

"""
Updates the radial dependent mass accretion rate to take into account the wind mass loss
self consistently. Initially it's initialized at a constant value of Mdot.
"""
function update_mdot_grid!(wind::WindStruct)
    accumulated_wind = 0.
    for line in reverse(wind.lines)
        r_0 = line.p.r_0
        width = line.p.line_width
        rmin_arg = get_index(wind.grids.disk_range, r_0 - width/2.)
        rmax_arg = get_index(wind.grids.disk_range, r_0 + width/2.)
        mw = compute_line_mdot(line, wind)
        mw_norm = mw / mass_accretion_rate(wind.bh)
        accumulated_wind += mw_norm
        mdot = wind.bh.mdot - accumulated_wind
        wind.grids.mdot[rmin_arg:rmax_arg] .= max(mdot , 0.)
    end
    rmin = wind.lines_initial_radius 
    rmin_arg = get_index(wind.grids.disk_range, rmin )
    wind.grids.mdot[1:rmin_arg] .= max(wind.bh.mdot - accumulated_wind, 0.)
    return wind.grids.mdot
end

"""
Novikov-Thorne relativistic factors for the AD spectrum.
If the passed radius is smaller than ISCO, it returns 0.
"""
function nt_rel_factors(r, astar, isco)
    if r <= isco
        return 0.
    end
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

"Simple X-Ray opacity as a function of ionization parameter xi"
function opacity_x(xi)
    if xi <= 1e5
        return 100
    else
        return 1
    end
end

"""
Computes the variation of optical depth inside a tree leaf.
We iterate multiple times to compute a consistent ionization parameter
and X-Ray optical depth
"""
function compute_taux_leaf(point, intersection, leaf, wind::WindStruct)
    deltad = distance2d(point, intersection) * wind.bh.R_g # cell size
    d = distance2d([0.,wind.z_0], point) * wind.bh.R_g # distance from the center
    density = leaf.data[1]
    xi0 = wind.radiation.xray_luminosity / (density * d^2)
    taux = 0.
    xi = xi0
    for i = 1:2 
        taux = density * opacity_x(xi) * deltad * SIGMA_T
        xi = xi0 * exp(-taux)
    end
    taux = density * opacity_x(xi) * deltad * SIGMA_T
    return taux
end

"""
Computes the X-Ray optical depth at a point (r,z).
Starting from the center, we compute the intersection to the next cell,
following the lightray direction (0,0) -> (r,z). At each cell,
we consistently compute tau_x and the ionization parameter.
"""
function compute_tau_x(r, z, wind::WindStruct ; return_coords = false)
    z = max(z, wind.z_0)
    @assert z >= 0
    point1 = [0.0, wind.z_0]
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

"Computes the ionization parameter at a point (r,z)"
function ionization_parameter(r, z, density, tau_x, wind::WindStruct)
    d2 = (r^2 + z^2) * wind.bh.R_g^2
    xi = wind.radiation.xray_luminosity * exp(-tau_x) / (density * d2) + 1e-20
    @assert xi >= 0
    return xi
end

"This is the sobolev optical depth parameter for the force multiplier"
function compute_tau_eff(density, dv_dr, v_th)
    if dv_dr == 0
        return 1.
    end
    @assert density >= 0
    t = density * SIGMA_T * abs(v_th / dv_dr)
    return t
end

"Fitting parameter for the fm"
function force_multiplier_k(xi)
    k = 0.03 + 0.385 * exp(-1.4 * xi^0.6)
    return k
end

"Fitting parameter for the fm"
function force_multiplier_eta(xi)
    if (log10(xi) < 0.5)
        aux = 6.9 * exp(0.16 * xi^0.4)
        eta_max = 10^aux
    else
        aux = 9.1 * exp(-7.96e-3 * xi)
        eta_max = 10^aux
    end
end

"""
Computes the analytical approximation for the force multiplier,
from Stevens and Kallman 1990. Note that we modify it slightly to avoid
numerical overflow.
"""
function force_multiplier(t, xi)
    @assert t>= 0
    @assert xi>= 0
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

"""
Computes the radiation force at a point (r,z) opacity (1+fm)sigma_t.
the option include_tau_uv decides wheter to consider the attenuation
of the UV flux along the disc annuli-gas patch los.
"""
function force_radiation(r, z, fm, wind::WindStruct ; include_tau_uv = false)
    @assert r >= 0
    @assert z >= 0
    if (wind.config["wind"]["gravity_only"])
        return [0.,0]
    end
    if (z < wind.config["radiation"]["constant_frad_height"])
        return [0.0, force_radiation(r, wind.config["radiation"]["constant_frad_height"], fm, wind, include_tau_uv = include_tau_uv)[2]]
    end
    int_values = integrate(r, z, wind, include_tau_uv=include_tau_uv)
    if wind.config["wind"]["nofm"]
        fm = 0.
    end
    force = wind.radiation.force_constant * (1. + fm) * int_values
    return force
end
