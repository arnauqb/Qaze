export tau_uv_disk_blob,
       integrate_kernel, integrate_notau_kernel, integrate,
       integrate_fromstreamline_kernel, integrate_fromstreamline,
       integrate_gauss_kernel_r,
       integrate_gauss_kernel_z,
       integrate_gauss_kernel,
       integrate_gauss,
       integrate_split,
       integrate_gauss_kernel_normalized,
       integrate_notau_gauss_kernel,
       precomputable,
       integrable,
       integrate_parallel,
       compute_delta2,
       integrate_for,
       findcommonparent,
       integral_gauss_kernel_r_phid,
       integral_gauss_kernel_r_rd,
       a
       
       #X_phi, W_phi, X_r, W_r
using Cubature
using QuadGK
#x, w = gauss(10, 0, pi)
#const X_phi = x
#const W_phi = w
#x, w = gauss(10, 6, 1600)
#const X_r = x
#const W_r = w
a=0

"Computes the UV optical depth inside a leaf, taking into account the opacity
boost by the force multiplier."
function compute_tauuv_leaf(point, intersection, leaf)
    deltad = distance2d(point, intersection) #* wind.bh.R_g
    cellheight = cell_width(leaf)
    density = leaf.data[3] / cellheight
    #println(leaf.data[3]/leaf.data[2])
    tauuv = density * deltad #* SIGMA_T 
    return tauuv
end

function compute_delta2(r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    return max(delta, 0.) # sometimes it overflows...
end

function allparents(cell::Cell)
    Channel() do c
        queue = [cell]
        while !isempty(queue)
            current = pop!(queue)
            p = parent(current)
            if ! (p === nothing)
                put!(c, p)
                push!(queue, p)
            end
        end
    end
end

function findcommonparent(leaf1, leaf2)
    for parent1 in allparents(leaf1)
        for parent2 in allparents(leaf2)
            if parent1 == parent2
                return parent1
            end
        end
    end
end

"""
Compute the UV optical depth from a disc patch located at (r_d, phi_d),
until a gas element at (r,z). 
"""
function tau_uv_disk_blob(r_d, r, z, quadtree, maxtau)
    r_d > r ? backwards = true : backwards = false
    point1 = [r_d, 0.0]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(quadtree, point2)
    if point1leaf == point2leaf
        tauuv = compute_tauuv_leaf(point1, point2, point1leaf)
        return tauuv
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    tauuv = compute_tauuv_leaf(point1, intersection, point1leaf)
    currentpoint = intersection
    backwards && (currentpoint[1] -= 1e-8)
    currentleaf = findleaf(quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf)
        tauuv >= maxtau && return tauuv
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf)
    return tauuv
end

function integrate_kernel(v, r_d, phi_d, r, z, wind, rgsigma, maxtau)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, r, z, wind.quadtree, maxtau) * rgsigma
    tau_uv = tau_uv / sqrt((r-r_d)^2 + z^2) * sqrt(delta2)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    v[1] = (r - r_d * cos(phi_d)) * aux
    v[2] = aux
end

function integrate_notau_kernel(v, r_d, phi_d, r, z, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  1. / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    v[1] = aux * (r - r_d * cos(phi_d)) 
    v[2] = aux
end

function integrate(r, z, wind::WindStruct; include_tau_uv=true, maxevals = 3000)
    global a = 0
    rgsigma = wind.bh.R_g * SIGMA_T
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    maxtau = 20 / rgsigma
    if include_tau_uv
        (val, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind, rgsigma, maxtau),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0, #minimum(abs.(grav))/100,
                maxevals=maxevals,
                )
    else
        (val, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0,
                maxevals=maxevals
                )
    end
    val .*= [2 * z, 2 * z^2]
    return val
end

function integrate_parallel_kernel(v, r_d, phi_d, args)
    r, z, wind, rgsigma, maxtau = args
    r_d_arg = get_index.(Ref(wind.grids.disk_range), r_d)
    delta2 = compute_delta2.(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob.(r_d, r, z, Ref(wind.quadtree), maxtau) * rgsigma
    tau_uv = @. tau_uv / sqrt((r-r_d)^2 + z^2) * sqrt(delta2)
    nt = nt_rel_factors.(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  @. exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = @. nt / r_d^2 * f_uv * mdot 
    aux = r_part .* phi_part
    v[1,:] = @. (r - r_d * cos(phi_d)) * aux
    v[2,:] = aux
end

function integrate_notau_parallel_kernel(v, r_d, phi_d, args)
    r, z, wind = args
    r_d_arg = get_index.(Ref(wind.grids.disk_range), r_d)
    delta2 = compute_delta2.(r_d, phi_d, r, z)
    nt = nt_rel_factors.(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  @. 1 / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = @. nt / r_d^2 * f_uv * mdot 
    aux = r_part .* phi_part
    v[1,:] = @. (r - r_d * cos(phi_d)) * aux
    v[2,:] = aux
end

function integrate_parallel(r, z, wind::WindStruct; include_tau_uv=true, maxevals=3000)
    rgsigma = wind.bh.R_g * SIGMA_T
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    maxtau = 20 / rgsigma
    if include_tau_uv
        (val, err) = hcubature_v(2, 
                (x,v) ->integrate_parallel_kernel(v, x[1,:], x[2,:], (r, z, wind, rgsigma, maxtau)),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0,
                maxevals=maxevals
                )
    else
        (val, err) = hcubature_v(2, 
                (x,v) ->integrate_notau_parallel_kernel(v, x[1,:], x[2,:], (r, z, wind)),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0,
                maxevals=maxevals
                )
    end
    val .*= [2 * z, 2 * z^2]
    return val
end

#### experimental ######



function integrate_fromstreamline_notau_kernel(v, p, psi, r, z, wind, r_max)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if r_d <= 6. || r_d > 1200.
        v[1] = 0.
        v[2] = 0.
        return nothing
    end
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    common_part = nt * mdot * f_uv * p / (p^2 + z^2)^2 / r_d^3
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline_kernel(v, p, psi, r, z, wind, r_max, rgsigma, maxtau)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if r_d <= 6. || r_d > 1200.
        v[1] = 0.
        v[2] = 0.
        return nothing
    end
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    delta2 = p^2 + z^2
    tauuv = tau_uv_disk_blob(r_d, r, z, wind.quadtree, maxtau)  * rgsigma
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * sqrt(delta2)
    common_part = nt * mdot * f_uv * p / delta2^2 / r_d^3 * exp(-tauuv)
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline(r, z, wind::WindStruct; include_tau_uv=true, maxevals=3000)
    rgsigma = wind.bh.R_g * SIGMA_T
    maxtau = 20 /rgsigma
    r_max = wind.config["disk"]["outer_radius"]
    xmin = (0., 0.)
    #xmax = (r+r_max, pi)
    xmax = (100, pi)
    if include_tau_uv
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_kernel(v, x[1], x[2], r, z, wind, r_max, rgsigma, maxtau),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol= 0,#abs(maximum(grav)/100),
            maxevals=maxevals,
        )
    else
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_notau_kernel(v, x[1], x[2], r, z, wind, r_max),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol=0.0,
            maxevals=maxevals,
        )
    end
    return [2*z, 2*z^2] .* val
end



# gauss thinks 

function integrate_kernel(v, r_d, phi_d, r, z, wind, rgsigma, maxtau)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, r, z, wind.quadtree, maxtau) * rgsigma
    tau_uv = tau_uv / sqrt((r-r_d)^2 + z^2) * sqrt(delta2)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    v[1] = (r - r_d * cos(phi_d)) * aux
    v[2] = aux
end

function integral_gauss_kernel_r_phid(phi_d, r_d, r, z, wind)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    phid_part =  (r - r_d * cos(phi_d)) / delta2^2
    return phid_part
end

function integral_gauss_kernel_r_rd(r_d, r, z, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    phi_part, _ = quadgk(x -> integral_gauss_kernel_r_phid(x, r_d, r, z, wind), 0, pi, rtol=1e-3, atol=0)
    println(phi_part)
    return r_part * phi_part
end