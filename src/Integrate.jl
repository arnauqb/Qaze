export tau_uv_disk_blob,
       tau_uv_disk_blob_fromstreamline,
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
a = 0

"Computes the UV optical depth inside a leaf, taking into account the opacity
boost by the force multiplier."
function compute_tauuv_leaf(point, intersection, leaf)
    deltad = distance2d(point, intersection) #* wind.bh.R_g
    density = leaf.data[1]
    tauuv = density * deltad #* SIGMA_T 
    return tauuv
end

function compute_delta2(r_d, phi_d, r, z)
    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
    return delta
end

"""
Compute the UV optical depth from a disc patch located at (r_d, phi_d),
until a gas element at (r,z). 
"""
function tau_uv_disk_blob(r_d, phi_d, r, z, delta, quadtree, rgsigma)
    r_d > r ? backwards = true : backwards = false
    point1 = [r_d, 0.0]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(quadtree, point2)
    maxtau = 14 / rgsigma
    if point1leaf == point2leaf
        tauuv = compute_tauuv_leaf(point1, point2, point1leaf)
        tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
        return tauuv * rgsigma
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    tauuv = compute_tauuv_leaf(point1, intersection, point1leaf)
    currentpoint = intersection
    backwards && (currentpoint[1] -= 1e-8)
    currentleaf = findleaf(quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf)
        if tauuv > maxtau
            return 14.0
        end
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf)
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta * rgsigma
    return tauuv
end

function integrate_kernel(v, r_d, phi_d, r, z, wind, rgsigma)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
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
    rgsigma = wind.bh.R_g * SIGMA_T
    xmin = (wind.config["disk"]["inner_radius"], 0.)
    xmax = (wind.config["disk"]["outer_radius"], pi)
    if include_tau_uv
        (val, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind, rgsigma),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=1e-15,
                maxevals=maxevals
                )
    else
        (val, err) = hcubature(2, 
                (x,v) ->integrate_notau_kernel(v, x[1], x[2], r, z, wind),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=1e-15,
                maxevals=maxevals
                )
    end
    val .*= [2 * z, 2 * z^2]
    return val
end

function integrate_parallel_kernel(v, r_d, phi_d, args)
    r, z, wind, rgsigma = args
    r_d_arg = get_index.(Ref(wind.grids.disk_range), r_d)
    delta2 = compute_delta2.(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob.(r_d, phi_d, r, z, sqrt.(delta2), Ref(wind.quadtree), rgsigma)
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
    if include_tau_uv
        (val, err) = hcubature_v(2, 
                (x,v) ->integrate_parallel_kernel(v, x[1,:], x[2,:], (r, z, wind, rgsigma)),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=1e-15,
                maxevals=maxevals
                )
    else
        (val, err) = hcubature_v(2, 
                (x,v) ->integrate_notau_parallel_kernel(v, x[1,:], x[2,:], (r, z, wind)),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=1e-15,
                maxevals=maxevals
                )
    end
    val .*= [2 * z, 2 * z^2]
    return val
end

#### experimental ######


function tau_uv_disk_blob_fromstreamline(p, psi, r, r_d, z, delta, quadtree, rgsigma)
    r_d > r ? backwards=true : backwards=false
    point1 = [r_d, 0]
    point1leaf = findleaf(quadtree, point1)
    point2 = [r,z]
    point2leaf = findleaf(quadtree, point2)
    maxtau = 14.0# * rgsigma
    if point1leaf == point2leaf
        tauuv = compute_tauuv_leaf(point1, point2, point1leaf)
        tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta
        return tauuv * rgsigma
    end
    intersection = compute_cell_intersection(point1, point1leaf, point1, point2)
    tauuv = compute_tauuv_leaf(point1, intersection, point1leaf)
    currentpoint = intersection
    backwards && (currentpoint[1] -= 1e-8)
    currentleaf = findleaf(quadtree, currentpoint)
    while currentleaf != point2leaf
        intersection = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        tauuv += compute_tauuv_leaf(currentpoint, intersection, currentleaf)
        if tauuv * rgsigma > maxtau
            return 14.0
        end
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf)
    tauuv = tauuv / sqrt((r-r_d)^2 + z^2) * delta * rgsigma 
    return tauuv
end

function integrate_fromstreamline_notau_kernel(v, p, psi, r, z, wind, r_max)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if r_d <= 6.
        v[1] = 0.
        v[2] = 0.
    end
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    common_part = nt * mdot * f_uv * p / (p^2 + z^2)^2 / r_d^3
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline_kernel(v, p, psi, r, z, wind, r_max, rgsigma)
    cosψ = cos(psi)
    r_d = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if r_d <= 6.
        v[1] = 0.
        v[2] = 0.
    end
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    delta = sqrt(p^2 + z^2)
    tauuv = tau_uv_disk_blob_fromstreamline(p, psi, r, r_d, z, delta, wind.quadtree, rgsigma)
    common_part = nt * mdot * f_uv * p / delta^4 / r_d^3 * exp(-tauuv)
    v[1] = common_part * p * cosψ
    v[2] = common_part
end

function integrate_fromstreamline(r, z, wind::WindStruct; include_tau_uv=true, maxevals=3000)
    rgsigma = wind.bh.R_g * SIGMA_T
    r_max = wind.config["disk"]["outer_radius"]
    xmin = (0., 0.)
    #xmax = (r + r_max, pi)
    if r < 100
        xmax = (r + r_max, pi)
    else
        xmax = (100, pi)
    end

    if include_tau_uv
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_kernel(v, x[1], x[2], r, z, wind, r_max, rgsigma),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol=1e-10,
            maxevals=maxevals,
        )
    else
        (val, err) = hcubature(2, 
            (x,v) ->integrate_fromstreamline_notau_kernel(v, x[1], x[2], r, z, wind, r_max),
            xmin,
            xmax,
            reltol = wind.config["radiation"]["integral_rtol"],
            abstol=1e-10,
            maxevals=maxevals,
        )
    end
    return [2*z, 2*z^2] .* val
end

function integrate_split(r, z, wind::WindStruct; include_tau_uv=true, eps=0)
    rgsigma = wind.bh.R_g * SIGMA_T
    if include_tau_uv
        xmin = (wind.config["disk"]["inner_radius"], 0.)
        xmax = (r-eps, pi)
        (val1, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind, rgsigma),
                xmin,
                xmax,
                reltol = wind.config["radiation"]["integral_rtol"],
                abstol=0.,
                )
        xmin = (r+eps, 0.)
        xmax = (wind.config["disk"]["outer_radius"], pi)
        (val2, err) = hcubature(2, 
                (x,v) ->integrate_kernel(v, x[1], x[2], r, z, wind, rgsigma),
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

# 1D integrals

function integrate_gauss_kernel_r(r_d, phi_d, r, z, wind, rgsigma)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    return (r - r_d * cos(phi_d)) * aux
end

function integrate_gauss_kernel_z(r_d, phi_d, r, z, wind, rgsigma)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    return aux
end

function integrate_gauss_kernel(r_d, phi_d, r, z, wind, rgsigma)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    return [(r - r_d * cos(phi_d)) * aux, aux]
end

function integrate_gauss_kernel_normalized(r_d, phi_d, r, z, wind, rgsigma)
    r_max = 1000
    r_min = 6
    r_d = r_d * (r_max- r_min)/2 + (r_max + r_min)/2
    dr = (r_max - r_min) / 2
    dphi = pi / 2.
    phi_d = pi/2 * phi_d + pi/2
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  exp(-tau_uv) / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part * dr * dphi
    return [(r - r_d * cos(phi_d)) * aux, aux]
end

function integrate_gauss(r, z, wind, X_r, X_phi, W_r, W_phi)
    rgsigma = wind.bh.R_g * SIGMA_T
    v = sum(integrate_gauss_kernel.(X_r, X_phi', r, z, Ref(wind)) .* (W_r * W_phi'), rgsigma)
    return 2 .* [z * v[1], z^2 * v[2]]
end

function integrate_notau_gauss_kernel(r_d, phi_d, r, z, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    delta2 = compute_delta2(r_d, phi_d, r, z)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    phi_part =  1 / delta2^2
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    aux = r_part * phi_part
    return [(r - r_d * cos(phi_d)) * aux, aux]
end

function precomputable(r_d, phi_d, wind)
    r_d_arg = get_index(wind.grids.disk_range, r_d)
    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
    f_uv = wind.grids.uv_fractions[r_d_arg]
    mdot = wind.grids.mdot[r_d_arg]
    r_part = nt / r_d^2 * f_uv * mdot 
    return r_part
end

function integrable(r_d, phi_d, r, z, wind)
    rgsigma = wind.bh.R_g * SIGMA_T

    delta2 = compute_delta2(r_d, phi_d, r, z)
    #tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, sqrt(delta2), wind.quadtree, rgsigma)
    abs = absorps.(r_d, r, z, sqrt.(delta2), Ref(wind), rgsigma)
    common = @. abs / delta2^2
    #common = exp(-tau_uv) / delta2^2
    #common =  1/ delta^2
    return common .* [r - r_d * cos(phi_d), 1.]
end

function absorps(r_d, r, z, delta, (quadtree, rgsigma))
    #factor = @. delta / sqrt((r-r_d)^2 + z^2)
    #taumax = 14 / rgsigma 
    #taus = @. factor * tau_uv_disk_blob(r_d, r, z, (quadtree, taumax))  * rgsigma
    #return exp.(-taus)
end

function integrate_for(x_r, x_phi, w_r, w_phi, precomputes, r, z, wind)
    result = [0.0, 0.0]
    sigmarg = wind.bh.R_g * SIGMA_T
    maxtau = 14 / sigmarg
    for (i, rd) in enumerate(x_r)
        wrd = w_r[i]
        linelength = sqrt((r-rd)^2 + z^2)
        tauuv = tau_uv_disk_blob(rd, r, z, (wind.quadtree, maxtau)) * sigmarg
        for (j, phi) in enumerate(x_phi)
            wphi = w_phi[j]
            delta2 = compute_delta2(rd, phi, r, z)
            factor = sqrt(delta2) / linelength
            abs = exp(-tauuv * factor)
            common = abs / delta2^2 * wrd * wphi * precomputes[i,j]
            result[1] += common * (r - rd * cos(phi))
            result[2] += common
        end
    end
    return 2 .* [z, z^2] .* result
end


function tau_uv_disk_blob(r_d, r, z, (quadtree, maxtau))
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
        if tauuv > maxtau
            return 14.0
        end
        currentpoint = intersection
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(quadtree, currentpoint)
    end
    tauuv += compute_tauuv_leaf(currentpoint, point2, currentleaf)
    return tauuv
end

#function precomputable(r_d, phi_d, wind)
#    r_d = r_d * (1600 - 6)/2 + (1600 + 6)/2
#    phi_d = pi/2 * phi_d + pi/2
#    r_d_arg = get_index(wind.grids.disk_range, r_d)
#    nt = nt_rel_factors(r_d, wind.bh.spin, wind.bh.isco)
#    f_uv = wind.grids.uv_fractions[r_d_arg]
#    mdot = wind.grids.mdot[r_d_arg]
#    r_part = nt / r_d^2 * f_uv * mdot 
#    return r_part
#end
#
#function integrable(r_d, phi_d, r, z, wind)
#    r_d = r_d * (1600 - 6)/2 + (1600 + 6)/2
#    phi_d = pi/2 * phi_d + pi/2
#    tau_uv = tau_uv_disk_blob(r_d, phi_d, r, z, wind)
#    delta = r^2 + z^2 + r_d^2 - 2. * r * r_d * cos(phi_d)
#    common = exp(-tau_uv) / delta^2
#    #common =  1/ delta^2
#    return common .* [r - r_d * cos(phi_d), 1.] .* (1600 - 6) / 2. .* pi
#end
