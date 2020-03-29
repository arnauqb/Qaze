using RegionTrees
using StaticArrays: SVector
using RegionTrees
using VoronoiDelaunay
using Distances
using Printf
export initialize_leaf,
       get_cell_coordinates,
       find_neighbour,
       quadtree_initialize,
       quadtree_fill_point,
       quadtree_fill_timestep,
       quadtree_fill_line,
       quadtree_fill_all_lines,
       quadtree_erase_line,
       compute_tau_cell,
       quadtree_deltatau,
       quadtree_effective_density,
       test_interp,
       fill_last_point,
       normalize_point,
       denormalize_point,
       create_tessellation,
       initialize_tessellation,
       tess_density



function initialize_leaf(cell, child_indices)
    #data = [0, Array{Float64}(undef, 2, 0), Float64[]] # line_id 0 -> vaccum
    #data = [0, Float64[], Float64[], 0] # line_id 0 -> vaccum
    data = CellData(Int[], Float64[], Float64[], Float64[])
    return data
end

function get_cell_coordinates(cell)
    vcs = vertices(cell)
    xmin = vcs[1][1]
    ymin = vcs[1][2]
    xmax = vcs[4][1]
    ymax = vcs[4][2]
    return xmin, ymin, xmax, ymax
end

function find_neighbour(currentleaf, nextpoint, direction)
    # direction 0 = South, 1 = West, 2 = North, 3 = East
    prnt = parent(currentleaf)
    if direction == N
        if currentleaf == prnt[SE]
            return findleaf(prnt[NE], nextpoint)
        elseif currentleaf == prnt[SW]
            return findleaf(prnt[NW], nextpoint)
        end
        node = find_neighbour(prnt, direction) 
        if isleaf(node)
            return node
        end
    elseif direction == E
        if currentleaf == prnt[NW]
            return findleaf(prnt[NE], nextpoint)
        elseif currentleaf == prnt[SW]
            return findleaf(prnt[SE], nextpoint)
        end
        node = find_neighbour(prnt, direction) 
        if isleaf(node)
            return node
        end
    elseif direction == W
        if currentleaf == prnt[NE]
            return findleaf(prnt[NW], nextpoint)
        elseif currentleaf == prnt[SE]
            return findleaf(prnt[SW], nextpoint)
        end
        node = find_neighbour(prnt, direction) 
        if isleaf(node)
            return node
        end
    elseif direction == S
        if currentleaf == prnt[NE]
            return findleaf(prnt[SE], nextpoint)
        elseif currentleaf == prnt[NW]
            return findleaf(prnt[SW], nextpoint)
        end
        node = find_neighbour(prnt, direction) 
        if isleaf(node)
            return node
        end
    end
end

function quadtree_initialize(wind)
    height = wind.config["grids"]["r_max"]
    wind.quadtree = Cell(SVector(0., 0.), 
                    SVector(2*height, 2*height),
                    CellData(Int[], Float64[], Float64[], Float64[])
                    #[0, Array{Float64}(undef, 2, 0), Float64[], 0] # line_id position density 
                    )
end

function refine_leaf(point, leaf, targetsize, wind)
    while cell_width(leaf) > targetsize 
        split!(leaf, initialize_leaf)
        leaf = findleaf(wind.quadtree, point)
    end
    return leaf
end

function compute_accumulative_taus(z_positions, densities, z0)
    taus = zero(densities)
    taus[1] = densities[1] * (z_positions[1] - z0)
    try
        @assert(taus[1]>=0)
    catch
        println("zposi: $z_positions")
        println("denss: $densities")
        println("z0 : $z0")
        @assert(taus[1]>=0)
    end
    tautotal = taus[1]
    for i in 2:length(taus)
        tautotal += densities[i] * (z_positions[i] - z_positions[i-1])
        taus[i] = tautotal
        @assert(taus[i]>=0)
    end
    return taus
end

function quadtree_fill_point(point, zstep, density, line_id, line, wind; needs_refinement=true)
    leaf = findleaf(wind.quadtree, point)
    z0 = leaf.boundary.origin[2]
    #if line_id == 0 #erasing line
    #    leaf.data.line_id = Int[]
    #    leaf.data.z_positions = Float64[]
    #    leaf.data.densities = Float64[]
    #    return leaf
    #end
    if length(leaf.data.line_id) == 0 # empty cell
        leaf = refine_leaf(point, leaf, zstep, wind)
        push!(leaf.data.line_id, line_id)
        push!(leaf.data.z_positions, point[2])
        push!(leaf.data.densities, density)
        push!(leaf.data.taus, abs(point[2] - z0) * density)
        return leaf
    end

    push!(leaf.data.line_id, line_id)
    push!(leaf.data.z_positions, point[2])
    push!(leaf.data.densities, density)
    sorting_idx = sortperm(leaf.data.z_positions)
    leaf.data.line_id = leaf.data.line_id[sorting_idx]
    leaf.data.z_positions = leaf.data.z_positions[sorting_idx]
    leaf.data.taus = compute_accumulative_taus(leaf.data.z_positions, leaf.data.densities, z0) 
    return leaf
end

function quadtree_fill_timestep(point1, point2, density, linewidthnorm, line_id, line, wind::WindStruct)
    r1, z1 = point1
    r2, z2 = point2
    if point1 == point2 
        return nothing
    end
    if r2 > wind.config["grids"]["r_max"] || z2 > wind.config["grids"]["z_max"]
        return nothing
    end
    if z2 < z1
        point1[2], point2[2] = point2[2], point1[2]#z1, z2 = z2, z1
    end
    r1 > r2 ? backwards = true : backwards = false
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    currentpoint = copy(point1)
    height = abs(point2[2] - point1[2])#point1[2] / 100 
    zstep = min(max(height, wind.config["grids"]["minimum_cell_size"]), linewidthnorm * point1[1])
    while true
        rleft = currentpoint[1] * (1 - linewidthnorm/2)
        rright = currentpoint[1] * (1 + linewidthnorm/2)
        d2 = currentpoint[1]^2 + currentpoint[2]^2
        quadtree_fill_horizontal(rleft, rright, currentpoint[2], zstep, density, line_id, line, wind, true)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        backwards && (currentpoint[1] -= 1e-8)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentpoint[2] > point2[2] || previousleaf == currentleaf
            break
        end
        previousleaf = currentleaf
    end
end

function quadtree_fill_horizontal(rleft, rright, z, zstep, density, line_id, line, wind, needs_refinement=true)
    rl = [rleft, z]
    rr = [rright, z]
    currentpoint = rl
    while true
        leaf = quadtree_fill_point(currentpoint, zstep, density, line_id, line, wind, needs_refinement=needs_refinement)
        currentpoint = compute_cell_intersection(currentpoint, leaf, rl, rr)
        if currentpoint[1] > rright
            break
        end
    end
end

function quadtree_fill_line(line, line_id, wind)
    for i in 1:length(line.p.u_hist[:,1])-1
        point1 = [line.p.u_hist[i,1], line.p.u_hist[i,2]]
        point2 = [line.p.u_hist[i+1,1], line.p.u_hist[i+1,2]]
        lw_norm = line.p.line_width / line.p.r_0 
        density = line.p.n_hist[i]
        quadtree_fill_timestep(point1, point2, density, lw_norm, line_id, line, wind)
    end
end

"Fills and refines data from all streamlines"
function quadtree_fill_all_lines(wind::WindStruct)
    for i in 1:length(wind.lines)
        println(@sprintf("Line %02d of %02d", i, length(wind.lines)))
        isassigned(wind.lines, i) || continue
        wind.lines[i] === nothing && continue
        line = wind.lines[i]
        quadtree_fill_line(line, i, wind)
    end
end

function quadtree_erase_line(line_id, wind)
    line = wind.lines[line_id]
    #line.p.n_hist .= wind.grids.n_vacuum 
    quadtree_fill_line(line, 0, wind)
end

function quadtree_deltatau(point1, point2, wind)
    leaf = findleaf(wind.quadtree, point1)
    line = wind.lines[leaf.data.line_id]
    compute_tau_cell(line, point1, point2, leaf, wind)
end

function quadtree_effective_density(point1, point2, wind)
    leaf = findleaf(wind.quadtree, point1)
    tau = compute_tau_cell(point1, point2, leaf, wind)
    deltad = evaluate(Euclidean(), point1, point2)
    n = tau / deltad
    return n
end

function compute_tau_cell(point1, point2, leaf, wind)
    if point1[2] == point2[2]
        return 0
    end
    @assert point2[2] > point1[2]
    deltad = evaluate(Euclidean(), point1, point2)
    if length(leaf.data.line_id) == 0
        return wind.grids.n_vacuum * deltad
    end
    z_idx_1 = get_index(leaf.data.z_positions, point1[2])
    z_idx_2 = get_index(leaf.data.z_positions, point2[2])
    deltaz = leaf.data.z_positions[z_idx_2] - leaf.data.z_positions[z_idx_1]
    if deltaz == 0
        tau = leaf.data.densities[z_idx_1] * deltad
        return tau
    end
    deltatau = leaf.data.taus[z_idx_2] - leaf.data.taus[z_idx_1]
    @assert deltatau >= 0
    tau = deltatau / deltaz * deltad
    return tau
     
    #idx = get_index(leaf.data.z_positions, point[2])
    #line_id = leaf.data.line_id[idx]
    #if point[2] > maximum(wind.lines[line_id].p.u_hist[:,2]) #leaf.data.z_max
    #    return wind.grids.n_vacuum
    #end
    #n = leaf.data.densities[idx]
    #return n
end

function fill_last_point(line)
    currentpoint = line.p.u_hist[end,1:2]
    previouspoint = line.p.u_hist[end-1,1:2]
    lwnorm = line.p.line_width / line.p.r_0
    density = line.p.n_hist[end-1]
    line_id = line.p.line_id
    wind = line.p.wind
    quadtree_fill_timestep(currentpoint, previouspoint, density, lwnorm, line_id, line, wind)
end

function initialize_tessellation(wind)
    wind.tessellation = DelaunayTessellation2D{CustomPoint}();
    #cpoints = [
    #    CustomPoint(min_coord, min_coord),
    #    CustomPoint(min_coord, max_coord),
    #    CustomPoint(max_coord, min_coord),
    #    CustomPoint(max_coord, max_coord),
    #]
    #push!(wind.tessellation, cpoints)
end

function tess_density(r, z, wind)
    rnorm, znorm = normalize_point(r, z, wind)
    t1 = locate(wind.tessellation, CustomPoint(rnorm, znorm))
    dmin = Inf
    drmin = Inf
    density = 0.0
    lw = 0.0
    vx = nothing
    for vertex in [t1._a, t1._b, t1._c]
        vdr, vdz = denormalize_point(vertex._x, vertex._y, wind)
        distance = sqrt((vdr - r)^2 + (vdz - z)^2)
        if distance < dmin
            dmin = distance
            vx = vertex
        end
    end
    rline = vx._m * z + vx._n
    dr = r - rline
    if vx._type == "left"
        if dr < 0 
            return 1e2
        else
            return vx._density
        end
    elseif vx._type == "right"
        if dr > 0
            return 1e2
        else
            return vx._density
        end
    else
        return vx._density
    end

    #println("vertex: $((t1._a._x - 1) * wind.r_max), $((t1._a._y -1) * wind.z_max), $(t1._a._density), $(t1._a._linewidth), m: $(t1._a._m), n: $(t1._a._n)")
    #println("vertex: $((t1._b._x - 1) * wind.r_max), $((t1._b._y -1) * wind.z_max), $(t1._b._density), $(t1._b._linewidth), m: $(t1._b._m), n: $(t1._b._n)")
    #println("vertex: $((t1._c._x - 1) * wind.r_max), $((t1._c._y -1) * wind.z_max), $(t1._c._density), $(t1._c._linewidth), m: $(t1._c._m), n: $(t1._c._n)")
    #println("rvx : $(denormalize_point(vx._x, 1, wind))")
    #println("rline: $rline")
    #println(vx)
    #println(dr)
    #println(vx._linewidth * rline/2)
    #if dr > vx._linewidth * rline / 2
    #    return 1e2
    #else
    #    return vx._density
    #end
end

function normalize_point(r, z, wind)
    r = min(max(r / (2 * wind.r_max) + 1, min_coord), max_coord)
    z = min(max(z / (2 * wind.z_max) + 1, min_coord), max_coord)
    return [r, z]
end

function denormalize_point(rnorm, znorm, wind)
    if !isnan(rnorm) && !isnan(znorm)
        @assert rnorm <= 2
        @assert rnorm >= 1
        @assert znorm <= 2
        @assert znorm >= 1
    end
    r = (rnorm - 1) * 2* wind.r_max
    z = (znorm - 1) * 2 * wind.z_max
    return [r, z]
end

function create_tessellation(wind)
    EPS = 2 * eps(Float64)
    points = Array{CustomPoint}(undef, 0)
    DZMAX = 3
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        lwnorm = line.p.line_width / line.p.r_0 
        for j in 3:length(line.p.u_hist[:,1])
            density = line.p.n_hist[j-1]
            r1, z1 = line.p.u_hist[j-1, 1:2]
            r2, z2 = line.p.u_hist[j, 1:2]
            if z1 == z2
                continue
            end
            #m = (r2 - r1) / (z2 - z1)
            #n = r2 - m * z2
            #rnorm, znorm = normalize_point(r1, z1, wind)
            #p = CustomPoint(rnorm, znorm, density, lwnorm, m, n)
            #push!(points, p)
            if (z2-z1) > DZMAX
                zp = z1
                m = (r2 - r1) / (z2 - z1)
                n = r2 - m * z2
                while true 
                    zp += DZMAX / 5
                    if zp >= z2
                        break
                    end
                    rp = m * zp + n
                    rnorm, znorm = normalize_point(rp, zp, wind)
                    p = CustomPoint(rnorm, znorm, density, lwnorm, "center", m, n)
                    push!(points, p)
                    # add borders
                    rm1 = rp - lwnorm * rp
                    rm2 = rp + lwnorm * rp
                    rm1norm, znorm = normalize_point(rm1, zp, wind)
                    rm2norm, znorm = normalize_point(rm2, zp, wind)
                    pmirror1 = CustomPoint(rm1norm, znorm, density, lwnorm, "left", m, n)
                    pmirror2 = CustomPoint(rm2norm, znorm, density, lwnorm, "right", m, n)
                    push!(points, pmirror1)
                    push!(points, pmirror2)
                end
            else
                m = (r2 - r1) / (z2 - z1)
                n = r2 - m * z2
                rnorm, znorm = normalize_point(r1, z1, wind)
                p = CustomPoint(rnorm, znorm, density, lwnorm, "center", m, n)
                push!(points, p)
                # add borders
                rm1 = r1 - lwnorm * r1
                rm2 = r1 + lwnorm * r1
                rm1norm, znorm = normalize_point(rm1, z1, wind)
                rm2norm, znorm = normalize_point(rm2, z1, wind)
                pmirror1 = CustomPoint(rm1norm, znorm, density, lwnorm, "left", m, n)
                pmirror2 = CustomPoint(rm2norm, znorm, density, lwnorm, "right", m, n)
                push!(points, pmirror1)
                push!(points, pmirror2)
            end
        end
    end
    points = unique(points)
    #for pt in points
        #println("inserting ... $pt")
        #flush(stdout)
        #push!(wind.tessellation, pt)
    #end
    push!(wind.tessellation, points)
end