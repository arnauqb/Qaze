using RegionTrees
using StaticArrays: SVector
using RegionTrees
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
       interpolate_density,
       quadtree_density,
       test_interp,
       fill_last_point


function initialize_leaf(cell, child_indices)
    #data = [0, Array{Float64}(undef, 2, 0), Float64[]] # line_id 0 -> vaccum
    #data = [0, Float64[], Float64[], 0] # line_id 0 -> vaccum
    data = CellData(0, Float64[], Float64[], 0, 0)
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
                    CellData(0, Float64[], Float64[], 0, 0)
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

function quadtree_fill_point(point, zstep, density, line_id, line, wind; needs_refinement=true)
    leaf = findleaf(wind.quadtree, point)
    if line_id == 0 #erasing line
        leaf.data.line_id = 0
        leaf.data.z_positions = Float64[]
        leaf.data.densities = Float64[]
        leaf.data.z_max = 0.0
        leaf.data.direction = 0
        return leaf
    end
    if leaf.data.line_id == 0 # empty cell
        leaf = refine_leaf(point, leaf, zstep, wind)
        leaf.data.line_id = line_id
        push!(leaf.data.z_positions, point[2])
        push!(leaf.data.densities, density)
        leaf.data.z_max = maximum(line.p.u_hist[:,2]) #point[2]
        leaf.data.direction = 0
        return leaf
    end
    den_cell = leaf.data.densities[end]
    if (density/den_cell) > 0.9 && (density / den_cell) < 1.1 # discard value if its within 10% of the others
        return leaf
    end
    if leaf.data.direction == 0
        if point[2] >= leaf.data.z_positions[1]
            leaf.data.direction = 1
        else
            leaf.data.direction = -1
        end
        push!(leaf.data.z_positions, point[2])
        push!(leaf.data.densities, density)
        leaf.data.z_max = max(point[2], leaf.data.z_max)
        return leaf
    end
    if ((leaf.data.direction==1) && (point[2]<leaf.data.z_positions[end]))
        return leaf
    end
    if ((leaf.data.direction==-1) && (point[2]>leaf.data.z_positions[end]))
        return leaf
    end
    push!(leaf.data.z_positions, point[2])
    push!(leaf.data.densities, density)
    leaf.data.z_max = max(point[2], leaf.data.z_max)
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
    d02 = point1[1]^2 + point1[2]^2
    zstep = min(max(height, wind.config["grids"]["minimum_cell_size"]), linewidthnorm * point1[1])
    interp_density(d2) = density #/ d2 * d02 
    while true
        rleft = max(currentpoint[1] * (1 - linewidthnorm/2), 6)
        rright = max(currentpoint[1] * (1 + linewidthnorm/2), 6)
        d2 = currentpoint[1]^2 + currentpoint[2]^2
        quadtree_fill_horizontal(rleft, rright, currentpoint[2], zstep, interp_density(d2), line_id, line, wind, true)
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

function quadtree_density(point, wind)
    leaf = findleaf(wind.quadtree, point)
    if leaf.data.line_id == 0
        return wind.grids.n_vacuum
    end
    line = wind.lines[leaf.data.line_id]
    interpolate_density(line, point, leaf, wind)
end

function interpolate_density(line, point, leaf, wind)
    if point[2] > leaf.data.z_max
        return wind.grids.n_vacuum
    end
    idx = get_index(leaf.data.z_positions, point[2])
    n = leaf.data.densities[idx]
    return n
    #X = leaf.data[2] #line.p.u_hist[:,1:2]
    #densities = leaf.data[3]
    #z_max = maximum(X[2,:])
    ##if point[2] > z_max
    ##    return wind.grids.n_vacuum
    ##end
    #point = reshape(point, 2, 1)
    #distances = pairwise(Euclidean(), point, X, dims=2)
    #closest = argmin(distances)[2]
    #distance = distances[closest]
    #if (abs(point[1] - X[1, closest])) > (X[1, closest] * line.p.line_width / line.p.r_0 / 2)
    #    return wind.grids.n_vacuum
    #end
    #n = densities[closest]
    #return n
end

function interpolate_density(line_id::Int, point, leaf, wind)
    if line_id == 0
        return wind.config["grids"]["n_vacuum"]
    else
        return interpolate_density(wind.lines[line_id], point, leaf, wind)
    end
end

function test_interp(points, line)
    points = reshape(points, 2, length(points))
    n_interp = interpolate_density.(Ref(line), points)
    return n_interp
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