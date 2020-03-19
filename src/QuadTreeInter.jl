using RegionTrees
using StaticArrays: SVector
using RegionTrees
using Distances
export initialize_leaf,
       quadtree_initialize,
       quadtree_fill_point,
       quadtree_fill_timestep,
       quadtree_fill_line,
       interpolate_density,
       test_interp

function initialize_leaf(cell, child_indices)
    data = [0] # line_id 0 -> vaccum
end

function quadtree_initialize(wind)
    height = wind.config["grids"]["r_max"]
    wind.quadtree = Cell(SVector(0., 0.), 
                    SVector(2000.0, 2000.0),
                    [0] # density, fm , cell optical thickness, lines that pass through this cell
                    )
end

function quadtree_fill_point(point, zstep, line_id, wind)
    leaf = findleaf(wind.quadtree, point)
    while cell_width(leaf) > zstep
        split!(leaf, initialize_leaf)
        leaf = findleaf(wind.quadtree, point)
    end
    leaf.data[1] = line_id
    return leaf
end

function quadtree_fill_timestep(point1, point2, linewidthnorm, line_id, wind::WindStruct)
    r1, z1 = point1
    r2, z2 = point2
    if point1 == point2 
        return nothing
    end
    if z2 < z1
        z1, z2 = z2, z1
    end
    zstep = wind.config["grids"]["minimum_cell_size"]
    height = z2 - z1
    m = (r2 - r1) / (z2 -z1)
    n = r2 - m * z2
    r_center(z) = m * z + n
    r_left(z) = max(r_center(z) * (1 - linewidthnorm/2.), 0)
    r_right(z) = max(r_center(z) * (1 + linewidthnorm/2.), 0)
    if height < zstep
        rl = [r_left(z1), z1]
        rr = [r_right(z1), z1]
        leaf = quadtree_fill_point(rl, zstep, line_id, wind)
        rrleaf = findleaf(wind.quadtree, rr)
        point = copy(rl)
        while true
            point = compute_cell_intersection(point, leaf, rl, rr)
            leaf = findleaf(wind.quadtree, point)
            leaf = quadtree_fill_point(point, zstep, line_id, wind)
            rrleaf = findleaf(wind.quadtree, rr)
            if leaf == rrleaf
                break
            end
        end
    else
        currentpoint = [r1, z1]
        z = z1
        while z < z2
            leftpoint = [r_left(z), z]
            rightpoint = [r_right(z), z]
            quadtree_fill_point(leftpoint, zstep, line_id, wind)
            quadtree_fill_point(rightpoint, zstep, line_id, wind)
            z += zstep
        end
        # fill inbetween
        z = z1
        while z < z2
            rl = [r_left(z), z]
            rr = [r_right(z), z]
            leaf = findleaf(wind.quadtree, rl)
            leaf.data[1] = line_id
            rrleaf = findleaf(wind.quadtree, rr)
            point = copy(rl)
            previousleaf = copy(leaf)
            while leaf != rrleaf 
                newpoint = compute_cell_intersection(point, leaf, rl, rr)
                previousleaf = copy(leaf)
                leaf = findleaf(wind.quadtree, newpoint)
                leaf.data[1] = line_id
            end
            z += zstep / 10
        end
    end

end

function quadtree_fill_line(line, line_id, wind)
    for i in 1:length(line.p.u_hist[:,1])-1
        point1 = [line.p.u_hist[i,1], line.p.u_hist[i,2]]
        point2 = [line.p.u_hist[i+1,1], line.p.u_hist[i+1,2]]
        lw_norm = line.p.line_width / point1[1]
        quadtree_fill_timestep(point1, point2, lw_norm, line_id, wind)
    end
end

function interpolate_density(line, point, wind)
    X = line.p.u_hist[:,1:2]
    z_max = maximum(X[:,2])
    if point[2] > z_max
        return wind.config["grids"]["n_vacuum"]
    end
    point = reshape(point, 2, 1)
    distances = pairwise(Euclidean(), point, X', dims=2)
    closest = argmin(distances)[2]
    n = line.p.n_hist[closest]
    return n
end

function interpolate_density(line_id::Int, point, wind)
    if line_id == 0
        return wind.config["grids"]["n_vacuum"]
    else
        return interpolate_density(wind.lines[line_id], point, wind)
    end
end

function test_interp(points, line)
    points = reshape(points, 2, length(points))
    n_interp = interpolate_density.(Ref(line), points)
    return n_interp
end