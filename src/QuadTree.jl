using RegionTrees
using StaticArrays: SVector
using RegionTrees
import Base: copy, ==
export copy,
       compute_cell_intersection,
       compute_cell_size,
       refine_leaf!,
       fill_and_refine_linewidth!,
       fill_and_refine!,
       fill_and_refine_line!,
       fill_and_refine_all_lines!,
       erase_line_from_tree!

"Computes the diagonal of the cell"
function compute_cell_size(cell::Cell)
    r1, z1 = vertices(cell.boundary)[1,1]
    r2, z2 = vertices(cell.boundary)[2,2]
    distance = sqrt((r2-r1)^2 + (z2-z1)^2)
    return distance
end

"Method of the copy function for Cell type"
function copy(cell::Cell)
    newcell = Cell(cell.boundary, cell.data)
    return newcell
end

"Method of the == function for the Cell type"
function ==(cell1::Cell, cell2::Cell)
    return cell1.boundary == cell2.boundary
end

"Reinitializes quadtree to one leaf and default values"
function reinitialize_tree!(wind::WindStruct)
    quadtree_max_radius = wind.config["grids"]["r_max"]
    quadtree_max_height = wind.config["grids"]["z_max"]
    n_vacuum = wind.grids.n_vacuum
    delta_tau_0 = n_vacuum * sqrt(quadtree_max_height^2 + quadtree_max_radius^2) * wind.bh.R_g * SIGMA_T
    wind.quadtree = Cell(SVector(0., 0.), SVector(2 * quadtree_max_radius, 2* quadtree_max_height), [n_vacuum, 0., delta_tau_0])
end

"Following the line initialpoint -> finalpoint, computes the intersection after the current point
to the next tree leaf"
function compute_cell_intersection(currentpoint, cell::Cell, initialpoint, finalpoint)
    r, z = currentpoint
    ri, zi = initialpoint
    rf, zf = finalpoint
    direction_r = sign(rf - ri)
    direction_z = sign(zf - zi)
    cell_vertices = vertices(cell)
    direction_r == 1 ? cell_r_boundary=cell_vertices[2,2][1] : cell_r_boundary=cell_vertices[1,1][1]
    direction_z == 1 ? cell_z_boundary=cell_vertices[2,2][2] : cell_z_boundary=cell_vertices[1,1][2]
    ri == rf ? lambda_r=Inf : lambda_r=(cell_r_boundary - r) / (rf - ri)
    zi == zf ? lambda_z=Inf : lambda_z=(cell_z_boundary - z) / (zf - zi)
    lambda = min(lambda_r, lambda_z)
    try
        @assert lambda >= 0.0
    catch
        println("current point : $currentpoint")
        println("cell: $cell")
        println("initial point: $initialpoint")
        println("final point: $finalpoint")
        throw(DomainError)
    end
    intersection = [r, z] + lambda .* [rf - ri, zf - zi]
    return intersection
end

"Initializes cell data after split"
function split_cell_initialization(cell, child_indices)
    newdata = [1e2, 0., 0.] 
    return newdata
end

"splits leaf until the max cell optical thickness is reached"
function refine_leaf!(currentpoint, leaf, n, fm, wind)
    optically_thin = false
    while !optically_thin
        optically_thin = true
        maxz = vertices(leaf.boundary)[2,2][2]
        refine_condition = (leaf.data[3] > wind.config["grids"]["cell_optical_thickness"]) && (maxz > wind.z_0) 
        if refine_condition
            optically_thin = false
            split!(leaf, split_cell_initialization)
            newleaf = findleaf(leaf, currentpoint)
            deltad = compute_cell_size(newleaf)
            deltatau = n * SIGMA_T * wind.bh.R_g * deltad
            newleaf.data = [n, fm, deltatau]
            leaf = newleaf
        end
    end
end

"Fills the line data to all leaves along the width of the streamline at point, and refines them if necessary"
function fill_and_refine_linewidth!(point, linewidth, density, fm, wind::WindStruct)
    rmin = 0.0
    rmax = 2 * wind.grids.r_max 
    point1 = [max(point[1] - linewidth / 2., rmin) , point[2]] 
    point2 = [min(point[1] + linewidth / 2., rmax) , point[2]] 
    sigma_rg = SIGMA_T * wind.bh.R_g
    currentpoint = copy(point1)
    currentleaf = findleaf(wind.quadtree, point1)
    point2leaf = findleaf(wind.quadtree, point2)
    n_fill = max(density, currentleaf.data[1])
    fm_fill = max(fm, currentleaf.data[2])
    deltad = compute_cell_size(currentleaf)
    deltatau = sigma_rg * deltad * n_fill
    currentleaf.data = [n_fill, fm_fill, deltatau]
    refine_leaf!(currentpoint, currentleaf, n_fill, fm_fill, wind)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    point2leaf = findleaf(wind.quadtree, point2)
    while (currentleaf != point2leaf)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        n_fill = max(density, currentleaf.data[1])
        fm_fill = max(fm, currentleaf.data[2])
        deltad = compute_cell_size(currentleaf)
        deltatau = sigma_rg * deltad * n_fill
        currentleaf.data = [n_fill, fm_fill, deltatau]
        refine_leaf!(currentpoint, currentleaf, n_fill, fm_fill, wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
    end
end

"Fills data to all leaves between point1 and point2 (including linewidth) and refines them if necessary. "
function fill_and_refine!(point1, point2, linewidth_normalized, n, fm, wind::WindStruct)
    lw = linewidth_normalized * point1[1]
    point1leaf = findleaf(wind.quadtree, point1)
    currentpoint = copy(point1)
    currentleaf = copy(point1leaf)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    fill_and_refine_linewidth!(currentpoint, lw, n, fm, wind)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    point2leaf = findleaf(wind.quadtree, point2)
    while(currentleaf != point2leaf)
        try
            currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        catch
            println("currentpoint : $currentpoint")
            println("currentleaf: $currentleaf")
            println("point1: $point1")
            println("point2: $point2")
            println("point1leaf: $point1leaf")
            println("point2leaf: $point2leaf")
            throw(DomainError)
        end
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf == previousleaf
            currentpoint .-= 1e-8
            currentleaf = findleaf(wind.quadtree, currentpoint)
        end
        fill_and_refine_linewidth!(currentpoint, lw, n, fm , wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
        previousleaf = copy(currentleaf) 
    end
end

"Fills and refines entire streamline"
function fill_and_refine_line!(line, wind::WindStruct)
    for i in 1:length(line.p.n_hist)-1
        point1 = line.p.u_hist[i, 1:2]
        point2 = line.p.u_hist[i+1, 1:2]
        linewidth_normalized = line.p.line_width / line.p.r_0
        n = line.p.n_hist[i]
        fm = line.p.fm_hist[i]
        fill_and_refine!(point1, point2, linewidth_normalized, n, fm, wind)
    end
end

"Fills and refines data from all streamlines"
function fill_and_refine_all_lines!(wind::WindStruct)
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        fill_and_refine_line!(line, wind)
    end
end

"Erases line density and fm from tree, setting the values to the initial ones."
function erase_line_from_tree!(line_id, wind::WindStruct)
    reinitialize_tree!(wind)
    line = wind.lines[line_id]
    n_vacuum = wind.grids.n_vacuum
    line.p.n_hist .= n_vacuum * ones(Float64, size(line.p.n_hist))
    line.p.fm_hist .= zeros(Float64, size(line.p.fm_hist))
    fill_and_refine_all_lines!(wind)
end

#=
function fill_linewidth!(point, linewidth, density, fm, wind::WindStruct)
    r_min = wind.grids.r_range[1]
    r_max = wind.grids.r_range[end]
    point1 = [max(point[1] - linewidth / 2., r_min) , point[2]] 
    point2 = [min(point[1] + linewidth / 2., r_max) , point[2]] 
    sigmarg = SIGMA_T * wind.bh.R_g
    currentpoint = copy(point1)
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    n_fill = max(density, currentleaf.data[1])
    fm_fill = max(fm, currentleaf.data[2])
    deltad = compute_cell_size(currentleaf)
    deltatau = sigmarg * deltad * n_fill
    currentleaf.data = [n_fill, fm_fill, deltatau]
    while (currentleaf != point2leaf)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf.boundary == previousleaf.boundary
            currentpoint .-= 1e-8
            currentleaf = findleaf(wind.quadtree, currentpoint)
        end
        n_fill = max(density, currentleaf.data[1])
        fm_fill = max(fm, currentleaf.data[2])
        deltad = compute_cell_size(currentleaf)
        deltatau = sigmarg * deltad * n_fill
        currentleaf.data = [n_fill, fm_fill, deltatau]
        previousleaf = copy(currentleaf)
    end
end

function fill_point!(point1, point2, linewidth_normalized, density, fm, wind::WindStruct)
    if point1[1] < point2[1]
        point1, point2 = point2, point1
    end
    currentleaf = findleaf(wind.quadtree, point1)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    linewidth = linewidth_normalized * point1[1]
    fill_linewidth(point1, linewidth, density, fm, wind)
    currentpoint = copy(point1)
    while currentleaf != point2leaf
        #println("point1 : $currentpoint, point2 : $point2, currentleaf : $currentleaf")
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        if currentleaf.boundary == previousleaf.boundary
            #println("currentleaf: $(currentleaf.boundary)")
            currentpoint .-= 1e-8
            #println(currentpoint)
            currentleaf = findleaf(wind.quadtree, currentpoint)
            #println("newleaf: $(currentleaf.boundary)")
        end
        linewidth = linewidth_normalized * currentpoint[1]
        fill_linewidth(currentpoint, linewidth, density, fm, wind)
        previousleaf = copy(currentleaf)
    end
end

function fill_line!(line, wind::WindStruct)
    if length(line.p.n_hist) < 3
        return nothing
    end
    linewidth_normalized = line.p.line_width / line.p.r_0 
    for i in 2:(length(line.p.n_hist) - 1)
        point1 = line.p.u_hist[i, 1:2]
        point2 = line.p.u_hist[i+1, 1:2]
        rho = line.p.n_hist[i]
        fm = line.p.fm_hist[i]
        fill_point(point1, point2, linewidth_normalized, rho, fm, wind)
    end
end

function fill_all_lines!(wind::WindStruct)
    #initialize them all
    for leaf in allleaves(wind.quadtree)
        leaf.data = [wind.grids.n_vacuum, 0., 0.01]
    end
    k=1
    for line in wind.lines
        #println("filling line ... $k")
        k += 1
        fill_line(line, wind)
    end
end
function refine_all(wind::WindStruct)
    fill_all_lines(wind)
    optically_thin = false
    while (!optically_thin)
        println("refining ...")
        optically_thin = true
        for leaf in allleaves(wind.quadtree)
            maxz = vertices(leaf.boundary)[2,2][2]
            refine_condition = (leaf.data[3] > 0.1) && (maxz > wind.z_0) 
            if refine_condition
                optically_thin = false
                split!(leaf)
            end
        end
        fill_all_lines(wind)
    end
end
=#