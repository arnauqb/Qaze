using RegionTrees
using StaticArrays: SVector
using RegionTrees
import Base: copy, ==
export copy,
       compute_cell_intersection,
       compute_cell_size,
       refine_leaf!,
       fill_cell!,
       fill_and_refine_leaf!,
       fill_and_refine_linewidth!,
       fill_and_refine!,
       fill_and_refine_line!,
       fill_and_refine_all_lines!,
       erase_line_from_tree!,
       reinitialize_tree!

"Computes the diagonal of the cell"
function compute_cell_size(cell::Cell)
    r1, z1 = vertices(cell.boundary)[1,1]
    r2, z2 = vertices(cell.boundary)[2,2]
    distance = sqrt((r2-r1)^2 + (z2-z1)^2)
    return distance
end

"Fills a cell with given density and force multiplier information."
function fill_cell!(density, fm, line_id, cell::Cell, wind::WindStruct)
    n_fill = max(density, cell.data[1])
    if wind.radiation.include_fm_tauuv
        fm_fill = max(fm, cell.data[2])
    else
        fm_fill = 0.
    end
    deltad = compute_cell_size(cell)
    deltatau = SIGMA_T * wind.bh.R_g * deltad * n_fill
    cell.data = [n_fill, fm_fill, deltatau, []]
    #if line_id in cell.data[4]
    #    cell.data = [n_fill, fm_fill, deltatau, cell.data[4]]
    #else
    #    cell.data = [n_fill, fm_fill, deltatau, push!(cell.data[4], line_id)]
    #end
end

function fill_all_lines!(wind::WindStruct)
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        for j in 1:size(line.p.u_hist)[1]
            r = line.p.u_hist[j, 1]
            z = line.p.u_hist[j, 2]
            n = line.p.n_hist[j]
            fm = line.p.fm_hist[j]
            cell = findleaf(wind.quadtree, [r,z])
            fill_cell!(n, fm, i, cell, wind)
        end
    end
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
    wind.quadtree = Cell(SVector(0., 0.), SVector(2 * quadtree_max_radius, 2* quadtree_max_height), [n_vacuum, 0., delta_tau_0, []])
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
    newdata = [cell.data[1], 0., cell.data[3] / 2.0, []]
    return newdata
end

"splits leaf until the max cell optical thickness is reached"
function refine_leaf!(currentpoint, leaf, density, fm, line_id, wind)
    optically_thin = false
    while !optically_thin
        optically_thin = true
        #maxz = vertices(leaf.boundary)[2,2][2]
        refine_condition = (leaf.data[3] > 0.1)#wind.config["grids"]["cell_optical_thickness"])# && (maxz > wind.z_0) 
        if refine_condition
            optically_thin = false
            f(cell, child_indices) = split_cell_initialization(cell, child_indices, line_id)
            split!(leaf, f)
            newleaf = findleaf(leaf, currentpoint)
            fill_cell!(density, fm, line_id, newleaf, wind)
            #fill_all_lines!(wind)
            leaf = newleaf
        end
    end
end

function fill_and_refine_leaf!(currentpoint, leaf, density, fm, wind)
    deltatau = SIGMA_T * wind.bh.R_g * density * compute_cell_size(leaf) 
    thicknessmax = wind.config["grids"]["cell_optical_thickness"]
    mincellsize = wind.config["grids"]["minimum_cell_size"] 
    refine_condition = (deltatau > thicknessmax) && (compute_cell_size(leaf) > mincellsize)
    if refine_condition
        N = ceil(log2(deltatau / thicknessmax))
        for i = 1:N
            split!(leaf, split_cell_initialization)
            leaf = findleaf(leaf, currentpoint)
            if compute_cell_size(leaf) < mincellsize
                break
            end
        end
    end
    n_fill = max(density, leaf.data[1])
    leaf.data = [n_fill, 0., deltatau, []]
end

"Fills the line data to all leaves along the width of the streamline at point, and refines them if necessary"
function fill_and_refine_linewidth!(point, linewidth, density, fm, line_id, wind::WindStruct)
    rmin = 0.0
    rmax = 2 * wind.grids.r_max 
    point1 = [max(point[1] - linewidth / 2., rmin) , point[2]] 
    point2 = [min(point[1] + linewidth / 2., rmax) , point[2]] 
    currentpoint = copy(point1)
    currentleaf = findleaf(wind.quadtree, point1)
    point2leaf = findleaf(wind.quadtree, point2)
    #fill_cell!(density, fm, line_id, currentleaf, wind)
    #refine_leaf!(currentpoint, currentleaf, density, fm, line_id, wind)
    fill_and_refine_leaf!(currentpoint, currentleaf, density, fm, wind)
    currentleaf = findleaf(wind.quadtree, currentpoint)
    point2leaf = findleaf(wind.quadtree, point2)
    while (currentleaf != point2leaf)
        currentpoint = compute_cell_intersection(currentpoint, currentleaf, point1, point2)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        #fill_cell!(density, fm, line_id, currentleaf, wind)
        #refine_leaf!(currentpoint, currentleaf, density, fm, line_id, wind)
        fill_and_refine_leaf!(currentpoint, currentleaf, density, fm, wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
    end
end

"Fills data to all leaves between point1 and point2 (including linewidth) and refines them if necessary. "
function fill_and_refine!(point1, point2, linewidth_normalized, n, fm, line_id, wind::WindStruct)
    lw = linewidth_normalized * point1[1]
    point1leaf = findleaf(wind.quadtree, point1)
    currentpoint = copy(point1)
    currentleaf = copy(point1leaf)
    previousleaf = copy(currentleaf)
    point2leaf = findleaf(wind.quadtree, point2)
    fill_and_refine_linewidth!(currentpoint, lw, n, fm, line_id, wind)
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
        fill_and_refine_linewidth!(currentpoint, lw, n, fm, line_id, wind)
        currentleaf = findleaf(wind.quadtree, currentpoint)
        point2leaf = findleaf(wind.quadtree, point2)
        previousleaf = copy(currentleaf) 
    end
end

"Fills and refines entire streamline"
function fill_and_refine_line!(line, line_id, wind::WindStruct)
    for i in 1:length(line.p.n_hist)-1
        point1 = line.p.u_hist[i, 1:2]
        point2 = line.p.u_hist[i+1, 1:2]
        linewidth_normalized = line.p.line_width / line.p.r_0
        n = line.p.n_hist[i]
        fm = line.p.fm_hist[i]
        fill_and_refine!(point1, point2, linewidth_normalized, n, fm, line_id, wind)
    end
end

"Fills and refines data from all streamlines"
function fill_and_refine_all_lines!(wind::WindStruct)
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        line = wind.lines[i]
        fill_and_refine_line!(line, i, wind)
    end
end

"Erases line density and fm from tree, setting the values to the initial ones."
function erase_line_from_tree!(line_id, wind::WindStruct)
    reinitialize_tree!(wind)
    line = wind.lines[line_id]
    n_vacuum = wind.grids.n_vacuum
    line.p.n_hist .= n_vacuum * ones(Float64, size(line.p.n_hist))
    line.p.fm_hist .= zeros(Float64, size(line.p.fm_hist))
    for i in 1:length(wind.lines)
        isassigned(wind.lines, i) || continue
        (i == line_id) && continue
        line = wind.lines[i]
        fill_and_refine_line!(line, i, wind)
    end
    #fill_and_refine_all_lines!(wind)
end