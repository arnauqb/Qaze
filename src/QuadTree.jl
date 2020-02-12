using RegionTrees: vertices, Cell, findleaf, children, split!
import Base: copy
using StaticArrays: SVector
import RegionTrees: needs_refinement, refine_data, AbstractRefinery, allleaves
export compute_cell_intersection, compute_cell_size, fill_linewidth, fill_point, fill_all_points, copy, needs_refinement, refine_data, refine, refine_mine


function compute_cell_size(cell::Cell)
    r1, z1 = vertices(cell.boundary)[1,1]
    r2, z2 = vertices(cell.boundary)[2,2]
    distance = sqrt((r2-r1)^2 + (z2-z1)^2)
    return distance
end

function copy(cell::Cell)
    newcell = Cell(cell.boundary, cell.data)
    return newcell
end

function compute_cell_intersection(current_point, cell::Cell, initial_point, final_point)
    r, z = current_point
    r_i, z_i = initial_point
    r_f, z_f = final_point
    direction_r = sign(r_f - r_i)
    direction_z = sign(z_f - z_i)
    cell_vertices = vertices(cell)
    if direction_r == 1
        cell_r_boundary = cell_vertices[2,2][1]
    else
        cell_r_boundary = cell_vertices[1,1][1]
    end
    if direction_z == 1
        cell_z_boundary = cell_vertices[2,2][2]
    else
        cell_z_boundary = cell_vertices[1,1][2]
    end
    if r_i == r_f
        lambda_r = Inf
    else
        lambda_r = (cell_r_boundary - r) / (r_f - r_i)
    end
    if z_i == z_f
        lambda_z = Inf
    else
        lambda_z = (cell_z_boundary - z) / (z_f - z_i)
    end
    lambda = min(lambda_r, lambda_z)
    @assert lambda >= 0.0
    intersection = [r, z] + lambda .* [r_f - r_i, z_f - z_i]
    return intersection
end

function fill_linewidth(point, linewidth, density, fm, wind::WindStruct)
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
    #println("pointleaf: $currentleaf")
    #println("deltad: $deltad, deltatau: $deltatau")
    while (currentleaf != point2leaf)
        #println("lwidth")
        #println("point1 : $point1, point2 : $point2")
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
        #println("deltad: $deltad, deltatau: $deltatau")
        currentleaf.data = [n_fill, fm_fill, deltatau]
        previousleaf = copy(currentleaf)
    end
end

function fill_point(point1, point2, linewidth_normalized, density, fm, wind::WindStruct)
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

function fill_all_points(wind::WindStruct)
    #initialize them all
    for leaf in allleaves(wind.quadtree)
        leaf.data = [wind.config["wind"]["n_shielding"], 0., 0.01]
    end
    k=1
    for line in wind.lines
        #println("filling line ... $k")
        k += 1
        for i in 2:(length(line.p.n_hist) - 1)
            point1 = line.p.u_hist[i, 1:2]
            point2 = line.p.u_hist[i+1, 1:2]
            rho = line.p.n_hist[i]
            fm = line.p.fm_hist[i]
            linewidth_normalized = line.p.line_width / line.p.r_0 
            fill_point(point1, point2, linewidth_normalized, rho, fm, wind)
        end
    end
end

function refine_mine(wind::WindStruct)
    fill_all_points(wind)
    optically_thin = false
    while (!optically_thin)
        println("refining ...")
        optically_thin = true
        for leaf in allleaves(wind.quadtree)
            if leaf.data[3] > 0.1
                optically_thin = false
                split!(leaf)
            end
        end
        fill_all_points(wind)
    end
end




function needs_refinement(r::GridRefinery, cell)
    cond = cell.data[3] > r.tau_tolerance
    return cond
end

function refine_data(r::GridRefinery, cell, indices)
    fill_all_points(r.wind)
    return cell.data
end

function adaptivesampling!(root::Cell, refinery::AbstractRefinery)
    refinement_queue = [root]
    refine_function =
        (cell, indices) -> refine_data(refinery, cell, indices)
    while !isempty(refinement_queue)
        cell = pop!(refinement_queue)
        println("cell: $cell")
        println("cell data : $(cell.data)")
        if needs_refinement(refinery, cell)
            println("needs refinement")
            split!(cell, refine_function)
            append!(refinement_queue, children(cell))
        end
    end
    root
end

function refine(cell::Cell, wind::WindStruct)
    r = GridRefinery(0.1, wind)
    adaptivesampling!(cell, r)
end